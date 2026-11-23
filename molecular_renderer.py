"""Antialiased, camera-fitted ball-and-stick rendering from RDKit 3D coordinates.

No Matplotlib axis or second window. This module uses orthographic 3D projection,
depth-sorted shaded atoms and cylindrical bonds. Drawing does not modify chemistry.
"""
import math
from PIL import Image, ImageDraw
import numpy as np

CPK = {
    "H": (220, 227, 236), "C": (66, 82, 105), "N": (42, 93, 224),
    "O": (231, 65, 74), "F": (39, 173, 111), "Cl": (33, 163, 96),
    "Br": (153, 74, 46), "I": (111, 58, 179), "S": (226, 184, 32),
    "P": (235, 127, 37), "B": (232, 147, 132), "Si": (205, 173, 152),
}

def orientation(points, yaw=35, pitch=22):
    """Center at geometric median of molecule, then rotate with a virtual trackball."""
    theta, phi = math.radians(yaw), math.radians(pitch)
    ry = np.array([[math.cos(theta), 0, math.sin(theta)],
                   [0, 1, 0],
                   [-math.sin(theta), 0, math.cos(theta)]])
    rx = np.array([[1, 0, 0],
                   [0, math.cos(phi), -math.sin(phi)],
                   [0, math.sin(phi), math.cos(phi)]])
    return np.asarray(points, dtype=float) @ (rx @ ry).T

def _blend(rgb, toward=(255, 255, 255), amount=0.25):
    return tuple(max(0, min(255, int(round(a*(1-amount)+b*amount)))) for a,b in zip(rgb, toward))

def _sphere(draw, cx, cy, radius, color):
    """Layered radial shading: a readable sphere without a wireframe 'blob'."""
    r = max(4, int(radius))
    # Solid dark rim and a subtle shadow at the lower right.
    draw.ellipse((cx-r+2, cy-r+3, cx+r+2, cy+r+3), fill=(192, 202, 216))
    draw.ellipse((cx-r, cy-r, cx+r, cy+r), fill=_blend(color, (12,19,33), .37))
    for fraction, lightness in ((.88,.05),(.72,.13),(.56,.24),(.40,.31),(.23,.42)):
        rr=r*fraction
        dx=-r*.20
        dy=-r*.24
        draw.ellipse((cx+dx-rr,cy+dy-rr,cx+dx+rr,cy+dy+rr),
                     fill=_blend(color, amount=lightness))
    hh=r*.16
    draw.ellipse((cx-r*.44-hh,cy-r*.49-hh,cx-r*.44+hh,cy-r*.49+hh),
                 fill=_blend(color, amount=.69))

def _bond(draw, p, q, r, color, width):
    """Outlined cylindrical bond with a slender lit highlight."""
    px,py=p
    qx,qy=q
    radius=max(3,int(round(width)))
    draw.line([(px+1,py+1),(qx+1,qy+1)],fill=(170,184,200),width=radius+5)
    draw.line([p,q], fill=_blend(color,(11,29,48),.16),width=radius+2)
    draw.line([(px-1,py-1),(qx-1,qy-1)],
              fill=_blend(color,amount=.43),width=max(1,radius//3))

def render_molecule(mol, width=720, height=450, yaw=35, pitch=22,
                    zoom=1.0, show_hydrogens=False, supersampling=2):
    """Return PIL RGB image; bonds, atoms, and depth ordered in one coordinate system.

    The viewing scale is derived from the XY projected extent of a rotated 3D
    conformer, not the original axes bounding box.
    """
    if mol.GetNumConformers()==0:
        raise ValueError("3D conformer required")
    W,H=max(180,int(width)), max(140,int(height))
    zoom=min(3.0,max(.45,float(zoom)))
    factor=max(1,min(3,int(supersampling)))
    w,h=W*factor,H*factor
    atoms=[a for a in mol.GetAtoms() if show_hydrogens or a.GetAtomicNum()!=1]
    if not atoms:
        raise ValueError("No visible atoms in conformer")
    ids=[a.GetIdx() for a in atoms]
    coords=mol.GetConformer().GetPositions().astype(float)
    # Center the molecule and give it an oblique view; no 2D projection flattened into a page.
    center=coords[ids].mean(axis=0)
    rotated=orientation(coords-center,yaw,pitch)
    pts=rotated[ids,:2]
    low,high=pts.min(axis=0),pts.max(axis=0)
    span=np.maximum(high-low,.75)
    # Reserve space for spheres, and do not let an extreme aspect ratio create huge atoms.
    fit=min(w*.79/(span[0]+1.4),h*.78/(span[1]+1.4))
    unit=fit*zoom
    screen={i:(w/2+rotated[i,0]*unit,h/2-rotated[i,1]*unit) for i in ids}
    image=Image.new("RGB",(w,h),(250,252,255))
    draw=ImageDraw.Draw(image)
    # Soft inset frame is outside the actual molecular geometry.
    draw.rounded_rectangle((4,4,w-5,h-5),radius=13*factor,outline=(227,234,243),width=2*factor)
    visible=set(ids)
    radius=min(21*factor,max(8*factor,unit*.19))
    bw=max(4*factor,min(12*factor,int(radius*.68)))
    # Draw back-to-front objects so that front atoms obscure back bonds naturally.
    primitives=[]
    for b in mol.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        if i not in visible or j not in visible:
            continue
        delta=np.array(screen[j])-np.array(screen[i])
        length=float(np.linalg.norm(delta))
        if length<1e-6:
            continue
        tangent=delta/length
        side=np.array([-tangent[1],tangent[0]])
        start=np.array(screen[i])+tangent*radius*.45
        finish=np.array(screen[j])-tangent*radius*.45
        order=2 if b.GetIsAromatic() else int(b.GetBondTypeAsDouble())
        n_lines=2 if order==2 else 3 if order>=3 else 1
        offsets=[0] if n_lines==1 else [-radius*.25,radius*.25] if n_lines==2 else [-radius*.35,0,radius*.35]
        depth=float((rotated[i,2]+rotated[j,2])/2)
        primitives.append((depth,0,("bond",i,j,start,finish,side,offsets)))
    for a in atoms:
        i=a.GetIdx()
        primitives.append((float(rotated[i,2]),1,("atom",i)))
    for _,_,item in sorted(primitives,key=lambda p:(p[0],p[1])):
        if item[0]=="atom":
            i=item[1]
            atom=mol.GetAtomWithIdx(i)
            _sphere(draw,*screen[i],radius,CPK.get(atom.GetSymbol(),(148,113,182)))
        else:
            _,i,j,start,finish,side,offsets=item
            ci=CPK.get(mol.GetAtomWithIdx(i).GetSymbol(),(134,143,161))
            cj=CPK.get(mol.GetAtomWithIdx(j).GetSymbol(),(134,143,161))
            mid=(start+finish)/2
            for offset in offsets:
                off=side*offset
                _bond(draw,tuple(start+off),tuple(mid+off),0,ci,max(3*factor,bw//len(offsets)))
                _bond(draw,tuple(mid+off),tuple(finish+off),0,cj,max(3*factor,bw//len(offsets)))
    if factor>1:
        image=image.resize((W,H),Image.Resampling.LANCZOS)
    return image
