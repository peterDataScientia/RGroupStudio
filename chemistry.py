"""RDKit molecular chemistry engine for scaffold R-group enumeration."""
from itertools import product
from io import BytesIO, StringIO
import csv
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors, QED

FRAGMENTS = {
    "Fluoro": "[*]F", "Chloro": "[*]Cl", "Bromo": "[*]Br",
    "Methyl": "[*]C", "Ethyl": "[*]CC", "Isopropyl": "[*]C(C)C",
    "Methoxy": "[*]OC", "Ethoxy": "[*]OCC", "Hydroxy": "[*]O",
    "Amino": "[*]N", "Cyano": "[*]C#N", "Trifluoromethyl": "[*]C(F)(F)F",
    "Carboxyl": "[*]C(=O)O", "Acetyl": "[*]C(C)=O",
    "Phenyl": "[*]c1ccccc1", "Pyridyl": "[*]c1ccncc1"
}

def parse_scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError("Invalid scaffold SMILES")
    return mol

def sites(mol):
    """Heavy atom indices with at least one implicit or explicit hydrogen."""
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.GetAtomicNum() > 1 and a.GetTotalNumHs() > 0]

def load_fragments(path):
    results = {}
    with open(path, encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            smiles = fields[0]
            name = " ".join(fields[1:]) if len(fields) > 1 else "Fragment_%04d" % number
            validate_fragment(smiles)
            results[name] = smiles
    return results

def validate_fragment(smiles):
    frag = Chem.MolFromSmiles(smiles)
    if frag is None:
        raise ValueError("Invalid fragment SMILES: " + smiles)
    atoms = [a for a in frag.GetAtoms() if a.GetAtomicNum() == 0]
    if len(atoms) != 1 or atoms[0].GetDegree() != 1:
        raise ValueError("Fragment must contain exactly one terminal '*' attachment point: " + smiles)
    return frag

def attach(mol, index, fragment_smiles):
    if index >= mol.GetNumAtoms() or index < 0:
        raise ValueError("Attachment site index out of range")
    atom = mol.GetAtomWithIdx(index)
    if atom.GetAtomicNum() == 1 or atom.GetTotalNumHs() < 1:
        raise ValueError("No available hydrogen at selected site " + str(index))
    frag = validate_fragment(fragment_smiles)
    dummy = next(a for a in frag.GetAtoms() if a.GetAtomicNum() == 0)
    neighbor = dummy.GetNeighbors()[0]
    bond = frag.GetBondBetweenAtoms(dummy.GetIdx(), neighbor.GetIdx())
    if bond.GetBondType() != Chem.BondType.SINGLE:
        raise ValueError("Attachment dummy must use a single bond")
    combo = Chem.CombineMols(mol, frag)
    edit = Chem.RWMol(combo)
    offset = mol.GetNumAtoms()
    edit.AddBond(index, offset + neighbor.GetIdx(), Chem.BondType.SINGLE)
    edit.RemoveAtom(offset + dummy.GetIdx())
    result = edit.GetMol()
    Chem.SanitizeMol(result)
    return result

def enumerate_library(scaffold_smiles, assignments, limit=2500):
    """assignments: {0:['[*]F', '[*]Cl'], 2:['[*]C']}."""
    scaffold = parse_scaffold(scaffold_smiles)
    allowed = set(sites(scaffold))
    for index, options in assignments.items():
        if index not in allowed or not options:
            raise ValueError("Invalid R-site or empty fragment selection: " + str(index))
        for fragment in options:
            validate_fragment(fragment)
    indexes = sorted(assignments)
    if not indexes:
        return [("Parent", scaffold)]
    requested = 1
    for i in indexes:
        requested *= len(assignments[i])
    if requested > limit:
        raise ValueError("Library has %d combinations; maximum is %d" % (requested, limit))
    unique = {}
    # Maintain original atom indices while modifying separate sites.
    # Substituent fragments are appended, so the scaffold indices remain stable.
    for n, options in enumerate(product(*(assignments[i] for i in indexes)), 1):
        mol = Chem.Mol(scaffold)
        for i, frag in zip(indexes, options):
            mol = attach(mol, i, frag)
        key = Chem.MolToSmiles(mol, canonical=True)
        unique.setdefault(key, ("RGS_%04d" % n, mol))
    return list(unique.values())

def molecular_properties(mol):
    return {
        "MW": round(Descriptors.MolWt(mol), 2),
        "cLogP": round(Crippen.MolLogP(mol), 2),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "RotB": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "QED": round(QED.qed(mol), 3),
    }

def conformer_3d(mol, seed=2024):
    mol = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) < 0:
        raise ValueError("3D conformer embedding failed")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=250)
    except Exception:
        pass
    return mol

def export_sdf(path, records, three_d=False):
    writer = Chem.SDWriter(str(path))
    errors = []
    for name, mol in records:
        try:
            m = conformer_3d(mol) if three_d else Chem.Mol(mol)
            if not three_d:
                AllChem.Compute2DCoords(m)
            m.SetProp("_Name", name)
            for k, v in molecular_properties(mol).items():
                m.SetProp(k, str(v))
            writer.write(m)
        except Exception as exc:
            errors.append(name + ": " + str(exc))
    writer.close()
    return errors

def export_csv(path, records):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ID", "SMILES", "MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "QED"])
        for name, mol in records:
            p = molecular_properties(mol)
            writer.writerow([name, Chem.MolToSmiles(mol)] + list(p.values()))

def export_smi(path, records):
    with open(path, "w", encoding="utf-8") as handle:
        for name, mol in records:
            handle.write(Chem.MolToSmiles(mol) + " " + name + "\n")
