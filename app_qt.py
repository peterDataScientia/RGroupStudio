"""RGroupStudio native Qt desktop application with embedded offline 3Dmol.js.

RDKit remains the source of chemical truth. SDF coordinates are exported to
the local Qt WebEngine GLViewer; no molecule data leave the computer.
"""
import io
import json
import os
import subprocess
import sys
import threading
from pathlib import Path

from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from PySide6.QtCore import Qt, QBuffer, QByteArray, QIODevice, QThread, Signal, QUrl
from PySide6.QtGui import QIcon, QPixmap, QImage, QColor
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFrame, QLabel, QPushButton, QLineEdit,
    QSpinBox, QCheckBox, QComboBox, QListWidget, QListWidgetItem, QAbstractItemView,
    QTableWidget, QTableWidgetItem, QSplitter, QScrollArea, QVBoxLayout, QHBoxLayout,
    QGridLayout, QTabWidget, QFileDialog, QMessageBox, QInputDialog, QHeaderView,
    QProgressBar
)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebEngineCore import QWebEngineSettings

from chemistry import (
    FRAGMENTS, parse_scaffold, sites, load_fragments, validate_fragment,
    enumerate_library, molecular_properties, conformer_3d,
    export_csv, export_smi, export_sdf
)
from updater import APP_VERSION, latest_update, get_auto_check, set_auto_check, download_installer

BLUE="#075AA8"
NAVY="#10273F"
BORDER="#B6C7D9"
BACKGROUND="#EAF1F8"

def resource(path):
    return Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent)) / path

def pixmap_from_pil(image):
    stream=io.BytesIO()
    image.save(stream,format="PNG")
    pixmap=QPixmap()
    if not pixmap.loadFromData(stream.getvalue(),"PNG"):
        raise RuntimeError("Could not decode RDKit structure image")
    return pixmap

class ScaffoldPicture(QLabel):
    siteClicked=Signal(int)
    def __init__(self):
        super().__init__("Load a scaffold")
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setMinimumSize(250,205)
        self.setStyleSheet("background:#fff;border:1px solid #9cb2cb;")
        self.points={}
        self._last_molecule=None
        self._current=0
        self._defined=set()
        self._show_h=False
        self._indices=True
        self.setCursor(Qt.CursorShape.CrossCursor)

    def draw(self, mol, current=-1, defined=(), indices=True, show_h=False):
        self._last_molecule=Chem.Mol(mol)
        self._current=current
        self._defined=set(defined)
        self._indices=indices
        self._show_h=show_h
        width=max(245,self.width()-8)
        height=max(175,self.height()-8)
        m=Chem.AddHs(self._last_molecule) if show_h else Chem.Mol(self._last_molecule)
        if indices:
            for atom in m.GetAtoms():
                atom.SetProp("atomNote",str(atom.GetIdx()))
        drawer=rdMolDraw2D.MolDraw2DCairo(width,height)
        marked=list(self._defined)
        if current>=0 and current<m.GetNumAtoms() and current not in marked:
            marked.append(current)
        colors={i:((.04,.67,.52) if i in self._defined else (.95,.61,.18)) for i in marked}
        drawer.DrawMolecule(m,highlightAtoms=marked,highlightAtomColors=colors)
        drawer.FinishDrawing()
        img=QPixmap()
        img.loadFromData(bytes(drawer.GetDrawingText()),"PNG")
        self.setPixmap(img)
        allowed=set(sites(self._last_molecule))
        self.points={}
        for a in m.GetAtoms():
            idx=a.GetIdx()
            parent=idx
            if a.GetSymbol()=="H" and a.GetDegree()==1:
                parent=a.GetNeighbors()[0].GetIdx()
            if parent in allowed:
                xy=drawer.GetDrawCoords(idx)
                self.points[idx]=(xy.x,xy.y,parent)
        self._image_dims=(width,height)

    def mousePressEvent(self,event):
        if self.points:
            # QLabel centers image when its native dimensions are smaller than the widget.
            off_x=(self.width()-self._image_dims[0])/2
            off_y=(self.height()-self._image_dims[1])/2
            x=float(event.position().x())-off_x
            y=float(event.position().y())-off_y
            near=min(self.points.values(),key=lambda p:(p[0]-x)**2+(p[1]-y)**2)
            if (near[0]-x)**2+(near[1]-y)**2 < 30**2:
                self.siteClicked.emit(near[2])
        super().mousePressEvent(event)

    def resizeEvent(self,event):
        super().resizeEvent(event)
        if self._last_molecule is not None:
            self.draw(self._last_molecule,self._current,self._defined,self._indices,self._show_h)

class Worker(QThread):
    completed=Signal(object)
    failed=Signal(str)
    def __init__(self,job,parent=None):
        super().__init__(parent)
        self.job=job
    def run(self):
        try:
            self.completed.emit(self.job())
        except Exception as exc:
            self.failed.emit(str(exc))

class Studio(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RGroupStudio v"+APP_VERSION+" | Molecular Library Designer")
        icon=resource("assets/approved_logo.png")
        if icon.exists():
            self.setWindowIcon(QIcon(str(icon)))
        self.resize(1400,900)
        self.setMinimumSize(1100,720)
        self.scaffold_mol=None
        self.loaded_scaffold=""
        self.allowed=[]
        self.fragments=dict(FRAGMENTS)
        self.sites_defined=set()
        self.assignments={}
        self.current_choices={}
        self.records=[]
        self.selected_index=None
        self.workers=[]
        self._web_ready=False
        self._pending_sdf=None
        self._active_mol3d=None
        self._fullscreen=False
        self._build()
        self.draw_scaffold()
        if get_auto_check() and not os.getenv("RGROUPSTUDIO_DISABLE_AUTO_CHECK"):
            from PySide6.QtCore import QTimer
            QTimer.singleShot(1800,lambda:self.check_update(silent=True))

    def _build(self):
        self.setStyleSheet("""
          QMainWindow,QWidget {font-family:'Segoe UI';font-size:10pt;color:#10273F;background:#EAF1F8}
          QLabel {background:transparent}
          QLineEdit,QSpinBox,QComboBox,QListWidget,QTableWidget {background:#fff;color:#10273F;
            border:1px solid #A7BBD0;padding:5px;selection-background-color:#075AA8;
            selection-color:white}
          QPushButton {background:#D4E4F3;color:#10273F;border:1px solid #95B2CE;
            border-radius:4px;padding:8px 11px;font-weight:600}
          QPushButton:hover {background:#BED7EE}
          QPushButton#primary {background:#075AA8;color:white;border:1px solid #064C8D}
          QPushButton#primary:hover {background:#06467F}
          QHeaderView::section {background:#C4D9EC;color:#173453;font-weight:700;
            border:1px solid #9FB5CD;padding:7px}
          QTabWidget::pane {border:1px solid #9BB5D0;background:#fff}
          QTabBar::tab {background:#D6E5F3;padding:9px 16px;border:1px solid #B0C5DA}
          QTabBar::tab:selected {background:#fff;color:#064C8D;font-weight:700}
        """)
        central=QWidget()
        self.setCentralWidget(central)
        root=QVBoxLayout(central)
        root.setContentsMargins(8,5,8,5)
        root.setSpacing(5)

        header=QFrame()
        header.setStyleSheet("QFrame{background:#183A61;border-radius:4px} QLabel{color:white}")
        h=QHBoxLayout(header)
        h.setContentsMargins(12,6,10,6)
        brand=QLabel("RGroupStudio")
        brand.setStyleSheet("font-size:15pt;font-weight:800;color:white")
        h.addWidget(brand)
        subtitle=QLabel("Molecular Library Designer")
        subtitle.setStyleSheet("color:#D9E8F5")
        h.addWidget(subtitle)
        h.addStretch(1)
        version=QLabel("v"+APP_VERSION)
        h.addWidget(version)
        self.auto_check=QCheckBox("Auto-check updates")
        self.auto_check.setChecked(get_auto_check())
        self.auto_check.setStyleSheet("color:white")
        self.auto_check.toggled.connect(set_auto_check)
        h.addWidget(self.auto_check)
        self.check_btn=self._button("Check for Updates",lambda:self.check_update(),True)
        h.addWidget(self.check_btn)
        root.addWidget(header)

        self.splitter=QSplitter(Qt.Orientation.Horizontal)
        self.splitter.setChildrenCollapsible(False)
        root.addWidget(self.splitter,1)
        left_scroll=QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left=QWidget()
        leftlayout=QVBoxLayout(left)
        leftlayout.setSpacing(8)
        self.splitter.addWidget(left_scroll)
        left_scroll.setWidget(left)
        left.setMinimumWidth(340)

        leftlayout.addWidget(self._heading("1  Scaffold"))
        smilesrow=QHBoxLayout()
        self.scaffold_input=QLineEdit("c1ccccc1")
        self.scaffold_input.returnPressed.connect(self.draw_scaffold)
        smilesrow.addWidget(self.scaffold_input,1)
        smilesrow.addWidget(self._button("Draw",self.draw_scaffold))
        leftlayout.addLayout(smilesrow)
        optrow=QHBoxLayout()
        self.show_indices=QCheckBox("Atom indices")
        self.show_indices.setChecked(True)
        self.show_indices.toggled.connect(self.refresh_scaffold_picture)
        optrow.addWidget(self.show_indices)
        self.show_h=QCheckBox("Show H")
        self.show_h.toggled.connect(self.refresh_scaffold_picture)
        optrow.addWidget(self.show_h)
        optrow.addStretch()
        leftlayout.addLayout(optrow)
        self.scaffold_pic=ScaffoldPicture()
        self.scaffold_pic.siteClicked.connect(self.select_atom)
        leftlayout.addWidget(self.scaffold_pic)
        helper=QLabel("Click a substitutable atom or use the site selector below")
        helper.setWordWrap(True)
        leftlayout.addWidget(helper)

        leftlayout.addWidget(self._heading("2  Attachment sites"))
        srow=QHBoxLayout()
        self.site_box=QComboBox()
        self.site_box.currentIndexChanged.connect(self.refresh_scaffold_picture)
        srow.addWidget(self.site_box,1)
        srow.addWidget(self._button("Define R-site",self.define_site,True))
        leftlayout.addLayout(srow)

        leftlayout.addWidget(self._heading("3  Fragments and assignments"))
        frow=QHBoxLayout()
        frow.addWidget(self._button("Import .smi",self.import_fragments))
        frow.addWidget(self._button("+ Fragment",self.add_fragment))
        frow.addWidget(self._button("Assign R-group",self.assign,True))
        leftlayout.addLayout(frow)
        self.fragment_list=QListWidget()
        self.fragment_list.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.fragment_list.setMinimumHeight(155)
        self.fragment_list.setMaximumHeight(235)
        leftlayout.addWidget(self.fragment_list)
        self.refresh_fragment_list()
        leftlayout.addWidget(self._heading("Defined R-groups"))
        self.assignment_list=QListWidget()
        self.assignment_list.setMinimumHeight(90)
        self.assignment_list.itemSelectionChanged.connect(self.describe_assignment)
        leftlayout.addWidget(self.assignment_list)
        self.assignment_detail=QLabel("Define R1/R2/R3 then select fragments.")
        self.assignment_detail.setWordWrap(True)
        leftlayout.addWidget(self.assignment_detail)
        leftlayout.addWidget(self._button("Remove Selected R-site",self.remove_assignment))
        leftlayout.addStretch(1)

        right=QWidget()
        r=QVBoxLayout(right)
        r.setContentsMargins(3,0,3,0)
        r.setSpacing(5)
        self.splitter.addWidget(right)
        self.splitter.setSizes([425,970])
        controls=QHBoxLayout()
        controls.setSpacing(7)
        self.enumerate_btn=self._button("ENUMERATE ALL COMBINATIONS",self.enumerate_all,True)
        self.current_btn=self._button("Enumerate Current",self.enumerate_current,True)
        controls.addWidget(self.enumerate_btn)
        controls.addWidget(self.current_btn)
        controls.addWidget(QLabel("Max:"))
        self.maximum=QSpinBox()
        self.maximum.setRange(1,50000)
        self.maximum.setValue(2500)
        controls.addWidget(self.maximum)
        controls.addStretch(1)
        r.addLayout(controls)
        counts=QHBoxLayout()
        counts.addWidget(self._heading("Enumerated molecules"))
        self.count_label=QLabel("0 compounds")
        counts.addStretch()
        counts.addWidget(self.count_label)
        r.addLayout(counts)
        self.table=QTableWidget(0,6)
        self.table.setHorizontalHeaderLabels(["ID","Canonical SMILES","MW","cLogP","TPSA","QED"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().hide()
        hdr=self.table.horizontalHeader()
        hdr.setSectionResizeMode(1,QHeaderView.ResizeMode.Stretch)
        for col in (0,2,3,4,5):
            hdr.setSectionResizeMode(col,QHeaderView.ResizeMode.ResizeToContents)
        self.table.itemSelectionChanged.connect(self.select_product)
        self.table.setMinimumHeight(135)
        self.table.setMaximumHeight(255)
        r.addWidget(self.table,0)
        exports=QHBoxLayout()
        self.selected_only=QCheckBox("Selected only")
        exports.addWidget(self.selected_only)
        for label,mode in [("CSV","csv"),("SMILES","smi"),("2D SDF","2d"),("3D SDF","3d")]:
            exports.addWidget(self._button("Export "+label,lambda checked=False,m=mode:self.export(m)))
        exports.addWidget(self._button("View 3D",lambda:self.tabs.setCurrentIndex(1),True))
        exports.addStretch()
        r.addLayout(exports)
        self.preview_heading=self._heading("Molecular preview and properties")
        r.addWidget(self.preview_heading)
        self.tabs=QTabWidget()
        r.addWidget(self.tabs,1)
        page2=QWidget()
        p2=QVBoxLayout(page2)
        self.preview2d=QLabel()
        self.preview2d.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.preview2d.setStyleSheet("background:#fff")
        self.preview2d.setMinimumHeight(160)
        p2.addWidget(self.preview2d,1)
        self.tabs.addTab(page2,"2D Structure")
        page3=QWidget()
        p3=QVBoxLayout(page3)
        p3.setContentsMargins(2,2,2,2)
        toolrow=QHBoxLayout()
        toolrow.addWidget(QLabel("WebGL 3D • Drag rotate • Wheel zoom"))
        toolrow.addStretch()
        self.style_box=QComboBox()
        self.style_box.addItems(["Ball and stick","Sticks","Space filling"])
        self.style_box.currentIndexChanged.connect(self.apply_style)
        toolrow.addWidget(self.style_box)
        self.labels_box=QCheckBox("Element labels")
        self.labels_box.toggled.connect(self.apply_labels)
        toolrow.addWidget(self.labels_box)
        toolrow.addWidget(self._button("Fit",self.fit_view))
        toolrow.addWidget(self._button("Fullscreen",self.toggle_fullscreen))
        p3.addLayout(toolrow)
        self.browser=QWebEngineView()
        self.browser.setMinimumHeight(180)
        settings=self.browser.settings()
        settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls,False)
        settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls,True)
        self.browser.loadFinished.connect(self._browser_loaded)
        p3.addWidget(self.browser,1)
        self.tabs.addTab(page3,"3D Molecule")
        self.tabs.currentChanged.connect(self._tab_changed)
        self.properties=QLabel("Load a scaffold to begin.")
        self.properties.setWordWrap(True)
        self.properties.setStyleSheet("background:white;padding:7px;font-weight:600")
        r.addWidget(self.properties)
        self.status=QLabel("Ready")
        self.status.setStyleSheet("padding:5px;color:#244766;background:#D7E6F4")
        root.addWidget(self.status)
        index=resource("viewer/index.html")
        if index.exists():
            self.browser.setUrl(QUrl.fromLocalFile(str(index.resolve())))
        else:
            self.status.setText("Offline molecular viewer assets missing: "+str(index))

    def _heading(self,text):
        label=QLabel(text)
        label.setStyleSheet("font-weight:750;color:#0B4B81;font-size:11pt;padding:2px 0")
        return label

    def _button(self,text,action,primary=False):
        b=QPushButton(text)
        if primary:
            b.setObjectName("primary")
        b.clicked.connect(action)
        return b

    def _message(self,title,message,warning=False):
        if warning:
            QMessageBox.warning(self,title,str(message))
        else:
            QMessageBox.information(self,title,str(message))

    def draw_scaffold(self):
        try:
            mol=parse_scaffold(self.scaffold_input.text())
            self.scaffold_mol=mol
            self.loaded_scaffold=self.scaffold_input.text().strip()
            self.allowed=sites(mol)
            self.site_box.blockSignals(True)
            self.site_box.clear()
            self.site_box.addItems([str(x) for x in self.allowed])
            self.site_box.blockSignals(False)
            self.sites_defined.clear()
            self.assignments.clear()
            self.current_choices.clear()
            self.refresh_assignments()
            self.clear_results()
            self.refresh_scaffold_picture()
            self.status.setText(f"Scaffold loaded — {len(self.allowed)} available sites")
        except Exception as exc:
            self._message("Scaffold error",exc,True)

    def refresh_scaffold_picture(self):
        if self.scaffold_mol is not None:
            self.scaffold_pic.draw(self.scaffold_mol,
                                   int(self.site_box.currentText()) if self.site_box.currentText().isdigit() else -1,
                                   self.sites_defined,self.show_indices.isChecked(),self.show_h.isChecked())

    def select_atom(self,index):
        self.site_box.setCurrentText(str(index))
        self.define_site()

    def define_site(self):
        try:
            site=int(self.site_box.currentText())
            if site not in self.allowed: raise ValueError("Selected atom is not substitutable")
            self.sites_defined.add(site)
            self.refresh_assignments()
            self.refresh_scaffold_picture()
            self.status.setText(f"R-site at atom {site} defined")
        except Exception as exc:
            self._message("R-site error",exc,True)

    def refresh_fragment_list(self):
        self.fragment_list.clear()
        for name,smiles in self.fragments.items():
            item=QListWidgetItem(name+"  |  "+smiles)
            item.setData(Qt.ItemDataRole.UserRole,name)
            self.fragment_list.addItem(item)

    def add_fragment(self):
        name,ok=QInputDialog.getText(self,"Fragment name","Name:")
        if not ok or not name.strip():return
        smiles,ok=QInputDialog.getText(self,"Fragment SMILES","SMILES with one [*] attachment point:")
        if not ok:return
        try:
            validate_fragment(smiles.strip())
            self.fragments[name.strip()]=smiles.strip()
            self.refresh_fragment_list()
            self.clear_results()
        except Exception as exc:
            self._message("Invalid fragment",exc,True)

    def import_fragments(self):
        path,_=QFileDialog.getOpenFileName(self,"Import fragment SMILES","","SMILES (*.smi *.smiles *.txt)")
        if not path:return
        try:
            new=load_fragments(path)
            self.fragments.update(new)
            self.refresh_fragment_list()
            self.clear_results()
            self.status.setText(f"Imported {len(new)} fragments")
        except Exception as exc:
            self._message("Fragment import failed",exc,True)

    def assign(self):
        try:
            site=int(self.site_box.currentText())
            if site not in self.allowed:raise ValueError("Select a valid atom first")
            selected=[x.data(Qt.ItemDataRole.UserRole) for x in self.fragment_list.selectedItems()]
            if not selected:raise ValueError("Select at least one fragment in the library")
            self.sites_defined.add(site)
            self.assignments[site]=selected
            self.current_choices[site]=selected[0]
            self.refresh_assignments()
            self.refresh_scaffold_picture()
            self.clear_results()
        except Exception as exc:
            self._message("R-group assignment failed",exc,True)

    def refresh_assignments(self):
        self.assignment_list.clear()
        for n,site in enumerate(sorted(self.sites_defined),1):
            count=len(self.assignments.get(site,[]))
            name="R"+str(n)+" · Atom "+str(site)
            text=name+" · "+(str(count)+" fragments" if count else "No fragments")
            item=QListWidgetItem(text)
            item.setData(Qt.ItemDataRole.UserRole,site)
            self.assignment_list.addItem(item)

    def describe_assignment(self):
        rows=self.assignment_list.selectedItems()
        if not rows:return
        site=rows[0].data(Qt.ItemDataRole.UserRole)
        names=self.assignments.get(site,[])
        self.assignment_detail.setText(
            "Atom "+str(site)+": "+(", ".join(names) if names else "No fragments assigned"))

    def remove_assignment(self):
        rows=self.assignment_list.selectedItems()
        if not rows:return
        site=rows[0].data(Qt.ItemDataRole.UserRole)
        self.sites_defined.discard(site)
        self.assignments.pop(site,None)
        self.current_choices.pop(site,None)
        self.refresh_assignments()
        self.refresh_scaffold_picture()
        self.clear_results()

    def _start(self,job,on_done,message):
        worker=Worker(job,self)
        self.workers.append(worker)
        self.enumerate_btn.setEnabled(False)
        self.current_btn.setEnabled(False)
        self.status.setText(message)
        def done(value):
            self.enumerate_btn.setEnabled(True)
            self.current_btn.setEnabled(True)
            on_done(value)
        def failed(detail):
            self.enumerate_btn.setEnabled(True)
            self.current_btn.setEnabled(True)
            self.status.setText("Failed: "+detail)
            self._message("Operation failed",detail,True)
        worker.completed.connect(done)
        worker.failed.connect(failed)
        worker.finished.connect(worker.deleteLater)
        worker.finished.connect(lambda w=worker:self.workers.remove(w) if w in self.workers else None)
        worker.start()

    def enumerate_all(self):
        options={site:[self.fragments[name] for name in names]
                 for site,names in self.assignments.items()}
        self._start(lambda:enumerate_library(self.loaded_scaffold,options,self.maximum.value()),
                    self.show_records,"Enumerating all R-group combinations...")

    def enumerate_current(self):
        options={site:[self.fragments[self.current_choices[site]]] for site in self.assignments}
        self._start(lambda:enumerate_library(self.loaded_scaffold,options,self.maximum.value()),
                    self.show_records,"Enumerating current R-group assignment...")

    def clear_results(self):
        self.records=[]
        self.selected_index=None
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        self.table.blockSignals(False)
        self.count_label.setText("0 compounds")
        self._active_mol3d=None
        self._pending_sdf=None
        if self._web_ready:
            self._js("if(window.viewerHealth){window.loadMolecule('', '');}")
        if self.scaffold_mol is not None:
            self.render_2d(self.scaffold_mol,"PARENT SCAFFOLD")
        else:
            self.preview2d.clear()

    def show_records(self,records):
        self.records=records
        self.table.blockSignals(True)
        self.table.setRowCount(len(records))
        for i,(name,mol) in enumerate(records):
            p=molecular_properties(mol)
            vals=[name,Chem.MolToSmiles(mol),str(p["MW"]),str(p["cLogP"]),str(p["TPSA"]),str(p["QED"])]
            for j,value in enumerate(vals):
                cell=QTableWidgetItem(value)
                self.table.setItem(i,j,cell)
        self.table.blockSignals(False)
        self.count_label.setText(str(len(records))+" compounds")
        if records:
            self.table.selectRow(0)
            self.select_product()
        self.status.setText("Enumerated "+str(len(records))+" unique molecules")

    def render_2d(self,mol,name):
        displayed=Chem.Mol(mol)
        rdDepictor.Compute2DCoords(displayed)
        pic=pixmap_from_pil(Draw.MolToImage(displayed,size=(750,440)))
        bounds=self.preview2d.size()
        self.preview2d.setPixmap(pic.scaled(max(250,bounds.width()-15),
                                            max(175,bounds.height()-10),
                                            Qt.AspectRatioMode.KeepAspectRatio,
                                            Qt.TransformationMode.SmoothTransformation))
        p=molecular_properties(mol)
        self.properties.setText(
            name+"    |    MW "+str(p["MW"])+"    cLogP "+str(p["cLogP"])+
            "    TPSA "+str(p["TPSA"])+"    QED "+str(p["QED"])+
            "    HBD "+str(p["HBD"])+"    HBA "+str(p["HBA"]))

    def select_product(self):
        indices=self.table.selectionModel().selectedRows()
        if not indices or not self.records:return
        idx=indices[0].row()
        if idx>=len(self.records):return
        self.selected_index=idx
        name,mol=self.records[idx]
        self.render_2d(mol,name)
        self._active_mol3d=None
        if self.tabs.currentIndex()==1:
            self.load_3d()

    def _tab_changed(self,index):
        if index==1:
            self.load_3d()

    def _browser_loaded(self,ok):
        self._web_ready=bool(ok)
        if ok:
            self.status.setText("Offline 3Dmol.js viewer loaded")
            if self._pending_sdf:
                info=self._pending_sdf
                self._pending_sdf=None
                self._show_sdf(*info)
        else:
            self.status.setText("3D viewer failed to load; check local viewer assets")

    def _js(self,expression,callback=None):
        if self._web_ready:
            self.browser.page().runJavaScript(expression,callback)

    def _show_sdf(self,name,sdf,optimization=""):
        if not self._web_ready:
            self._pending_sdf=(name,sdf,optimization)
            return
        self._js("window.loadMolecule("+json.dumps(sdf)+","+json.dumps(name)+")",
                 lambda ok:self.status.setText(
                     ("3D model loaded: "+name+" | "+optimization) if ok else "3D renderer error; check viewer panel"))
        self.apply_style()
        self.apply_labels()

    def load_3d(self):
        if self.selected_index is None or not self.records:
            self.status.setText("Enumerate and select a molecule to inspect its 3D geometry")
            return
        idx=self.selected_index
        if self._active_mol3d==idx:return
        self._active_mol3d=idx
        name,mol=self.records[idx]
        def make_sdf():
            conformer=conformer_3d(mol)
            return (name,Chem.MolToMolBlock(conformer),conformer.GetProp("_RGS_OptimizationStatus"))
        def finished(info):
            if self.selected_index!=idx:return
            self._show_sdf(*info)
        self._start(make_sdf,finished,"Generating 3D conformer with RDKit ETKDG...")
    def apply_style(self):
        mode=["ballstick","stick","sphere"][self.style_box.currentIndex()]
        self._js("window.setDisplayStyle("+json.dumps(mode)+")")
    def apply_labels(self):
        self._js("window.setAtomLabels("+("true" if self.labels_box.isChecked() else "false")+")")
    def fit_view(self):
        self._js("window.fitMolecule()")
    def toggle_fullscreen(self):
        self._fullscreen=not self._fullscreen
        if self._fullscreen:
            self._stored_split=self.splitter.sizes()
            self.splitter.widget(0).hide()
            self.table.hide()
            self.tabs.setCurrentIndex(1)
        else:
            self.splitter.widget(0).show()
            self.table.show()
            self.splitter.setSizes(getattr(self,"_stored_split",[425,970]))
        self.fit_view()

    def export(self,kind):
        if not self.records:
            self._message("No molecules","Enumerate products before exporting")
            return
        records=self.records
        if self.selected_only.isChecked():
            if self.selected_index is None:
                self._message("Select molecule","Select a molecule in results table")
                return
            records=[self.records[self.selected_index]]
        ext="sdf" if kind in ("2d","3d") else kind
        path,_=QFileDialog.getSaveFileName(self,"Export molecular library",
                                           "RGroupStudio."+ext,"*."+ext)
        if not path:return
        def output():
            if kind=="csv":export_csv(path,records);return []
            if kind=="smi":export_smi(path,records);return []
            return export_sdf(path,records,three_d=(kind=="3d"))
        def finished(errors):
            count=len(records)-len(errors)
            self.status.setText(f"Exported {count}/{len(records)} molecules to "+path)
            if errors:self._message("Partial export","\n".join(errors[:8]),True)
        self._start(output,finished,"Exporting "+str(len(records))+" structures...")

    def check_update(self,silent=False):
        self.check_btn.setEnabled(False)
        def done(info):
            self.check_btn.setEnabled(True)
            if info is None:
                self.status.setText("RGroupStudio v"+APP_VERSION+" is up to date")
                if not silent:self._message("Updates","You have the latest version.")
                return
            choice=QMessageBox.question(self,"Update available",
                "RGroupStudio v"+info["version"]+" is available.\nDownload the verified installer?")
            if choice==QMessageBox.StandardButton.Yes:
                self._download_update(info)
        worker=Worker(latest_update,self)
        self.workers.append(worker)
        worker.completed.connect(done)
        worker.failed.connect(lambda e:(self.check_btn.setEnabled(True),self.status.setText("Update check failed: "+e)))
        worker.finished.connect(worker.deleteLater)
        worker.start()

    def _download_update(self,info):
        def done(path):
            choice=QMessageBox.question(self,"Install update",
                "Installer checksum verified. Close RGroupStudio and install v"+
                info["version"]+" now?")
            if choice==QMessageBox.StandardButton.Yes:
                subprocess.Popen([str(path),"/SILENT","/SUPPRESSMSGBOXES","/NORESTART",
                                   "/CLOSEAPPLICATIONS","/RESTARTAPPLICATIONS"],
                                  cwd=str(Path(path).parent))
                QApplication.quit()
        job=Worker(lambda:download_installer(info),self)
        self.workers.append(job)
        job.completed.connect(done)
        job.failed.connect(lambda e:self._message("Update failed",e,True))
        job.finished.connect(job.deleteLater)
        job.start()

def self_test():
    """Chemistry and local offline viewer contract, no GUI required."""
    from chemistry import attach
    records=enumerate_library("c1ccccc1",{0:["[*]F"],2:["[*]OC"]})
    assert len(records)==1
    mol=conformer_3d(records[0][1])
    sdf=Chem.MolToMolBlock(mol)
    assert "M  END" in sdf and mol.GetNumConformers()==1
    html=resource("viewer/index.html")
    js=resource("viewer/3Dmol-min.js")
    assert html.is_file(),f"Offline viewer HTML missing: {html}"
    assert js.is_file() and js.stat().st_size>100000,f"Offline 3Dmol bundle missing: {js}"
    return 0

def main():
    if "--self-test" in sys.argv:
        return self_test()
    QApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
    app=QApplication(sys.argv)
    window=Studio()
    window.show()
    return app.exec()

if __name__=="__main__":
    sys.exit(main())
