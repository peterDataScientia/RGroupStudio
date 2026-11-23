"""Native Qt desktop regression coverage; no network or installer needed."""
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import time

os.environ["RGROUPSTUDIO_DISABLE_AUTO_CHECK"]="1"
os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS","--enable-unsafe-swiftshader")

from rdkit import Chem
from PySide6.QtCore import Qt, QEventLoop, QTimer
from PySide6.QtWidgets import QApplication
from PySide6.QtWebEngineWidgets import QWebEngineView
from chemistry import enumerate_library
from app_qt import Studio, resource

class QtDesktopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.application=QApplication.instance() or QApplication([])
        cls.gui=Studio()
        cls.gui.show()
        cls.application.processEvents()

    @classmethod
    def tearDownClass(cls):
        cls.gui.close()
        cls.application.processEvents()

    def setUp(self):
        gui=self.gui
        gui.scaffold_input.setText("c1ccccc1")
        gui.draw_scaffold()
        self.application.processEvents()

    def test_embedded_offline_viewer_and_visible_export(self):
        gui=self.gui
        self.assertIsInstance(gui.browser,QWebEngineView)
        self.assertTrue(gui.browser.url().toString().startswith("file:"))
        self.assertTrue(gui.tabs.isVisible())
        self.assertTrue(gui.table.isVisible())
        self.assertGreater(gui.browser.minimumHeight(),100)
        self.assertGreater(gui.enumerate_btn.width(),70)
        self.assertTrue((resource("viewer")/"3Dmol-min.js").is_file())

    def test_scaffold_atom_click_assignment_enumeration(self):
        gui=self.gui
        self.assertEqual(len(gui.allowed),6)
        gui.select_atom(2)
        self.assertIn(2,gui.sites_defined)
        gui.site_box.setCurrentText("2")
        first=gui.fragment_list.item(0)
        first.setSelected(True)
        gui.assign()
        self.assertIn(2,gui.assignments)
        records=enumerate_library(gui.loaded_scaffold,{2:["[*]F"]})
        gui.show_records(records)
        self.assertEqual(gui.table.rowCount(),1)
        self.assertEqual(gui.count_label.text(),"1 compounds")
        self.assertIn("MW",gui.properties.text())
        self.assertIsNotNone(gui.preview2d.pixmap())

    def test_actual_export_buttons_create_correct_formats(self):
        gui=self.gui
        gui.show_records(enumerate_library("c1ccccc1",{0:["[*]F"]}))
        with tempfile.TemporaryDirectory() as folder:
            for name,extension in (("Export CSV","csv"),("Export SMILES","smi"),
                                   ("Export 2D SDF","sdf"),("Export 3D SDF","sdf")):
                path=Path(folder)/(name.replace(" ","_")+"."+extension)
                buttons=[b for b in gui.findChildren(__import__(
                    "PySide6.QtWidgets",fromlist=["QPushButton"]).QPushButton)
                    if b.text()==name]
                self.assertEqual(len(buttons),1)
                with patch("app_qt.QFileDialog.getSaveFileName",
                           return_value=(str(path),extension)):
                    buttons[0].click()
                    for _ in range(120):
                        self.application.processEvents()
                        if path.exists() and not any(w.isRunning() for w in gui.workers):
                            break
                        time.sleep(.025)
                self.assertTrue(path.exists(),name+" did not create file")
                self.assertGreater(path.stat().st_size,10)
                if extension=="sdf":
                    mols=[m for m in Chem.SDMolSupplier(str(path)) if m]
                    self.assertEqual(len(mols),1)

    def test_site_remove_and_3d_selection_contract(self):
        gui=self.gui
        gui.select_atom(3)
        self.assertEqual(gui.assignment_list.count(),1)
        gui.assignment_list.setCurrentRow(0)
        gui.remove_assignment()
        self.assertEqual(gui.assignment_list.count(),0)
        self.assertEqual(gui.assignments,{})
        self.assertEqual(gui.selected_index,None)

    def test_qt_desktop_export_all(self):
        gui=self.gui
        rec=enumerate_library("c1ccccc1",{0:["[*]F","[*]Cl"]})
        gui.show_records(rec)
        with tempfile.TemporaryDirectory() as folder:
            from chemistry import export_sdf
            path=Path(folder)/"two.sdf"
            errors=export_sdf(path,gui.records,three_d=True)
            self.assertFalse(errors)
            self.assertEqual(len([m for m in Chem.SDMolSupplier(str(path)) if m]),2)

    def test_zz_real_webgl_sdf_loading(self):
        """Prove JS library parses a real RDKit 3D SDF inside Qt WebEngine."""
        gui=self.gui
        gui.tabs.setCurrentIndex(1)
        for _ in range(160):
            self.application.processEvents()
            if gui._web_ready:break
            time.sleep(.05)
        self.assertTrue(gui._web_ready,"Offline Qt WebEngine page did not load")
        from chemistry import conformer_3d
        mol=conformer_3d(Chem.MolFromSmiles("CC(C)c1ccc(F)cc1"))
        sdf=Chem.MolToMolBlock(mol)
        import json
        result=[]
        gui.browser.page().runJavaScript(
            "window.loadMolecule("+json.dumps(sdf)+",'RGS_Test')",
            lambda value:result.append(value))
        for _ in range(160):
            self.application.processEvents()
            if result:break
            time.sleep(.05)
        self.assertEqual(result,[True],"3Dmol JS failed to load RDKit SDF")
        health=[]
        gui.browser.page().runJavaScript(
            "window.viewerHealth()",lambda value:health.append(value))
        for _ in range(160):
            self.application.processEvents()
            if health:break
            time.sleep(.05)
        self.assertTrue(health and health[0]["library"])
        self.assertTrue(health[0]["model"])
        self.assertGreater(health[0]["loadedSdfBytes"],100)
        gui.tabs.setCurrentIndex(0)

    def test_local_html_3dmol_contract(self):
        html=(resource("viewer")/"index.html").read_text(encoding="utf-8")
        for text in ('src="3Dmol-min.js"',"window.loadMolecule","viewer.zoomTo","window.viewerHealth"):
            self.assertIn(text,html)
        self.assertNotIn('src="https://',html)

if __name__=="__main__":
    unittest.main()
