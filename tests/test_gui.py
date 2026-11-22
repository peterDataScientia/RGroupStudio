"""Windows desktop GUI end-to-end exercise with real Tk widgets and real file outputs."""
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from rdkit import Chem
from app import Studio

class GUIIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gui = Studio()
        cls.gui.update()

    @classmethod
    def tearDownClass(cls):
        cls.gui.destroy()

    def setUp(self):
        self.gui.scaffold.set("c1ccccc1")
        self.gui.refresh_scaffold()
        self.gui.update()

    def test_scaffold_change_resets_results_and_assignments(self):
        gui=self.gui
        gui.site_choice.set("0")
        gui.fragment_list.selection_clear(0,"end")
        gui.fragment_list.selection_set(0)
        gui.assign()
        gui.enumerate()
        self.assertGreater(len(gui.records),0)
        gui.show_hydrogens.set(True)\n        gui.draw_scaffold()\n        self.assertIsNotNone(gui.scaffold_display.image)\n        gui.scaffold.set("CCO")\n        gui.refresh_scaffold()
        self.assertEqual(len(gui.records),0)
        self.assertIsNone(gui.selected)
        self.assertEqual(len(gui.assigned),0)
        self.assertEqual(gui.results_count.cget("text"),"0 compounds")

    def test_assign_multiple_sites_enumerate_select_view_export(self):
        gui=self.gui
        # Assign 2 options at position 0 and one at position 2.
        gui.site_choice.set("0")
        gui.fragment_list.selection_clear(0,"end")
        gui.fragment_list.selection_set(0,1)
        gui.assign()
        gui.site_choice.set("2")
        gui.fragment_list.selection_clear(0,"end")
        gui.fragment_list.selection_set(3)
        gui.assign()
        self.assertEqual(len(gui.assigned),2)
        gui.enumerate_current()
        self.assertEqual(len(gui.records),1)
        gui.enumerate()
        self.assertEqual(len(gui.records),2)
        self.assertEqual(len(gui.table.get_children()),2)
        self.assertIsNotNone(gui.preview_image.image)
        self.assertIn("MW",gui.properties.cget("text"))
        gui.table.selection_set("1")
        gui.choose_record()
        self.assertEqual(gui.selected,1)
        with tempfile.TemporaryDirectory() as tmp:
            for kind,extension in (("csv","csv"),("smi","smi"),("2d","sdf"),("3d","sdf")):
                path=Path(tmp)/("out_"+kind+"."+extension)
                with patch("app.filedialog.asksaveasfilename",return_value=str(path)):
                    gui.export(kind)
                self.assertTrue(path.exists(),f"{kind} export missing")
                self.assertGreater(path.stat().st_size,30)
                self.assertEqual(gui.last_export_errors,[])
                if extension=="sdf":
                    result=[m for m in Chem.SDMolSupplier(str(path)) if m]
                    self.assertEqual(len(result),2)
                    self.assertTrue(all(m.GetNumConformers()==1 for m in result))
                    if kind=="3d":
                        self.assertTrue(any(abs(m.GetConformer().GetAtomPosition(0).z)>1e-5 for m in result))
            # Verify selected-only export writes exactly one product.\n            gui.export_selection_only.set(True)\n            selected_path=Path(tmp)/"selected.sdf"\n            with patch("app.filedialog.asksaveasfilename",return_value=str(selected_path)):\n                gui.export("3d")\n            selected_mols=[m for m in Chem.SDMolSupplier(str(selected_path)) if m]\n            self.assertEqual(len(selected_mols),1)\n            gui.export_selection_only.set(False)\n            # Import another fragment library using the actual handler.
            infile=Path(tmp)/"fragments.smi"
            infile.write_text("[*]N CustomAmino\n[*]OC CustomMethoxy\n",encoding="utf-8")
            with patch("app.filedialog.askopenfilename",return_value=str(infile)):
                gui.import_fragments()
            self.assertIn("CustomAmino",gui.fragments)
            self.assertGreaterEqual(gui.fragment_list.size(),2)
            self.assertFalse(gui.records) # prevent stale products
        # 3D preview opens without an exception and can be destroyed.
        existing=set(gui.winfo_children())
        with patch("app.messagebox.showerror") as errors:
            gui.enumerate()
            gui.view3d()
            gui.update_idletasks()
            self.assertFalse(errors.called, str(errors.call_args))
        new_windows=[w for w in gui.winfo_children() if w not in existing and w.winfo_class()=="Toplevel"]
        self.assertTrue(new_windows,"View 3D failed to create a window")
        for w in new_windows:
            w.destroy()

    def test_remove_assignment_and_invalid_input(self):
        gui=self.gui
        gui.site_choice.set("1")
        gui.fragment_list.selection_clear(0,"end")
        gui.fragment_list.selection_set(1)
        gui.assign()
        gui.assignment_list.selection_set(0)
        gui.remove_assignment()
        self.assertEqual(gui.assigned,{})
        # Changing text without clicking Draw cannot silently shift atom indices.
        gui.site_choice.set("0")
        gui.fragment_list.selection_clear(0,"end")
        gui.fragment_list.selection_set(0)
        gui.assign()
        gui.scaffold.set("CO")
        gui.enumerate()
        self.assertEqual(len(gui.records),1)
        self.assertEqual(Chem.MolToSmiles(gui.records[0][1]),"Fc1ccccc1")

if __name__=="__main__":
    unittest.main()
