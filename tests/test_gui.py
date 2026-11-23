"""Windows desktop GUI end-to-end exercise with real Tk widgets and real file outputs."""
import os
os.environ["RGROUPSTUDIO_DISABLE_AUTO_CHECK"] = "1"
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

    def test_v030_atom_click_and_complete_preview(self):
        from types import SimpleNamespace
        gui = self.gui
        gui.geometry("1100x690")
        gui.update()
        self.assertTrue(gui.export_bar.winfo_viewable())
        self.assertGreater(gui.export_bar.winfo_width(), 150)
        self.assertTrue(gui.preview_image.winfo_viewable())
        self.assertGreater(gui.preview_image.winfo_height(), 90)
        point = next((p for p in gui._scaffold_hit_points.values() if p[2] == 2), None)
        self.assertIsNotNone(point)
        gui.select_scaffold_atom(SimpleNamespace(x=round(point[0]), y=round(point[1])))
        self.assertEqual(gui.site_choice.get(), "2")
        self.assertIn(2, gui.r_sites)
        self.assertGreaterEqual(gui.assignment_list.size(), 1)
        gui.fragment_list.selection_clear(0, "end")
        gui.fragment_list.selection_set(0)
        gui.assign()
        gui.enumerate()
        gui.update()
        self.assertGreater(len(gui.records), 0)
        self.assertIsNotNone(gui.preview_image.image)
        self.assertLessEqual(gui.preview_image.image.width(), gui.preview_image.winfo_width())
        self.assertLessEqual(gui.preview_image.image.height(), gui.preview_image.winfo_height())
        self.assertTrue(gui.export_bar.winfo_viewable())
        gui.geometry("1380x820")
        gui.update()

    def test_updater_controls_and_manual_check(self):
        from unittest.mock import patch
        gui = self.gui
        self.assertTrue(gui.check_update_button.winfo_viewable())
        self.assertIn("Check for Updates", gui.check_update_button.cget("text"))
        self.assertTrue(gui.auto_check_updates.get())
        with patch("app.latest_update", return_value=None):
            gui.check_for_updates(silent=True)
            gui.update()
            import time
            for _ in range(50):
                gui.update()
                if not gui.update_busy:
                    break
                time.sleep(0.03)
            self.assertFalse(gui.update_busy)
            self.assertIn("up to date",gui.status.get())
        with patch("app.set_auto_check") as setting:
            gui.auto_check_updates.set(False)
            gui._set_auto_update_preference()
            setting.assert_called_once_with(False)
            gui.auto_check_updates.set(True)

    def test_high_contrast_theme_and_primary_buttons(self):
        from tkinter import ttk
        gui = self.gui
        theme = ttk.Style(gui)
        self.assertEqual(theme.theme_use(), "clam")
        self.assertEqual(theme.lookup("Primary.TButton", "background"), "#0759A7")
        self.assertEqual(theme.lookup("Primary.TButton", "foreground"), "#FFFFFF")
        self.assertEqual(theme.lookup("Treeview", "foreground"), "#10243B")
        self.assertEqual(gui.enumerate_all_button.cget("style"), "Primary.TButton")
        self.assertEqual(gui.enumerate_current_button.cget("style"), "Primary.TButton")
        self.assertEqual(gui.fragment_list.cget("fg"), "#10243B")
        self.assertEqual(gui.assignment_list.cget("selectbackground"), "#0864B1")

    def test_editor_scrollbar_and_mousewheel(self):
        gui = self.gui
        gui.geometry("1040x670")
        gui.update()
        canvas = gui.editor_canvas
        self.assertTrue(canvas.winfo_viewable())
        self.assertTrue(gui.editor_scrollbar.winfo_viewable())
        self.assertIsNotNone(canvas.cget("yscrollcommand"))
        self.assertTrue(gui.assignment_list.winfo_exists())
        # Canvas must have valid content and a functioning scroll region.
        bbox = canvas.bbox("all")
        self.assertIsNotNone(bbox)
        self.assertGreater(bbox[3], 0)
        canvas.yview_moveto(1.0)
        gui.update()
        bottom = canvas.yview()
        self.assertGreaterEqual(bottom[1], 0.99)
        canvas.yview_moveto(0.0)
        gui.update()
        gui.geometry("1280x760")
        gui.update()

    def test_primary_controls_visible_and_initial_preview(self):
        gui = self.gui
        gui.geometry("1040x670")
        gui.update()
        self.assertIsNotNone(gui.preview_image.image)
        self.assertIn("PARENT SCAFFOLD", gui.properties.cget("text"))
        for button in (gui.enumerate_all_button, gui.enumerate_current_button):
            self.assertEqual(button.winfo_viewable(), 1)
            self.assertGreater(button.winfo_width(), 70)
            self.assertGreaterEqual(button.winfo_rootx(), gui.winfo_rootx())
            self.assertLessEqual(button.winfo_rootx()+button.winfo_width(),
                                 gui.winfo_rootx()+gui.winfo_width())
            self.assertLessEqual(button.winfo_rooty()+button.winfo_height(),
                                 gui.winfo_rooty()+gui.winfo_height())
        gui.geometry("1280x760")
        gui.update()

    def test_scaffold_change_resets_results_and_assignments(self):
        gui=self.gui
        gui.site_choice.set("0")
        gui.fragment_list.selection_clear(0,"end")
        gui.fragment_list.selection_set(0)
        gui.assign()
        gui.enumerate()
        self.assertGreater(len(gui.records),0)
        gui.show_hydrogens.set(True)
        gui.draw_scaffold()
        self.assertIsNotNone(gui.scaffold_display.image)
        gui.scaffold.set("CCO")
        gui.refresh_scaffold()
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
            # Verify selected-only export writes exactly one product.
            gui.export_selection_only.set(True)
            selected_path=Path(tmp)/"selected.sdf"
            with patch("app.filedialog.asksaveasfilename",return_value=str(selected_path)):
                gui.export("3d")
            selected_mols=[m for m in Chem.SDMolSupplier(str(selected_path)) if m]
            self.assertEqual(len(selected_mols),1)
            gui.export_selection_only.set(False)
            # Import another fragment library using the actual handler.
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
