"""RGroupStudio: standalone Tkinter desktop application."""
import os
import sys
import logging
import traceback
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, simpledialog
from PIL import ImageTk
from rdkit import Chem
from rdkit.Chem import Draw
from chemistry import FRAGMENTS, parse_scaffold, sites, validate_fragment, load_fragments, enumerate_library, molecular_properties, export_sdf, export_smi, export_csv, conformer_3d

LOGDIR = Path(os.getenv("LOCALAPPDATA") or Path.home()) / "RGroupStudio" / "logs"
LOGDIR.mkdir(parents=True, exist_ok=True)
LOGFILE = LOGDIR / "RGroupStudio.log"
logging.basicConfig(filename=LOGFILE, level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s", encoding="utf-8")

class Studio(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("RGroupStudio | Molecular Library Designer")
        try:
            asset_root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
            self.logo_image = tk.PhotoImage(file=str(asset_root / "assets" / "approved_logo.png"))
            self.iconphoto(True, self.logo_image)
        except Exception:
            logging.exception("Application icon could not be loaded")
        self.geometry("1350x810")
        self.minsize(1060, 650)
        self.configure(bg="#f1f5fa")
        self.scaffold = tk.StringVar(value="c1ccccc1")
        self.site_choice = tk.StringVar()
        self.fragments = dict(FRAGMENTS)
        self.assigned = {}
        self.records = []
        self.selected = None
        self.images = []
        self.last_export_errors = []
        self.limit = tk.IntVar(value=2500)
        self.show_indices = tk.BooleanVar(value=True)
        self.show_hydrogens = tk.BooleanVar(value=False)
        self.export_selection_only = tk.BooleanVar(value=False)
        self.current_fragments = {}
        self.status = tk.StringVar(value="Ready")
        self._build()
        self.refresh_scaffold()

    def _build(self):
        style = ttk.Style(self)
        try:
            style.theme_use("vista")
        except tk.TclError:
            style.theme_use("clam")
        style.configure("TButton", padding=(8, 4), font=("Segoe UI", 9))
        style.configure("TLabel", font=("Segoe UI", 9))
        style.configure("TLabelframe.Label", font=("Segoe UI", 9, "bold"))
        style.configure("Treeview", rowheight=25, font=("Segoe UI", 9))
        style.configure("Treeview.Heading", font=("Segoe UI", 9, "bold"))
        self.geometry("1280x760")
        self.minsize(1040, 670)
        # Compact toolbar: no oversized branding or unused banner.
        header = tk.Frame(self, bg="#17365d", height=44)
        header.pack(fill="x")
        header.pack_propagate(False)
        tk.Label(header, text="RGroupStudio", font=("Segoe UI", 14, "bold"),
                 fg="white", bg="#17365d", padx=13).pack(side="left", fill="y")
        tk.Label(header, text="Molecular Library Designer", font=("Segoe UI", 9),
                 fg="#d9e7f8", bg="#17365d").pack(side="left")
        main = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        main.pack(fill="both", expand=True, padx=9, pady=7)
        left = ttk.Frame(main, width=362)
        right = ttk.Frame(main)
        main.add(left, weight=1)
        main.add(right, weight=3)

        # Left: fit the full editing workflow within the visible window.
        a = ttk.LabelFrame(left, text="1  Scaffold", padding=7)
        a.pack(fill="x", pady=(0, 5))
        inputbar = ttk.Frame(a)
        inputbar.pack(fill="x")
        ttk.Entry(inputbar, textvariable=self.scaffold).pack(side="left", fill="x", expand=True)
        ttk.Button(inputbar, text="Draw", command=self.refresh_scaffold).pack(side="left", padx=(5, 0))
        optionrow = ttk.Frame(a)
        optionrow.pack(fill="x", pady=(3,0))
        ttk.Checkbutton(optionrow, text="Atom indices", variable=self.show_indices,
                        command=self.draw_scaffold).pack(side="left")
        ttk.Checkbutton(optionrow, text="Show H", variable=self.show_hydrogens,
                        command=self.draw_scaffold).pack(side="left", padx=12)
        self.scaffold_display = ttk.Label(a, anchor="center")
        self.scaffold_display.pack(fill="x", pady=(0, 1))

        b = ttk.LabelFrame(left, text="2  Attachment site", padding=7)
        b.pack(fill="x", pady=(0, 5))
        self.site_combo = ttk.Combobox(b, state="readonly", textvariable=self.site_choice)
        self.site_combo.pack(fill="x")

        c = ttk.LabelFrame(left, text="3  Fragments and assignments", padding=7)
        c.pack(fill="both", expand=True, pady=(0, 5))
        fragmentbar = ttk.Frame(c)
        fragmentbar.pack(fill="x", pady=(0, 4))
        ttk.Button(fragmentbar, text="Import .smi", command=self.import_fragments).pack(side="left")
        ttk.Button(fragmentbar, text="+ Fragment", command=self.add_fragment).pack(side="left", padx=4)
        ttk.Button(fragmentbar, text="Add R-group", command=self.assign).pack(side="right")
        listframe = ttk.Frame(c)
        listframe.pack(fill="both", expand=True)
        self.fragment_list = tk.Listbox(listframe, height=6, selectmode="extended",
                                        exportselection=False, activestyle="dotbox",
                                        font=("Segoe UI", 9), relief="solid", bd=1)
        self.fragment_list.pack(side="left", fill="both", expand=True)
        fragscroll = ttk.Scrollbar(listframe, orient="vertical", command=self.fragment_list.yview)
        fragscroll.pack(side="right", fill="y")
        self.fragment_list.configure(yscrollcommand=fragscroll.set)
        self.refresh_fragment_list()
        ttk.Label(c, text="R-group assignments").pack(anchor="w", pady=(5, 1))
        assignmentrow = ttk.Frame(c)
        assignmentrow.pack(fill="x")
        self.assignment_list = tk.Listbox(assignmentrow, height=3, font=("Segoe UI", 9),
                                         relief="solid", bd=1)
        self.assignment_list.pack(side="left", fill="x", expand=True)
        ttk.Button(assignmentrow, text="Remove", command=self.remove_assignment).pack(side="left", padx=(5, 0))

        controls = ttk.Frame(left)
        controls.pack(fill="x", pady=(0, 3))
        ttk.Label(controls, text="Max products:").pack(side="left")
        ttk.Entry(controls, textvariable=self.limit, width=7).pack(side="left", padx=(5, 8))
        ttk.Button(controls, text="Current", command=self.enumerate_current).pack(side="left", padx=(0,3))
        ttk.Button(controls, text="All combinations", command=self.enumerate).pack(side="left", fill="x", expand=True)

        # Right: a compact results table with room for an actual molecule preview.
        resulthead = ttk.Frame(right)
        resulthead.pack(fill="x", pady=(0, 5))
        ttk.Label(resulthead, text="Enumerated molecules", font=("Segoe UI", 10, "bold")).pack(side="left")
        self.results_count = ttk.Label(resulthead, text="0 compounds")
        self.results_count.pack(side="right")
        results = ttk.Frame(right)
        results.pack(fill="x")
        cols = ("id", "smiles", "mw", "logp", "tpsa", "qed")
        self.table = ttk.Treeview(results, columns=cols, show="headings", height=9)
        labels = {"id":"ID", "smiles":"Canonical SMILES", "mw":"MW", "logp":"cLogP", "tpsa":"TPSA", "qed":"QED"}
        for col in cols:
            self.table.heading(col, text=labels[col])
            self.table.column(col, width=90 if col!="smiles" else 300, minwidth=65,
                              stretch=(col=="smiles"))
        self.table.pack(side="left", fill="x", expand=True)
        scrollbar=ttk.Scrollbar(results, orient="vertical", command=self.table.yview)
        scrollbar.pack(side="right", fill="y")
        self.table.configure(yscrollcommand=scrollbar.set)
        self.table.bind("<<TreeviewSelect>>", self.choose_record)

        preview = ttk.LabelFrame(right, text="Molecular preview and properties", padding=7)
        preview.pack(fill="both", expand=True, pady=(7, 5))
        self.preview_image=ttk.Label(preview, anchor="center")
        self.preview_image.pack(side="left", fill="both", expand=True)
        self.properties=ttk.Label(preview, text="Load a scaffold and enumerate to view products.\n\nSelect any product to inspect properties.",
                                  font=("Consolas", 10), justify="left")
        self.properties.pack(side="left", anchor="n", padx=(8, 22), pady=12)
        actions = ttk.Frame(right)
        actions.pack(fill="x", pady=(0, 2))
        ttk.Checkbutton(actions, text="Selected only", variable=self.export_selection_only).pack(side="left", padx=(0, 6))
        ttk.Button(actions, text="CSV", command=lambda:self.export("csv")).pack(side="left", padx=(0, 4))
        ttk.Button(actions, text="SMILES", command=lambda:self.export("smi")).pack(side="left", padx=4)
        ttk.Button(actions, text="2D SDF", command=lambda:self.export("2d")).pack(side="left", padx=4)
        ttk.Button(actions, text="3D SDF", command=lambda:self.export("3d")).pack(side="left", padx=4)
        ttk.Button(actions, text="View 3D", command=self.view3d).pack(side="left", padx=4)
        tk.Label(self, textvariable=self.status, bg="#e2eaf4", fg="#14325a",
                 anchor="w", padx=9, pady=3, font=("Segoe UI", 9)).pack(fill="x")

    def error(self, title, exc):
        logging.exception("%s: %s", title, exc)
        messagebox.showerror(title, str(exc) + "\n\nLog: " + str(LOGFILE))

    def reset_results(self):
        self.records = []
        self.selected = None
        self.table.delete(*self.table.get_children())
        self.results_count.configure(text="0 compounds")
        self.preview_image.configure(image="")
        self.preview_image.image = None
        self.properties.configure(text="Enumerate a library to inspect molecules.")

    def refresh_scaffold(self):
        try:
            candidate = parse_scaffold(self.scaffold.get())
            self.mol = candidate
            self.loaded_scaffold = self.scaffold.get()
            allowed = sites(self.mol)
            self.site_combo["values"]=[str(i) for i in allowed]
            self.site_choice.set(str(allowed[0]) if allowed else "")
            self.assigned.clear()
            self.current_fragments.clear()
            self.refresh_assignments()
            self.reset_results()
            self.draw_scaffold()
            self.status.set("Scaffold loaded: %d substitutable positions" % len(allowed))
        except Exception as exc:
            self.error("Scaffold error", exc)

    def draw_scaffold(self):
        if not hasattr(self, "mol"):
            return
        m = Chem.AddHs(self.mol) if self.show_hydrogens.get() else Chem.Mol(self.mol)
        if self.show_indices.get():
            for atom in m.GetAtoms():
                atom.SetProp("atomNote", str(atom.GetIdx()))
        img = Draw.MolToImage(m, size=(310, 145))
        photo = ImageTk.PhotoImage(img)
        self.scaffold_display.configure(image=photo)
        self.scaffold_display.image = photo

    def refresh_fragment_list(self):
        self.fragment_names = list(self.fragments)
        self.fragment_list.delete(0, tk.END)
        for name in self.fragment_names:
            self.fragment_list.insert(tk.END, "%s  |  %s" % (name, self.fragments[name]))

    def import_fragments(self):
        path=filedialog.askopenfilename(filetypes=[("SMILES library", "*.smi *.smiles *.txt"), ("All files","*.*")])
        if not path:
            return
        try:
            imported=load_fragments(path)
            self.fragments.update(imported)
            self.refresh_fragment_list()
            self.reset_results()
            self.status.set("Imported %d fragments" % len(imported))
        except Exception as exc:
            self.error("Fragment import failed",exc)

    def add_fragment(self):
        name=simpledialog.askstring("Fragment name","Name of new fragment:",parent=self)
        if not name:
            return
        smiles=simpledialog.askstring("Fragment SMILES","SMILES with one terminal [*] attachment point:",parent=self)
        if not smiles:
            return
        try:
            validate_fragment(smiles)
            self.fragments[name]=smiles
            self.refresh_fragment_list()
            self.reset_results()
        except Exception as exc:
            self.error("Invalid fragment",exc)

    def refresh_assignments(self):
        self.assignment_list.delete(0,tk.END)
        for n,(site, items) in enumerate(sorted(self.assigned.items()),1):
            self.assignment_list.insert(tk.END,"R%d @ atom %d: %s"%(n,site,", ".join(items)))

    def assign(self):
        try:
            site=int(self.site_choice.get())
            selected=[self.fragment_names[i] for i in self.fragment_list.curselection()]
            if not selected:
                raise ValueError("Select at least one fragment")
            self.assigned[site]=selected
            self.current_fragments[site]=selected[0]
            self.refresh_assignments()
            self.reset_results()
            self.status.set("Assigned %d fragments to R @ %d"%(len(selected),site))
        except Exception as exc:
            self.error("Assignment failed",exc)

    def remove_assignment(self):
        selected=self.assignment_list.curselection()
        if not selected:
            return
        site=sorted(self.assigned)[selected[0]]
        self.assigned.pop(site)
        self.current_fragments.pop(site, None)
        self.refresh_assignments()
        self.reset_results()

    def enumerate_current(self):
        try:
            self._enumerate({site: [self.fragments[self.current_fragments[site]]]
                             for site in self.assigned}, "current assignment")
        except Exception as exc:
            self.error("Current assignment failed", exc)

    def enumerate(self):
        try:
            self._enumerate({site:[self.fragments[name] for name in names]
                             for site,names in self.assigned.items()}, "all combinations")
        except Exception as exc:
            self.error("Enumeration failed", exc)

    def _enumerate(self, assignments, label):
        try:
            self.status.set("Enumerating molecules...")
            self.update_idletasks()
            self.records=enumerate_library(self.loaded_scaffold,assignments,self.limit.get())
            self.table.delete(*self.table.get_children())
            self.results_count.configure(text=f"{len(self.records)} compounds")
            for n,(name,mol) in enumerate(self.records):
                p=molecular_properties(mol)
                self.table.insert("",tk.END,iid=str(n),
                                  values=(name,Chem.MolToSmiles(mol),p["MW"],p["cLogP"],p["TPSA"],p["QED"]))
            if self.records:
                self.table.selection_set("0")
                self.choose_record()
            else:
                self.selected = None
            self.status.set("Generated %d unique molecules (%s)" % (len(self.records), label))
        except Exception as exc:
            self.error("Enumeration failed",exc)

    def choose_record(self,event=None):
        selected=self.table.selection()
        if not selected:
            return
        self.selected=int(selected[0])
        name,mol=self.records[self.selected]
        pic=ImageTk.PhotoImage(Draw.MolToImage(mol,size=(360,245)))
        self.preview_image.configure(image=pic)
        self.preview_image.image=pic
        p=molecular_properties(mol)
        self.properties.configure(text="\n".join([name, ""] + [f"{k:>7}: {v}" for k,v in p.items()]))

    def export(self,kind):
        if not self.records:
            messagebox.showinfo("No molecules","Enumerate a library first")
            return
        records = self.records
        if self.export_selection_only.get():
            ids = self.table.selection()
            if not ids:
                messagebox.showinfo("No selection","Select at least one molecule first")
                return
            records = [self.records[int(i)] for i in ids]
        ext="sdf" if kind in ("2d","3d") else kind
        path=filedialog.asksaveasfilename(defaultextension="."+ext,
                   initialfile="RGroupStudio_%s.%s"%(kind,ext))
        if not path:
            return
        try:
            errors = []
            if kind=="csv":
                export_csv(path,records)
            elif kind=="smi":
                export_smi(path,records)
            else:
                errors=export_sdf(path,records,three_d=(kind=="3d"))
            self.last_export_errors = errors
            if errors:
                messagebox.showwarning("Partial export", "%d molecules failed. First errors: %s" % (len(errors), errors[:3]))
            completed = len(records) - len(errors)
            self.status.set("Exported %d / %d molecules to %s" % (completed,len(records),path))
        except Exception as exc:
            self.error("Export failed",exc)

    def view3d(self):
        if self.selected is None:
            messagebox.showinfo("Select molecule","Select a generated molecule first")
            return
        try:
            import matplotlib
            matplotlib.use("TkAgg")
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            from mpl_toolkits.mplot3d import Axes3D
            m=conformer_3d(self.records[self.selected][1])
            positions=m.GetConformer().GetPositions()
            window=tk.Toplevel(self)
            window.title("3D conformer — " + self.records[self.selected][0])
            window.geometry("650x610")
            fig=Figure(figsize=(6,5),dpi=100)
            ax=fig.add_subplot(111,projection="3d")
            colors={"C":"#3e4f61","N":"#2f6de0","O":"#e34349","F":"#3ca86b",
                    "Cl":"#3ca86b","Br":"#a44d32","S":"#e0ba32","H":"#c6cbd1"}
            for b in m.GetBonds():
                a1,a2=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
                ax.plot(positions[[a1,a2],0],positions[[a1,a2],1],positions[[a1,a2],2],color="#8c98a8",lw=2)
            for atom in m.GetAtoms():
                if atom.GetSymbol()=="H":
                    continue
                p=positions[atom.GetIdx()]
                ax.scatter([p[0]],[p[1]],[p[2]],s=125,c=colors.get(atom.GetSymbol(),"#8268ad"),depthshade=True)
            ax.set_axis_off()
            fig.tight_layout()
            canvas=FigureCanvasTkAgg(fig,master=window)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both",expand=True)
        except Exception as exc:
            self.error("3D visualization failed",exc)

def self_test():
    from chemistry import enumerate_library
    import tempfile
    from chemistry import export_sdf, export_csv, export_smi
    records=enumerate_library("c1ccccc1",{0:["[*]F","[*]Cl"],2:["[*]C"]})
    assert len(records)==2
    assert all(m is not None for _,m in records)
    with tempfile.TemporaryDirectory() as folder:
        for extension, export in (("csv",export_csv),("smi",export_smi)):
            dest = str(Path(folder) / ("library." + extension))
            export(dest,records)
            assert Path(dest).stat().st_size > 0
        for dimension in (False,True):
            dest = str(Path(folder) / ("library_%s.sdf" % ("3d" if dimension else "2d")))
            assert not export_sdf(dest,records,three_d=dimension)
            got = [m for m in Chem.SDMolSupplier(dest) if m]
            assert len(got) == len(records)
            assert all(m.GetNumConformers() == 1 for m in got)
    return 0

def main():
    try:
        if "--self-test" in sys.argv:
            return self_test()
        app=Studio()
        app.mainloop()
        return 0
    except Exception:
        logging.exception("Uncaught fatal application error")
        try:
            root=tk.Tk()
            root.withdraw()
            messagebox.showerror("RGroupStudio startup error",
                                 "The application failed to start. Log: "+str(LOGFILE)+"\n\n"+traceback.format_exc()[-1800:])
            root.destroy()
        except Exception:
            pass
        return 1

if __name__=="__main__":
    sys.exit(main())
