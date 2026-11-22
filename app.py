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
        self.limit = tk.IntVar(value=2500)
        self.show_indices = tk.BooleanVar(value=True)
        self.status = tk.StringVar(value="Ready")
        self._build()
        self.refresh_scaffold()

    def _build(self):
        title = tk.Frame(self, bg="#17365d")
        title.pack(fill="x")
        tk.Label(title, text="RGroupStudio", font=("Segoe UI", 23, "bold"),
                 fg="white", bg="#17365d", padx=20, pady=12).pack(side="left")
        tk.Label(title, text="Scaffold-driven 2D / 3D ligand design",
                 font=("Segoe UI", 12), fg="#d7e8ff", bg="#17365d").pack(side="left")
        main = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        main.pack(fill="both", expand=True, padx=12, pady=10)
        left = ttk.Frame(main, width=410)
        right = ttk.Frame(main)
        main.add(left, weight=1)
        main.add(right, weight=3)

        a = ttk.LabelFrame(left, text="1. Scaffold SMILES", padding=9)
        a.pack(fill="x", pady=5)
        ttk.Entry(a, textvariable=self.scaffold, width=42).pack(fill="x")
        ttk.Button(a, text="Load / Draw scaffold", command=self.refresh_scaffold).pack(anchor="w", pady=4)
        ttk.Checkbutton(a, text="Show atom indices", variable=self.show_indices,
                        command=self.draw_scaffold).pack(anchor="w")
        self.scaffold_display = ttk.Label(a)
        self.scaffold_display.pack(pady=4)

        b = ttk.LabelFrame(left, text="2. R-group positions", padding=9)
        b.pack(fill="x", pady=5)
        self.site_combo = ttk.Combobox(b, state="readonly", textvariable=self.site_choice)
        self.site_combo.pack(fill="x")
        ttk.Label(b, text="Heavy-atom indices with replaceable H atoms").pack(anchor="w", pady=4)

        c = ttk.LabelFrame(left, text="3. Fragment library", padding=9)
        c.pack(fill="both", expand=True, pady=5)
        bar = ttk.Frame(c)
        bar.pack(fill="x")
        ttk.Button(bar, text="Import .smi", command=self.import_fragments).pack(side="left")
        ttk.Button(bar, text="Add SMILES", command=self.add_fragment).pack(side="left", padx=5)
        self.fragment_list = tk.Listbox(c, height=10, selectmode="extended", exportselection=False)
        self.fragment_list.pack(fill="both", expand=True, pady=6)
        self.refresh_fragment_list()
        ttk.Button(c, text="Assign selected fragments to R-site",
                   command=self.assign).pack(fill="x", pady=3)
        ttk.Label(c, text="Assignments").pack(anchor="w")
        self.assignment_list = tk.Listbox(c, height=5)
        self.assignment_list.pack(fill="x", pady=3)
        ttk.Button(c, text="Remove selected assignment",
                   command=self.remove_assignment).pack(anchor="w")

        controls = ttk.Frame(left)
        controls.pack(fill="x", pady=8)
        ttk.Label(controls, text="Max products").pack(side="left")
        ttk.Entry(controls, textvariable=self.limit, width=8).pack(side="left", padx=5)
        ttk.Button(controls, text="Enumerate all combinations",
                   command=self.enumerate).pack(side="left", padx=5)

        results = ttk.LabelFrame(right, text="Enumerated molecules", padding=8)
        results.pack(fill="both", expand=True)
        cols = ("id", "smiles", "mw", "logp", "tpsa", "qed")
        self.table = ttk.Treeview(results, columns=cols, show="headings", height=14)
        labels = {"id":"ID", "smiles":"Canonical SMILES", "mw":"MW", "logp":"cLogP", "tpsa":"TPSA", "qed":"QED"}
        for col in cols:
            self.table.heading(col, text=labels[col])
            self.table.column(col, width=105 if col!="smiles" else 350, stretch=(col=="smiles"))
        self.table.pack(side="left", fill="both", expand=True)
        scrollbar=ttk.Scrollbar(results, orient="vertical", command=self.table.yview)
        scrollbar.pack(side="right", fill="y")
        self.table.configure(yscrollcommand=scrollbar.set)
        self.table.bind("<<TreeviewSelect>>", self.choose_record)

        preview = ttk.LabelFrame(right, text="Selected molecule — 2D preview and properties", padding=8)
        preview.pack(fill="both", expand=True, pady=8)
        self.preview_image=ttk.Label(preview)
        self.preview_image.pack(side="left", padx=15)
        self.properties=ttk.Label(preview, text="Select a molecule in the results table",
                                  font=("Consolas", 11), justify="left")
        self.properties.pack(side="left", anchor="n", padx=15, pady=12)

        actions=ttk.Frame(right)
        actions.pack(fill="x", pady=5)
        ttk.Button(actions, text="Export CSV", command=lambda:self.export("csv")).pack(side="left", padx=3)
        ttk.Button(actions, text="Export SMILES", command=lambda:self.export("smi")).pack(side="left", padx=3)
        ttk.Button(actions, text="Export 2D SDF", command=lambda:self.export("2d")).pack(side="left", padx=3)
        ttk.Button(actions, text="Export 3D SDF", command=lambda:self.export("3d")).pack(side="left", padx=3)
        ttk.Button(actions, text="View 3D", command=self.view3d).pack(side="left", padx=3)
        tk.Label(self, textvariable=self.status, bg="#e2eaf4", fg="#14325a",
                 anchor="w", padx=12, pady=5).pack(fill="x")

    def error(self, title, exc):
        logging.exception("%s: %s", title, exc)
        messagebox.showerror(title, str(exc) + "\n\nLog: " + str(LOGFILE))

    def refresh_scaffold(self):
        try:
            self.mol = parse_scaffold(self.scaffold.get())
            allowed = sites(self.mol)
            self.site_combo["values"]=[str(i) for i in allowed]
            self.site_choice.set(str(allowed[0]) if allowed else "")
            self.assigned.clear()
            self.refresh_assignments()
            self.draw_scaffold()
            self.status.set("Scaffold loaded: %d substitutable positions" % len(allowed))
        except Exception as exc:
            self.error("Scaffold error", exc)

    def draw_scaffold(self):
        if not hasattr(self, "mol"):
            return
        m = Chem.Mol(self.mol)
        if self.show_indices.get():
            for atom in m.GetAtoms():
                atom.SetProp("atomNote", str(atom.GetIdx()))
        img = Draw.MolToImage(m, size=(355, 190))
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
        except Exception as exc:
            self.error("Invalid fragment",exc)

    def refresh_assignments(self):
        self.assignment_list.delete(0,tk.END)
        for site, items in sorted(self.assigned.items()):
            self.assignment_list.insert(tk.END,"R @ %d: %s"%(site,", ".join(items)))

    def assign(self):
        try:
            site=int(self.site_choice.get())
            selected=[self.fragment_names[i] for i in self.fragment_list.curselection()]
            if not selected:
                raise ValueError("Select at least one fragment")
            self.assigned[site]=selected
            self.refresh_assignments()
            self.status.set("Assigned %d fragments to R @ %d"%(len(selected),site))
        except Exception as exc:
            self.error("Assignment failed",exc)

    def remove_assignment(self):
        selected=self.assignment_list.curselection()
        if not selected:
            return
        site=sorted(self.assigned)[selected[0]]
        self.assigned.pop(site)
        self.refresh_assignments()

    def enumerate(self):
        try:
            self.status.set("Enumerating molecules...")
            self.update_idletasks()
            assignments={site:[self.fragments[name] for name in names]
                         for site,names in self.assigned.items()}
            self.records=enumerate_library(self.scaffold.get(),assignments,self.limit.get())
            self.table.delete(*self.table.get_children())
            for n,(name,mol) in enumerate(self.records):
                p=molecular_properties(mol)
                self.table.insert("",tk.END,iid=str(n),
                                  values=(name,Chem.MolToSmiles(mol),p["MW"],p["cLogP"],p["TPSA"],p["QED"]))
            if self.records:
                self.table.selection_set("0")
                self.choose_record()
            self.status.set("Generated %d unique molecules" % len(self.records))
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
        ext="sdf" if kind in ("2d","3d") else kind
        path=filedialog.asksaveasfilename(defaultextension="."+ext,
                   initialfile="RGroupStudio_%s.%s"%(kind,ext))
        if not path:
            return
        try:
            if kind=="csv":
                export_csv(path,self.records)
            elif kind=="smi":
                export_smi(path,self.records)
            else:
                errors=export_sdf(path,self.records,three_d=(kind=="3d"))
                if errors:
                    messagebox.showwarning("Partial 3D export", "%d structures failed 3D generation. Details: %s"%(len(errors),errors[:3]))
            self.status.set("Exported %d molecules to %s"%(len(self.records),path))
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
    records=enumerate_library("c1ccccc1",{0:["[*]F","[*]Cl"],2:["[*]C"]})
    assert len(records)==2
    assert all(m is not None for _,m in records)
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
