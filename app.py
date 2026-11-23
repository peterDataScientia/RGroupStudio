"""RGroupStudio: standalone Tkinter desktop application."""
import os
import sys
import logging
import traceback
import io
import math
import queue
import threading
import subprocess
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, simpledialog
from PIL import ImageTk, Image
from rdkit.Chem.Draw import rdMolDraw2D
from updater import APP_VERSION, latest_update, download_installer, get_auto_check, set_auto_check
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
        self.title(f"RGroupStudio v{APP_VERSION} | Molecular Library Designer")
        try:
            asset_root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
            self.logo_image = tk.PhotoImage(file=str(asset_root / "assets" / "approved_logo.png"))
            self.iconphoto(True, self.logo_image)
        except Exception:
            logging.exception("Application icon could not be loaded")
        self.geometry("1350x810")
        self.minsize(1060, 650)
        self.configure(bg="#E9EEF5")
        self.scaffold = tk.StringVar(value="c1ccccc1")
        self.site_choice = tk.StringVar()
        self.fragments = dict(FRAGMENTS)
        self.assigned = {}
        self.r_sites = set()
        self._scaffold_hit_points = {}
        self._preview_mol = None
        self._preview_render_size = None
        self._preview_refresh_id = None
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
        self.auto_check_updates = tk.BooleanVar(value=get_auto_check())
        self.update_events = queue.Queue()
        self.update_busy = False
        self._update_info = None
        self._build()
        self.refresh_scaffold()
        if self.auto_check_updates.get() and not os.getenv("RGROUPSTUDIO_DISABLE_AUTO_CHECK"):
            self.after(2200, lambda: self.check_for_updates(silent=True))

    def _build(self):
        # CLAM allows consistent, high-contrast palette on Windows.
        # Native Vista ignores many ttk background and text color properties.
        style = ttk.Style(self)
        style.theme_use("clam")
        self.palette = {
            "background": "#E9EEF5", "surface": "#FFFFFF", "text": "#10243B",
            "muted": "#334E69", "border": "#A8BBCF", "primary": "#0759A7",
            "primary_hover": "#034580", "accent": "#0F766E",
            "selection": "#0864B1", "selection_text": "#FFFFFF",
        }
        P = self.palette
        style.configure(".", background=P["background"], foreground=P["text"],
                        font=("Segoe UI", 10))
        style.configure("TFrame", background=P["background"])
        style.configure("TLabel", background=P["background"], foreground=P["text"],
                        font=("Segoe UI", 10, "bold"))
        style.configure("TButton", padding=(10, 6), font=("Segoe UI", 10, "bold"),
                        background="#D7E4F1", foreground=P["text"], bordercolor=P["border"],
                        lightcolor=P["border"], darkcolor=P["border"], relief="raised")
        style.map("TButton", background=[("active", "#C1D6EB"), ("pressed", "#AFCAE3")],
                  foreground=[("disabled", "#53697F")])
        style.configure("Primary.TButton", padding=(13, 7), font=("Segoe UI", 10, "bold"),
                        foreground="#FFFFFF", background=P["primary"],
                        bordercolor=P["primary"], relief="raised")
        style.map("Primary.TButton",
                  background=[("active", P["primary_hover"]), ("pressed", "#02376A")],
                  foreground=[("active", "#FFFFFF"), ("disabled", "#DAE1E8")])
        style.configure("TLabelframe", background=P["surface"], borderwidth=2,
                        bordercolor=P["border"], relief="solid")
        style.configure("TLabelframe.Label", background=P["background"],
                        foreground="#0A426F", font=("Segoe UI", 10, "bold"))
        style.configure("TEntry", fieldbackground=P["surface"], foreground=P["text"],
                        insertcolor=P["text"], bordercolor=P["border"], padding=5)
        style.configure("TCombobox", fieldbackground=P["surface"], foreground=P["text"],
                        bordercolor=P["border"], padding=4)
        style.map("TCombobox", fieldbackground=[("readonly", P["surface"])],
                  foreground=[("readonly", P["text"])])
        style.configure("TCheckbutton", background=P["background"], foreground=P["text"],
                        font=("Segoe UI", 10, "bold"))
        style.map("TCheckbutton", background=[("active", P["background"])])
        style.configure("Treeview", rowheight=29, background=P["surface"],
                        fieldbackground=P["surface"], foreground=P["text"],
                        font=("Segoe UI", 10), borderwidth=1)
        style.map("Treeview", background=[("selected", P["selection"])],
                  foreground=[("selected", P["selection_text"])])
        style.configure("Treeview.Heading", background="#C7D9ED", foreground="#0B3259",
                        font=("Segoe UI", 10, "bold"), padding=7, relief="raised")
        style.map("Treeview.Heading", background=[("active", "#AFC9E4")])
        style.configure("Vertical.TScrollbar", background="#A3BED8", troughcolor="#E5ECF4",
                        arrowcolor="#143C5F", bordercolor="#7A99B6")
        self.geometry("1380x820")
        self.minsize(1100, 690)
        # Compact toolbar: no oversized branding or unused banner.
        header = tk.Frame(self, bg="#17365d", height=44)
        header.pack(fill="x")
        header.pack_propagate(False)
        tk.Label(header, text="RGroupStudio", font=("Segoe UI", 14, "bold"),
                 fg="white", bg="#17365d", padx=13).pack(side="left", fill="y")
        tk.Label(header, text="Molecular Library Designer", font=("Segoe UI", 9),
                 fg="#d9e7f8", bg="#17365d").pack(side="left")
        self.check_update_button = tk.Button(
            header, text="Check for Updates", command=self.check_for_updates,
            font=("Segoe UI", 9, "bold"), fg="#FFFFFF", bg="#0759A7",
            activebackground="#1274C5", activeforeground="#FFFFFF",
            relief="flat", cursor="hand2", padx=11, pady=4)
        self.check_update_button.pack(side="right", padx=(8, 10), pady=6)
        tk.Checkbutton(
            header, text="Auto-check at startup", variable=self.auto_check_updates,
            command=self._set_auto_update_preference, bg="#17365d",
            fg="#FFFFFF", activebackground="#17365d", activeforeground="#FFFFFF",
            selectcolor="#17365d", font=("Segoe UI", 9, "bold")
        ).pack(side="right", padx=(0, 6))
        tk.Label(header, text="v" + APP_VERSION, fg="#F0F6FF",
                 bg="#17365d", font=("Segoe UI", 9, "bold")).pack(side="right", padx=(0, 8))
        main = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        main.pack(fill="both", expand=True, padx=9, pady=7)
        # The complete editor (including assignments) is scrollable.
        left_host = ttk.Frame(main, width=405)
        right = ttk.Frame(main)
        main.add(left_host, weight=2)
        main.add(right, weight=4)
        self.editor_canvas = tk.Canvas(left_host, borderwidth=0, highlightthickness=0,
                                       background="#E9EEF5")
        self.editor_scrollbar = ttk.Scrollbar(left_host, orient="vertical",
                                               command=self.editor_canvas.yview)
        self.editor_canvas.configure(yscrollcommand=self.editor_scrollbar.set)
        self.editor_scrollbar.pack(side="right", fill="y")
        self.editor_canvas.pack(side="left", fill="both", expand=True)
        left = ttk.Frame(self.editor_canvas)
        self.editor_window = self.editor_canvas.create_window((0, 0), window=left,
                                                               anchor="nw")
        left.bind("<Configure>", lambda event: self.editor_canvas.configure(
            scrollregion=self.editor_canvas.bbox("all")))
        self.editor_canvas.bind("<Configure>", lambda event: self.editor_canvas.itemconfigure(
            self.editor_window, width=event.width))
        # Mouse-wheel events are routed only while the pointer is inside this pane.
        self.bind_all("<MouseWheel>", self._editor_mousewheel, add="+")
        self.bind_all("<Button-4>", self._editor_mousewheel, add="+")
        self.bind_all("<Button-5>", self._editor_mousewheel, add="+")

        # Left editor: scaffold, sites, fragments, assignments.
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
        self.scaffold_display = tk.Canvas(a, bg="white", height=185,
                                          highlightthickness=1, highlightbackground="#A8BBCF",
                                          cursor="crosshair")
        self.scaffold_display.pack(fill="x", pady=(0, 1))
        self.scaffold_display.bind("<Button-1>", self.select_scaffold_atom)
        self.scaffold_display.bind("<Configure>", lambda event: self.draw_scaffold())
        ttk.Label(a, text="Click a highlighted atom or choose its index below.").pack(anchor="w")

        b = ttk.LabelFrame(left, text="2  Attachment site", padding=7)
        b.pack(fill="x", pady=(0, 5))
        self.site_combo = ttk.Combobox(b, state="readonly", textvariable=self.site_choice)
        self.site_combo.pack(side="left", fill="x", expand=True)
        self.site_combo.bind("<<ComboboxSelected>>", lambda event: self.draw_scaffold())
        ttk.Button(b, text="Define R-site", command=self.define_site,
                   style="Primary.TButton").pack(side="left", padx=(6,0))

        c = ttk.LabelFrame(left, text="3  Fragments and assignments", padding=7)
        c.pack(fill="both", expand=True, pady=(0, 5))
        fragmentbar = ttk.Frame(c)
        fragmentbar.pack(fill="x", pady=(0, 4))
        ttk.Button(fragmentbar, text="Import .smi", command=self.import_fragments).pack(side="left")
        ttk.Button(fragmentbar, text="+ Fragment", command=self.add_fragment).pack(side="left", padx=4)
        ttk.Button(fragmentbar, text="Add R-group", command=self.assign, style="Primary.TButton").pack(side="right")
        listframe = ttk.Frame(c)
        listframe.pack(fill="both", expand=True)
        self.fragment_list = tk.Listbox(listframe, height=6, selectmode="extended",
                                        exportselection=False, activestyle="dotbox",
                                        font=("Segoe UI", 10, "bold"), bg="#FFFFFF", fg="#10243B",
                                        selectbackground="#0864B1", selectforeground="#FFFFFF",
                                        highlightbackground="#A8BBCF", highlightthickness=1, relief="solid", bd=1)
        self.fragment_list.pack(side="left", fill="both", expand=True)
        fragscroll = ttk.Scrollbar(listframe, orient="vertical", command=self.fragment_list.yview)
        fragscroll.pack(side="right", fill="y")
        self.fragment_list.configure(yscrollcommand=fragscroll.set)
        self.refresh_fragment_list()
        ttk.Label(c, text="R-group assignments").pack(anchor="w", pady=(5, 1))
        assignmentrow = ttk.Frame(c)
        assignmentrow.pack(fill="x")
        self.assignment_list = tk.Listbox(assignmentrow, height=3, font=("Segoe UI", 10, "bold"), bg="#FFFFFF", fg="#10243B",
                                         selectbackground="#0864B1", selectforeground="#FFFFFF",
                                         highlightbackground="#A8BBCF", highlightthickness=1,
                                         relief="solid", bd=1)
        self.assignment_list.pack(side="left", fill="x", expand=True)
        self.assignment_list.bind("<<ListboxSelect>>", self.show_assignment_details)
        ttk.Button(assignmentrow, text="Remove", command=self.remove_assignment).pack(side="left", padx=(5, 0))
        self.assignment_details = ttk.Label(c, text="Choose an R-site and assign fragments.",
                                            wraplength=360, justify="left")
        self.assignment_details.pack(fill="x", pady=(4, 2))

        # Right: a compact results table with room for an actual molecule preview.
        # Primary actions are pinned to the TOP of the right workspace.
        # They remain visible even when Windows display scaling clips the left pane.
        controls = ttk.Frame(right)
        controls.pack(fill="x", pady=(0, 8))
        self.enumerate_all_button = ttk.Button(controls, text="ENUMERATE ALL COMBINATIONS",
                    command=self.enumerate, style="Primary.TButton")
        self.enumerate_all_button.pack(side="left", padx=(0, 6))
        self.enumerate_current_button = ttk.Button(controls, text="Enumerate Current",
                    command=self.enumerate_current, style="Primary.TButton")
        self.enumerate_current_button.pack(side="left", padx=(0, 10))
        ttk.Label(controls, text="Max:").pack(side="left")
        ttk.Entry(controls, textvariable=self.limit, width=7).pack(side="left", padx=(4, 0))
        ttk.Label(controls, text="Select fragments, Add R-group, then enumerate.",
                  foreground="#244766").pack(side="right")
        resulthead = ttk.Frame(right)
        resulthead.pack(fill="x", pady=(0, 5))
        ttk.Label(resulthead, text="Enumerated molecules", font=("Segoe UI", 10, "bold")).pack(side="left")
        self.results_count = ttk.Label(resulthead, text="0 compounds")
        self.results_count.pack(side="right")
        results = ttk.Frame(right)
        results.pack(fill="x")
        cols = ("id", "smiles", "mw", "logp", "tpsa", "qed")
        self.table = ttk.Treeview(results, columns=cols, show="headings", height=7)
        labels = {"id":"ID", "smiles":"Canonical SMILES", "mw":"MW", "logp":"cLogP", "tpsa":"TPSA", "qed":"QED"}
        for col in cols:
            self.table.heading(col, text=labels[col])
            self.table.column(col, width=98 if col!="smiles" else 260,
                              minwidth=62, stretch=(col=="smiles"))
        self.table.pack(side="left", fill="x", expand=True)
        scrollbar=ttk.Scrollbar(results, orient="vertical", command=self.table.yview)
        scrollbar.pack(side="right", fill="y")
        self.table.configure(yscrollcommand=scrollbar.set)
        self.table.bind("<<TreeviewSelect>>", self.choose_record)

        # Critical export actions stay in view regardless of preview content size.
        actions = ttk.Frame(right)
        actions.pack(fill="x", pady=(7, 4))
        ttk.Checkbutton(actions, text="Selected only", variable=self.export_selection_only).pack(side="left", padx=(0, 6))
        ttk.Button(actions, text="Export CSV", command=lambda:self.export("csv")).pack(side="left", padx=3)
        ttk.Button(actions, text="SMILES", command=lambda:self.export("smi")).pack(side="left", padx=3)
        ttk.Button(actions, text="2D SDF", command=lambda:self.export("2d")).pack(side="left", padx=3)
        ttk.Button(actions, text="3D SDF", command=lambda:self.export("3d")).pack(side="left", padx=3)
        ttk.Button(actions, text="View 3D", command=self.view3d,
                   style="Primary.TButton").pack(side="left", padx=3)

        preview = ttk.LabelFrame(right, text="Molecular preview — complete structure", padding=6)
        preview.pack(fill="both", expand=True, pady=(3, 2))
        self.properties=ttk.Label(preview, text="Load a scaffold to begin.",
                                  font=("Consolas", 10, "bold"), justify="left")
        self.properties.pack(side="bottom", fill="x", anchor="w", padx=8, pady=(2,4))
        self.preview_image=tk.Canvas(preview, bg="#FFFFFF", highlightthickness=0,
                                     borderwidth=0, height=190)
        self.preview_image.pack(side="top", fill="both", expand=True)
        self.preview_image.bind("<Configure>", self.on_preview_resize)
        tk.Label(self, textvariable=self.status, bg="#e2eaf4", fg="#10243B",
                 anchor="w", padx=9, pady=3, font=("Segoe UI", 9)).pack(fill="x")

    def _set_auto_update_preference(self):
        try:
            set_auto_check(self.auto_check_updates.get())
            self.status.set("Automatic update checks " +
                            ("enabled" if self.auto_check_updates.get() else "disabled"))
        except Exception as exc:
            self.error("Update preference failed", exc)

    def check_for_updates(self, silent=False):
        """Check without freezing the Tk interface; show no dialog if already current."""
        if self.update_busy:
            if not silent:
                self.status.set("An update operation is already running.")
            return
        self.update_busy = True
        self.check_update_button.configure(state="disabled")
        self.status.set("Checking GitHub Releases for updates...")
        def worker():
            try:
                self.update_events.put(("check", silent, latest_update()))
            except Exception as exc:
                self.update_events.put(("check_error", silent, str(exc)))
        threading.Thread(target=worker, daemon=True).start()
        self.after(120, self._poll_updates)

    def _poll_updates(self):
        try:
            kind, *payload = self.update_events.get_nowait()
        except queue.Empty:
            if self.update_busy:
                self.after(120, self._poll_updates)
            return
        if kind == "check":
            silent, info = payload
            if info is None:
                self.status.set(f"RGroupStudio v{APP_VERSION} is up to date.")
                if not silent:
                    messagebox.showinfo("No updates", f"You have the latest version: v{APP_VERSION}.", parent=self)
                self._update_finished()
            else:
                self._update_info = info
                self.status.set("RGroupStudio v" + info["version"] + " is available.")
                self._update_finished()
                if messagebox.askyesno("RGroupStudio Update Available",
                        "Installed: v" + APP_VERSION + "\nAvailable: v" + info["version"] +
                        "\n\nDownload the verified installer and update now?\n"
                        "Your app will close when installation starts.", parent=self):
                    self._download_update(info)
        elif kind == "check_error":
            silent, message = payload
            self.status.set("Update check unavailable: " + message)
            self._update_finished()
            if not silent:
                messagebox.showwarning("Update check failed", message +
                                       "\n\nCheck your internet connection and try again.", parent=self)
        elif kind == "progress":
            received, total = payload
            self.status.set("Downloading verified update: %d%%" % (100 * received // max(total,1)))
        elif kind == "downloaded":
            path = payload[0]
            self._update_finished()
            if not messagebox.askyesno("Install Verified Update",
                "Installer downloaded and SHA-256 verified.\n\n"
                "Close RGroupStudio and run the official installer now?", parent=self):
                self.status.set("Verified installer saved at " + str(path))
                return
            try:
                if sys.platform != "win32":
                    raise RuntimeError("Automatic installer launch is supported on Windows only.")
                flags = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
                subprocess.Popen([str(path), "/SILENT", "/SUPPRESSMSGBOXES",
                                  "/NORESTART", "/CLOSEAPPLICATIONS", "/RESTARTAPPLICATIONS"],
                                 cwd=str(path.parent), creationflags=flags, close_fds=True)
                self.status.set("Installer launched; closing RGroupStudio...")
                self.after(400, self.destroy)
            except Exception as exc:
                self.error("Update installation failed", exc)
        elif kind == "download_error":
            message = payload[0]
            self.status.set("Update download failed: " + message)
            self._update_finished()
            messagebox.showerror("Update download failed", message, parent=self)
        if self.update_busy:
            self.after(120, self._poll_updates)

    def _update_finished(self):
        self.update_busy = False
        self.check_update_button.configure(state="normal")

    def _download_update(self, info):
        self.update_busy = True
        self.check_update_button.configure(state="disabled")
        self.status.set("Downloading and verifying official installer...")
        def worker():
            try:
                def progress(received,total):
                    self.update_events.put(("progress", received, total))
                path = download_installer(info, progress=progress)
                self.update_events.put(("downloaded", path))
            except Exception as exc:
                self.update_events.put(("download_error", str(exc)))
        threading.Thread(target=worker, daemon=True).start()
        self.after(120, self._poll_updates)

    def _editor_mousewheel(self, event):
        """Scroll the entire assignment editor when the pointer is inside it."""
        canvas = getattr(self, "editor_canvas", None)
        if canvas is None:
            return
        px, py = self.winfo_pointerxy()
        cx, cy = canvas.winfo_rootx(), canvas.winfo_rooty()
        if not (cx <= px < cx + canvas.winfo_width() and
                cy <= py < cy + canvas.winfo_height()):
            return
        if getattr(event, "num", None) == 4:
            steps = -1
        elif getattr(event, "num", None) == 5:
            steps = 1
        else:
            delta = getattr(event, "delta", 0)
            steps = -max(1, abs(delta)//120) if delta > 0 else max(1, abs(delta)//120)
            if delta == 0:
                return
        canvas.yview_scroll(steps, "units")
        return "break"

    def error(self, title, exc):
        logging.exception("%s: %s", title, exc)
        messagebox.showerror(title, str(exc) + "\n\nLog: " + str(LOGFILE))

    def reset_results(self):
        self.records = []
        self.selected = None
        self.table.delete(*self.table.get_children())
        self.results_count.configure(text="0 compounds")
        # Show the loaded scaffold immediately instead of an empty preview.
        if hasattr(self, "mol"):
            image = ImageTk.PhotoImage(Draw.MolToImage(self.mol, size=(360, 245)))
            self.preview_image.configure(image=image)
            self.preview_image.image = image
            props = molecular_properties(self.mol)
            self.properties.configure(text="PARENT SCAFFOLD\n\n" +
                "\n".join(f"{k:>7}: {v}" for k, v in props.items()))
        else:
            self.preview_image.configure(image="")
            self.preview_image.image = None
            self.properties.configure(text="Load a scaffold to begin.")

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
