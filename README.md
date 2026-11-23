# RGroupStudio — Molecular Library Designer

RGroupStudio is a **native Windows desktop** RDKit molecular library enumeration application.

## v0.4.0 desktop architecture

- **PySide6 / Qt** for scaffold and fragment editing, tabular enumeration, exports, and update controls.
- **3Dmol.js 2.5.5 with Qt WebEngine** for 3D molecules rendered **inside the desktop window**. WebGL offers rotate, zoom, fit, atom labels, and stick/ball-and-stick/space-filling representations.
- **RDKit** for SMILES parsing, R-group attachment, enumerated products, descriptors, 2D SDF, and ETKDGv3 3D SDF.
- The browser bundle is packaged **offline inside the installer**. No molecular data are uploaded.
- The update checker contacts official GitHub Releases only if enabled and verifies installer SHA-256 before installation.

## Windows installation

Download the official `RGroupStudio_Setup_v0.4.0_Windows_x64.exe` from the GitHub Releases page. Python is not needed.

### Workflow

1. Enter a scaffold SMILES and click **Draw**.
2. Click a highlighted atom or choose an atom index, then click **Define R-site**.
3. Select one or more fragments and click **Assign R-group**. Repeat for multiple sites.
4. Click **Enumerate Current** or **ENUMERATE ALL COMBINATIONS**.
5. Select a product in the table. Use **2D Structure** or **3D Molecule** tabs.
6. In 3D, rotate by dragging, zoom with the wheel, change representation, click **Fit**, or enlarge using **Fullscreen**.
7. Export selected or all products as CSV, SMILES, 2D SDF, or 3D SDF.

## Build from source

This repository keeps the previous Tk application in `app.py` for legacy comparison.
The **new installed desktop entry point is `app_qt.py`**.

Install Python 3.11 and run:

```powershell
python -m pip install -r requirements.txt
npm install --prefix .3dmol_npm --ignore-scripts --no-audit --no-fund --no-save 3dmol@2.5.5
Copy-Item .3dmol_npm/node_modules/3dmol/build/3Dmol-min.js viewer/3Dmol-min.js
python app_qt.py
```

The Windows CI workflow performs these steps and bundles the local assets with PyInstaller and Inno Setup.

## Graphics and scientific interpretation

3Dmol.js uses WebGL and therefore requires a working graphics driver and WebGL support. On machines without GPU support, the viewer may report that WebGL cannot initialize; install appropriate graphics drivers. The app does not use unsafe software-rendering flags by default.

RDKit's ETKDGv3 conformer is computationally generated, not an experimental crystal structure or guaranteed global energy minimum. MMFF/UFF optimization outcomes are recorded in the molecular metadata. Conformer quality should be checked before scientific conclusions.

Third-party license acknowledgments are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
