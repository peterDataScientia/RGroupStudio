# Third-party components — RGroupStudio v0.4.0

RGroupStudio bundles independent scientific and UI software. The full applicable upstream license texts and source distribution requirements must be honored when redistributing the Windows installer.

- **3Dmol.js 2.5.5** — WebGL molecular viewer, BSD-3-Clause. The Windows build copies the pinned package's LICENSE file to `viewer/3DMOL-LICENSE.txt` when present. Upstream: https://github.com/3dmol/3Dmol.js
- **RDKit** — chemistry, ETKDG conformers, molecular descriptors and SDF handling; BSD-3-Clause. https://github.com/rdkit/rdkit
- **PySide6 / Qt WebEngine** — native UI and embedded Chromium-based viewer. Qt for Python components are available under LGPL and other applicable license options. See Qt's license documentation and bundled package license texts. https://doc.qt.io/qtforpython-6/licenses.html
- **Pillow** — RDKit 2D image handling; HPND. https://python-pillow.github.io/
- **NumPy** — scientific array operations; BSD-3-Clause. https://numpy.org/

RGroupStudio does **not** fetch a 3D library or upload structures while rendering. Its 3Dmol.js script is downloaded once from npm during the controlled Windows build and bundled locally.

All LGPL-related redistribution obligations, notices and replacement/library linking rights must be considered before further commercial redistribution. RGroupStudio is not a replacement for an instrument-validated or experimental molecular structure determination.
