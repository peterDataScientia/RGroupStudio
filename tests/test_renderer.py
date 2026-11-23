"""Image-level scientific 3D viewer regression tests (no GUI required)."""
import unittest
import numpy as np
from rdkit import Chem
from chemistry import conformer_3d
from molecular_renderer import render_molecule, orientation

class MolecularRendererTests(unittest.TestCase):
    def setUp(self):
        # Match the user's substituted aromatic example, plus heteroatom ligands.
        self.smiles = ["CC(C)c1ccc(O)cc1", "COc1ccc(Cl)cc1",
                       "CC(=O)Nc1ccc(F)cc1", "c1ccccc1"]

    def test_molecules_have_well_spread_visible_3d_bonds(self):
        for smiles in self.smiles:
            with self.subTest(smiles=smiles):
                molecule = conformer_3d(Chem.MolFromSmiles(smiles))
                image = render_molecule(molecule,width=720,height=450)
                pixels = np.asarray(image)
                self.assertEqual(pixels.shape,(450,720,3))
                # Molecular colors must occupy much more than a tiny central dot.
                colored = (pixels.min(axis=2)<145) & (pixels.max(axis=2)-pixels.min(axis=2)>8)
                yy, xx = np.where(colored)
                self.assertGreater(len(xx),500)
                self.assertGreater(xx.max()-xx.min(),200)
                self.assertGreater(yy.max()-yy.min(),60)

    def test_rotation_and_zoom_change_real_image(self):
        molecule = conformer_3d(Chem.MolFromSmiles("CC(C)c1ccc(O)cc1"))
        first = np.asarray(render_molecule(molecule,width=500,height=330,yaw=35,pitch=22))
        turned = np.asarray(render_molecule(molecule,width=500,height=330,yaw=105,pitch=12))
        enlarged = np.asarray(render_molecule(molecule,width=500,height=330,yaw=35,pitch=22,zoom=1.4))
        self.assertGreater(np.mean(np.abs(first.astype(float)-turned.astype(float))),1.5)
        self.assertGreater(np.mean(np.abs(first.astype(float)-enlarged.astype(float))),1.5)

    def test_render_does_not_change_chemical_structure(self):
        mol=conformer_3d(Chem.MolFromSmiles("CC(O)c1ccccc1"))
        smiles_before=Chem.MolToSmiles(Chem.RemoveHs(mol))
        xyz=mol.GetConformer().GetPositions().copy()
        render_molecule(mol,show_hydrogens=False)
        render_molecule(mol,show_hydrogens=True,width=460,height=290)
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(mol)),smiles_before)
        np.testing.assert_allclose(mol.GetConformer().GetPositions(),xyz)

    def test_projection_is_valid(self):
        xyz=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        projected=orientation(xyz,35,22)
        np.testing.assert_allclose(np.linalg.norm(projected,axis=1),np.ones(3),atol=1e-8)

if __name__=="__main__":
    unittest.main()
