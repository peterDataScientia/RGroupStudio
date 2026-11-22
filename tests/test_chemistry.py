import os
import tempfile
import unittest
from rdkit import Chem
from chemistry import (FRAGMENTS, parse_scaffold, sites, validate_fragment,
                       enumerate_library, molecular_properties, export_sdf,
                       export_csv, export_smi, conformer_3d)

class ChemistryTests(unittest.TestCase):
    def test_scaffold(self):
        self.assertEqual(parse_scaffold("c1ccccc1").GetNumAtoms(), 6)
    def test_sites(self):
        self.assertEqual(len(sites(parse_scaffold("c1ccccc1"))),6)
    def test_fragments(self):
        self.assertIsNotNone(validate_fragment("[*]OC"))
        with self.assertRaises(ValueError):
            validate_fragment("CC")
    def test_enumeration(self):
        records=enumerate_library("c1ccccc1",{0:["[*]F","[*]Cl"],2:["[*]C"]})
        self.assertEqual(len(records),2)
        for _,mol in records:
            Chem.SanitizeMol(mol)
    def test_duplicates(self):
        records=enumerate_library("c1ccccc1",{0:["[*]F","[*]F"]})
        self.assertEqual(len(records),1)
    def test_properties(self):
        p=molecular_properties(parse_scaffold("c1ccccc1"))
        self.assertAlmostEqual(p["MW"],78.11,places=1)
    def test_3d(self):
        mol=conformer_3d(parse_scaffold("CCO"))
        self.assertTrue(mol.GetNumConformers()>0)
    def test_exports(self):
        records=enumerate_library("c1ccccc1",{0:["[*]Cl"]})
        with tempfile.TemporaryDirectory() as d:
            for fmt,fn in (("csv",export_csv),("smi",export_smi)):
                path=os.path.join(d,"test."+fmt)
                fn(path,records)
                self.assertTrue(os.path.getsize(path)>0)
            path=os.path.join(d,"test.sdf")
            self.assertEqual(export_sdf(path,records),[])
            self.assertEqual(len([m for m in Chem.SDMolSupplier(path) if m]),1)
if __name__=="__main__":
    unittest.main()
