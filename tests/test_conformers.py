"""Real RDKit geometry checks plus force-field failure handling."""
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem

from app.conformers import calculate_conformers


class ConformerTests(unittest.TestCase):
    def test_real_ethanol_geometry_energy_order_and_sdf_are_consistent(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "conformers.sdf"
            result = calculate_conformers("CCO", path, num_confs=4)
            self.assertEqual(result["status"], "completed")
            self.assertGreater(result["converged"], 0)
            self.assertEqual(result["converged"], len(result["items"]))
            self.assertEqual(len(result["atoms"]), 9)
            self.assertEqual(result["items"][0]["relative_energy_kcal_mol"], 0)
            self.assertTrue(all(row["converged"] for row in result["items"]))
            for row in result["items"]:
                self.assertEqual(len(row["coordinates"]), len(result["atoms"]))
                self.assertTrue(all(len(point) == 3 and all(math.isfinite(value) for value in point) for point in row["coordinates"]))
                self.assertGreaterEqual(row["relative_energy_kcal_mol"], 0)
            molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=False) if mol is not None]
            self.assertEqual(len(molecules), len(result["items"]))
            self.assertEqual(molecules[0].GetNumAtoms(), 9)
            self.assertAlmostEqual(float(molecules[0].GetProp("MMFF94_energy_kcal_mol")), result["items"][0]["energy_kcal_mol"])
            self.assertGreater(result["shape_descriptors"]["radius_of_gyration"], 0)
            self.assertTrue(any("不是量子" in warning for warning in result["warnings"]))

    def test_large_or_highly_flexible_molecule_is_skipped_before_embedding(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch("app.conformers.AllChem.EmbedMultipleConfs") as embed:
                result = calculate_conformers("C" * 61, Path(directory) / "large.sdf")
            embed.assert_not_called()
            self.assertEqual(result["status"], "skipped")
            self.assertEqual(result["items"], [])

    def test_absent_force_field_parameters_do_not_become_fake_geometry(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch("app.conformers.AllChem.MMFFHasAllMoleculeParams", return_value=False):
                result = calculate_conformers("CCO", Path(directory) / "missing.sdf")
            self.assertEqual(result["status"], "skipped")
            self.assertEqual(result["reason"], "missing_mmff94_parameters")

    def test_no_converged_conformer_has_no_lowest_energy_or_download(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "unconverged.sdf"
            def unconverged(mol, **kwargs):
                return [(1, -100.0) for _ in mol.GetConformers()]
            with patch("app.conformers.AllChem.MMFFOptimizeMoleculeConfs", side_effect=unconverged):
                result = calculate_conformers("CCO", path)
            self.assertEqual(result["status"], "failed")
            self.assertEqual(result["reason"], "no_converged_conformers")
            self.assertIsNone(result["lowest_energy_id"])
            self.assertEqual(result["shape_descriptors"], {})
            self.assertEqual(result["items"], [])
            self.assertGreater(result["unconverged_count"], 0)
            self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main()
