"""Audit real release artifacts and held-out evidence independently of training.

These tests skip only while one or more trained artifacts are absent. If all
artifacts exist, missing/inconsistent reports are failures, not skipped audits.
"""

from __future__ import annotations

import itertools
import hashlib
import json
import math
import pickle
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from threadpoolctl import threadpool_limits

from app.chemistry import FEATURE_NAMES as DOMAIN_FEATURES, calculate_descriptors, identity_smiles, parse_smiles
from app.features import feature_vector, fingerprint
from app.service import PredictionService, TARGET_IDS

ROOT = Path(__file__).resolve().parents[1]
MISSING_MODELS = [target for target in TARGET_IDS if not (ROOT / "models" / f"{target}.pkl").is_file()]


@unittest.skipIf(bool(MISSING_MODELS), "Trained artifacts are absent: " + ", ".join(MISSING_MODELS))
class TrainedArtifactAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts = {}
        for target in TARGET_IDS:
            with (ROOT / "models" / f"{target}.pkl").open("rb") as handle:
                bundle = pickle.load(handle)
            report_dir = ROOT / "reports" / target
            manifest = pd.read_csv(report_dir / "split_manifest.csv", keep_default_na=False)
            predictions = pd.read_csv(report_dir / "test_predictions.csv", keep_default_na=False)
            metadata = json.loads((report_dir / "metrics.json").read_text(encoding="utf-8"))
            cls.artifacts[target] = {"bundle": bundle, "manifest": manifest, "predictions": predictions, "metadata": metadata}
        cls.molecules = {}
        cls.recomputed_predictions = {}

    @classmethod
    def molecule(cls, smiles):
        if smiles not in cls.molecules:
            cls.molecules[smiles] = parse_smiles(smiles)
        return cls.molecules[smiles]

    @classmethod
    def predict_rows(cls, target, smiles):
        bundle = cls.artifacts[target]["bundle"]
        matrix = np.stack([feature_vector(cls.molecule(value), bundle["feature_mode"]) for value in smiles])
        with threadpool_limits(limits=2):
            return np.asarray(bundle["model"].predict(matrix), dtype=float)

    def test_fit_calibration_and_test_are_disjoint_by_identity_and_scaffold(self):
        for target, data in self.artifacts.items():
            with self.subTest(target=target):
                manifest = data["manifest"]
                self.assertEqual(set(manifest["split"]), {"train", "validation", "calibration", "test"})
                self.assertFalse(manifest["canonical_smiles"].eq("").any())
                for row in manifest.itertuples(index=False):
                    mol = self.molecule(row.canonical_smiles)
                    self.assertEqual(row.identity_smiles, identity_smiles(mol))
                    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
                    self.assertEqual(row.scaffold, scaffold)
                partitions = {
                    "fit": manifest[manifest["split"].isin(["train", "validation"])],
                    "calibration": manifest[manifest["split"] == "calibration"],
                    "test": manifest[manifest["split"] == "test"],
                }
                for partition in partitions.values():
                    self.assertGreater(len(partition), 0)
                for left, right in itertools.combinations(partitions, 2):
                    self.assertFalse(set(partitions[left]["identity_smiles"]) & set(partitions[right]["identity_smiles"]), (target, left, right, "identity overlap"))
                    self.assertFalse(set(partitions[left]["scaffold"]) & set(partitions[right]["scaffold"]), (target, left, right, "scaffold overlap"))
                # Selection validation is also isolated from initial fitting.
                train = manifest[manifest["split"] == "train"]
                validation = manifest[manifest["split"] == "validation"]
                self.assertFalse(set(train["identity_smiles"]) & set(validation["identity_smiles"]))
                self.assertFalse(set(train["scaffold"]) & set(validation["scaffold"]))

    def test_artifact_similarity_records_contain_fit_molecules_only(self):
        for target, data in self.artifacts.items():
            with self.subTest(target=target):
                bundle, manifest = data["bundle"], data["manifest"]
                PredictionService._validate_bundle(target, bundle)
                fit = manifest[manifest["split"].isin(["train", "validation"])]
                self.assertEqual(set(bundle["training_smiles"]), set(fit["canonical_smiles"]))
                self.assertEqual(set(bundle["training_identities"]), set(fit["identity_smiles"]))
                self.assertEqual(len(bundle["training_smiles"]), len(fit))
                for smiles, identity, encoded in zip(bundle["training_smiles"], bundle["training_identities"], bundle["training_fingerprints"]):
                    mol = self.molecule(smiles)
                    self.assertEqual(identity, identity_smiles(mol))
                    self.assertEqual(encoded, DataStructs.BitVectToBinaryText(fingerprint(mol)))

    def test_report_hashes_and_candidate_selection_match_saved_evidence(self):
        source_hash = hashlib.sha256((ROOT / "data" / "curated_observations.csv").read_bytes()).hexdigest()
        for target, data in self.artifacts.items():
            with self.subTest(target=target):
                metadata = data["metadata"]
                report_dir = ROOT / "reports" / target
                self.assertEqual(metadata["data_sha256"], source_hash)
                self.assertEqual(metadata["split_sha256"], hashlib.sha256((report_dir / "split_manifest.csv").read_bytes()).hexdigest())
                self.assertEqual(metadata["model_sha256"], hashlib.sha256((ROOT / "models" / f"{target}.pkl").read_bytes()).hexdigest())
                candidates = pd.read_csv(report_dir / "candidate_validation.csv")
                self.assertEqual(len(candidates), 7)
                selection = json.loads((report_dir / "ensemble_selection.json").read_text(encoding="utf-8"))
                self.assertEqual(metadata["ensemble"]["member_names"], selection["member_names"])
                self.assertEqual(metadata["ensemble"]["weights"], selection["weights"])
                self.assertEqual(metadata["feature_mode"], "combined")
                for metric in ("rmse", "mae", "r2"):
                    self.assertAlmostEqual(metadata["metrics"]["validation_selected"][metric], selection["metrics"][metric], places=9)

    def test_domain_bounds_contain_every_raw_training_descriptor(self):
        for target, data in self.artifacts.items():
            with self.subTest(target=target):
                bundle = data["bundle"]
                raw = [calculate_descriptors(self.molecule(smiles)) for smiles in bundle["training_smiles"]]
                for name in DOMAIN_FEATURES:
                    values = np.asarray([record[name] for record in raw], dtype=float)
                    lower, upper = bundle["descriptor_bounds"][name]
                    self.assertTrue(np.all(values >= lower), (target, name, "raw lower bound exclusion"))
                    self.assertTrue(np.all(values <= upper), (target, name, "raw upper bound exclusion"))
                    self.assertEqual(lower, float(values.min()))
                    self.assertEqual(upper, float(values.max()))

    def test_frozen_logS_test_preserves_all_original_224_molecules_and_values(self):
        original = pd.read_csv(ROOT / "reports" / "v1_test_predictions.csv", keep_default_na=False)
        current = self.artifacts["logS"]["predictions"]
        self.assertEqual(len(original), 224)
        self.assertEqual(len(current), 224)
        self.assertTrue(original["canonical_smiles"].is_unique)
        self.assertTrue(current["canonical_smiles"].is_unique)
        self.assertEqual(set(original["canonical_smiles"]), set(current["canonical_smiles"]))
        original = original.set_index("canonical_smiles").sort_index()
        current = current.set_index("canonical_smiles").sort_index()
        np.testing.assert_allclose(current["y_true"].to_numpy(), original["observed"].to_numpy(), rtol=0, atol=1e-12)

    def test_all_saved_test_predictions_recompute_from_the_release_artifact(self):
        for target, data in self.artifacts.items():
            with self.subTest(target=target):
                predictions, manifest = data["predictions"], data["manifest"]
                test = manifest[manifest["split"] == "test"]
                self.assertEqual(len(test), len(predictions))
                self.assertEqual(set(test["canonical_smiles"]), set(predictions["canonical_smiles"]))
                self.assertTrue(predictions["canonical_smiles"].is_unique)
                recomputed = self.predict_rows(target, predictions["canonical_smiles"])
                self.recomputed_predictions[target] = recomputed
                np.testing.assert_allclose(predictions["y_pred"].to_numpy(dtype=float), recomputed, rtol=0, atol=1e-9)
                indexed_test = test.set_index("canonical_smiles")
                np.testing.assert_allclose(predictions["y_true"].to_numpy(dtype=float), indexed_test.loc[predictions["canonical_smiles"], "value"].to_numpy(dtype=float), rtol=0, atol=1e-12)
                radius = float(data["bundle"]["conformal_radius"])
                np.testing.assert_allclose(predictions["lower90"].to_numpy(dtype=float), recomputed - radius, rtol=0, atol=1e-9)
                np.testing.assert_allclose(predictions["upper90"].to_numpy(dtype=float), recomputed + radius, rtol=0, atol=1e-9)

    def test_heldout_metrics_baseline_and_empirical_coverage_recompute(self):
        for target, data in self.artifacts.items():
            with self.subTest(target=target):
                frame, metadata = data["predictions"], data["metadata"]
                truth = frame["y_true"].to_numpy(dtype=float)
                predicted = self.recomputed_predictions.get(target)
                if predicted is None:
                    predicted = self.predict_rows(target, frame["canonical_smiles"])
                residual = truth - predicted
                expected = {
                    "rmse": float(np.sqrt(np.mean(residual**2))),
                    "mae": float(np.mean(np.abs(residual))),
                    "r2": float(1 - np.sum(residual**2) / np.sum((truth - np.mean(truth))**2)),
                }
                for name, value in expected.items():
                    self.assertAlmostEqual(metadata["metrics"]["test"][name], value, places=9)
                    self.assertAlmostEqual(data["bundle"]["metadata"]["metrics"]["test"][name], value, places=9)
                radius = float(data["bundle"]["conformal_radius"])
                coverage = float(np.mean((truth >= predicted - radius) & (truth <= predicted + radius)))
                self.assertAlmostEqual(metadata["calibration"]["test_coverage"], coverage, places=12)
                self.assertAlmostEqual(metadata["calibration"]["radius"], radius, places=12)
                self.assertEqual(metadata["calibration"]["nominal_coverage"], 0.9)
                fit_values = data["manifest"].loc[data["manifest"]["split"].isin(["train", "validation"]), "value"].to_numpy(dtype=float)
                baseline_errors = truth - np.mean(fit_values)
                baseline = {"rmse": float(np.sqrt(np.mean(baseline_errors**2))), "mae": float(np.mean(np.abs(baseline_errors))), "r2": float(1 - np.sum(baseline_errors**2) / np.sum((truth - np.mean(truth))**2))}
                for name, value in baseline.items():
                    self.assertAlmostEqual(metadata["metrics"]["baseline"][name], value, places=9)

    def test_calibration_radius_uses_only_the_independent_calibration_rows(self):
        for target, data in self.artifacts.items():
            with self.subTest(target=target):
                calibration = data["manifest"].loc[data["manifest"]["split"] == "calibration"]
                self.assertEqual(data["metadata"]["calibration"]["n_calibration"], len(calibration))
                predicted = self.predict_rows(target, calibration["canonical_smiles"])
                residuals = np.sort(np.abs(calibration["value"].to_numpy(dtype=float) - predicted))
                rank = min(len(residuals), math.ceil((len(residuals) + 1) * 0.9))
                self.assertAlmostEqual(float(data["bundle"]["conformal_radius"]), float(residuals[rank - 1]), places=9)


if __name__ == "__main__":
    unittest.main()
