"""V3 representation and validation-only ensemble release evidence checks."""
import csv
import hashlib
import inspect
import itertools
import json
import pickle
import runpy
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from rdkit.Chem import Descriptors
from sklearn.dummy import DummyRegressor
from sklearn.impute import SimpleImputer
from app.chemistry import calculate_descriptors, parse_smiles
from app.ensemble import FeatureSubset, ValidationEnsemble
from app.features import FEATURE_NAMES, FEATURE_VERSION, LEGACY_FEATURE_NAMES, _safe_descriptor, feature_vector, model_feature_names

ROOT = Path(__file__).resolve().parents[1]
TRAINING = runpy.run_path(str(ROOT / "scripts/train_models.py"))
TARGETS = ("logS", "logD74", "hydration_free_energy")


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


class V3FeatureTests(unittest.TestCase):
    def test_fixed_order_unique_names_and_core_values(self):
        self.assertEqual(len(FEATURE_NAMES), 219)
        self.assertEqual(len(set(FEATURE_NAMES)), 219)
        self.assertEqual(FEATURE_NAMES[:28], LEGACY_FEATURE_NAMES)
        self.assertEqual(len(model_feature_names("combined")), 2267)
        mol = parse_smiles("CCO")
        vector = feature_vector(mol, "combined")
        self.assertEqual(vector.dtype, np.float32)
        for i, (name, value) in enumerate(calculate_descriptors(mol).items()):
            self.assertAlmostEqual(float(vector[FEATURE_NAMES.index(name)]), value, places=5)
        self.assertTrue(np.isin(vector[219:], [0, 1]).all())

    def test_equivalent_smiles_have_identical_vectors(self):
        np.testing.assert_allclose(feature_vector(parse_smiles("CCO"), "combined"),
                                   feature_vector(parse_smiles("OCC"), "combined"), equal_nan=True)

    def test_nonfinite_or_failed_extra_descriptor_is_missing(self):
        mol = parse_smiles("CCO")
        for value in [float("nan"), float("inf"), -float("inf"), 1e100]:
            self.assertTrue(np.isnan(_safe_descriptor(lambda m: value, mol)))
        self.assertTrue(np.isnan(_safe_descriptor(lambda m: 1 / 0, mol)))
        self.assertEqual(_safe_descriptor(lambda m: 2.5, mol), 2.5)

    def test_imputation_uses_fitting_rows_and_retains_empty_columns(self):
        original = np.array([[1., np.nan], [3., np.nan]], dtype=np.float32)
        selector = FeatureSubset((0, 1)).fit(original)
        selected = selector.transform(original)
        imputer = SimpleImputer(strategy="median", keep_empty_features=True, copy=False).fit(selected)
        validation = np.array([[np.nan, 100.], [1000., np.nan]], dtype=np.float32)
        result = imputer.transform(selector.transform(validation))
        self.assertEqual(result[0, 0], 2.)
        self.assertEqual(result[1, 1], 0.)
        self.assertTrue(np.isnan(original[:, 1]).all())
        self.assertTrue(np.isnan(validation[0, 0]))

    def test_ensemble_weighted_member_predictions_and_invalid_weights(self):
        x = np.ones((3, 4), dtype=np.float32)
        y = np.array([1., 2., 3.])
        members = [("low", DummyRegressor(strategy="constant", constant=1)),
                   ("high", DummyRegressor(strategy="constant", constant=3))]
        model = ValidationEnsemble(members, (.25, .75)).fit(x, y)
        self.assertEqual(model.predict_members(x).shape, (3, 2))
        np.testing.assert_allclose(model.predict(x), 2.5)
        for weights in [(-.1, 1.1), (.2, .2), (float("nan"), 1.)]:
            with self.assertRaises(ValueError):
                ValidationEnsemble(members, weights).fit(x, y)

    def test_selection_api_accepts_only_validation_inputs(self):
        signature = inspect.signature(TRAINING["select_validation"])
        self.assertEqual(list(signature.parameters), ["y_validation", "predictions"])
        candidates = TRAINING["candidates"]()
        self.assertEqual(len(candidates), 7)
        self.assertEqual(len({row["name"] for row in candidates}), 7)


class V3ReleaseEvidenceTests(unittest.TestCase):
    def test_frozen_v2_partitions_are_byte_identical(self):
        for target in TARGETS:
            directory = ROOT / "reports" / target
            previous = read_json(directory / "v2_baseline_metrics.json")
            frozen = (directory / "frozen_split_manifest.csv").read_bytes()
            self.assertEqual(frozen, (directory / "split_manifest.csv").read_bytes())
            self.assertEqual(hashlib.sha256(frozen).hexdigest(), previous["split_sha256"])

    def test_validation_selection_recomputes_without_test_predictions(self):
        for target in TARGETS:
            with self.subTest(target=target):
                directory = ROOT / "reports" / target
                frame = read_csv(directory / "validation_predictions.csv")
                manifest = read_csv(directory / "split_manifest.csv")
                validation = {row["canonical_smiles"]: float(row["value"]) for row in manifest if row["split"] == "validation"}
                self.assertEqual({row["canonical_smiles"] for row in frame}, set(validation))
                self.assertEqual(len(frame), len(validation))
                truth = np.array([float(row["y_true"]) for row in frame])
                np.testing.assert_allclose(truth, [validation[row["canonical_smiles"]] for row in frame], atol=1e-12)
                singles = read_csv(directory / "candidate_validation.csv")
                predictions = {row["model_name"]: np.array([float(item[row["model_name"]]) for item in frame]) for row in singles}
                rmse = {name: float(np.sqrt(np.mean((values - truth) ** 2))) for name, values in predictions.items()}
                for candidate in singles:
                    self.assertAlmostEqual(float(candidate["rmse"]), rmse[candidate["model_name"]], places=10)
                ranked = sorted(rmse, key=lambda name: (rmse[name], name))
                pair_options = []
                for left, right in itertools.combinations(ranked[:3], 2):
                    difference = predictions[left] - predictions[right]
                    denom = float(np.dot(difference, difference))
                    weight = float(np.clip(np.dot(truth - predictions[right], difference) / denom, .1, .9)) if denom else .5
                    score = float(np.sqrt(np.mean((weight * predictions[left] + (1 - weight) * predictions[right] - truth) ** 2)))
                    pair_options.append((score, [left, right], [weight, 1 - weight]))
                score, names, weights = min(pair_options, key=lambda item: (item[0], item[1]))
                if score > .99 * rmse[ranked[0]]:
                    names, weights = [ranked[0]], [1.]
                selection = read_json(directory / "ensemble_selection.json")
                self.assertEqual(selection["member_names"], names)
                np.testing.assert_allclose(selection["weights"], weights, atol=1e-12)
                combined = sum(weight * predictions[name] for weight, name in zip(weights, names))
                np.testing.assert_allclose(combined, [float(row["selected_prediction"]) for row in frame], atol=1e-12)

    def test_protocol_and_validation_hashes_match_release(self):
        for target in TARGETS:
            directory = ROOT / "reports" / target
            metadata = read_json(directory / "metrics.json")
            protocol = read_json(directory / "training_protocol.json")
            self.assertEqual(protocol["feature_version"], FEATURE_VERSION)
            self.assertEqual(len(protocol["candidates"]), 7)
            self.assertEqual(metadata["training_protocol_sha256"], hashlib.sha256((directory / "training_protocol.json").read_bytes()).hexdigest())
            self.assertEqual(metadata["selection_evidence_sha256"], hashlib.sha256((directory / "validation_predictions.csv").read_bytes()).hexdigest())
            self.assertEqual(metadata["candidate_count"], 7)
            self.assertEqual(metadata["feature_count"], 2267)

    def test_release_members_are_nonnegative_normalized_and_match_metadata(self):
        for target in TARGETS:
            with self.subTest(target=target):
                with (ROOT / "models" / f"{target}.pkl").open("rb") as handle:
                    bundle = pickle.load(handle)
                model = bundle["model"]
                self.assertIsInstance(model, ValidationEnsemble)
                self.assertGreaterEqual(len(model.member_names_), 1)
                self.assertLessEqual(len(model.member_names_), 2)
                self.assertTrue((model.member_weights_ >= 0).all())
                self.assertAlmostEqual(float(model.member_weights_.sum()), 1., places=12)
                self.assertEqual(model.member_names_, bundle["metadata"]["ensemble"]["member_names"])
                np.testing.assert_allclose(model.member_weights_, bundle["metadata"]["ensemble"]["weights"])
                x = np.stack([feature_vector(parse_smiles(value), "combined") for value in bundle["training_smiles"][:3]])
                member_values = model.predict_members(x)
                self.assertEqual(member_values.shape, (3, len(model.member_names_)))
                np.testing.assert_allclose(model.predict(x), member_values @ model.member_weights_, atol=1e-12)
                del bundle, model


if __name__ == "__main__":
    unittest.main()
