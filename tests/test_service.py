"""Service fault/domain tests use tiny fitted fixtures, never the release data."""

from __future__ import annotations

import copy
import math
import pickle
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import sklearn
from rdkit import DataStructs, rdBase
from sklearn.dummy import DummyRegressor

from app.chemistry import identity_smiles, parse_smiles
from app.features import FEATURE_VERSION, feature_vector, fingerprint, model_feature_names
from app.service import ModelUnavailable, PredictionService, TARGET_IDS, TARGETS


class FakeDatabase:
    """An explicit observation fixture distinguishable from model output."""

    def __init__(self, path=None):
        self.path = path

    def stats(self):
        return {"molecules": 1, "observations": 1, "datasets": [], "targets": {"logS": {"observations": 1, "molecules": 1, "unit": "log10(mol/L)"}}}

    def search(self, query="", target=None, limit=20, offset=0):
        matches = query in {"", "ethanol", "CCO"} and target in {None, "logS"}
        item = {"id": 1, "name": "ethanol", "canonical_smiles": "CCO", "formula": "C2H6O", "mw": 46.069, "targets": ["logS"]}
        return {"total": int(matches), "items": [item] if matches and offset == 0 else [], "limit": limit, "offset": offset}

    def lookup(self, canonical_smiles):
        if canonical_smiles != "CCO":
            return []
        return [{"dataset": "test-observation-fixture", "source_id": "row-1", "target": "logS", "value": -9.99, "unit": "log10(mol/L)", "quality": "fixture", "source_url": "https://example.invalid/fixture"}]


def make_bundle(target="logS", mode="descriptors", value=-2.125):
    smiles = ["CCO", "c1ccccc1"]
    molecules = [parse_smiles(s) for s in smiles]
    x = np.stack([feature_vector(mol, mode) for mol in molecules])
    estimator = DummyRegressor(strategy="constant", constant=value).fit(x, [value, value])
    return {
        "schema_version": 2,
        "feature_mode": mode,
        "feature_names": model_feature_names(mode),
        "model": estimator,
        "metadata": {"target": target, **TARGETS[target], "model_name": "Test fixture DummyRegressor", "version": "test-only-v2", "sklearn_version": sklearn.__version__, "rdkit_version": rdBase.rdkitVersion, "feature_version": FEATURE_VERSION},
        "training_smiles": smiles,
        "training_identities": [identity_smiles(mol) for mol in molecules],
        "training_names": ["Ethanol", "Benzene"],
        "training_fingerprints": [DataStructs.BitVectToBinaryText(fingerprint(mol)) for mol in molecules],
        "descriptor_bounds": {"MW": [0.0, 1000.0], "LogP": [-50.0, 50.0], "HBD": [0.0, 100.0], "HBA": [0.0, 100.0], "TPSA": [0.0, 2000.0], "RotatableBonds": [0.0, 200.0], "RingCount": [0.0, 100.0], "AromaticFraction": [0.0, 1.0], "HeavyAtoms": [0.0, 200.0]},
        "conformal_radius": 0.625,
    }


def write_bundle(directory, target="logS", **kwargs):
    bundle = make_bundle(target, **kwargs)
    path = Path(directory) / f"{target}.pkl"
    with path.open("wb") as handle:
        pickle.dump(bundle, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return bundle


class PredictionServiceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.model_dir = Path(self.temp.name)
        self.database_patch = patch("app.service.MoleculeDatabase", FakeDatabase)
        self.database_patch.start()
        self.addCleanup(self.database_patch.stop)

    def service(self, targets=TARGET_IDS, mode="descriptors"):
        for target in targets:
            write_bundle(self.model_dir, target, mode=mode)
        service = PredictionService(self.model_dir)
        service.load()
        return service

    def test_models_load_only_once(self):
        for target in TARGET_IDS:
            write_bundle(self.model_dir, target)
        with patch("app.service.pickle.load", wraps=pickle.load) as load:
            service = PredictionService(self.model_dir)
            self.assertTrue(service.load())
            self.assertTrue(service.load())
            self.assertEqual(load.call_count, 3)
        self.assertEqual(service.models_loaded, list(TARGET_IDS))

    def test_invalid_artifact_versions_schema_and_calibration_fail_closed(self):
        mutations = [
            lambda b: b.update(schema_version=1),
            lambda b: b["metadata"].update(sklearn_version="0.invalid"),
            lambda b: b["metadata"].update(rdkit_version="0.invalid"),
            lambda b: b["metadata"].update(feature_version="wrong-fingerprint-settings"),
            lambda b: b.update(feature_names=list(reversed(b["feature_names"]))),
            lambda b: b.update(conformal_radius=float("nan")),
            lambda b: b.update(training_fingerprints=[b"bad", b"bad"]),
        ]
        for mutate in mutations:
            bundle = make_bundle()
            mutate(bundle)
            with self.subTest(mutation=mutate):
                with self.assertRaises(ValueError):
                    PredictionService._validate_bundle("logS", bundle)

    def test_partial_missing_models_are_unavailable_without_fabricated_values(self):
        service = self.service(targets=["logS"])
        result = service.predict("CCO")
        self.assertEqual(result["predictions"]["logS"]["status"], "ok")
        for target in ("logD74", "hydration_free_energy"):
            self.assertEqual(result["predictions"][target]["status"], "unavailable")
            self.assertIsNone(result["predictions"][target]["value"])
            self.assertIsNone(result["predictions"][target]["interval90"])

    def test_no_models_has_no_fallback_prediction(self):
        service = PredictionService(self.model_dir)
        self.assertFalse(service.load())
        with self.assertRaises(ModelUnavailable):
            service.predict("CCO")

    def test_both_feature_modes_match_the_estimator_and_keep_precision(self):
        for mode in ("descriptors", "combined"):
            with self.subTest(mode=mode):
                service = self.service(mode=mode)
                result = service.predict("CCO", ["logS"])
                prediction = result["predictions"]["logS"]
                bundle = service._models["logS"].bundle
                expected = float(bundle["model"].predict(feature_vector(parse_smiles("CCO"), mode).reshape(1, -1))[0])
                self.assertEqual(prediction["value"], expected)
                self.assertEqual(prediction["interval90"], [expected - 0.625, expected + 0.625])
                self.assertAlmostEqual(math.log10(result["solubility_mol_L"]), result["logS_prediction"], places=12)

    def test_large_benign_peg_withholds_values_and_does_not_invoke_estimator(self):
        service = self.service()
        estimator = service._models["logS"].bundle["model"]
        with patch.object(estimator, "predict", side_effect=AssertionError("OOD inference must not run")) as predict:
            result = service.predict("O" + "CCO" * 35, ["logS"])
            predict.assert_not_called()
        target = result["predictions"]["logS"]
        self.assertEqual(target["status"], "out_of_domain")
        self.assertIsNone(target["value"])
        self.assertIsNone(target["interval90"])
        self.assertIsNone(result["logS_prediction"])
        self.assertIsNone(result["solubility_mol_L"])
        self.assertTrue(target["domain_violations"])

    def test_three_range_violations_withhold_but_one_only_cautions(self):
        service = self.service()
        bounds = service._models["logS"].bundle["descriptor_bounds"]
        bounds["TPSA"] = [0.0, 1.0]
        one = service.predict("CCO", ["logS"])["predictions"]["logS"]
        self.assertEqual(one["status"], "caution")
        self.assertIsNotNone(one["value"])
        bounds["HBD"] = [0.0, 0.0]
        bounds["HBA"] = [0.0, 0.0]
        three = service.predict("CCO", ["logS"])["predictions"]["logS"]
        self.assertEqual(three["status"], "out_of_domain")
        self.assertIsNone(three["value"])

    def test_similarity_is_training_only_and_membership_uses_connectivity(self):
        service = self.service()
        result = service.predict("OCC", ["logS"])["predictions"]["logS"]
        self.assertTrue(result["training_match"])
        self.assertEqual(result["nearest_similarity"], 1.0)
        self.assertEqual(result["nearest_neighbors"][0]["canonical_smiles"], "CCO")
        candidate = "C[C@H](O)Cl"
        bundle = copy.deepcopy(make_bundle())
        bundle["training_smiles"][0] = candidate
        bundle["training_identities"][0] = identity_smiles(parse_smiles(candidate))
        bundle["training_fingerprints"][0] = DataStructs.BitVectToBinaryText(fingerprint(parse_smiles(candidate)))
        bundle["training_names"][0] = "Stereochemistry test fixture"
        service._models["logS"] = PredictionService._validate_bundle("logS", bundle)
        reverse_stereo = service.predict("C[C@@H](O)Cl", ["logS"])["predictions"]["logS"]
        self.assertTrue(reverse_stereo["training_match"])

    def test_low_similarity_returns_a_caution_not_false_certainty(self):
        service = self.service()
        target = service.predict("Cn1c(=O)c2c(ncn2C)n(C)c1=O", ["logS"])["predictions"]["logS"]
        self.assertLess(target["nearest_similarity"], 0.3)
        self.assertEqual(target["status"], "caution")
        self.assertIsNotNone(target["value"])

    def test_observation_is_never_substituted_for_prediction(self):
        result = self.service().predict("CCO")
        self.assertEqual(result["logS_prediction"], -2.125)
        self.assertEqual(result["measurements"][0]["value"], -9.99)
        self.assertEqual(result["measurements"][0]["source_id"], "row-1")
        self.assertIn("source_url", result["measurements"][0])

    def test_one_nonfinite_model_output_does_not_destroy_other_properties(self):
        service = self.service()
        with patch.object(service._models["logS"].bundle["model"], "predict", return_value=np.array([np.nan])):
            with self.assertLogs("app.service", level="ERROR"):
                result = service.predict("CCO")
        self.assertEqual(result["predictions"]["logS"]["status"], "unavailable")
        self.assertIsNone(result["logS_prediction"])
        self.assertIsNotNone(result["predictions"]["logD74"]["value"])

    def test_charged_input_has_an_explicit_ionization_warning(self):
        result = self.service().predict("C[NH3+]", ["logD74"])
        self.assertTrue(any("离子化" in warning for warning in result["warnings"]))


if __name__ == "__main__":
    unittest.main()
