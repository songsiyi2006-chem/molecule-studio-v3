"""HTTP contract, independent batch failures, request bounds, and provenance."""

from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fastapi.testclient import TestClient

from app.main import MAX_REQUEST_BYTES, create_app
from app.service import TARGET_IDS

try:
    from .test_service import FakeDatabase, write_bundle
except ImportError:
    from test_service import FakeDatabase, write_bundle


class APITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.model_dir = Path(cls.temp.name)
        for target in TARGET_IDS:
            write_bundle(cls.model_dir, target)
        cls.database_patch = patch("app.service.MoleculeDatabase", FakeDatabase)
        cls.database_patch.start()
        cls.app = create_app(cls.model_dir)
        cls.client_context = TestClient(cls.app)
        cls.client = cls.client_context.__enter__()

    @classmethod
    def tearDownClass(cls):
        cls.client_context.__exit__(None, None, None)
        cls.database_patch.stop()
        cls.temp.cleanup()

    def test_health_and_catalog_report_all_properties(self):
        health = self.client.get("/health").json()
        self.assertEqual(health["status"], "ok")
        self.assertTrue(health["database_ready"])
        self.assertEqual(health["models_loaded"], list(TARGET_IDS))
        catalog = self.client.get("/model-info").json()
        self.assertEqual(catalog["version"], "3.0.0")
        self.assertEqual({model["id"] for model in catalog["models"]}, set(TARGET_IDS))

    def test_default_prediction_and_legacy_fields_are_consistent(self):
        response = self.client.post("/predict", json={"smiles": "CCO"})
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result["valid"])
        self.assertEqual(set(result["predictions"]), set(TARGET_IDS))
        self.assertEqual(result["logS_prediction"], result["predictions"]["logS"]["value"])
        self.assertAlmostEqual(math.log10(result["solubility_mol_L"]), result["logS_prediction"], places=12)
        self.assertAlmostEqual(result["descriptor"]["MW"], 46.069, places=3)
        self.assertEqual(result["descriptor"]["HBD"], 1)
        self.assertEqual(result["formula"], "C2H6O")
        self.assertIn("<svg", result["structure_svg"])
        self.assertNotIn("<script", result["structure_svg"].lower())
        self.assertNotEqual(result["measurements"][0]["value"], result["logS_prediction"])

    def test_selected_properties_do_not_fabricate_legacy_logS(self):
        result = self.client.post("/predict", json={"smiles": "CCO", "properties": ["logD74"]}).json()
        self.assertEqual(list(result["predictions"]), ["logD74"])
        self.assertIsNone(result["logS_prediction"])
        self.assertIsNone(result["solubility_mol_L"])
        self.assertEqual(result["measurements"][0]["target"], "logS")

    def test_invalid_requests_are_structured_422(self):
        payloads = [{}, {"smiles": None}, {"smiles": 42}, {"smiles": True}, {"smiles": []}, {"smiles": ""}, {"smiles": "C" * 2001}, {"smiles": "CCO", "properties": []}, {"smiles": "CCO", "properties": ["unknown"]}, {"smiles": "CCO", "properties": ["logS", "logS"]}, {"smiles": "CCO", "model_path": "other.pkl"}]
        for payload in payloads:
            with self.subTest(payload=str(payload)[:80]):
                response = self.client.post("/predict", json=payload)
                self.assertEqual(response.status_code, 422)
                self.assertFalse(response.json()["valid"])
                self.assertIn("message", response.json()["error"])
        malformed = self.client.post("/predict", content=b'{"smiles":', headers={"Content-Type": "application/json"})
        self.assertEqual(malformed.status_code, 422)

    def test_invalid_and_unsupported_smiles_remain_distinct(self):
        for smiles, valid in [("broken(", False), ("<script>alert(1)</script>", False), ("CCO.[Na+]", True), ("[Cu]", True), ("*C", True), ("[CH3]", True), ("C" * 201, True)]:
            with self.subTest(smiles=smiles[:40]):
                response = self.client.post("/predict", json={"smiles": smiles})
                self.assertEqual(response.status_code, 422)
                self.assertEqual(response.json()["valid"], valid)
                self.assertNotIn("predictions", response.json())

    def test_large_peg_never_exposes_numeric_predictions(self):
        result = self.client.post("/predict", json={"smiles": "O" + "CCO" * 35}).json()
        self.assertTrue(result["valid"])
        for prediction in result["predictions"].values():
            self.assertEqual(prediction["status"], "out_of_domain")
            self.assertIsNone(prediction["value"])
            self.assertIsNone(prediction["interval90"])
        self.assertIsNone(result["logS_prediction"])

    def test_batch_mixed_rows_continue_and_keep_input_order(self):
        inputs = ["CCO", "broken(", "CCO.[Na+]", "", "C" * 2001, "c1ccccc1"]
        response = self.client.post("/predict/batch", json={"smiles": inputs, "properties": ["logS"]})
        self.assertEqual(response.status_code, 200)
        items = response.json()["items"]
        self.assertEqual([item["input"] for item in items], inputs)
        self.assertEqual([item["index"] for item in items], list(range(len(inputs))))
        self.assertIn("result", items[0])
        self.assertIn("result", items[-1])
        for item in items[1:-1]:
            self.assertIn("error", item)
            self.assertNotIn("result", item)
        self.assertFalse(items[1]["error"]["valid"])
        self.assertTrue(items[2]["error"]["valid"])

    def test_batch_bounds_prevent_unbounded_work(self):
        for smiles in ([], ["CCO"] * 51, ["CCO", 123]):
            with self.subTest(length=len(smiles)):
                with patch.object(self.app.state.prediction_service, "predict") as predict:
                    response = self.client.post("/predict/batch", json={"smiles": smiles})
                    self.assertEqual(response.status_code, 422)
                    predict.assert_not_called()
        accepted = self.client.post("/predict/batch", json={"smiles": ["CCO"] * 50, "properties": ["logS"]})
        self.assertEqual(accepted.status_code, 200)
        self.assertEqual(len(accepted.json()["items"]), 50)

    def test_request_body_limit_rejects_before_prediction(self):
        with patch.object(self.app.state.prediction_service, "predict") as predict:
            response = self.client.post("/predict", content=b"x" * (MAX_REQUEST_BYTES + 1), headers={"Content-Type": "application/json"})
            self.assertEqual(response.status_code, 413)
            self.assertEqual(response.json()["error"]["code"], "REQUEST_TOO_LARGE")
            predict.assert_not_called()

    def test_database_queries_are_validated_and_observations_keep_provenance(self):
        stats = self.client.get("/database/stats")
        self.assertEqual(stats.status_code, 200)
        self.assertEqual(stats.json()["observations"], 1)
        search = self.client.get("/database/search", params={"q": "ethanol", "target": "logS", "limit": 1})
        self.assertEqual(search.status_code, 200)
        self.assertEqual(search.json()["items"][0]["targets"], ["logS"])
        for params in ({"target": "unknown"}, {"limit": 101}, {"limit": 0}, {"offset": -1}, {"offset": 1_000_001}, {"q": "x" * 2001}):
            self.assertEqual(self.client.get("/database/search", params=params).status_code, 422)

    def test_home_and_static_assets_are_served(self):
        self.assertEqual(self.client.get("/").status_code, 200)
        self.assertEqual(self.client.get("/static/app.js").status_code, 200)
        self.assertEqual(self.client.get("/static/styles.css").status_code, 200)


class UnavailableAPITests(unittest.TestCase):
    def test_missing_artifacts_returns_503_but_health_still_works(self):
        with tempfile.TemporaryDirectory() as directory:
            with TestClient(create_app(directory, Path(directory) / "missing.sqlite")) as client:
                self.assertFalse(client.get("/health").json()["model_loaded"])
                self.assertFalse(client.get("/health").json()["database_ready"])
                for path in ("/predict", "/predict/batch"):
                    payload = {"smiles": "CCO" if path == "/predict" else ["CCO"]}
                    response = client.post(path, json=payload)
                    self.assertEqual(response.status_code, 503)
                    self.assertNotIn("Traceback", response.text)
                self.assertEqual(client.get("/model-info").status_code, 503)
                self.assertEqual(client.get("/database/stats").status_code, 503)
                self.assertEqual(client.get("/database/search").status_code, 503)

    def test_corrupted_artifact_is_unavailable_instead_of_a_dummy_prediction(self):
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory) / "logS.pkl").write_bytes(b"not a valid pickle artifact")
            with self.assertLogs("app.service", level="WARNING"):
                with TestClient(create_app(directory, Path(directory) / "missing.sqlite")) as client:
                    response = client.post("/predict", json={"smiles": "CCO"})
                    self.assertEqual(response.status_code, 503)
                    self.assertNotIn(directory, response.text)


if __name__ == "__main__":
    unittest.main()
