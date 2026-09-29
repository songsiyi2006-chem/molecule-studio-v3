"""V3 asynchronous HTTP and process-lifecycle tests; fixtures are explicit."""
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from fastapi.testclient import TestClient

from app.jobs import AnalysisManager, QueueFull
from app.main import create_app
from app.service import TARGET_IDS
from app.features import all_descriptors, feature_vector

try:
    from .test_service import FakeDatabase, write_bundle
except ImportError:
    from test_service import FakeDatabase, write_bundle


def await_condition(condition, timeout=8):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = condition()
        if value:
            return value
        time.sleep(0.02)
    raise AssertionError("Timed out waiting for test condition")


class TinyService:
    def predict(self, smiles, properties=None):
        return {"valid": True, "canonical_smiles": smiles, "predictions": {}, "timing": {"prediction_seconds": 0.001}, "evidence": {}}


class AnalysisManagerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "workspace.sqlite"

    def manager(self, **kwargs):
        manager = AnalysisManager(TinyService(), self.path, **kwargs)
        manager.start()
        self.addCleanup(manager.close)
        return manager

    def sleeping_child(self):
        original = subprocess.Popen
        processes = []
        def spawn(*args, **kwargs):
            process = original([sys.executable, "-c", "import time; time.sleep(60)"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            processes.append(process)
            return process
        return processes, spawn

    def test_deep_timeout_terminates_real_child_and_preserves_model_result(self):
        manager = self.manager(timeout_seconds=0.15)
        processes, spawn = self.sleeping_child()
        with patch("app.jobs.subprocess.Popen", side_effect=spawn):
            job = manager.submit("CCO", mode="deep")
            result = await_condition(lambda: (item if (item := manager.get(job["id"]))["status"] == "completed" else None))
        self.assertEqual(result["result"]["conformers"]["status"], "timeout")
        self.assertEqual(len(processes), 1)
        self.assertIsNotNone(processes[0].poll())
        self.assertIsNone(manager._active_process)
        with self.assertRaises(KeyError):
            manager.sdf_path(job["id"])

    def test_running_cancel_reaps_real_subprocess(self):
        manager = self.manager()
        processes, spawn = self.sleeping_child()
        with patch("app.jobs.subprocess.Popen", side_effect=spawn):
            job = manager.submit("CCO", mode="deep")
            await_condition(lambda: manager._active_process is not None)
            cancelled = manager.cancel(job["id"])
            self.assertEqual(cancelled["status"], "cancelled")
            self.assertIsNotNone(processes[0].poll())
            await_condition(lambda: manager._active_process is None)
        self.assertIsNone(manager.get(job["id"])["result"])

    def test_service_close_reaps_real_subprocess_and_marks_history(self):
        manager = self.manager()
        processes, spawn = self.sleeping_child()
        with patch("app.jobs.subprocess.Popen", side_effect=spawn):
            job = manager.submit("CCO", mode="deep")
            await_condition(lambda: manager._active_process is not None)
            manager.close()
        self.assertIsNotNone(processes[0].poll())
        self.assertFalse(manager._thread.is_alive())
        self.assertEqual(manager.get(job["id"])["status"], "cancelled")

    def test_finite_queue_and_queued_cancellation_do_not_run_cancelled_row(self):
        manager = self.manager(max_queued=1)
        entered, release = threading.Event(), threading.Event()
        calls = []
        def slow(smiles, properties=None):
            calls.append(smiles)
            entered.set()
            release.wait(timeout=5)
            return TinyService().predict(smiles)
        with patch.object(manager.service, "predict", side_effect=slow):
            first = manager.submit("CCO")
            self.assertTrue(entered.wait(2))
            queued = manager.submit("CCC")
            with self.assertRaises(QueueFull):
                manager.submit("CCCC")
            self.assertEqual(manager.cancel(queued["id"])["status"], "cancelled")
            release.set()
            await_condition(lambda: manager.get(first["id"])["status"] == "completed")
            await_condition(lambda: not manager._events)
        self.assertEqual(calls, ["CCO"])

    def test_history_persists_across_restart_and_retention_is_bounded(self):
        manager = self.manager(max_history=2)
        ids = []
        for _ in range(4):
            job = manager.submit("CCO")
            ids.append(job["id"])
            await_condition(lambda: manager.get(job["id"])["status"] == "completed")
            await_condition(lambda: not manager._events)
        self.assertEqual(manager.history()["total"], 2)
        manager.close()
        restarted = self.manager(max_history=2)
        self.assertEqual(restarted.get(ids[-1])["status"], "completed")
        self.assertNotIn("result", restarted.history()["items"][0])
        with self.assertRaises(KeyError):
            restarted.get(ids[0])

    def test_locked_old_artifact_does_not_kill_worker_and_cleanup_retries(self):
        manager = self.manager(max_history=1)
        first = manager.submit("CCO")
        await_condition(lambda: manager.get(first["id"])["status"] == "completed")
        await_condition(lambda: not manager._events)
        old_directory = manager.artifact_root / first["id"]
        old_directory.mkdir()
        (old_directory / "conformers.sdf").write_text("locked download fixture", encoding="utf-8")
        with patch("app.jobs.shutil.rmtree", side_effect=PermissionError("simulated Windows open-file lock")):
            with self.assertLogs("app.jobs", level="WARNING"):
                second = manager.submit("CCC")
                await_condition(lambda: manager.get(second["id"])["status"] == "completed")
                await_condition(lambda: not manager._events)
                manager._prune()
        self.assertEqual(manager.history()["total"], 1)
        self.assertTrue(old_directory.exists())
        self.assertTrue(manager._thread.is_alive())
        with self.assertRaises(KeyError):
            manager.get(first["id"])
        third = manager.submit("CCCC")
        await_condition(lambda: manager.get(third["id"])["status"] == "completed")
        await_condition(lambda: not manager._events)
        manager._prune()
        self.assertFalse(old_directory.exists())
        self.assertTrue(manager._thread.is_alive())


class V3APITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.directory = Path(cls.temp.name)
        for target in TARGET_IDS:
            write_bundle(cls.directory, target)
        cls.mock_database = patch("app.service.MoleculeDatabase", FakeDatabase)
        cls.mock_database.start()
        cls.application = create_app(models_dir=cls.directory, workspace_path=cls.directory / "workspace.sqlite")
        cls.context = TestClient(cls.application)
        cls.client = cls.context.__enter__()

    @classmethod
    def tearDownClass(cls):
        cls.context.__exit__(None, None, None)
        cls.mock_database.stop()
        cls.temp.cleanup()

    def wait_job(self, job_id, timeout=10):
        return await_condition(lambda: (job if (job := self.client.get(f"/analyses/{job_id}").json())["status"] in {"completed", "failed", "cancelled"} else None), timeout)

    def test_standard_analysis_is_asynchronous_and_skips_3d(self):
        response = self.client.post("/analyses", json={"smiles": "CCO", "mode": "standard"})
        self.assertEqual(response.status_code, 202)
        self.assertEqual(response.json()["status"], "queued")
        result = self.wait_job(response.json()["id"])
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["result"]["conformers"]["status"], "skipped")
        self.assertEqual(result["result"]["conformers"]["reason"], "standard_mode")
        self.assertFalse(result["result"]["evidence"]["quantum_calculation_executed"])
        self.assertEqual(self.client.get("/model-info").json()["version"], "3.0.0")

    def test_real_deep_analysis_exposes_3d_and_downloadable_sdf(self):
        response = self.client.post("/analyses", json={"smiles": "CCO", "mode": "deep", "properties": ["logS"]})
        self.assertEqual(response.status_code, 202)
        job = self.wait_job(response.json()["id"], timeout=22)
        self.assertEqual(job["status"], "completed")
        result = job["result"]["conformers"]
        self.assertEqual(result["status"], "completed")
        self.assertGreater(result["converged"], 0)
        sdf = self.client.get(result["sdf_url"])
        self.assertEqual(sdf.status_code, 200)
        self.assertIn("$$$$", sdf.text)
        self.assertIn("MMFF94_energy_kcal_mol", sdf.text)

    def test_health_remains_responsive_during_worker_prediction(self):
        entered, release = threading.Event(), threading.Event()
        original = self.application.state.prediction_service.predict
        def blocked(*args, **kwargs):
            entered.set()
            release.wait(timeout=3)
            return original(*args, **kwargs)
        with patch.object(self.application.state.prediction_service, "predict", side_effect=blocked):
            response = self.client.post("/analyses", json={"smiles": "CCO"})
            self.assertTrue(entered.wait(2))
            started = time.monotonic()
            self.assertEqual(self.client.get("/health").status_code, 200)
            self.assertLess(time.monotonic() - started, 0.5)
            release.set()
            self.wait_job(response.json()["id"])

    def test_analysis_validation_history_and_missing_routes(self):
        for payload in ({"smiles": "CCO", "mode": "DFT"}, {"smiles": "broken("}, {"smiles": "CCO", "properties": []}):
            self.assertEqual(self.client.post("/analyses", json=payload).status_code, 422)
        self.assertEqual(self.client.get("/analyses?limit=101").status_code, 422)
        self.assertEqual(self.client.get("/analyses/unknown").status_code, 404)
        self.assertEqual(self.client.delete("/analyses/unknown").status_code, 404)
        self.assertEqual(self.client.get("/analyses/unknown/conformers.sdf").status_code, 404)
        self.assertIn("items", self.client.get("/analyses?limit=5").json())

    def test_assignment_endpoint_reproduces_original_aspirin_and_withholds_large_ood(self):
        result = self.client.post("/assignment/predict", json={"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"})
        self.assertEqual(result.status_code, 200)
        self.assertAlmostEqual(result.json()["logS_prediction"], -1.9999771, places=6)
        self.assertEqual(result.json()["model"]["purpose"], "assignment_compatibility")
        ood = self.client.post("/assignment/predict", json={"smiles": "O" + "CCO" * 35}).json()
        self.assertIsNone(ood["logS_prediction"])
        self.assertEqual(ood["applicability_domain"]["status"], "out_of_domain")

    def test_ensemble_values_and_spread_are_consistent_and_suppressed_for_ood(self):
        estimator = self.application.state.prediction_service._models["logS"].bundle["model"]
        with patch.object(estimator, "predict_members", return_value=np.array([[-1.0, -3.0]]), create=True) as members, patch.object(estimator, "member_names_", ["A", "B"], create=True), patch.object(estimator, "member_weights_", np.array([0.25, 0.75]), create=True), patch.object(estimator, "predict", return_value=np.array([-2.5])):
            response = self.client.post("/predict", json={"smiles": "CCO", "properties": ["logS"]}).json()
            ensemble = response["predictions"]["logS"]["ensemble"]
            self.assertEqual(ensemble["weights"], [0.25, 0.75])
            self.assertEqual([row["value"] for row in ensemble["member_values"]], [-1.0, -3.0])
            self.assertAlmostEqual(ensemble["spread"], np.sqrt(0.75))
            members.reset_mock()
            ood = self.client.post("/predict", json={"smiles": "O" + "CCO" * 35, "properties": ["logS"]}).json()["predictions"]["logS"]
            members.assert_not_called()
            self.assertIsNone(ood["value"])
            self.assertEqual(ood["ensemble"]["member_values"], [])
            self.assertIsNone(ood["ensemble"]["spread"])

    def test_missing_extended_descriptor_is_null_and_not_fabricated_zero(self):
        def missing_descriptors(mol):
            values = all_descriptors(mol)
            values["BertzCT"] = float("nan")
            return values
        def missing_vector(mol, mode):
            values = feature_vector(mol, mode)
            values[-1] = float("nan")
            return values
        with patch("app.service.all_descriptors", side_effect=missing_descriptors), patch("app.service.feature_vector", side_effect=missing_vector):
            response = self.client.post("/predict", json={"smiles": "CCO", "properties": ["logS"]})
        self.assertEqual(response.status_code, 200)
        self.assertIsNone(response.json()["descriptor"]["BertzCT"])
        self.assertIn("BertzCT", response.json()["missing_descriptors"])
        self.assertIsNotNone(response.json()["logS_prediction"])


if __name__ == "__main__":
    unittest.main()
