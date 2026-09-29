"""One analysis worker, bounded queue, cancellable 3D children, SQLite history."""
from __future__ import annotations

import json
import logging
import queue
import shutil
import sqlite3
import subprocess
import sys
import threading
import time
import uuid
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

from .chemistry import canonical_smiles, parse_smiles
from .conformers import empty_conformers

LOGGER = logging.getLogger(__name__)
TERMINAL = {"completed", "failed", "cancelled"}


class QueueFull(RuntimeError):
    pass


def now():
    return datetime.now(timezone.utc).isoformat()


def terminate_child(process):
    if process is None:
        return
    if process.poll() is None:
        try:
            process.terminate()
        except OSError:
            pass
        try:
            process.wait(timeout=1)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=2)


class AnalysisManager:
    def __init__(self, service, path=None, *, timeout_seconds=20.0, max_queued=8, max_history=100):
        self.service = service
        self.path = Path(path) if path is not None else Path(__file__).resolve().parents[1] / "data" / "workspace.sqlite"
        self.artifact_root = self.path.parent / (self.path.stem + "_artifacts")
        self.timeout_seconds = min(20.0, max(0.05, timeout_seconds))
        self.max_history = max_history
        self._queue = queue.Queue(maxsize=max_queued)
        self._lock = threading.RLock()
        self._stop = threading.Event()
        self._events = {}
        self._active_process = None
        self._active_id = None
        self._thread = None

    def _connect(self):
        connection = sqlite3.connect(self.path, timeout=10)
        connection.row_factory = sqlite3.Row
        return connection

    def start(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.artifact_root.mkdir(parents=True, exist_ok=True)
        with closing(self._connect()) as connection, connection:
            connection.execute("PRAGMA journal_mode=WAL")
            connection.execute("CREATE TABLE IF NOT EXISTS analyses(id TEXT PRIMARY KEY,smiles TEXT NOT NULL,properties TEXT,mode TEXT NOT NULL,status TEXT NOT NULL,stage TEXT NOT NULL,created_at TEXT NOT NULL,started_at TEXT,finished_at TEXT,result TEXT,error TEXT)")
            connection.execute("UPDATE analyses SET status='failed',stage='interrupted',finished_at=?,error=? WHERE status IN ('queued','running')", (now(), json.dumps({"code": "INTERRUPTED", "message": "上次服务结束前任务未完成，请重新提交。"}, ensure_ascii=False)))
        self._prune()
        self._thread = threading.Thread(target=self._work, name="molecule-analysis-worker", daemon=True)
        self._thread.start()

    def _update(self, job_id, **values):
        with self._lock, closing(self._connect()) as connection, connection:
            connection.execute("UPDATE analyses SET " + ",".join(f"{key}=?" for key in values) + " WHERE id=?", [*values.values(), job_id])

    @staticmethod
    def _public(row, include_result=True):
        item = dict(row)
        result = item.pop("result", None)
        item["properties"] = json.loads(item["properties"]) if item.get("properties") else None
        item["error"] = json.loads(item["error"]) if item.get("error") else None
        end = datetime.fromisoformat(item["finished_at"]) if item.get("finished_at") else datetime.now(timezone.utc)
        start = datetime.fromisoformat(item.get("started_at") or item["created_at"])
        item["elapsed_seconds"] = max(0.0, (end - start).total_seconds())
        item["poll_url"] = f"/analyses/{item['id']}"
        if include_result:
            item["result"] = json.loads(result) if result else None
        return item

    def get(self, job_id):
        with self._lock, closing(self._connect()) as connection:
            row = connection.execute("SELECT * FROM analyses WHERE id=?", (job_id,)).fetchone()
        if row is None:
            raise KeyError(job_id)
        return self._public(row)

    def history(self, limit=20):
        with self._lock, closing(self._connect()) as connection:
            total = connection.execute("SELECT COUNT(*) FROM analyses").fetchone()[0]
            rows = connection.execute("SELECT * FROM analyses ORDER BY created_at DESC,id DESC LIMIT ?", (limit,)).fetchall()
        return {"items": [self._public(row, False) for row in rows], "total": total, "limit": limit}

    def submit(self, smiles, properties=None, mode="standard"):
        canonical = canonical_smiles(parse_smiles(smiles))
        with self._lock:
            if self._stop.is_set() or self._queue.full():
                raise QueueFull("计算队列已满，请等候当前任务完成后重试。")
            job_id = str(uuid.uuid4())
            with closing(self._connect()) as connection, connection:
                connection.execute("INSERT INTO analyses(id,smiles,properties,mode,status,stage,created_at) VALUES(?,?,?,?,?,?,?)", (job_id, canonical, json.dumps(properties) if properties else None, mode, "queued", "queued", now()))
            self._events[job_id] = threading.Event()
            self._queue.put_nowait(job_id)
            response = self.get(job_id)
        self._prune()
        return response

    def cancel(self, job_id):
        with self._lock:
            job = self.get(job_id)
            if job["status"] in TERMINAL:
                return job
            event = self._events.get(job_id)
            if event:
                event.set()
            self._update(job_id, status="cancelled", stage="cancelled", finished_at=now())
            process = self._active_process if self._active_id == job_id else None
        terminate_child(process)
        return self.get(job_id)

    def _conformers(self, job_id, smiles, event):
        directory = self.artifact_root / job_id
        directory.mkdir(parents=True, exist_ok=True)
        request = directory / "request.json"
        output = directory / "result.json"
        sdf = directory / "conformers.sdf"
        request.write_text(json.dumps({"smiles": smiles}), encoding="utf-8")
        started = time.monotonic()
        with self._lock:
            if event.is_set() or self._stop.is_set():
                return empty_conformers("skipped", "cancelled")
            process = subprocess.Popen([sys.executable, "-m", "app.conformers", str(request), str(output), str(sdf)], cwd=Path(__file__).resolve().parents[1], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            self._active_process, self._active_id = process, job_id
        try:
            while process.poll() is None:
                if event.wait(0.04) or self._stop.is_set():
                    terminate_child(process)
                    return empty_conformers("skipped", "cancelled")
                if time.monotonic() - started >= self.timeout_seconds:
                    terminate_child(process)
                    return empty_conformers("timeout", "wall_time_limit_20_seconds")
            if event.is_set():
                return empty_conformers("skipped", "cancelled")
            if process.returncode != 0 or not output.is_file():
                return empty_conformers("failed", "worker_failed")
            result = json.loads(output.read_text(encoding="utf-8"))
            if result["status"] == "completed" and sdf.is_file():
                result["sdf_url"] = f"/analyses/{job_id}/conformers.sdf"
            return result
        finally:
            terminate_child(process)
            with self._lock:
                if self._active_process is process:
                    self._active_process, self._active_id = None, None

    def _work(self):
        while not self._stop.is_set():
            try:
                job_id = self._queue.get(timeout=0.1)
            except queue.Empty:
                continue
            event = self._events[job_id]
            try:
                with self._lock:
                    if event.is_set() or self._stop.is_set():
                        continue
                    job = self.get(job_id)
                    self._update(job_id, status="running", stage="predicting", started_at=now())
                started = time.monotonic()
                result = self.service.predict(job["smiles"], job["properties"])
                if event.is_set() or self._stop.is_set():
                    continue
                conformer_started = time.monotonic()
                if job["mode"] == "deep":
                    self._update(job_id, stage="conformers")
                    result["conformers"] = self._conformers(job_id, job["smiles"], event)
                else:
                    result["conformers"] = empty_conformers()
                if event.is_set() or self._stop.is_set():
                    continue
                result.setdefault("timing", {})["conformer_seconds"] = time.monotonic() - conformer_started if job["mode"] == "deep" else 0.0
                result["timing"]["total_seconds"] = time.monotonic() - started
                result.setdefault("evidence", {})["conformers"] = "molecular_mechanics" if result["conformers"]["status"] == "completed" else "not_computed"
                with self._lock:
                    if not event.is_set():
                        self._update(job_id, status="completed", stage="completed", finished_at=now(), result=json.dumps(result, ensure_ascii=False, allow_nan=False))
            except Exception:
                LOGGER.exception("Analysis %s failed", job_id)
                with self._lock:
                    if not event.is_set():
                        self._update(job_id, status="failed", stage="failed", finished_at=now(), error=json.dumps({"code": "ANALYSIS_FAILED", "message": "分析未完成，请查看服务日志；没有生成伪造结果。"}, ensure_ascii=False))
            finally:
                self._queue.task_done()
                with self._lock:
                    self._events.pop(job_id, None)
                self._prune()

    def sdf_path(self, job_id):
        job = self.get(job_id)
        if job["status"] != "completed" or (job["result"] or {}).get("conformers", {}).get("status") != "completed":
            raise KeyError(job_id)
        path = self.artifact_root / job_id / "conformers.sdf"
        if not path.is_file():
            raise KeyError(job_id)
        return path

    def _prune(self):
        # History retention commits independently of best-effort disk cleanup.
        # Windows may hold an SDF open while it downloads. Keep that orphaned
        # artifact directory and retry it on the next prune/start, never kill the
        # worker or retain unbounded database rows because a file is locked.
        with self._lock:
            try:
                with closing(self._connect()) as connection, connection:
                    rows = connection.execute("SELECT id FROM analyses WHERE status IN ('completed','failed','cancelled') ORDER BY created_at DESC,id DESC LIMIT -1 OFFSET ?", (self.max_history,)).fetchall()
                    for row in rows:
                        if row["id"] not in self._events:
                            connection.execute("DELETE FROM analyses WHERE id=?", (row["id"],))
                    retained = {row[0] for row in connection.execute("SELECT id FROM analyses")}
                retained.update(self._events)
            except (sqlite3.Error, OSError):
                LOGGER.exception("History retention could not run; artifact cleanup deferred")
                return
            try:
                directories = list(self.artifact_root.iterdir())
            except OSError:
                LOGGER.exception("Analysis artifact directory could not be listed")
                return
            for candidate in directories:
                if candidate.name in retained:
                    continue
                try:
                    if candidate.name != str(uuid.UUID(candidate.name)):
                        continue
                    directory = candidate.resolve()
                    if directory.parent != self.artifact_root.resolve() or not directory.is_dir():
                        continue
                    shutil.rmtree(directory)
                except ValueError:
                    continue  # Non-job files are never cleanup targets.
                except OSError:
                    LOGGER.warning("Deferred cleanup of analysis artifacts %s", candidate.name, exc_info=True)

    def close(self):
        self._stop.set()
        with self._lock:
            for event in self._events.values():
                event.set()
            process = self._active_process
            with closing(self._connect()) as connection, connection:
                connection.execute("UPDATE analyses SET status='cancelled',stage='cancelled',finished_at=? WHERE status IN ('queued','running')", (now(),))
        terminate_child(process)
        if self._thread:
            self._thread.join(timeout=10)
