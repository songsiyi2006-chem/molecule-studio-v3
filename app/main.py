"""FastAPI serving, catalog search, and bounded batch prediction."""

from __future__ import annotations

import logging
import sqlite3
from contextlib import closing
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Annotated, Any, Literal

from fastapi import FastAPI, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, ConfigDict, Field, StrictStr, field_validator

from .chemistry import MoleculeError
from .service import AssignmentPredictionService, DatabaseUnavailable, ModelUnavailable, PredictionService
from .jobs import AnalysisManager, QueueFull

LOGGER = logging.getLogger(__name__)
STATIC_DIR = Path(__file__).resolve().parent / "static"
TargetId = Literal["logS", "logD74", "hydration_free_energy"]
MAX_REQUEST_BYTES = 256 * 1024


class PredictionBodyLimit:
    """Bound input buffering before JSON decoding or chemistry calculation."""

    def __init__(self, app: Any, max_bytes: int = MAX_REQUEST_BYTES) -> None:
        self.app, self.max_bytes = app, max_bytes

    async def __call__(self, scope: dict[str, Any], receive: Any, send: Any) -> None:
        if scope["type"] != "http" or scope.get("path") not in {"/predict", "/predict/batch", "/analyses", "/assignment/predict"} or scope.get("method") != "POST":
            await self.app(scope, receive, send)
            return
        headers = dict(scope.get("headers", []))
        try:
            declared_size = int(headers.get(b"content-length", b"0"))
        except ValueError:
            declared_size = 0
        if declared_size > self.max_bytes:
            await error_response(413, "REQUEST_TOO_LARGE", "预测请求体不能超过 256 KiB，请减少输入长度或拆分批次。")(scope, receive, send)
            return
        chunks: list[bytes] = []
        size = 0
        while True:
            message = await receive()
            if message["type"] == "http.disconnect":
                return
            chunk = message.get("body", b"")
            size += len(chunk)
            if size > self.max_bytes:
                await error_response(413, "REQUEST_TOO_LARGE", "预测请求体不能超过 256 KiB，请减少输入长度或拆分批次。")(scope, receive, send)
                return
            chunks.append(chunk)
            if not message.get("more_body", False):
                break
        body = b"".join(chunks)
        replayed = False

        async def replay_receive() -> dict[str, Any]:
            nonlocal replayed
            if not replayed:
                replayed = True
                return {"type": "http.request", "body": body, "more_body": False}
            return await receive()

        await self.app(scope, replay_receive, send)


class PropertySelection(BaseModel):
    model_config = ConfigDict(extra="forbid")
    properties: list[TargetId] | None = Field(default=None, min_length=1, max_length=3)

    @field_validator("properties")
    @classmethod
    def unique_properties(cls, values: list[str] | None) -> list[str] | None:
        if values is not None and len(values) != len(set(values)):
            raise ValueError("Properties must be unique")
        return values


class PredictRequest(PropertySelection):
    smiles: StrictStr = Field(..., min_length=1, max_length=2000, examples=["CCO"])


class AnalysisRequest(PredictRequest):
    mode: Literal["standard", "deep"] = "standard"


class AssignmentRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    smiles: StrictStr = Field(..., min_length=1, max_length=2000)


class BatchRequest(PropertySelection):
    # Individual empty/oversized strings are rejected per row by the chemistry
    # parser so one bad molecule does not discard valid batch members.
    smiles: list[StrictStr] = Field(..., min_length=1, max_length=50)


def error_response(status_code: int, code: str, message: str, *, valid: bool = False) -> JSONResponse:
    return JSONResponse(status_code=status_code, content={"valid": valid, "error": {"code": code, "message": message}})


def create_app(
    models_dir: Path | str | None = None,
    database_path: Path | str | None = None,
    workspace_path: Path | str | None = None,
    assignment_model_path: Path | str | None = None,
    analysis_timeout_seconds: float = 20.0,
) -> FastAPI:
    service = PredictionService(models_dir, database_path)
    assignment = AssignmentPredictionService(assignment_model_path)
    # An explicitly supplied model directory is commonly an isolated test or
    # alternate deployment. Keep its writable history separate by default.
    selected_workspace = workspace_path if workspace_path is not None else Path(models_dir) / "workspace.sqlite" if models_dir is not None else None
    manager = AnalysisManager(service, selected_workspace, timeout_seconds=analysis_timeout_seconds)

    @asynccontextmanager
    async def lifespan(application: FastAPI):
        service.load()
        assignment.load()
        quantum_path = Path(__file__).resolve().parents[1] / "data" / "quantum.sqlite"
        try:
            with closing(sqlite3.connect(quantum_path.resolve().as_uri() + "?mode=ro", uri=True, timeout=2)) as connection:
                application.state.quantum_ready = connection.execute("SELECT 1 FROM quantum LIMIT 1").fetchone() is not None
        except (sqlite3.Error, OSError):
            application.state.quantum_ready = False
            LOGGER.warning("Quantum reference catalog is unavailable at startup")
        manager.start()
        try:
            yield
        finally:
            manager.close()

    application = FastAPI(
        title="Molecule Studio V3",
        version="3.0.0",
        description="Local molecular models, cancellable asynchronous ETKDGv3/MMFF94 analysis, and separate experimental and quantum-computed reference catalogs.",
        lifespan=lifespan,
    )
    application.state.prediction_service = service
    application.state.analysis_manager = manager
    application.state.assignment_service = assignment
    application.state.quantum_ready = False
    application.add_middleware(PredictionBodyLimit)
    from .catalog_api import router as catalog_router
    application.include_router(catalog_router)

    @application.exception_handler(RequestValidationError)
    async def request_validation_handler(request: Request, exc: RequestValidationError) -> JSONResponse:
        return error_response(422, "INVALID_REQUEST", "请求格式不正确：单分子 SMILES 限 1–2000 字符；批量为 1–50 个字符串；性质必须是不重复的受支持 ID；搜索分页或筛选参数必须在允许范围内。")

    @application.exception_handler(MoleculeError)
    async def molecule_error_handler(request: Request, exc: MoleculeError) -> JSONResponse:
        return error_response(422, exc.code, exc.message, valid=exc.valid)

    @application.exception_handler(ModelUnavailable)
    async def model_error_handler(request: Request, exc: ModelUnavailable) -> JSONResponse:
        return error_response(503, "MODEL_UNAVAILABLE", str(exc), valid=True)

    @application.exception_handler(DatabaseUnavailable)
    async def database_error_handler(request: Request, exc: DatabaseUnavailable) -> JSONResponse:
        return error_response(503, "DATABASE_UNAVAILABLE", str(exc))

    @application.exception_handler(QueueFull)
    async def queue_error_handler(request: Request, exc: QueueFull) -> JSONResponse:
        return error_response(429, "ANALYSIS_QUEUE_FULL", str(exc))

    @application.exception_handler(Exception)
    async def unexpected_error_handler(request: Request, exc: Exception) -> JSONResponse:
        LOGGER.exception("Unexpected API error", exc_info=exc)
        return error_response(500, "INTERNAL_ERROR", "服务暂时无法处理此请求，请查看服务端日志。")

    @application.get("/health")
    def health() -> dict[str, Any]:
        return {**service.health(), "assignment_model_loaded": assignment.bundle is not None,
                "quantum_ready": application.state.quantum_ready}

    @application.get("/model-info")
    def model_info() -> dict[str, Any]:
        return service.model_info()

    @application.get("/database/stats")
    def database_stats() -> dict[str, Any]:
        return service.database_stats()

    @application.get("/database/search")
    def database_search(
        q: Annotated[str, Query(max_length=2000)] = "",
        target: TargetId | None = None,
        limit: Annotated[int, Query(ge=1, le=100)] = 20,
        offset: Annotated[int, Query(ge=0, le=1_000_000)] = 0,
    ) -> dict[str, Any]:
        return service.database_search(q, target, limit, offset)

    @application.post("/predict")
    def predict(payload: PredictRequest) -> dict[str, Any]:
        return service.predict(payload.smiles, payload.properties)

    @application.post("/assignment/predict")
    def assignment_predict(payload: AssignmentRequest) -> dict[str, Any]:
        return assignment.predict(payload.smiles)

    @application.post("/analyses", status_code=202)
    def submit_analysis(payload: AnalysisRequest) -> dict[str, Any]:
        if not service.model_loaded:
            raise ModelUnavailable("没有可用主模型，无法提交分析任务。")
        return manager.submit(payload.smiles, payload.properties, payload.mode)

    @application.get("/analyses")
    def analysis_history(limit: Annotated[int, Query(ge=1, le=100)] = 20) -> dict[str, Any]:
        return manager.history(limit)

    @application.get("/analyses/{job_id}")
    def analysis_status(job_id: str):
        try:
            return manager.get(job_id)
        except KeyError:
            return error_response(404, "ANALYSIS_NOT_FOUND", "未找到该分析任务。")

    @application.delete("/analyses/{job_id}")
    def cancel_analysis(job_id: str):
        try:
            return manager.cancel(job_id)
        except KeyError:
            return error_response(404, "ANALYSIS_NOT_FOUND", "未找到该分析任务。")

    @application.get("/analyses/{job_id}/conformers.sdf")
    def download_conformers(job_id: str):
        try:
            path = manager.sdf_path(job_id)
        except KeyError:
            return error_response(404, "CONFORMERS_NOT_AVAILABLE", "此任务没有可下载的已完成构象文件。")
        return FileResponse(path, media_type="chemical/x-mdl-sdfile", filename=f"{job_id}-MMFF94.sdf")

    @application.post("/predict/batch")
    def batch_predict(payload: BatchRequest) -> dict[str, Any]:
        if not service.model_loaded:
            raise ModelUnavailable("没有可用模型，请先运行训练脚本并重新启动服务。")
        items = []
        for index, smiles in enumerate(payload.smiles):
            item: dict[str, Any] = {"index": index, "input": smiles}
            try:
                item["result"] = service.predict(smiles, payload.properties)
            except MoleculeError as exc:
                item["error"] = {"code": exc.code, "message": exc.message, "valid": exc.valid, "status_code": 422}
            except ModelUnavailable as exc:
                item["error"] = {"code": "MODEL_UNAVAILABLE", "message": str(exc), "valid": True, "status_code": 503}
            except Exception:
                LOGGER.exception("Batch row %d failed", index)
                item["error"] = {"code": "INTERNAL_ERROR", "message": "此行暂时无法处理，请查看服务端日志。", "valid": False, "status_code": 500}
            items.append(item)
        return {"items": items, "total": len(items)}

    @application.get("/", include_in_schema=False)
    def index() -> FileResponse:
        return FileResponse(STATIC_DIR / "index.html", media_type="text/html")

    application.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
    return application


app = create_app()
