"""Trusted local model catalog with explicit domain and evidence boundaries.

Pickles are executable data. Only the offline project trainer writes these
artifacts; no HTTP endpoint accepts a model path or an uploaded artifact.
"""

from __future__ import annotations

import logging
import csv
import math
import pickle
import time
from dataclasses import dataclass
from numbers import Real
from pathlib import Path
from typing import Any

import numpy as np
import sklearn
from rdkit import DataStructs, rdBase
from sklearn.utils.validation import check_is_fitted

from .chemistry import (
    FEATURE_NAMES as DOMAIN_FEATURES,
    MoleculeError,
    calculate_descriptors,
    canonical_smiles,
    identity_smiles,
    known_name,
    molecule_formula,
    molecule_svg,
    parse_smiles,
)
from .database import MoleculeDatabase
from .conformers import empty_conformers
from .features import FEATURE_VERSION, all_descriptors, feature_vector, fingerprint, model_feature_names

LOGGER = logging.getLogger(__name__)
PROJECT_ROOT = Path(__file__).resolve().parents[1]
TARGETS = {
    "logS": {"label": "水溶解度 logS", "unit": "log10(mol/L)"},
    "logD74": {"label": "分配系数 logD（pH 7.4）", "unit": "log10(D), pH 7.4"},
    "hydration_free_energy": {"label": "水合自由能", "unit": "kcal/mol"},
}
TARGET_IDS = tuple(TARGETS)
INTERVAL_WARNING = "区间为独立校准集残差得到的名义 90% 经验区间；不保证新分子或分布外分子的覆盖率。"


class ModelUnavailable(RuntimeError):
    """No usable local model is available."""


class DatabaseUnavailable(RuntimeError):
    """The local reference-data catalog cannot be queried."""


@dataclass(frozen=True)
class LoadedModel:
    bundle: dict[str, Any]
    fingerprints: tuple[Any, ...]
    training_identities: frozenset[str]


def _finite_real(value: Any) -> bool:
    return isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(float(value))


class PredictionService:
    """Load each fixed artifact once and keep scientific observations separate."""

    def __init__(
        self,
        models_dir: Path | str | None = None,
        database_path: Path | str | None = None,
    ) -> None:
        self.models_dir = Path(models_dir) if models_dir is not None else PROJECT_ROOT / "models"
        self.database_path = Path(database_path) if database_path is not None else None
        self._models: dict[str, LoadedModel] = {}
        self._database: MoleculeDatabase | None = None
        self._load_attempted = False

    @property
    def model_loaded(self) -> bool:
        return bool(self._models)

    @property
    def models_loaded(self) -> list[str]:
        return [target for target in TARGET_IDS if target in self._models]

    @property
    def database_ready(self) -> bool:
        return self._database is not None

    def load(self) -> bool:
        if self._load_attempted:
            return self.model_loaded
        self._load_attempted = True
        for target in TARGET_IDS:
            path = self.models_dir / f"{target}.pkl"
            try:
                with path.open("rb") as handle:
                    bundle = pickle.load(handle)
                self._models[target] = self._validate_bundle(target, bundle)
                LOGGER.info("Loaded %s from %s", target, path)
            except FileNotFoundError:
                LOGGER.warning("Model artifact is missing: %s", path)
            except Exception:
                LOGGER.exception("Could not load model artifact %s", path)
        try:
            database = MoleculeDatabase(self.database_path)
            database.stats()
            self._database = database
        except FileNotFoundError:
            LOGGER.warning("Local molecular observation database is missing")
        except Exception:
            LOGGER.exception("Could not initialize the local molecular observation database")
        return self.model_loaded

    @staticmethod
    def _validate_bundle(target: str, bundle: Any) -> LoadedModel:
        if not isinstance(bundle, dict) or type(bundle.get("schema_version")) is not int or bundle["schema_version"] != 2:
            raise ValueError("Unsupported model artifact schema")
        mode = bundle.get("feature_mode")
        if mode not in {"descriptors", "combined"}:
            raise ValueError("Unknown model feature mode")
        names = model_feature_names(mode)
        if bundle.get("feature_names") != names:
            raise ValueError("Artifact feature names do not match the shared feature implementation")
        metadata = bundle.get("metadata")
        if not isinstance(metadata, dict) or metadata.get("target") != target:
            raise ValueError("Artifact metadata target does not match its file")
        for field in ("model_name", "version", "label", "unit", "sklearn_version", "rdkit_version"):
            if not isinstance(metadata.get(field), str) or not metadata[field].strip():
                raise ValueError(f"Missing model metadata: {field}")
        if metadata["sklearn_version"] != sklearn.__version__ or metadata["rdkit_version"] != rdBase.rdkitVersion:
            raise ValueError("Model requires the recorded scikit-learn and RDKit versions")
        if metadata.get("feature_version") != FEATURE_VERSION:
            raise ValueError("Model requires the recorded shared feature implementation version")
        estimator = bundle.get("model")
        if not callable(getattr(estimator, "predict", None)):
            raise ValueError("Model does not expose predict")
        check_is_fitted(estimator)
        if getattr(estimator, "n_features_in_", len(names)) != len(names):
            raise ValueError("Estimator feature count does not match the artifact")
        fitted_names = getattr(estimator, "feature_names_in_", None)
        if fitted_names is not None and list(fitted_names) != names:
            raise ValueError("Estimator feature order does not match the artifact")
        if not _finite_real(bundle.get("conformal_radius")) or bundle["conformal_radius"] < 0:
            raise ValueError("Invalid calibration radius")
        bounds = bundle.get("descriptor_bounds")
        if not isinstance(bounds, dict):
            raise ValueError("Missing descriptor domain bounds")
        for name in DOMAIN_FEATURES:
            interval = bounds.get(name)
            if (
                not isinstance(interval, (list, tuple))
                or len(interval) != 2
                or not all(_finite_real(value) for value in interval)
                or interval[0] > interval[1]
            ):
                raise ValueError(f"Invalid descriptor bounds for {name}")
        if bounds["MW"][1] <= 0:
            raise ValueError("Training maximum molecular weight must be positive")
        training_smiles = bundle.get("training_smiles")
        identities = bundle.get("training_identities")
        training_names = bundle.get("training_names")
        encoded_fps = bundle.get("training_fingerprints")
        if not isinstance(training_smiles, list) or not training_smiles or not all(isinstance(s, str) and s for s in training_smiles):
            raise ValueError("Missing training SMILES")
        count = len(training_smiles)
        if not isinstance(identities, list) or len(identities) != count or not all(isinstance(s, str) and s for s in identities):
            raise ValueError("Training identities do not align with training molecules")
        if not isinstance(training_names, list) or len(training_names) != count or not all(name is None or isinstance(name, str) for name in training_names):
            raise ValueError("Training names do not align with training molecules")
        if not isinstance(encoded_fps, list) or len(encoded_fps) != count:
            raise ValueError("Training fingerprints do not align with training molecules")
        fingerprints = []
        for encoded in encoded_fps:
            if not isinstance(encoded, bytes) or len(encoded) != 256:
                raise ValueError("Training fingerprint must contain exactly 2048 bits")
            fp = DataStructs.CreateFromBinaryText(encoded)
            if fp.GetNumBits() != 2048:
                raise ValueError("Unexpected training fingerprint size")
            fingerprints.append(fp)
        return LoadedModel(bundle, tuple(fingerprints), frozenset(identities))

    def health(self) -> dict[str, Any]:
        complete = len(self._models) == len(TARGETS) and self.database_ready
        return {
            "status": "ok" if complete else "degraded",
            "model_loaded": self.model_loaded,
            "models_loaded": self.models_loaded,
            "database_ready": self.database_ready,
        }

    def model_info(self) -> dict[str, Any]:
        if not self._models:
            raise ModelUnavailable("没有可用模型，请先运行训练脚本并重新启动服务。")
        catalog = []
        for target, defaults in TARGETS.items():
            loaded = self._models.get(target)
            if loaded is None:
                catalog.append({"id": target, **defaults, "available": False, "status": "unavailable"})
            else:
                item = {"id": target, **loaded.bundle["metadata"], "available": True}
                path = self.models_dir.parent / "reports" / target / "candidate_validation.csv"
                if path.is_file():
                    try:
                        with path.open(encoding="utf-8-sig", newline="") as handle:
                            candidates = list(csv.DictReader(handle))
                        for row in candidates:
                            for name in ("rmse", "mae", "r2", "fit_seconds"):
                                if name in row:
                                    row[name] = float(row[name])
                                    if not math.isfinite(row[name]):
                                        raise ValueError("Candidate evidence contains a non-finite metric")
                            if "feature_count" in row:
                                row["feature_count"] = int(row["feature_count"])
                        item["candidate_validation"] = candidates
                    except (OSError, ValueError, TypeError, csv.Error):
                        LOGGER.exception("Could not read candidate validation evidence for %s", target)
                catalog.append(item)
        return {
            "version": "3.0.0",
            "models": catalog,
            "models_loaded": self.models_loaded,
            "database_ready": self.database_ready,
            "interval_note": INTERVAL_WARNING,
        }

    def database_stats(self) -> dict[str, Any]:
        if self._database is None:
            raise DatabaseUnavailable("本地实验观测数据库不可用，请先运行数据库构建脚本。")
        try:
            return self._database.stats()
        except Exception as exc:
            LOGGER.exception("Local observation database statistics failed")
            raise DatabaseUnavailable("本地实验观测数据库暂时无法读取。") from exc

    def database_search(self, query: str = "", target: str | None = None, limit: int = 20, offset: int = 0) -> dict[str, Any]:
        if self._database is None:
            raise DatabaseUnavailable("本地实验观测数据库不可用，请先运行数据库构建脚本。")
        try:
            return self._database.search(query=query, target=target, limit=limit, offset=offset)
        except ValueError as exc:
            raise MoleculeError("INVALID_DATABASE_QUERY", str(exc)) from exc
        except Exception as exc:
            LOGGER.exception("Local observation database search failed")
            raise DatabaseUnavailable("本地实验观测数据库暂时无法查询。") from exc

    @staticmethod
    def _unavailable_prediction(target: str, message: str, loaded: LoadedModel | None = None) -> dict[str, Any]:
        metadata = loaded.bundle["metadata"] if loaded else TARGETS[target]
        return {
            "label": metadata["label"], "unit": metadata["unit"],
            "status": "unavailable", "value": None, "interval90": None,
            "model_name": metadata.get("model_name"), "training_match": False,
            "nearest_similarity": None, "nearest_neighbors": [],
            "domain_violations": [], "warnings": [message],
            "ensemble": {"member_values": [], "spread": None, "spread_kind": "weighted_population_std", "weights": []},
        }

    @staticmethod
    def _target_prediction(
        target: str,
        loaded: LoadedModel,
        mol: Any,
        identity: str,
        descriptors: dict[str, float],
        query_fingerprint: Any,
        vectors: dict[str, np.ndarray],
    ) -> dict[str, Any]:
        bundle = loaded.bundle
        metadata = bundle["metadata"]
        bounds = bundle["descriptor_bounds"]
        range_violations = [name for name in DOMAIN_FEATURES if not bounds[name][0] <= descriptors[name] <= bounds[name][1]]
        molecular_weight_limit = min(1200.0, float(bounds["MW"][1]) * 1.25)
        violations = [f"{name} 超出训练集范围" for name in range_violations]
        too_heavy = descriptors["MW"] > molecular_weight_limit
        too_many_atoms = descriptors["HeavyAtoms"] > 80
        if too_heavy:
            violations.append(f"MW 超过适用域上限 {molecular_weight_limit:.2f} g/mol")
        if too_many_atoms:
            violations.append("重原子数超过 80")
        severe = too_heavy or too_many_atoms or len(range_violations) >= 3
        similarities = DataStructs.BulkTanimotoSimilarity(query_fingerprint, loaded.fingerprints)
        nearest_indices = sorted(range(len(similarities)), key=similarities.__getitem__, reverse=True)[:3]
        nearest_similarity = float(similarities[nearest_indices[0]])
        nearest = [{
            "similarity": float(similarities[index]),
            "canonical_smiles": bundle["training_smiles"][index],
            "name": bundle["training_names"][index] or bundle["training_smiles"][index],
        } for index in nearest_indices]
        training_match = identity in loaded.training_identities
        warnings = [INTERVAL_WARNING]
        if training_match:
            warnings.append("模型训练集合包含该分子的连接结构（可能包括立体异构体）；此结果不是独立外部验证。")
        if nearest_similarity < 0.3:
            warnings.append("与训练集中最相似分子的 Morgan 指纹 Tanimoto 相似度低于 0.30，外推可靠性可能降低。")
        if range_violations:
            warnings.append("存在超出训练集范围的描述符；适用域筛选是启发式规则，不是可靠性的保证。")
        status = "out_of_domain" if severe else "caution" if range_violations or nearest_similarity < 0.3 else "ok"
        prediction: dict[str, Any] = {
            "label": metadata["label"], "unit": metadata["unit"],
            "status": status, "value": None, "interval90": None,
            "model_name": metadata["model_name"], "training_match": training_match,
            "nearest_similarity": nearest_similarity, "nearest_neighbors": nearest,
            "domain_violations": violations, "warnings": warnings,
            "ensemble": {"member_values": [], "spread": None, "spread_kind": "weighted_population_std", "weights": []},
        }
        if severe:
            warnings.append("分子明显超出模型适用域，已停止数值预测并隐藏区间。")
            return prediction
        mode = bundle["feature_mode"]
        if mode not in vectors:
            vectors[mode] = np.asarray(feature_vector(mol, mode), dtype=np.float32).reshape(1, -1)
        vector: Any = vectors[mode]
        if vector.shape != (1, len(bundle["feature_names"])) or np.isinf(vector).any():
            raise ValueError("Input feature vector does not match the trained estimator")
        estimator = bundle["model"]
        if getattr(estimator, "feature_names_in_", None) is not None:
            import pandas as pd

            vector = pd.DataFrame(vector, columns=bundle["feature_names"])
        values = np.asarray(estimator.predict(vector), dtype=float).reshape(-1)
        if values.size != 1 or not math.isfinite(float(values[0])):
            raise ValueError("Estimator returned a non-finite or non-scalar prediction")
        value = float(values[0])
        radius = float(bundle["conformal_radius"])
        interval = [value - radius, value + radius]
        if not all(math.isfinite(endpoint) for endpoint in interval):
            raise ValueError("Calibrated interval is not finite")
        prediction.update(value=value, interval90=interval)
        if callable(getattr(estimator, "predict_members", None)):
            member_values = np.asarray(estimator.predict_members(vector), dtype=float)
            names = list(estimator.member_names_)
            weights = np.asarray(estimator.member_weights_, dtype=float)
            if member_values.shape != (1, len(names)) or weights.shape != (len(names),) or not np.isfinite(member_values).all() or not np.isfinite(weights).all() or np.any(weights < 0) or not np.isclose(weights.sum(), 1):
                raise ValueError("Invalid ensemble member evidence")
            if not np.isclose(float(member_values[0] @ weights), value, rtol=1e-8, atol=1e-8):
                raise ValueError("Ensemble evidence does not reproduce its reported prediction")
            spread = float(np.sqrt(np.sum(weights * (member_values[0] - value) ** 2)))
            prediction["ensemble"] = {"member_values": [{"name": name, "value": float(member)} for name, member in zip(names, member_values[0])], "spread": spread, "spread_kind": "weighted_population_std", "weights": weights.tolist()}
            warnings.append("集成成员差异仅表示模型间分歧，不是实验误差或校准置信区间。")
        return prediction

    def predict(self, smiles: str, properties: list[str] | None = None) -> dict[str, Any]:
        started = time.perf_counter()
        requested = list(TARGET_IDS) if properties is None else list(properties)
        if not requested or len(requested) != len(set(requested)) or any(target not in TARGETS for target in requested):
            raise MoleculeError("INVALID_PROPERTIES", "请选择互不重复的已支持性质：logS、logD74、hydration_free_energy。")
        mol = parse_smiles(smiles)
        if not self._models:
            raise ModelUnavailable("没有可用模型，请先运行训练脚本并重新启动服务。")
        canonical = canonical_smiles(mol)
        identity = identity_smiles(mol)
        try:
            descriptors = {name: float(value) for name, value in all_descriptors(mol).items()}
        except (ValueError, RuntimeError) as exc:
            raise MoleculeError("INVALID_DESCRIPTORS", "此分子无法生成有限的模型描述符。", valid=True) from exc
        if not all(math.isfinite(descriptors[name]) for name in DOMAIN_FEATURES):
            raise MoleculeError("INVALID_DESCRIPTORS", "此分子无法生成有限的模型描述符。", valid=True)
        missing_descriptors = [name for name, value in descriptors.items() if not math.isfinite(value)]
        query_fingerprint = fingerprint(mol)
        vectors: dict[str, np.ndarray] = {}
        predictions = {}
        for target in requested:
            loaded = self._models.get(target)
            if loaded is None:
                predictions[target] = self._unavailable_prediction(target, "此性质的模型尚未加载，未生成预测值。")
                continue
            try:
                predictions[target] = self._target_prediction(target, loaded, mol, identity, descriptors, query_fingerprint, vectors)
            except Exception:
                LOGGER.exception("Prediction failed for target %s", target)
                predictions[target] = self._unavailable_prediction(target, "此性质的模型无法完成本次预测，请检查服务端日志。", loaded)
        warnings = ["模型估计与数据库原始观测分别展示；匹配的数据库记录不代表该结构存在唯一、无条件的实验真值。"]
        if missing_descriptors:
            warnings.append("部分扩展描述符不可计算，已显示为空值；模型仅使用其训练阶段拟合的缺失值填补规则。")
        if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
            warnings.append("输入含带电原子；模型未执行酸碱微态枚举、pKa 预测或离子化校正，性质可能依赖 pH 和实验条件。")
        measurements: list[dict[str, Any]] = []
        if self._database is None:
            warnings.append("本地观测数据库不可用，无法查询对应实验或文献记录。")
        else:
            try:
                measurements = self._database.lookup(canonical)
            except Exception:
                LOGGER.exception("Observation lookup failed")
                warnings.append("数据库记录查询失败；下方仅包含模型估计。")
        log_s = predictions.get("logS", {}).get("value")
        concentration = None
        if log_s is not None:
            try:
                concentration = float(10.0**log_s)
                if not math.isfinite(concentration) or concentration <= 0:
                    concentration = None
            except OverflowError:
                concentration = None
        return {
            "valid": True,
            "canonical_smiles": canonical,
            "molecule": known_name(mol) or canonical,
            "formula": molecule_formula(mol),
            "descriptor": {name: value if math.isfinite(value) else None for name, value in descriptors.items()},
            "missing_descriptors": missing_descriptors,
            "structure_svg": molecule_svg(mol),
            "predictions": predictions,
            "measurements": measurements,
            "warnings": warnings,
            "logS_prediction": log_s,
            "solubility_mol_L": concentration,
            "conformers": empty_conformers(),
            "timing": {"prediction_seconds": time.perf_counter() - started, "conformer_seconds": 0.0, "total_seconds": time.perf_counter() - started},
            "evidence": {"predictions": "machine_learning", "measurements": "source_observations", "conformers": "not_computed", "quantum_calculation_executed": False},
        }


class AssignmentPredictionService:
    """The original ESOL model is a separate, explicitly labelled compatibility endpoint."""

    def __init__(self, path=None):
        self.path = Path(path) if path is not None else PROJECT_ROOT / "model.pkl"
        self.bundle = None

    def load(self):
        try:
            with self.path.open("rb") as handle:
                bundle = pickle.load(handle)
            if bundle.get("schema_version") != 1 or bundle.get("feature_names") != list(DOMAIN_FEATURES):
                raise ValueError("Original assignment feature schema does not match")
            if bundle["metadata"]["sklearn_version"] != sklearn.__version__ or bundle["metadata"]["rdkit_version"] != rdBase.rdkitVersion:
                raise ValueError("Original assignment model runtime version does not match")
            check_is_fitted(bundle["model"])
            self.bundle = bundle
        except Exception:
            LOGGER.exception("Could not load the original ESOL assignment artifact")

    def predict(self, smiles):
        mol = parse_smiles(smiles)
        if self.bundle is None:
            raise ModelUnavailable("题目兼容模型 model.pkl 不可用。")
        bundle = self.bundle
        canonical = canonical_smiles(mol)
        descriptors = calculate_descriptors(mol)
        bounds = bundle["descriptor_bounds"]
        violations = [name for name in DOMAIN_FEATURES if not bounds[name][0] <= descriptors[name] <= bounds[name][1]]
        severe = descriptors["MW"] > min(1200, bounds["MW"][1] * 1.25) or descriptors["HeavyAtoms"] > 80 or len(violations) >= 3
        value = None
        if not severe:
            value = float(bundle["model"].predict(np.asarray([[descriptors[name] for name in DOMAIN_FEATURES]], dtype=float))[0])
            if not math.isfinite(value):
                raise ModelUnavailable("题目兼容模型未返回有限预测值。")
        warnings = ["这是题目兼容接口，使用原始纯 ESOL 随机森林 model.pkl；与当前多数据集主模型分开。", "该模型输出是估计，不是实验测量。"]
        if severe:
            warnings.append("输入明显超出原模型适用域，已隐藏数值预测。")
        return {"valid": True, "molecule": known_name(mol) or canonical, "canonical_smiles": canonical,
                "formula": molecule_formula(mol), "descriptor": descriptors, "structure_svg": molecule_svg(mol),
                "logS_prediction": value, "solubility_mol_L": 10.0**value if value is not None else None,
                "model": {"name": bundle["metadata"]["model_name"], "version": bundle["metadata"]["model_version"], "training_match": canonical in bundle["training_smiles"], "purpose": "assignment_compatibility"},
                "applicability_domain": {"status": "out_of_domain" if severe else "caution" if violations else "ok", "violations": violations}, "warnings": warnings}
