"""V3: seven prespecified candidates, validation-only blending, frozen V2 splits.

CPU is limited to two threads. Float32 features are disk cached. Candidate
estimators are fitted sequentially and immediately released after validation.
Only the final one or two selected members are retained/refitted.
"""
from __future__ import annotations
import argparse
import csv
import gc
import hashlib
import itertools
import json
import math
import os
import pickle
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "2"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import sklearn
from rdkit import DataStructs, rdBase
from sklearn.base import clone
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from threadpoolctl import threadpool_limits
from app.chemistry import FEATURE_NAMES as DOMAIN_FEATURES, calculate_descriptors, parse_smiles
from app.features import FEATURE_NAMES, FEATURE_VERSION, MISSING_VALUE_POLICY, feature_vector, fingerprint, model_feature_names
from app.ensemble import FeatureSubset, ValidationEnsemble

TARGETS = {
    "logS": ("水溶解度 logS", "log10(mol/L)"),
    "logD74": ("分配系数 logD（pH 7.4）", "log10(D), pH 7.4"),
    "hydration_free_energy": ("水合自由能", "kcal/mol"),
}
SEED = 42
RULE = ("Fit seven prespecified candidates on train only. Rank by validation RMSE (name tie-break). "
        "Among the top three evaluate all pairs with least-squares simplex weights clipped to [0.1,0.9]. "
        "Use the lowest-RMSE pair only if it improves best-single validation RMSE by at least 1%; "
        "otherwise retain the best single. Freeze weights; refit selected members on train+validation. "
        "Calibration/test labels never enter candidate selection or weight fitting.")


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields=None):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metrics(y, prediction):
    return {"rmse": float(math.sqrt(mean_squared_error(y, prediction))),
            "mae": float(mean_absolute_error(y, prediction)), "r2": float(r2_score(y, prediction))}


def candidates():
    """Predeclared before any V3 test evaluation; no automatic test tuning."""
    n = len(FEATURE_NAMES)
    sets = {
        "legacy28": tuple(range(28)),
        "legacy28_morgan": tuple(range(28)) + tuple(range(n, n + 2048)),
        "expanded2d": tuple(range(n)),
        "expanded2d_morgan": tuple(range(n + 2048)),
    }
    specifications = [
        ("RandomForest_28D", "legacy28", RandomForestRegressor(n_estimators=160, max_depth=18,
             min_samples_leaf=2, max_features=.85, random_state=SEED, n_jobs=2), None),
        ("ExtraTrees_28D_Morgan2048", "legacy28_morgan", ExtraTreesRegressor(n_estimators=192, max_depth=18,
             min_samples_leaf=2, max_features=.25, random_state=SEED, n_jobs=2), None),
        ("ExtraTrees_219D", "expanded2d", ExtraTreesRegressor(n_estimators=180, max_depth=18,
             min_samples_leaf=2, max_features=.8, random_state=SEED, n_jobs=2), None),
        ("HistGradientBoosting_219D", "expanded2d", HistGradientBoostingRegressor(max_iter=350,
             max_leaf_nodes=15, min_samples_leaf=15, learning_rate=.05, l2_regularization=2.5,
             early_stopping=False, random_state=SEED), None),
        ("Ridge_219D_Morgan2048", "expanded2d_morgan", Ridge(alpha=30., solver="lsqr", max_iter=800,
             tol=.001, copy_X=False), False),
        ("SVR_219D", "expanded2d", SVR(C=8., epsilon=.1, gamma="scale", cache_size=96), True),
        ("ExtraTrees_219D_Morgan2048", "expanded2d_morgan", ExtraTreesRegressor(n_estimators=160, max_depth=18,
             min_samples_leaf=2, max_features=.25, random_state=SEED, n_jobs=2), None),
    ]
    result = []
    for name, representation, estimator, center in specifications:
        steps = [("representation", FeatureSubset(sets[representation])),
                 ("imputer", SimpleImputer(strategy="median", keep_empty_features=True, copy=False))]
        if center is not None:
            steps.append(("scaler", StandardScaler(with_mean=center, copy=False)))
        steps.append(("regressor", estimator))
        result.append({"name": name, "representation": representation, "feature_count": len(sets[representation]),
                       "parameters": estimator.get_params(), "scaling": center,
                       "estimator": Pipeline(steps)})
    return result


def select_validation(y_validation, predictions):
    """Pure selection function deliberately accepts no calibration/test data."""
    scores = {name: metrics(y_validation, values) for name, values in predictions.items()}
    ranked = sorted(scores, key=lambda name: (scores[name]["rmse"], name))
    best = ranked[0]
    pairs = []
    for left, right in itertools.combinations(ranked[:3], 2):
        a, b = predictions[left], predictions[right]
        diff = a - b
        denom = float(np.dot(diff, diff))
        weight = float(np.clip(np.dot(y_validation - b, diff) / denom, .1, .9)) if denom > 0 else .5
        score = metrics(y_validation, weight * a + (1 - weight) * b)
        pairs.append({"member_names": [left, right], "weights": [weight, 1 - weight], "metrics": score})
    pair = min(pairs, key=lambda item: (item["metrics"]["rmse"], item["member_names"]))
    if pair["metrics"]["rmse"] <= .99 * scores[best]["rmse"]:
        selected = pair
    else:
        selected = {"member_names": [best], "weights": [1.], "metrics": scores[best]}
    return {**selected, "rule": RULE, "best_single": best, "ranked_candidates": ranked,
            "relative_improvement_required": .01, "pairs_considered": pairs,
            "selection_label_partition": "validation", "candidate_fit_partition": "train"}


def frozen_partitions(report_dir):
    """Snapshot V2 evidence once, then use exactly those rows/labels/splits."""
    baseline_path = report_dir / "v2_baseline_metrics.json"
    if not baseline_path.exists():
        source = json.loads((report_dir / "metrics.json").read_text(encoding="utf-8"))
        if not source["version"].startswith("physchem-v2-"):
            raise ValueError("Missing original V2 baseline; cannot silently invent fresh partitions")
        shutil.copyfile(report_dir / "metrics.json", baseline_path)
        shutil.copyfile(report_dir / "test_predictions.csv", report_dir / "v2_test_predictions.csv")
        shutil.copyfile(report_dir / "split_manifest.csv", report_dir / "frozen_split_manifest.csv")
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    frozen = report_dir / "frozen_split_manifest.csv"
    if sha256(frozen) != baseline["split_sha256"]:
        raise ValueError("Frozen V2 split hash changed")
    if sha256(ROOT / "data/curated_observations.csv") != baseline["data_sha256"]:
        raise ValueError("Input observations changed; requires a separately documented experiment")
    shutil.copyfile(frozen, report_dir / "split_manifest.csv")
    rows = read_csv(frozen)
    for row in rows:
        row["value"] = float(row["value"])
    indices = {key: np.array([i for i, row in enumerate(rows) if row["split"] == key], dtype=int)
               for key in ("train", "validation", "calibration", "test")}
    for left, right in itertools.combinations(indices, 2):
        for key in ("identity_smiles", "scaffold"):
            if {rows[i][key] for i in indices[left]} & {rows[i][key] for i in indices[right]}:
                raise ValueError(f"Frozen partition leakage: {left}/{right}/{key}")
    return rows, indices, baseline


def prepare_features(rows, target, cache_dir, split_hash):
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256((FEATURE_VERSION + split_hash).encode()).hexdigest()[:16]
    path = cache_dir / f"{target}-{key}.npy"
    audit_path = cache_dir / f"{target}-{key}.json"
    if not path.exists() or not audit_path.exists():
        matrix = np.lib.format.open_memmap(path, mode="w+", dtype=np.float32,
                                          shape=(len(rows), len(model_feature_names("combined"))))
        missing = {split: np.zeros(len(FEATURE_NAMES), dtype=int) for split in ("train", "validation", "calibration", "test")}
        for i, row in enumerate(rows):
            matrix[i] = feature_vector(parse_smiles(row["canonical_smiles"]), "combined")
            missing[row["split"]] += np.isnan(matrix[i, :len(FEATURE_NAMES)])
            if (i + 1) % 1000 == 0:
                print(f"{target}: feature cache {i + 1}/{len(rows)}", flush=True)
        matrix.flush()
        del matrix
        audit = {"feature_version": FEATURE_VERSION, "feature_names": model_feature_names("combined"),
                 "dtype": "float32", "shape": [len(rows), len(model_feature_names("combined"))],
                 "missing_value_policy": MISSING_VALUE_POLICY,
                 "missing_counts": {split: {name: int(count) for name, count in zip(FEATURE_NAMES, counts) if count}
                                    for split, counts in missing.items()}}
        write_json(audit_path, audit)
    matrix = np.load(path, mmap_mode="r")
    if matrix.shape != (len(rows), len(model_feature_names("combined"))) or np.isinf(matrix).any():
        raise ValueError("Invalid feature cache")
    return matrix, json.loads(audit_path.read_text(encoding="utf-8"))


def train_target(target, cache_dir):
    started = time.perf_counter()
    report_dir = ROOT / "reports" / target
    rows, indices, previous = frozen_partitions(report_dir)
    configurations = candidates()
    protocol = {"version": "v3-prespecified-7-candidates", "created_utc": datetime.now(timezone.utc).isoformat(),
                "selection_rule": RULE, "seed": SEED, "max_threads": 2,
                "candidates": [{key: value for key, value in spec.items() if key != "estimator"} for spec in configurations],
                "split_sha256": sha256(report_dir / "split_manifest.csv"),
                "data_sha256": sha256(ROOT / "data/curated_observations.csv"),
                "feature_version": FEATURE_VERSION, "missing_value_policy": MISSING_VALUE_POLICY}
    write_json(report_dir / "training_protocol.json", protocol)
    print(f"{target}: frozen partitions { {k:len(v) for k,v in indices.items()} }; protocol written before fitting", flush=True)
    matrix, audit = prepare_features(rows, target, cache_dir, protocol["split_sha256"])
    write_json(report_dir / "feature_audit.json", audit)
    y = np.asarray([row["value"] for row in rows], dtype=float)
    train_x, validation_x = matrix[indices["train"]], matrix[indices["validation"]]
    predictions, leaderboard = {}, []
    for spec in configurations:
        tick = time.perf_counter()
        print(f"{target}: fitting {spec['name']} ({spec['feature_count']} selected inputs)", flush=True)
        fitted = clone(spec["estimator"]).fit(train_x, y[indices["train"]])
        predictions[spec["name"]] = fitted.predict(validation_x)
        score = metrics(y[indices["validation"]], predictions[spec["name"]])
        leaderboard.append({"model_name": spec["name"], "feature_mode": "combined", "representation": spec["representation"],
                            "feature_count": spec["feature_count"], **score, "fit_seconds": time.perf_counter() - tick})
        print(f"{target}: {spec['name']} validation RMSE={score['rmse']:.6f}; {leaderboard[-1]['fit_seconds']:.1f}s", flush=True)
        del fitted
        gc.collect()
    del train_x, validation_x
    selection = select_validation(y[indices["validation"]], predictions)
    selected_prediction = sum(weight * predictions[name] for name, weight in zip(selection["member_names"], selection["weights"]))
    write_csv(report_dir / "candidate_validation.csv", leaderboard)
    write_json(report_dir / "ensemble_selection.json", selection)
    validation_rows = [{"canonical_smiles": rows[i]["canonical_smiles"], "identity_smiles": rows[i]["identity_smiles"],
                        "y_true": float(y[i]), **{name: float(values[j]) for name, values in predictions.items()},
                        "selected_prediction": float(selected_prediction[j])} for j, i in enumerate(indices["validation"])]
    write_csv(report_dir / "validation_predictions.csv", validation_rows)
    by_name = {spec["name"]: spec for spec in configurations}
    members = [(name, by_name[name]["estimator"]) for name in selection["member_names"]]
    fit_indices = np.concatenate([indices["train"], indices["validation"]])
    fit_rows = [rows[i] for i in fit_indices]
    print(f"{target}: selected {selection['member_names']} weights={selection['weights']}; refitting train+validation", flush=True)
    fit_x = matrix[fit_indices]
    model = ValidationEnsemble(members, tuple(selection["weights"])).fit(fit_x, y[fit_indices])
    del fit_x
    gc.collect()
    # First use of calibration labels and test predictions occurs after selection
    # has been saved and final fitting has finished.
    calibration_prediction = model.predict(matrix[indices["calibration"]])
    errors = np.abs(calibration_prediction - y[indices["calibration"]])
    rank = min(len(errors), math.ceil((len(errors) + 1) * .9))
    radius = float(np.sort(errors)[rank - 1])
    write_csv(report_dir / "calibration_predictions.csv", [
        {"canonical_smiles": rows[i]["canonical_smiles"], "y_true": float(y[i]), "y_pred": float(pred), "absolute_error": float(error)}
        for i, pred, error in zip(indices["calibration"], calibration_prediction, errors)])
    test_members = model.predict_members(matrix[indices["test"]])
    predicted = test_members @ model.member_weights_
    test_y = y[indices["test"]]
    lower, upper = predicted - radius, predicted + radius
    test_metrics = metrics(test_y, predicted)
    coverage = float(np.mean((test_y >= lower) & (test_y <= upper)))
    write_csv(report_dir / "test_predictions.csv", [
        {"canonical_smiles": rows[i]["canonical_smiles"], "identity_smiles": rows[i]["identity_smiles"],
         "y_true": float(y[i]), "y_pred": float(pred), "lower90": float(lo), "upper90": float(hi)}
        for i, pred, lo, hi in zip(indices["test"], predicted, lower, upper)])
    write_csv(report_dir / "test_member_predictions.csv", [
        {"canonical_smiles": rows[i]["canonical_smiles"], **{name: float(test_members[j, k]) for k, name in enumerate(model.member_names_)},
         "weighted_population_std": float(np.sqrt(np.dot(model.member_weights_, (test_members[j] - predicted[j]) ** 2)))}
        for j, i in enumerate(indices["test"])])
    name = selection["member_names"][0] if len(members) == 1 else "ValidationEnsemble[" + "+".join(selection["member_names"]) + "]"
    label, unit = TARGETS[target]
    metadata = {
        "target": target, "label": label, "unit": unit, "version": f"physchem-v3-{target}",
        "created_utc": datetime.now(timezone.utc).isoformat(), "model_name": name, "feature_mode": "combined",
        "feature_version": FEATURE_VERSION, "feature_count": len(model_feature_names("combined")),
        "descriptor_count": len(FEATURE_NAMES), "fingerprint_bits": 2048, "candidate_count": len(configurations),
        "sklearn_version": sklearn.__version__, "rdkit_version": rdBase.rdkitVersion, "numpy_version": np.__version__,
        "seed": SEED, "parameters": {member: by_name[member]["parameters"] for member in selection["member_names"]},
        "ensemble": {"member_names": selection["member_names"], "weights": selection["weights"], "selection_rule": RULE},
        "missing_value_policy": MISSING_VALUE_POLICY,
        "split_counts": {key: len(value) for key, value in indices.items()}, "fit_count": len(fit_rows),
        "unique_identity_count_before_quarantine": previous["unique_identity_count_before_quarantine"],
        "quarantined_from_fit": previous["quarantined_from_fit"], "split_method": previous["split_method"] + "; frozen from V2 byte-for-byte",
        "identity_overlap": 0, "scaffold_overlap": 0, "selection": RULE, "aggregation": previous["aggregation"],
        "metrics": {"test": test_metrics, "baseline": metrics(test_y, np.repeat(y[fit_indices].mean(), len(test_y))),
                    "validation_selected": selection["metrics"]},
        "calibration": {"radius": radius, "nominal_coverage": .9, "test_coverage": coverage,
                        "n_calibration": len(errors), "quantile_rank": rank, "mean_width": 2 * radius},
        "evaluation_note": previous["evaluation_note"],
        "data_sha256": protocol["data_sha256"], "split_sha256": protocol["split_sha256"],
        "training_protocol_sha256": sha256(report_dir / "training_protocol.json"),
        "selection_evidence_sha256": sha256(report_dir / "validation_predictions.csv"),
        "limitations": previous["limitations"] + ["Only seven fixed candidates; no external validation or guarantee of improvement.",
            "Member disagreement describes the fitted models; it is not calibrated uncertainty.",
            "V2 and V3 share a previously inspected test set; comparison is descriptive and not a fresh prospective test."],
        "v2_same_test_comparison": {"test_count": len(test_y), "v2_metrics": previous["metrics"]["test"], "v3_metrics": test_metrics,
                                    "rmse_change": test_metrics["rmse"] - previous["metrics"]["test"]["rmse"],
                                    "note": "Identical V2 rows, labels, training/validation/calibration/test assignment; no test-driven retries."},
    }
    if target == "logS":
        old = json.loads((ROOT / "reports/v1_baseline_metrics.json").read_text(encoding="utf-8"))
        metadata["v1_same_test_comparison"] = {"test_count": len(test_y), "v1_rmse": old["metrics"]["test"]["rmse"],
            "v3_rmse": test_metrics["rmse"], "rmse_change": test_metrics["rmse"] - old["metrics"]["test"]["rmse"],
            "note": "Original 224 ESOL labels; matching test identities/scaffolds excluded from every other partition."}
    domain_values, fps = [], []
    for row in fit_rows:
        mol = parse_smiles(row["canonical_smiles"])
        domain_values.append(calculate_descriptors(mol))
        fps.append(DataStructs.BitVectToBinaryText(fingerprint(mol)))
    bounds = {key: [float(min(values[key] for values in domain_values)), float(max(values[key] for values in domain_values))]
              for key in DOMAIN_FEATURES}
    bundle = {"schema_version": 2, "feature_mode": "combined", "feature_names": model_feature_names("combined"),
              "model": model, "metadata": metadata, "conformal_radius": radius, "descriptor_bounds": bounds,
              "training_smiles": [row["canonical_smiles"] for row in fit_rows],
              "training_identities": [row["identity_smiles"] for row in fit_rows],
              "training_names": [row["name"] or None for row in fit_rows], "training_fingerprints": fps}
    model_path = ROOT / "models" / f"{target}.pkl"
    with model_path.open("wb") as handle:
        pickle.dump(bundle, handle, protocol=pickle.HIGHEST_PROTOCOL)
    report = {**metadata, "model_sha256": sha256(model_path), "elapsed_seconds": time.perf_counter() - started}
    write_json(report_dir / "metrics.json", report)
    print(f"{target}: COMPLETE test RMSE={test_metrics['rmse']:.6f}, MAE={test_metrics['mae']:.6f}, coverage={coverage:.3%}, radius={radius:.6f}", flush=True)
    del matrix, bundle, model, domain_values, predictions, configurations
    gc.collect()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=list(TARGETS), action="append")
    parser.add_argument("--cache-dir", type=Path, default=ROOT.parent.parent / "work/v3-model-training")
    args = parser.parse_args()
    with threadpool_limits(limits=2):
        for target in args.target or list(TARGETS):
            train_target(target, args.cache_dir)
    catalog = {"version": "3.0.0", "feature_version": FEATURE_VERSION, "models": []}
    for target in TARGETS:
        path = ROOT / "reports" / target / "metrics.json"
        if path.exists():
            catalog["models"].append(json.loads(path.read_text(encoding="utf-8")))
    write_json(ROOT / "reports/model_catalog.json", catalog)


if __name__ == "__main__":
    main()
