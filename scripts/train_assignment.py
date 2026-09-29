"""Train a reproducible ESOL baseline; hold out whole molecular scaffolds.

Run from any directory: python path/to/scripts/train.py
No selection/tuning is performed on the held-out test set.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import pickle
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import rdkit
from rdkit.Chem.Scaffolds import MurckoScaffold
import sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from app.chemistry import (ALLOWED_ELEMENTS, FEATURE_NAMES, MoleculeError,
                           calculate_descriptors, canonical_smiles, parse_smiles)

TARGET = "measured log solubility in mols per litre"
EXPECTED_SHA256 = "8c06a76f0c6487d29ab0f903e6a7a7139f189ab3c1178f159c8be8964602f189"


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def save_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def score(y, pred):
    return {"rmse": float(np.sqrt(mean_squared_error(y, pred))),
            "mae": float(mean_absolute_error(y, pred)), "r2": float(r2_score(y, pred))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data/raw/delaney-processed.csv")
    parser.add_argument("--output", type=Path, default=ROOT)
    args = parser.parse_args()
    if not args.data.exists():
        raise SystemExit("ESOL CSV missing. Restore the included data/raw snapshot.")
    data_hash = hashlib.sha256(args.data.read_bytes()).hexdigest()
    if data_hash != EXPECTED_SHA256:
        raise SystemExit("CSV hash differs from the recorded ESOL source. Verify provenance before training.")
    with args.data.open(encoding="utf-8-sig", newline="") as f:
        raw = list(csv.DictReader(f))
    grouped = defaultdict(list)
    excluded = []
    for i, row in enumerate(raw):
        try:
            mol = parse_smiles(row["smiles"])
            target = float(row[TARGET])
            if not np.isfinite(target):
                raise ValueError("Non-finite measured label")
            descriptors = calculate_descriptors(mol)
            smi = canonical_smiles(mol)
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
            grouped[smi].append((target, descriptors, scaffold, i + 2, row["Compound ID"]))
        except (MoleculeError, ValueError) as exc:
            excluded.append({"csv_line": i + 2, "smiles": row["smiles"], "reason": str(exc)})
    # Aggregate repeated measurements of the same canonical isomeric structure.
    # All replicates stay together and contribute exactly one supervised example.
    records = []
    duplicate_rows = []
    for smi, entries in sorted(grouped.items()):
        values = [e[0] for e in entries]
        records.append({"canonical_smiles": smi, "scaffold": entries[0][2],
                        "logS": float(np.mean(values)), "replicate_count": len(entries),
                        "features": entries[0][1]})
        if len(entries) > 1:
            duplicate_rows.append({"canonical_smiles": smi, "replicates": len(entries),
                                   "values": values, "mean": float(np.mean(values)),
                                   "csv_lines": [e[3] for e in entries]})
    scaffold_groups = defaultdict(list)
    for idx, record in enumerate(records):
        scaffold_groups[record["scaffold"]].append(idx)
    # Deterministic largest-first allocation. All acyclic molecules share the
    # empty scaffold group. This is a deliberately strict, asymmetric holdout,
    # not a random split nor the claim of an official MoleculeNet benchmark.
    train_idx, test_idx = [], []
    target_train = int(len(records) * 0.8)
    for scaffold, ids in sorted(scaffold_groups.items(), key=lambda item: (-len(item[1]), item[0])):
        (train_idx if len(train_idx) + len(ids) <= target_train else test_idx).extend(ids)
    if len(test_idx) < 2 or len(train_idx) < 2:
        raise SystemExit("Insufficient independent scaffold groups to train/evaluate.")
    X = np.asarray([[r["features"][f] for f in FEATURE_NAMES] for r in records], dtype=float)
    y = np.asarray([r["logS"] for r in records])
    parameters = {"n_estimators": 400, "min_samples_leaf": 2, "max_features": 1.0,
                  "random_state": 42, "n_jobs": 2}
    model = RandomForestRegressor(**parameters)
    model.fit(X[train_idx], y[train_idx])
    test_pred = model.predict(X[test_idx])
    baseline = np.full(len(test_idx), float(np.mean(y[train_idx])))
    metrics = {"train": score(y[train_idx], model.predict(X[train_idx])),
               "test": score(y[test_idx], test_pred), "baseline": score(y[test_idx], baseline)}
    training_smiles = [records[i]["canonical_smiles"] for i in train_idx]
    train_scaffolds = {records[i]["scaffold"] for i in train_idx}
    test_scaffolds = {records[i]["scaffold"] for i in test_idx}
    assert train_scaffolds.isdisjoint(test_scaffolds)
    assert set(training_smiles).isdisjoint(records[i]["canonical_smiles"] for i in test_idx)
    metadata = {
        "model_name": "RandomForestRegressor", "model_version": "esol-rf-scaffold-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(), "sklearn_version": sklearn.__version__,
        "rdkit_version": rdkit.__version__, "numpy_version": np.__version__,
        "target": TARGET, "unit": "log10(S / [mol/L])", "dataset": "Delaney ESOL",
        "data_sha256": data_hash, "raw_rows": len(raw), "excluded_rows": len(excluded),
        "unique_molecules": len(records),
        "replicate_rows_aggregated": sum(len(e) - 1 for e in grouped.values()),
        "feature_names": FEATURE_NAMES, "parameters": parameters,
        "split": {"method": "deterministic largest-first Bemis-Murcko scaffold 80/20",
                  "train": len(train_idx), "test": len(test_idx),
                  "train_scaffolds": len(train_scaffolds), "test_scaffolds": len(test_scaffolds),
                  "acyclic_policy": "one shared empty scaffold group",
                  "canonical_overlap": 0, "scaffold_overlap": 0},
        "metrics": metrics,
        "evaluation_note": "One fixed internal scaffold holdout; no hyperparameter selection; served model fits train only. Not external validation.",
        "limitations": ["未建模 pH、温度、晶型及具体实验条件。", "描述符相似不代表预测可靠；未校准逐分子置信区间。",
                        "训练数据内分子的演示结果不能作为泛化能力证据。", "二维描述符不区分所有立体异构体。"],
    }
    bundle = {"schema_version": 1, "model": model, "feature_names": FEATURE_NAMES,
              "metadata": metadata, "training_smiles": training_smiles,
              "descriptor_bounds": {f: [float(X[train_idx, i].min()), float(X[train_idx, i].max())]
                                    for i, f in enumerate(FEATURE_NAMES)},
              "allowed_elements": ALLOWED_ELEMENTS}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    reports = output / "reports" / "assignment"
    reports.mkdir(parents=True, exist_ok=True)
    with (output / "model.pkl").open("wb") as f:
        pickle.dump(bundle, f, protocol=pickle.HIGHEST_PROTOCOL)
    metadata["model_sha256"] = hashlib.sha256((output / "model.pkl").read_bytes()).hexdigest()
    save_json(reports / "metrics.json", metadata)
    save_json(reports / "data_audit.json", {"raw_rows": len(raw), "excluded": excluded, "duplicates": duplicate_rows})
    train_ids = set(train_idx)
    save_csv(reports / "split_manifest.csv", [
        {k: r[k] for k in ["canonical_smiles", "scaffold", "logS", "replicate_count"]} |
        {"split": "train" if i in train_ids else "test"} for i, r in enumerate(records)
    ], ["canonical_smiles", "scaffold", "split", "logS", "replicate_count"])
    save_csv(reports / "test_predictions.csv", [
        {"canonical_smiles": records[i]["canonical_smiles"], "observed": float(y[i]), "predicted": float(p)}
        for i, p in zip(test_idx, test_pred)
    ], ["canonical_smiles", "observed", "predicted"])
    print(json.dumps({"rows": len(raw), "unique": len(records), "excluded": len(excluded),
                      "split": metadata["split"], "metrics": metrics}, indent=2))


if __name__ == "__main__":
    main()
