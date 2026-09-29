"""Bounded ETKDGv3/MMFF94 subprocess worker; no quantum calculation is implied."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski, rdMolDescriptors

from .chemistry import parse_smiles

METHOD = "ETKDGv3 + MMFF94"
WARNINGS = [
    "此步骤是 ETKDGv3 构象生成与 MMFF94 分子力场优化，不是量子化学或 DFT 计算。",
    "能量仅用于同一分子的构象比较；不能跨分子比较，也不能作为溶解度、水合自由能或溶液相热力学结果。",
    "构象采样有限；最低采样能量不保证是全局最低能构象。坐标单位为 Å，能量单位为 kcal/mol。",
]


def empty_conformers(status="skipped", reason="standard_mode"):
    return {"status": status, "reason": reason, "method": METHOD, "requested": 16,
            "embedded": 0, "converged": 0, "atoms": [], "bonds": [], "items": [],
            "unconverged_count": 0, "unconverged": [],
            "lowest_energy_id": None, "shape_descriptors": {}, "sdf_url": None,
            "warnings": list(WARNINGS)}


def calculate_conformers(smiles: str, sdf_path: Path | str, num_confs: int = 16):
    """Execute only inside an externally timed child process in production."""
    mol = parse_smiles(smiles)
    if mol.GetNumHeavyAtoms() > 60 or Lipinski.NumRotatableBonds(mol) > 12:
        return empty_conformers("skipped", "too_large_or_flexible")
    if mol.GetNumHeavyAtoms() < 2:
        return empty_conformers("skipped", "too_few_heavy_atoms")
    mol = Chem.AddHs(mol)
    if mol.GetNumAtoms() > 200:
        return empty_conformers("skipped", "too_many_atoms_with_hydrogens")
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        return empty_conformers("skipped", "missing_mmff94_parameters")
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 42
    parameters.numThreads = 1
    parameters.pruneRmsThresh = 0.35
    parameters.maxIterations = 300
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=max(1, min(16, int(num_confs))), params=parameters))
    if not ids:
        return empty_conformers("failed", "embedding_failed")
    optimized = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=500, mmffVariant="MMFF94")
    rows = []
    for conf_id, (status, energy) in zip(ids, optimized):
        if math.isfinite(float(energy)):
            coordinates = [[float(x) for x in point] for point in mol.GetConformer(conf_id).GetPositions()]
            if all(math.isfinite(x) for point in coordinates for x in point):
                rows.append({"id": int(conf_id), "energy_kcal_mol": float(energy), "converged": status == 0, "coordinates": coordinates})
    if not rows:
        return empty_conformers("failed", "optimization_failed")
    unconverged = [{"id": row["id"], "converged": False} for row in rows if not row["converged"]]
    rows = [row for row in rows if row["converged"]]
    if not rows:
        result = empty_conformers("failed", "no_converged_conformers")
        result.update(embedded=len(ids), unconverged_count=len(unconverged), unconverged=unconverged)
        result["warnings"].append("没有构象在迭代上限内达到收敛；未报告最低能构象、形状性质或优化完成。")
        return result
    rows.sort(key=lambda row: row["energy_kcal_mol"])
    minimum = rows[0]["energy_kcal_mol"]
    for row in rows:
        row["relative_energy_kcal_mol"] = row["energy_kcal_mol"] - minimum
    lowest = rows[0]["id"]
    shape_functions = {
        "radius_of_gyration": rdMolDescriptors.CalcRadiusOfGyration,
        "asphericity": rdMolDescriptors.CalcAsphericity,
        "eccentricity": rdMolDescriptors.CalcEccentricity,
        "inertial_shape_factor": rdMolDescriptors.CalcInertialShapeFactor,
        "spherocity_index": rdMolDescriptors.CalcSpherocityIndex,
    }
    shape = {}
    for name, function in shape_functions.items():
        value = float(function(mol, confId=lowest))
        shape[name] = value if math.isfinite(value) else None
    path = Path(sdf_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(path)) as writer:
        for row in rows:
            mol.SetProp("Conformer_ID", str(row["id"]))
            mol.SetProp("MMFF94_energy_kcal_mol", str(row["energy_kcal_mol"]))
            mol.SetProp("Relative_MMFF94_energy_kcal_mol", str(row["relative_energy_kcal_mol"]))
            mol.SetProp("Optimization_converged", str(row["converged"]))
            mol.SetProp("Method", METHOD)
            writer.write(mol, confId=row["id"])
    warnings = list(WARNINGS)
    converged = sum(row["converged"] for row in rows)
    if unconverged:
        warnings.append("部分构象在 500 次迭代内未收敛，已从能量排序、形状性质和 SDF 下载中排除。")
    return {"status": "completed", "reason": None, "method": METHOD, "requested": num_confs,
            "embedded": len(ids), "converged": converged,
            "unconverged_count": len(unconverged), "unconverged": unconverged,
            "atoms": [{"index": atom.GetIdx(), "element": atom.GetSymbol(), "atomic_number": atom.GetAtomicNum()} for atom in mol.GetAtoms()],
            "bonds": [{"begin": bond.GetBeginAtomIdx(), "end": bond.GetEndAtomIdx(), "order": float(bond.GetBondTypeAsDouble())} for bond in mol.GetBonds()],
            "items": rows, "lowest_energy_id": lowest, "shape_descriptors": shape,
            "sdf_url": None, "warnings": warnings}


def worker_main():
    request_path, output_path, sdf_path = map(Path, sys.argv[1:4])
    request = json.loads(request_path.read_text(encoding="utf-8"))
    try:
        result = calculate_conformers(request["smiles"], sdf_path)
    except Exception as exc:
        result = empty_conformers("failed", "calculation_error")
        result["warnings"].append("分子力场步骤未完成：" + type(exc).__name__)
    output_path.write_text(json.dumps(result, ensure_ascii=False, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    worker_main()
