"""One shared, deterministic RDKit feature pipeline for training and serving."""
from __future__ import annotations

import math
from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

FEATURE_NAMES = ["MW", "LogP", "HBD", "HBA", "TPSA", "RotatableBonds",
                 "RingCount", "AromaticFraction", "HeavyAtoms"]
ALLOWED_ELEMENTS = ["H", "B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"]
MAX_SMILES_LENGTH = 2000
MAX_ATOMS = 200


class MoleculeError(ValueError):
    def __init__(self, code: str, message: str, valid: bool = False):
        super().__init__(message)
        self.code, self.message, self.valid = code, message, valid


def parse_smiles(text: str) -> Chem.Mol:
    if not isinstance(text, str) or not text.strip():
        raise MoleculeError("EMPTY_SMILES", "请输入非空的 SMILES 字符串。")
    text = text.strip()
    if len(text) > MAX_SMILES_LENGTH:
        raise MoleculeError("SMILES_TOO_LONG", "SMILES 最多为 2000 个字符。")
    params = Chem.SmilesParserParams()
    params.parseName = False
    params.allowCXSMILES = False
    with rdBase.BlockLogs():
        try:
            mol = Chem.MolFromSmiles(text, params)
        except (ValueError, RuntimeError):
            mol = None
    if mol is None or mol.GetNumAtoms() == 0:
        raise MoleculeError("INVALID_SMILES", "无法解析此 SMILES，请检查原子、括号和环编号。")
    if mol.GetNumAtoms() > MAX_ATOMS:
        raise MoleculeError("MOLECULE_TOO_LARGE", "本演示最多支持 200 个原子的分子。", True)
    if len(Chem.GetMolFrags(mol)) != 1:
        raise MoleculeError("MULTIPLE_FRAGMENTS", "请输入单一分子；本模型不支持盐或多组分混合物。", True)
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise MoleculeError("WILDCARD_ATOM", "含通配原子的结构无法用于此性质模型。", True)
    unsupported = sorted({a.GetSymbol() for a in mol.GetAtoms()} - set(ALLOWED_ELEMENTS))
    if unsupported:
        raise MoleculeError("UNSUPPORTED_ELEMENTS", "模型不支持这些元素：" + ", ".join(unsupported), True)
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise MoleculeError("RADICAL_MOLECULE", "本模型未验证自由基结构，请输入闭壳层分子。", True)
    # Atom maps describe a reaction, not a chemical identity.
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    # Use the same ordering for equivalent SMILES in every code path.
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True), params)


def canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def identity_smiles(mol: Chem.Mol) -> str:
    """Connectivity-level grouping prevents stereoisomer leakage between splits."""
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def calculate_descriptors(mol: Chem.Mol) -> dict:
    heavy = mol.GetNumHeavyAtoms()
    result = {
        "MW": float(Descriptors.MolWt(mol)),
        "LogP": float(Crippen.MolLogP(mol)),
        "HBD": int(Lipinski.NumHDonors(mol)),
        "HBA": int(Lipinski.NumHAcceptors(mol)),
        "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
        "RotatableBonds": int(Lipinski.NumRotatableBonds(mol)),
        "RingCount": int(rdMolDescriptors.CalcNumRings(mol)),
        "AromaticFraction": float(sum(a.GetIsAromatic() for a in mol.GetAtoms()) / max(heavy, 1)),
        "HeavyAtoms": int(heavy),
    }
    if not all(math.isfinite(v) for v in result.values()):
        raise MoleculeError("NONFINITE_DESCRIPTOR", "此结构无法生成有限的分子描述符。", True)
    return result


def molecule_formula(mol: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(mol)


def molecule_svg(mol: Chem.Mol) -> str:
    drawer = rdMolDraw2D.MolDraw2DSVG(520, 300)
    drawer.drawOptions().clearBackground = False
    drawer.drawOptions().padding = 0.12
    drawer.drawOptions().bondLineWidth = 2.0
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, Chem.Mol(mol))
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg[svg.index("<svg"):]


_NAMES = {
    canonical_smiles(parse_smiles(s)): name for s, name in [
        ("CC(=O)OC1=CC=CC=C1C(=O)O", "Aspirin"),
        ("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "Caffeine"),
        ("CCO", "Ethanol"),
        ("CC(=O)Nc1ccc(O)cc1", "Paracetamol"),
        ("c1ccccc1", "Benzene"),
    ]
}


def known_name(mol: Chem.Mol) -> str | None:
    """Small explicit offline lookup; never invent names for arbitrary inputs."""
    return _NAMES.get(canonical_smiles(mol))
