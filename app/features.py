"""Frozen V3 two-dimensional RDKit features; no fitted preprocessing here.

The first 28 positions preserve the V2 representation for explicit controls.
Descriptor exceptions, nonfinite values and float32 overflow become NaN.
Each estimator fits its own median imputer using only its fitting partition.
QED is omitted (composite score), Ipc is omitted (unbounded magnitudes; AvgIpc
is retained), and radical count is omitted (radicals are rejected by parser).
"""
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator
from .chemistry import calculate_descriptors

LEGACY_FEATURE_NAMES = ["MW", "LogP", "HBD", "HBA", "TPSA", "RotatableBonds", "RingCount",
    "AromaticFraction", "HeavyAtoms", "FractionCSP3", "MolMR", "Heteroatoms", "FormalCharge",
    "AromaticRings", "AliphaticRings", "LabuteASA", "ExactMW", "BertzCT", "BalabanJ",
    "Kappa1", "Kappa2", "Kappa3", "Chi0v", "Chi1v", "NCount", "OCount", "HalogenCount", "SaturatedRings"]

# Explicit order, not runtime registry iteration: package upgrades cannot silently
# add, remove or reorder model inputs.
RDKIT_EXTRA_NAMES = """
MaxAbsEStateIndex MaxEStateIndex MinAbsEStateIndex MinEStateIndex SPS HeavyAtomMolWt
NumValenceElectrons MaxPartialCharge MinPartialCharge MaxAbsPartialCharge MinAbsPartialCharge
FpDensityMorgan1 FpDensityMorgan2 FpDensityMorgan3 BCUT2D_MWHI BCUT2D_MWLOW BCUT2D_CHGHI
BCUT2D_CHGLO BCUT2D_LOGPHI BCUT2D_LOGPLOW BCUT2D_MRHI BCUT2D_MRLOW AvgIpc Chi0 Chi0n
Chi1 Chi1n Chi2n Chi2v Chi3n Chi3v Chi4n Chi4v HallKierAlpha
PEOE_VSA1 PEOE_VSA10 PEOE_VSA11 PEOE_VSA12 PEOE_VSA13 PEOE_VSA14 PEOE_VSA2 PEOE_VSA3
PEOE_VSA4 PEOE_VSA5 PEOE_VSA6 PEOE_VSA7 PEOE_VSA8 PEOE_VSA9 SMR_VSA1 SMR_VSA10 SMR_VSA2
SMR_VSA3 SMR_VSA4 SMR_VSA5 SMR_VSA6 SMR_VSA7 SMR_VSA8 SMR_VSA9 SlogP_VSA1 SlogP_VSA10
SlogP_VSA11 SlogP_VSA12 SlogP_VSA2 SlogP_VSA3 SlogP_VSA4 SlogP_VSA5 SlogP_VSA6 SlogP_VSA7
SlogP_VSA8 SlogP_VSA9 EState_VSA1 EState_VSA10 EState_VSA11 EState_VSA2 EState_VSA3 EState_VSA4
EState_VSA5 EState_VSA6 EState_VSA7 EState_VSA8 EState_VSA9 VSA_EState1 VSA_EState10 VSA_EState2
VSA_EState3 VSA_EState4 VSA_EState5 VSA_EState6 VSA_EState7 VSA_EState8 VSA_EState9 NHOHCount
NOCount NumAliphaticCarbocycles NumAliphaticHeterocycles NumAmideBonds NumAromaticCarbocycles
NumAromaticHeterocycles NumAtomStereoCenters NumBridgeheadAtoms NumHeterocycles
NumSaturatedCarbocycles NumSaturatedHeterocycles NumSpiroAtoms NumUnspecifiedAtomStereoCenters Phi
fr_Al_COO fr_Al_OH fr_Al_OH_noTert fr_ArN fr_Ar_COO fr_Ar_N fr_Ar_NH fr_Ar_OH fr_COO fr_COO2
fr_C_O fr_C_O_noCOO fr_C_S fr_HOCCN fr_Imine fr_NH0 fr_NH1 fr_NH2 fr_N_O fr_Ndealkylation1
fr_Ndealkylation2 fr_Nhpyrrole fr_SH fr_aldehyde fr_alkyl_carbamate fr_alkyl_halide fr_allylic_oxid
fr_amide fr_amidine fr_aniline fr_aryl_methyl fr_azide fr_azo fr_barbitur fr_benzene fr_benzodiazepine
fr_bicyclic fr_diazo fr_dihydropyridine fr_epoxide fr_ester fr_ether fr_furan fr_guanido fr_halogen
fr_hdrzine fr_hdrzone fr_imidazole fr_imide fr_isocyan fr_isothiocyan fr_ketone fr_ketone_Topliss
fr_lactam fr_lactone fr_methoxy fr_morpholine fr_nitrile fr_nitro fr_nitro_arom fr_nitro_arom_nonortho
fr_nitroso fr_oxazole fr_oxime fr_para_hydroxylation fr_phenol fr_phenol_noOrthoHbond fr_phos_acid
fr_phos_ester fr_piperdine fr_piperzine fr_priamide fr_prisulfonamd fr_pyridine fr_quatN fr_sulfide
fr_sulfonamd fr_sulfone fr_term_acetylene fr_tetrazole fr_thiazole fr_thiocyan fr_thiophene
fr_unbrch_alkane fr_urea
""".split()
FEATURE_NAMES = LEGACY_FEATURE_NAMES + RDKIT_EXTRA_NAMES
FEATURE_VERSION = "rdkit-2d219-morgan2048-radius2-v3"
MISSING_VALUE_POLICY = "Per-estimator fit-only median; all-missing fitting column retained and filled with zero."
_EXTRA_FUNCTIONS = [(name, getattr(Descriptors, name)) for name in RDKIT_EXTRA_NAMES]
_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)


def fingerprint(mol):
    return _GENERATOR.GetFingerprint(mol)


def _safe_descriptor(function, mol):
    try:
        value = float(function(mol))
        return value if np.isfinite(value) and abs(value) < float(np.finfo(np.float32).max) else float("nan")
    except Exception:
        return float("nan")


def all_descriptors(mol):
    result = calculate_descriptors(mol)
    functions = {
        "FractionCSP3": rdMolDescriptors.CalcFractionCSP3, "MolMR": Descriptors.MolMR,
        "Heteroatoms": rdMolDescriptors.CalcNumHeteroatoms, "FormalCharge": Chem.GetFormalCharge,
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings,
        "AliphaticRings": rdMolDescriptors.CalcNumAliphaticRings, "LabuteASA": Descriptors.LabuteASA,
        "ExactMW": Descriptors.ExactMolWt, "BertzCT": Descriptors.BertzCT, "BalabanJ": Descriptors.BalabanJ,
        "Kappa1": Descriptors.Kappa1, "Kappa2": Descriptors.Kappa2, "Kappa3": Descriptors.Kappa3,
        "Chi0v": Descriptors.Chi0v, "Chi1v": Descriptors.Chi1v,
        "NCount": lambda m: sum(a.GetAtomicNum() == 7 for a in m.GetAtoms()),
        "OCount": lambda m: sum(a.GetAtomicNum() == 8 for a in m.GetAtoms()),
        "HalogenCount": lambda m: sum(a.GetAtomicNum() in {9, 17, 35, 53} for a in m.GetAtoms()),
        "SaturatedRings": rdMolDescriptors.CalcNumSaturatedRings,
    }
    result.update({name: _safe_descriptor(function, mol) for name, function in functions.items()})
    result.update({name: _safe_descriptor(function, mol) for name, function in _EXTRA_FUNCTIONS})
    return result


def feature_vector(mol, mode="descriptors"):
    descriptors = all_descriptors(mol)
    values = np.asarray([descriptors[n] for n in FEATURE_NAMES], dtype=np.float32)
    if mode == "combined":
        bits = np.zeros(2048, dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fingerprint(mol), bits)
        return np.concatenate([values, bits])
    if mode != "descriptors":
        raise ValueError("Unknown feature mode")
    return values


def model_feature_names(mode="descriptors"):
    if mode not in {"descriptors", "combined"}:
        raise ValueError("Unknown feature mode")
    return FEATURE_NAMES + ([f"Morgan_{i}" for i in range(2048)] if mode == "combined" else [])
