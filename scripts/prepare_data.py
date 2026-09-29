"""Build an auditable SQLite store from four public physicochemical datasets.

Raw labels are retained separately. This step performs no averaging, model
training, feature selection, or selection using a test set.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import closing
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import sqlite3
import sys
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rdkit.Chem.Scaffolds import MurckoScaffold
from app.chemistry import (MoleculeError, calculate_descriptors, canonical_smiles,
                           identity_smiles, known_name, molecule_formula, parse_smiles)

LOCAL_DATA = Path.home() / ".codex" / "tools" / "chem-ai4s" / "data"
MOLNET = "https://deepchem.readthedocs.io/en/latest/api_reference/moleculenet.html"
SOURCES = [
    dict(id="ESOL", name="Delaney ESOL (MoleculeNet)", filename="delaney-processed.csv", local_key="esol",
         source_url="https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv",
         reference_url=MOLNET, doi="10.1021/ci034243x", raw_count=1128,
         sha256="8c06a76f0c6487d29ab0f903e6a7a7139f189ab3c1178f159c8be8964602f189",
         smiles_column="smiles", label_column="measured log solubility in mols per litre",
         name_column="Compound ID", target="logS", unit="log10(mol/L)",
         license_note="MoleculeNet-distributed benchmark; original Delaney source terms and attribution apply. No standalone dataset license was verified in the CSV."),
    dict(id="AqSolDB", name="AqSolDB author-curated dataset", filename="aqsoldb.csv", local_key=None,
         source_url="https://raw.githubusercontent.com/mcsorkun/AqSolDB/master/results/data_curated.csv",
         reference_url="https://github.com/mcsorkun/AqSolDB", doi="10.1038/s41597-019-0151-1", raw_count=9982,
         sha256="363a42a4e4ea8e6039758146d34efdd95cde2da81c665316cacc1869e0a805ed",
         smiles_column="SMILES", label_column="Solubility", name_column="Name", target="logS", unit="log10(mol/L)",
         license_note="Author GitHub repository includes MIT license (raw/AqSolDB-LICENSE.txt); original source datasets retain their terms. The separate Dataverse data license was not independently retrieved."),
    dict(id="Lipophilicity", name="Lipophilicity / ChEMBL AZ (MoleculeNet)", filename="Lipophilicity.csv", local_key="lipophilicity",
         source_url="https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/Lipophilicity.csv",
         reference_url=MOLNET, doi="10.6019/chembl3301361", raw_count=4200,
         sha256="aed41590cb30609d51d8e08ad3ff06495a76e80e211358801f596b10da69bacd",
         smiles_column="smiles", label_column="exp", name_column="CMPD_CHEMBLID", target="logD74", unit="log10(D), pH 7.4",
         license_note="MoleculeNet distribution of ChEMBL AZ deposited data; preserve ChEMBL/source attribution and terms. This project does not assign a new license to source data."),
    dict(id="FreeSolv", name="FreeSolv / SAMPL (MoleculeNet)", filename="SAMPL.csv", local_key="freesolv",
         source_url="https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/SAMPL.csv",
         reference_url="https://github.com/MobleyLab/FreeSolv", doi="10.1007/s10822-014-9747-x", raw_count=642,
         sha256="ab5895d914ee87cb563bd7b9611e869527bba45bec6b014d34dc495a0f9dcb72",
         smiles_column="smiles", label_column="expt", name_column="iupac", target="hydration_free_energy", unit="kcal/mol",
         license_note="FreeSolv data repository specifies CC-BY-4.0 to the extent the authors can license the underlying data, with original-source caveats; code license is separate."),
]
EXPORT_FIELDS = ["dataset", "source_id", "target", "canonical_smiles", "identity_smiles", "scaffold", "name", "value", "unit", "quality"]


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire(data_dir, local_data, offline):
    raw_dir = data_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    old_path = data_dir / "sources.json"
    prior = {row["id"]: row for row in json.loads(old_path.read_text(encoding="utf-8")).get("datasets", [])} if old_path.exists() else {}
    manifest_path = local_data / "manifest.json"
    local_manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    records = []
    for source in SOURCES:
        record = dict(source)
        destination = raw_dir / source["filename"]
        original = local_manifest.get(source["local_key"], {})
        local_file = local_data / source["filename"]
        if destination.exists():
            previous = prior.get(source["id"], {})
            record["acquisition"] = previous.get("acquisition", "project_raw_snapshot")
            record["retrieved_utc"] = previous.get("retrieved_utc", datetime.fromtimestamp(destination.stat().st_mtime, timezone.utc).isoformat())
            if previous.get("reused_from"):
                record["reused_from"] = previous["reused_from"]
        elif source["local_key"] and local_file.exists():
            if sha256(local_file) != source["sha256"] or original.get("sha256") != source["sha256"]:
                raise ValueError(f"Unverified local dataset: {source['id']}")
            shutil.copyfile(local_file, destination)
            record["acquisition"] = "reused_verified_local_snapshot"
            record["reused_from"] = str(local_file)
            record["retrieved_utc"] = original.get("retrieved_utc")
        else:
            if offline:
                raise FileNotFoundError(f"Missing verified raw data in offline mode: {destination}")
            request = Request(source["source_url"], headers={"User-Agent": "MoleculeLab-data-preparation/2.0"})
            with urlopen(request, timeout=60) as response:
                content = response.read(20_000_001)
            if len(content) > 20_000_000 or hashlib.sha256(content).hexdigest() != source["sha256"]:
                raise ValueError(f"Downloaded bytes differ from the pinned snapshot: {source['id']}")
            destination.write_bytes(content)
            record["acquisition"] = "downloaded_from_recorded_source"
            record["retrieved_utc"] = utc_now()
        actual_hash = sha256(destination)
        if actual_hash != source["sha256"]:
            raise ValueError(f"Raw file hash differs from pinned source: {source['id']}")
        record.update(local_file=f"raw/{source['filename']}", verified_utc=utc_now(), bytes=destination.stat().st_size)
        records.append(record)
    return records


def quality_label(source, row):
    if source["id"] == "AqSolDB":
        return f"author_curated_experimental;Group={row['Group']};SD={row['SD']};Occurrences={row['Occurrences']}"
    if source["id"] == "Lipophilicity":
        return "experimental_reference;pH=7.4;per_row_conditions_not_in_CSV"
    if source["id"] == "FreeSolv":
        return "experimental_reference;expt_column;per_row_uncertainty_not_in_CSV"
    return "experimental_reference;measured_column;per_row_conditions_not_in_CSV"


def source_identifier(source, row, row_number):
    if source["id"] == "AqSolDB":
        return row["ID"]
    if source["id"] == "Lipophilicity":
        return row["CMPD_CHEMBLID"]
    # These distributed CSVs have names, but no stable unique upstream row ID.
    return f"csv-row-{row_number}"


def prepare(data_dir, reports_dir, sources):
    observations, molecules, audits = [], {}, []
    for source in sources:
        with (data_dir / source["local_file"]).open(encoding="utf-8-sig", newline="") as stream:
            raw_rows = list(csv.DictReader(stream))
        if len(raw_rows) != source["raw_count"]:
            raise ValueError(f"Unexpected row count for {source['id']}: {len(raw_rows)}")
        rejected = []
        accepted = 0
        for line_number, row in enumerate(raw_rows, start=2):
            source_id = source_identifier(source, row, line_number)
            try:
                value = float(row[source["label_column"]])
                if not math.isfinite(value):
                    raise ValueError("nonfinite_label")
                mol = parse_smiles(row[source["smiles_column"]])
                canonical = canonical_smiles(mol)
                identity = identity_smiles(mol)
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
                name = row[source["name_column"]].strip() or canonical
                observations.append(dict(dataset=source["id"], source_id=source_id, target=source["target"],
                    canonical_smiles=canonical, identity_smiles=identity, scaffold=scaffold, name=name,
                    value=value, unit=source["unit"], quality=quality_label(source, row)))
                if canonical not in molecules:
                    molecules[canonical] = dict(canonical_smiles=canonical, identity_smiles=identity,
                        scaffold=scaffold, name=known_name(mol) or name, formula=molecule_formula(mol),
                        mw=calculate_descriptors(mol)["MW"])
                accepted += 1
            except (MoleculeError, ValueError, KeyError) as exc:
                rejected.append(dict(csv_line=line_number, source_id=source_id,
                    smiles=row.get(source["smiles_column"], ""), reason=getattr(exc, "code", str(exc))))
        source["accepted_count"] = accepted
        source["excluded_count"] = len(rejected)
        audits.append(dict(dataset=source["id"], raw_count=len(raw_rows), accepted_count=accepted,
            excluded_count=len(rejected), exclusion_reasons=dict(Counter(item["reason"] for item in rejected)), excluded=rejected))
        print(f"{source['id']}: raw={len(raw_rows)}, accepted={accepted}, excluded={len(rejected)}", flush=True)
    with (data_dir / "curated_observations.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=EXPORT_FIELDS)
        writer.writeheader()
        writer.writerows(observations)
    db_path = data_dir / "chemistry.sqlite"
    building = db_path.with_suffix(".sqlite.building")
    if building.exists():
        building.unlink()
    with closing(sqlite3.connect(building)) as connection:
        with connection:
            connection.executescript("""
                PRAGMA foreign_keys = ON;
                CREATE TABLE datasets(id TEXT PRIMARY KEY,name TEXT NOT NULL,source_url TEXT NOT NULL,
                    sha256 TEXT NOT NULL,raw_count INTEGER NOT NULL,accepted_count INTEGER NOT NULL,license_note TEXT NOT NULL);
                CREATE TABLE molecules(id INTEGER PRIMARY KEY,canonical_smiles TEXT NOT NULL UNIQUE,
                    identity_smiles TEXT NOT NULL,scaffold TEXT NOT NULL,name TEXT NOT NULL,formula TEXT NOT NULL,mw REAL NOT NULL);
                CREATE TABLE observations(id INTEGER PRIMARY KEY,molecule_id INTEGER NOT NULL REFERENCES molecules(id),
                    dataset TEXT NOT NULL REFERENCES datasets(id),source_id TEXT NOT NULL,target TEXT NOT NULL,
                    value REAL NOT NULL,unit TEXT NOT NULL,quality TEXT NOT NULL,UNIQUE(dataset,source_id,target));
                CREATE INDEX idx_molecule_identity ON molecules(identity_smiles);
                CREATE INDEX idx_observation_molecule ON observations(molecule_id);
                CREATE INDEX idx_observation_target ON observations(target,molecule_id);
                CREATE INDEX idx_observation_source ON observations(source_id);
            """)
            connection.executemany("INSERT INTO datasets VALUES(?,?,?,?,?,?,?)", [
                tuple(source[key] for key in ("id","name","source_url","sha256","raw_count","accepted_count","license_note")) for source in sources])
            sorted_molecules = sorted(molecules.items())
            molecule_ids = {canonical: idx for idx, (canonical, _) in enumerate(sorted_molecules, start=1)}
            connection.executemany("INSERT INTO molecules VALUES(?,?,?,?,?,?,?)", [
                (molecule_ids[canonical], canonical, row["identity_smiles"], row["scaffold"], row["name"], row["formula"], row["mw"])
                for canonical, row in sorted_molecules])
            connection.executemany("INSERT INTO observations(molecule_id,dataset,source_id,target,value,unit,quality) VALUES(?,?,?,?,?,?,?)", [
                (molecule_ids[row["canonical_smiles"]], row["dataset"], row["source_id"], row["target"], row["value"], row["unit"], row["quality"])
                for row in observations])
            assert connection.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
            assert not connection.execute("PRAGMA foreign_key_check").fetchall()
    building.replace(db_path)
    audit = dict(created_utc=utc_now(), schema_version=1, raw_count=sum(s["raw_count"] for s in sources),
        observations=len(observations), molecules=len(molecules), identity_count=len({r["identity_smiles"] for r in observations}),
        targets=dict(Counter(r["target"] for r in observations)), datasets=audits,
        note="No averaging at database stage. Original row labels remain separate; AqSolDB already includes literature sources overlapping ESOL.",
        database_sha256=sha256(db_path), curated_csv_sha256=sha256(data_dir / "curated_observations.csv"))
    reports_dir.mkdir(parents=True, exist_ok=True)
    write_json(reports_dir / "data_audit.json", audit)
    write_json(data_dir / "sources.json", dict(schema_version=1, prepared_utc=utc_now(), datasets=sources))
    print(json.dumps({key: audit[key] for key in ("raw_count","observations","molecules","identity_count","targets")}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--reports-dir", type=Path, default=ROOT / "reports")
    parser.add_argument("--local-data", type=Path, default=LOCAL_DATA)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    sources = acquire(args.data_dir, args.local_data, args.offline)
    prepare(args.data_dir, args.reports_dir, sources)


if __name__ == "__main__":
    main()
