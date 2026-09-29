"""Check stored labels, provenance, SQL safety, and observable database behavior."""
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

from app.database import MoleculeDatabase
from app.chemistry import canonical_smiles, parse_smiles

ROOT = Path(__file__).resolve().parents[1]


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


class DatabaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = MoleculeDatabase()
        cls.sources = json.loads((ROOT / "data/sources.json").read_text(encoding="utf-8"))["datasets"]
        cls.audit = json.loads((ROOT / "reports/data_audit.json").read_text(encoding="utf-8"))
        cls.export = read_csv(ROOT / "data/curated_observations.csv")

    def test_stats_match_all_rows_and_audit(self):
        stats = self.db.stats()
        self.assertEqual(stats["dataset_count"], 4)
        self.assertEqual(stats["molecules"], self.audit["molecules"])
        self.assertEqual(stats["observations"], len(self.export))
        self.assertEqual(stats["observations"], self.audit["observations"])
        self.assertEqual(sum(row["accepted_count"] for row in stats["datasets"]), stats["observations"])
        self.assertEqual(sum(stats["target_counts"].values()), stats["observations"])
        self.assertEqual(set(stats["targets"]), {"logS", "logD74", "hydration_free_energy"})
        self.assertEqual(sum(row["raw_count"] for row in stats["datasets"]), 15952)
        for row in self.audit["datasets"]:
            self.assertEqual(row["raw_count"], row["accepted_count"] + row["excluded_count"])
            self.assertEqual(row["excluded_count"], len(row["excluded"]))
            self.assertEqual(row["excluded_count"], sum(row["exclusion_reasons"].values()))

    def test_all_raw_snapshot_hashes_and_counts_match_provenance(self):
        for source in self.sources:
            with self.subTest(dataset=source["id"]):
                path = ROOT / "data" / source["local_file"]
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), source["sha256"])
                self.assertEqual(len(read_csv(path)), source["raw_count"])
                self.assertTrue(source["license_note"])
                self.assertTrue(source["source_url"].startswith("https://"))

    def test_export_preserves_experimental_columns_and_units(self):
        mapping = {
            "ESOL": ("delaney-processed.csv", "measured log solubility in mols per litre", None, "logS", "log10(mol/L)"),
            "AqSolDB": ("aqsoldb.csv", "Solubility", "ID", "logS", "log10(mol/L)"),
            "Lipophilicity": ("Lipophilicity.csv", "exp", "CMPD_CHEMBLID", "logD74", "log10(D), pH 7.4"),
            "FreeSolv": ("SAMPL.csv", "expt", None, "hydration_free_energy", "kcal/mol"),
        }
        originals = {}
        for dataset, (filename, label, id_column, _, _) in mapping.items():
            originals[dataset] = {
                row[id_column] if id_column else f"csv-row-{number}": row
                for number, row in enumerate(read_csv(ROOT / "data/raw" / filename), start=2)
            }
        keys = set()
        for row in self.export:
            dataset = row["dataset"]
            _, label, _, target, unit = mapping[dataset]
            self.assertEqual(row["target"], target)
            self.assertEqual(row["unit"], unit)
            self.assertEqual(float(row["value"]), float(originals[dataset][row["source_id"]][label]))
            key = (dataset, row["source_id"], target)
            self.assertNotIn(key, keys)
            keys.add(key)
            if dataset == "AqSolDB":
                self.assertIn("Group=" + originals[dataset][row["source_id"]]["Group"], row["quality"])

    def test_search_is_paged_and_returns_supported_contract(self):
        first = self.db.search(limit=10, offset=0)
        second = self.db.search(limit=10, offset=10)
        self.assertEqual(first["total"], self.db.stats()["molecules"])
        self.assertEqual(len(first["items"]), 10)
        self.assertFalse({row["id"] for row in first["items"]} & {row["id"] for row in second["items"]})
        for row in first["items"]:
            self.assertEqual(set(row), {"id", "name", "canonical_smiles", "formula", "mw", "targets"})
            self.assertGreater(row["mw"], 0)
            self.assertIsInstance(row["targets"], list)
            self.assertTrue(row["targets"])

    def test_target_filter_and_chembl_identifier_search(self):
        for target in ("logS", "logD74", "hydration_free_energy"):
            result = self.db.search(target=target)
            self.assertEqual(result["total"], self.db.stats()["targets"][target]["molecules"])
            self.assertTrue(all(target in row["targets"] for row in result["items"]))
        result = self.db.search(query="CHEMBL596271", target="logD74")
        self.assertEqual(result["total"], 1)
        measurements = self.db.lookup(result["items"][0]["canonical_smiles"])
        self.assertTrue(any(row["source_id"] == "CHEMBL596271" and row["value"] == 3.54 for row in measurements))

    def test_ethanol_lookup_exposes_separate_sources_and_property_units(self):
        rows = self.db.lookup("CCO")
        self.assertTrue(rows)
        self.assertTrue(any(row["target"] == "logS" for row in rows))
        self.assertTrue(any(row["target"] == "hydration_free_energy" for row in rows))
        for row in rows:
            self.assertEqual(set(row), {"dataset", "source_id", "target", "value", "unit", "quality", "source_url"})
            self.assertTrue(row["source_url"].startswith("https://"))

    def test_known_chinese_aliases_find_the_exact_molecule(self):
        for alias, smiles in [("乙醇", "CCO"), ("咖啡因", "Cn1c(=O)c2c(ncn2C)n(C)c1=O")]:
            with self.subTest(alias=alias):
                canonical = canonical_smiles(parse_smiles(smiles))
                result = self.db.search(query=alias)
                self.assertEqual(result["total"], 1)
                self.assertEqual([item["canonical_smiles"] for item in result["items"]], [canonical])
                self.assertEqual(self.db.search(query="  " + alias + "  ")["items"], result["items"])

    def test_chinese_alias_does_not_become_a_smiles_substring_or_ignore_filters(self):
        exact = self.db.search(query="乙醇")
        substring = self.db.search(query="CCO")
        self.assertEqual(exact["total"], 1)
        self.assertGreater(substring["total"], exact["total"])
        self.assertEqual(self.db.search(query="乙醇", target="hydration_free_energy")["total"], 1)
        self.assertEqual(self.db.search(query="乙醇", target="logD74")["total"], 0)
        self.assertEqual(self.db.search(query="乙醇", offset=1)["items"], [])
        self.assertEqual(self.db.search(query="乙醇溶液")["total"], 0)
        self.assertEqual(self.db.search(query="咖啡")["total"], 0)

    def test_exact_english_name_and_smiles_rank_before_partial_matches(self):
        for query in ("ethanol", "ETHANOL", "CCO"):
            with self.subTest(query=query):
                first = self.db.search(query=query, limit=1)
                second = self.db.search(query=query, limit=1, offset=1)
                self.assertGreater(first["total"], 1)
                self.assertEqual(first["items"][0]["canonical_smiles"], "CCO")
                self.assertEqual(first["total"], second["total"])
                self.assertNotEqual(first["items"][0]["id"], second["items"][0]["id"])

    def test_parameterized_queries_reject_injection_and_literal_wildcards(self):
        before = self.db.stats()
        self.assertEqual(self.db.search("' OR 1=1 --")["total"], 0)
        self.assertEqual(self.db.lookup("' OR 1=1 --"), [])
        self.assertLess(self.db.search("%")["total"], before["molecules"])
        self.assertEqual(self.db.stats(), before)

    def test_query_bounds_are_explicit(self):
        for kwargs in [{"limit": 0}, {"limit": 101}, {"limit": True}, {"offset": -1},
                       {"offset": 1000001}, {"query": None}, {"query": "x" * 2001}, {"target": "LogP"}]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.db.search(**kwargs)

    def test_connections_are_read_only_and_schema_has_no_orphan_observations(self):
        connection = self.db._connect()
        try:
            with self.assertRaises(sqlite3.OperationalError):
                connection.execute("CREATE TABLE accidental_mutation(id INTEGER)")
            self.assertEqual(connection.execute("PRAGMA integrity_check").fetchone()[0], "ok")
            self.assertEqual(connection.execute("PRAGMA foreign_key_check").fetchall(), [])
            self.assertEqual(connection.execute("SELECT COUNT(*)-COUNT(DISTINCT canonical_smiles) FROM molecules").fetchone()[0], 0)
        finally:
            connection.close()

    def test_absent_database_fails_without_creating_an_empty_file(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "missing.sqlite"
            with self.assertRaises(FileNotFoundError):
                MoleculeDatabase(path)
            self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main()
