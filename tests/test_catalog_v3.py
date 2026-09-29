"""Real-catalog tests: identity, filter semantics, units and provenance."""
import sqlite3,unittest,tempfile
from contextlib import closing
from pathlib import Path
from app.catalog import Catalog,QuantumCatalog,ROOT

class CatalogTests(unittest.TestCase):
    def test_missing_fingerprint_index_is_service_unavailable(self):
        from app.catalog_api import checked
        from fastapi import HTTPException
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'incomplete.sqlite'
            with closing(sqlite3.connect(path)) as c:
                c.execute('CREATE TABLE molecules(id INTEGER PRIMARY KEY)');c.commit()
            with self.assertRaises(HTTPException) as caught:
                checked(lambda:Catalog(path).search('CCO','similarity'))
            self.assertEqual(caught.exception.status_code,503)
    def test_equivalent_smiles_exact(self):
        left=Catalog().search('OCC','exact');right=Catalog().search('乙醇','exact')
        self.assertEqual([r['id'] for r in left['items']],[r['id'] for r in right['items']])
        self.assertEqual(left['total'],1)
    def test_similarity_self_first(self):
        result=Catalog().search('CCO','similarity',threshold=.9)
        self.assertEqual(result['items'][0]['canonical_smiles'],'CCO')
        self.assertAlmostEqual(result['items'][0]['similarity'],1.)
        self.assertTrue(all(i['similarity']>=.9 for i in result['items']))
    def test_filters_apply_to_same_observation(self):
        result=Catalog().search(source='FreeSolv',target='logS')
        self.assertEqual(result['total'],0)
        result=Catalog().search(source='FreeSolv',target='hydration_free_energy')
        self.assertGreater(result['total'],500)
        self.assertTrue(all('FreeSolv' in r['sources'] for r in result['items']))
    def test_stable_pages(self):
        a=Catalog().search(limit=3);b=Catalog().search(limit=3,offset=3)
        self.assertEqual(a['total'],b['total'])
        self.assertFalse({r['id'] for r in a['items']} & {r['id'] for r in b['items']})
    def test_parameterized_literal_query(self):
        self.assertEqual(Catalog().search("' OR 1=1 --")['total'],0)
        self.assertEqual(Catalog().search('%_') ['total'],0)
    def test_invalid_structure_is_not_empty_success(self):
        with self.assertRaises(ValueError):Catalog().search('not a molecule','similarity')
        with self.assertRaises(ValueError):Catalog().search('CCO',source='invented')
    def test_record_keeps_source_observations(self):
        item=Catalog().search('CCO','exact')['items'][0]
        detail=Catalog().detail(item['id'])
        self.assertGreaterEqual(len(detail['observations']),3)
        self.assertTrue(all(r['source_url'].startswith('https://') for r in detail['observations']))
        self.assertIn('<svg',detail['structure_svg'])
    def test_database_integrity(self):
        for name in ['chemistry.sqlite','quantum.sqlite']:
            with closing(sqlite3.connect(ROOT/'data'/name)) as c:
                self.assertEqual(c.execute('PRAGMA integrity_check').fetchone()[0],'ok')
                self.assertEqual(c.execute('PRAGMA foreign_key_check').fetchall(),[])
    def test_qm9_known_methane_units(self):
        result=QuantumCatalog().search('gdb_1')
        self.assertEqual(result['total'],1)
        row=result['items'][0];self.assertEqual(row['formula'],'CH4')
        self.assertAlmostEqual(row['properties']['homo'],-.3877*27.211386245988)
        self.assertEqual(row['properties']['dipole_moment'],0)
        self.assertEqual(result['units']['gap'],'eV')
        self.assertEqual(row['kind'],'computed_reference')
    def test_qm9_official_bad_geometry_excluded(self):
        self.assertEqual(QuantumCatalog().search('gdb_58')['total'],0)
    def test_qm9_miss_is_not_zero_property(self):
        result=QuantumCatalog().search('nonexistent-name-453')
        self.assertEqual(result['items'],[])
    def test_qm9_formula_and_alias(self):
        self.assertGreaterEqual(QuantumCatalog().search('C2H6O')['total'],2)
        self.assertEqual(QuantumCatalog().search('乙醇')['items'][0]['canonical_smiles'],'CCO')

if __name__=='__main__':unittest.main()
