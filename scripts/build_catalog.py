"""Build V3 indexed catalog and separate QM9 reference database, entirely offline.

Input snapshots and exclusions are included in data/raw; no labels here are
added to the three physicochemical prediction training targets.
"""
from pathlib import Path
import csv, gzip, hashlib, json, math, re, sqlite3, sys, time
from contextlib import closing
from datetime import datetime, timezone
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rdkit import DataStructs, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator
from app.chemistry import parse_smiles, canonical_smiles

HARTREE_EV = 27.211386245988
EXPECTED_QM9_SHA = '3e668f8c34e4bc392a90d417a50a5eed3b64b842a817a633024bdc054c68ccb4'

def build():
    started = time.perf_counter()
    raw = ROOT/'data/raw'
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    with closing(sqlite3.connect(ROOT/'data/chemistry.sqlite')) as con:
        con.execute('CREATE TABLE IF NOT EXISTS fingerprints(molecule_id INTEGER PRIMARY KEY REFERENCES molecules(id), morgan2048 BLOB NOT NULL)')
        con.execute('CREATE INDEX IF NOT EXISTS idx_observation_dataset_target ON observations(dataset,target,molecule_id)')
        con.execute('CREATE INDEX IF NOT EXISTS idx_molecule_formula ON molecules(formula)')
        con.execute('CREATE INDEX IF NOT EXISTS idx_molecule_name ON molecules(name COLLATE NOCASE)')
        for mol_id,smiles in con.execute('SELECT id,canonical_smiles FROM molecules').fetchall():
            mol=parse_smiles(smiles)
            con.execute('INSERT OR REPLACE INTO fingerprints VALUES(?,?)',(mol_id,DataStructs.BitVectToBinaryText(gen.GetFingerprint(mol))))
        con.commit()
        assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
        assert not con.execute('PRAGMA foreign_key_check').fetchall()
        n_fps=con.execute('SELECT COUNT(*) FROM fingerprints').fetchone()[0]
    skip={int(m.group(1)) for line in (raw/'qm9-uncharacterized.txt').read_text().splitlines() if (m:=re.match(r'^\s+(\d+)\s+',line))}
    assert len(skip)==3054, f'Unexpected exclusion list: {len(skip)}'
    with gzip.open(raw/'qm9.csv.gz','rb') as f:
        digest=hashlib.file_digest(f,'sha256').hexdigest()
    assert digest==EXPECTED_QM9_SHA, 'QM9 snapshot hash mismatch'
    dest=ROOT/'data/quantum.sqlite'
    temporary=dest.with_suffix('.building.sqlite')
    if temporary.exists(): temporary.unlink()
    counts={'raw':0,'geometry_excluded':0,'unsupported':0,'accepted':0}
    rejected=[]
    with closing(sqlite3.connect(temporary)) as con:
        con.execute('CREATE TABLE quantum(id TEXT PRIMARY KEY,canonical_smiles TEXT NOT NULL,formula TEXT NOT NULL,mw REAL NOT NULL,dipole_moment REAL NOT NULL,polarizability REAL NOT NULL,homo REAL NOT NULL,lumo REAL NOT NULL,gap REAL NOT NULL)')
        with gzip.open(raw/'qm9.csv.gz','rt',encoding='utf-8',newline='') as f:
            for row in csv.DictReader(f):
                counts['raw']+=1
                source_id=row['mol_id']; idx=int(source_id.split('_')[1])
                if idx in skip: counts['geometry_excluded']+=1; continue
                try:
                    mol=parse_smiles(row['smiles'])
                    numbers=[float(row[key]) for key in ['mu','alpha','homo','lumo','gap']]
                    assert all(math.isfinite(v) for v in numbers)
                    assert abs(numbers[3]-numbers[2]-numbers[4])<0.00021
                    values=(source_id,canonical_smiles(mol),rdMolDescriptors.CalcMolFormula(mol),float(Descriptors.MolWt(mol)),numbers[0],numbers[1],*[v*HARTREE_EV for v in numbers[2:]])
                    con.execute('INSERT INTO quantum VALUES(?,?,?,?,?,?,?,?,?)',values)
                    counts['accepted']+=1
                except Exception as exc:
                    counts['unsupported']+=1;rejected.append({'id':source_id,'reason':getattr(exc,'code',type(exc).__name__)})
                if counts['raw']%20000==0:
                    con.commit();print(f"QM9 {counts['raw']} rows processed",flush=True)
        con.executescript('CREATE INDEX idx_quantum_smiles ON quantum(canonical_smiles); CREATE INDEX idx_quantum_formula ON quantum(formula);')
        counts['unique_structures']=con.execute('SELECT COUNT(DISTINCT canonical_smiles) FROM quantum').fetchone()[0]
        con.commit()
        assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    temporary.replace(dest)
    audit={'version':'3.0.0','created_utc':datetime.now(timezone.utc).isoformat(),'rdkit_version':rdBase.rdkitVersion,
        'experimental_fingerprint_count':n_fps,'qm9':counts,'excluded_other':rejected,
        'qm9_snapshot_sha256':digest,'geometry_exclusions_sha256':hashlib.sha256((raw/'qm9-uncharacterized.txt').read_bytes()).hexdigest(),
        'quantum_database_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'hartree_to_ev':HARTREE_EV,
        'method':'B3LYP/6-31G(2df,p)','kind':'computed_reference','training_use':'none',
        'source_url':'https://doi.org/10.1038/sdata.2014.22','data_url':'https://doi.org/10.6084/m9.figshare.c.978904',
        'note':'QM9 calculations are reference lookups, not experimental measurements, newly performed DFT, or predictions for arbitrary molecules. Official 3054 geometry inconsistencies excluded; remaining parse/finite-value failures reported.',
        'elapsed_seconds':time.perf_counter()-started}
    (ROOT/'reports/catalog_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(audit,ensure_ascii=True),flush=True)

if __name__=='__main__':build()
