"""Indexed reference retrieval; experimental and quantum evidence never merged."""
from contextlib import closing
from functools import lru_cache
from pathlib import Path
import json, sqlite3
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
from .chemistry import parse_smiles,canonical_smiles, _NAMES, molecule_svg
from .database import MoleculeDatabase, _CHINESE_NAME_ALIASES, TARGETS

ROOT=Path(__file__).resolve().parents[1]
ALIASES={v.casefold():k for k,v in _NAMES.items()}
ALIASES.update({zh:ALIASES[en.casefold()] for zh,en in _CHINESE_NAME_ALIASES.items()})
GENERATOR=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
QUANTUM_FIELDS=[{'id':'dipole_moment','label':'偶极矩','unit':'D'},
    {'id':'polarizability','label':'各向同性极化率','unit':'a0³'},
    {'id':'homo','label':'HOMO 能量','unit':'eV'}, {'id':'lumo','label':'LUMO 能量','unit':'eV'},
    {'id':'gap','label':'HOMO–LUMO 能隙','unit':'eV'}]
UNITS={r['id']:r['unit'] for r in QUANTUM_FIELDS}
QUANTUM_NOTE='文献 B3LYP/6-31G(2df,p) 计算参考值；不是实验测量，也不是本次运行的量子计算。未命中表示库中没有对应结构，不代表该性质为零。'
CONVERGENCE_CAUTION={f'gdb_{n}' for n in [21725,87037,59827,117523,128113,129053,129152,129158,130535,6620,59818]}

def connect(path):
    if not Path(path).is_file():raise FileNotFoundError('参考数据库尚未就绪。')
    con=sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro',uri=True,timeout=10)
    con.row_factory=sqlite3.Row;con.execute('PRAGMA query_only=ON');return con

def structure(query):
    return canonical_smiles(parse_smiles(ALIASES.get(query.strip().casefold(),query)))

def _filters(source,target):
    clauses=[];params=[]
    # One observation must meet both filters, rather than two unrelated records.
    if source:clauses.append('f.dataset=?');params.append(source)
    if target:clauses.append('f.target=?');params.append(target)
    if clauses:return 'EXISTS(SELECT 1 FROM observations f WHERE f.molecule_id=m.id AND '+' AND '.join(clauses)+')',params
    return '1=1',params

@lru_cache(maxsize=2)
def _fingerprints(db_path,modified_ns):
    with closing(connect(db_path)) as con:
        return tuple((r[0],DataStructs.CreateFromBinaryText(r[1])) for r in con.execute('SELECT molecule_id,morgan2048 FROM fingerprints ORDER BY molecule_id'))

class Catalog:
    def __init__(self,path=None):self.path=Path(path) if path else ROOT/'data/chemistry.sqlite'
    def search(self,query='',mode='text',source=None,target=None,limit=20,offset=0,threshold=.2):
        if mode not in {'text','exact','similarity'}:raise ValueError('未知检索模式。')
        if target and target not in TARGETS:raise ValueError('未知性质。')
        if not 1<=limit<=100 or not 0<=offset<=1000000 or not 0<=threshold<=1:raise ValueError('分页或相似度阈值无效。')
        if not isinstance(query,str) or len(query)>2000:raise ValueError('查询超过长度上限。')
        query=query.strip();where,params=_filters(source,target);similarities={}
        with closing(connect(self.path)) as con:
            if source and not con.execute('SELECT 1 FROM datasets WHERE id=?',(source,)).fetchone():raise ValueError('未知数据来源。')
            if mode=='similarity':
                if not query:raise ValueError('相似结构检索需要 SMILES 或已收录的常见名称。')
                fp=GENERATOR.GetFingerprint(parse_smiles(structure(query)))
                indexed=_fingerprints(str(self.path.resolve()),self.path.stat().st_mtime_ns)
                scores=DataStructs.BulkTanimotoSimilarity(fp,[r[1] for r in indexed])
                eligible={r[0] for r in con.execute('SELECT m.id FROM molecules m WHERE '+where,params)}
                ranked=sorted(((i,float(s)) for (i,_),s in zip(indexed,scores) if i in eligible and s>=threshold),key=lambda x:(-x[1],x[0]))
                total=len(ranked);selected=ranked[offset:offset+limit];similarities=dict(selected)
                rows=[con.execute('SELECT * FROM molecules WHERE id=?',(i,)).fetchone() for i,_ in selected]
            else:
                ordering='m.id';orderparams=[]
                if query and mode=='exact':where+=' AND m.canonical_smiles=?';params.append(structure(query))
                elif query:
                    escaped=query.replace('!','!!').replace('%','!%').replace('_','!_');pattern='%'+escaped+'%'
                    alias=ALIASES.get(query.casefold(),query)
                    where+=" AND (m.name LIKE ? ESCAPE '!' OR m.canonical_smiles LIKE ? ESCAPE '!' OR m.formula LIKE ? ESCAPE '!' OR m.canonical_smiles=? OR EXISTS(SELECT 1 FROM observations s WHERE s.molecule_id=m.id AND s.source_id LIKE ? ESCAPE '!'))"
                    params += [pattern,pattern,pattern,alias,pattern]
                    ordering='CASE WHEN m.canonical_smiles=? THEN 0 WHEN m.name COLLATE NOCASE=? THEN 1 ELSE 2 END,m.id';orderparams=[alias,query]
                elif mode=='exact':raise ValueError('精确结构检索需要输入结构。')
                total=con.execute('SELECT COUNT(*) FROM molecules m WHERE '+where,params).fetchone()[0]
                rows=con.execute('SELECT m.* FROM molecules m WHERE '+where+' ORDER BY '+ordering+' LIMIT ? OFFSET ?',params+orderparams+[limit,offset]).fetchall()
            items=[]
            for r in rows:
                item={k:r[k] for k in ['id','canonical_smiles','name','formula','mw']}
                obs=con.execute('SELECT DISTINCT target,dataset FROM observations WHERE molecule_id=?',(r['id'],)).fetchall()
                item.update(targets=sorted({o[0] for o in obs}),sources=sorted({o[1] for o in obs}),kind='experimental_reference')
                if mode=='similarity':item['similarity']=similarities[r['id']]
                items.append(item)
        return dict(total=total,items=items,limit=limit,offset=offset,mode=mode,similarity_method='Morgan radius 2 / 2048 bits / Tanimoto / chirality excluded' if mode=='similarity' else None)
    def detail(self,molecule_id):
        with closing(connect(self.path)) as con:r=con.execute('SELECT * FROM molecules WHERE id=?',(molecule_id,)).fetchone()
        if r is None:return None
        observations=MoleculeDatabase(self.path).lookup(r['canonical_smiles'])
        summary={}
        for target in sorted({o['target'] for o in observations}):
            vals=[o['value'] for o in observations if o['target']==target]
            summary[target]={'count':len(vals),'min':min(vals),'max':max(vals),'range':max(vals)-min(vals)}
        return {**dict(r),'structure_svg':molecule_svg(parse_smiles(r['canonical_smiles'])),'observations':observations,'summary':summary,'note':'保留各来源原值；数值差异可能来自条件、来源或测量方法，不自动判定某条记录错误。'}

class QuantumCatalog:
    def __init__(self,path=None):self.path=Path(path) if path else ROOT/'data/quantum.sqlite'
    def stats(self):
        with closing(connect(self.path)) as con:
            records=con.execute('SELECT COUNT(*) FROM quantum').fetchone()[0]
            unique=con.execute('SELECT COUNT(DISTINCT canonical_smiles) FROM quantum').fetchone()[0]
        return dict(records=records,unique_structures=unique,dataset='QM9',method='B3LYP/6-31G(2df,p)',kind='computed_reference',fields=QUANTUM_FIELDS,source_url='https://doi.org/10.1038/sdata.2014.22',note=QUANTUM_NOTE)
    def search(self,query='',limit=20,offset=0):
        if not 1<=limit<=100 or not 0<=offset<=1000000 or len(query)>2000:raise ValueError('分页或查询无效。')
        query=query.strip();params=[];where='';canonical=ALIASES.get(query.casefold())
        if query:
            if canonical is None:
                try:canonical=structure(query)
                except ValueError:canonical=None
            # Formula, source ID, and normalized structure exact match. No
            # broad SMILES substring matching of this 130k reference corpus.
            where=' WHERE id=? OR formula=? OR canonical_smiles=?';params=[query,query,canonical or '']
        with closing(connect(self.path)) as con:
            total=con.execute('SELECT COUNT(*) FROM quantum'+where,params).fetchone()[0]
            rows=con.execute('SELECT * FROM quantum'+where+' ORDER BY rowid LIMIT ? OFFSET ?',params+[limit,offset]).fetchall()
        items=[{**{k:r[k] for k in ['id','canonical_smiles','formula','mw']},'properties':{k:r[k] for k in UNITS},'kind':'computed_reference','quality_note':'原作者标注此条目存在收敛困难、放宽阈值或低频鞍点，谨慎使用。' if r['id'] in CONVERGENCE_CAUTION else None} for r in rows]
        return dict(total=total,items=items,limit=limit,offset=offset,method='B3LYP/6-31G(2df,p)',units=UNITS,note=QUANTUM_NOTE)
