"""Package source, data and trained models; verify CRC, hashes and extracted API.

Uses the current installed Python dependencies for the extracted smoke test.
The virtual environment is deliberately excluded from the portable archive.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
IGNORED = {'.venv','__pycache__','.pytest_cache','node_modules','.git','workspace_artifacts','training_cache','.training_cache','cache'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT.parent/'Molecule_Studio_V3.zip')
    parser.add_argument('--scratch',type=Path,default=Path(tempfile.gettempdir()))
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_relative_to(ROOT):
        raise SystemExit('Place the archive outside the project directory.')
    files = [path for path in sorted(ROOT.rglob('*')) if path.is_file()
             and not IGNORED.intersection(path.relative_to(ROOT).parts)
             and path.suffix not in {'.pyc','.pyo','.log','.npy','.npz'} and path.name != 'PACKAGE_MANIFEST.json'
             and not path.name.startswith('workspace.sqlite') and '.building' not in path.name]
    manifest = {'version':'3.0.0', 'hash':'sha256',
                'note':'Manifest excludes itself; environments, caches and personal analysis history excluded.',
                'files':{path.relative_to(ROOT).as_posix():digest(path) for path in files}}
    manifest_path = ROOT/'PACKAGE_MANIFEST.json'
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    files.append(manifest_path)
    output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for path in files:
            archive.write(path,Path('Molecule_Studio_V3')/path.relative_to(ROOT))
    args.scratch.mkdir(parents=True,exist_ok=True)
    extracted_parent = Path(tempfile.mkdtemp(prefix='v3-release-',dir=args.scratch))
    with zipfile.ZipFile(output) as archive:
        assert archive.testzip() is None, 'ZIP CRC failure'
        archive.extractall(extracted_parent)
    extracted = extracted_parent/'Molecule_Studio_V3'
    for name, expected in manifest['files'].items():
        assert digest(extracted/name) == expected, name
    smoke_code = '''
import json,time
from pathlib import Path
from fastapi.testclient import TestClient
import app
from app.main import create_app
assert Path(app.__file__).resolve().is_relative_to(Path.cwd())
with TestClient(create_app()) as client:
    health = client.get('/health').json()
    assert health['status'] == 'ok'
    assert client.get('/database/stats').json()['observations'] == 14685
    assert client.get('/quantum/stats').json()['records'] == 130831
    assert client.get('/model-info').json()['version'] == '3.0.0'
    single = client.post('/predict',json={'smiles':'CCO'}).json()
    assert set(single['predictions']) == {'logS','logD74','hydration_free_energy'}
    assert all(p['value'] is not None for p in single['predictions'].values())
    assert len(single['measurements']) == 3
    assert client.post('/predict',json={'smiles':'C1('}).status_code == 422
    results = client.post('/predict/batch',json={'smiles':['CCO','C1(','O'+'CCO'*40]}).json()['items']
    assert 'error' in results[1]
    assert all(p['value'] is None and p['status']=='out_of_domain' for p in results[2]['result']['predictions'].values())
    assert client.get('/database/search',params={'q':'乙醇'}).json()['items'][0]['canonical_smiles'] == 'CCO'
    assert client.get('/database/advanced',params={'q':'OCC','mode':'exact'}).json()['items'][0]['canonical_smiles'] == 'CCO'
    assert client.get('/quantum/search',params={'q':'gdb_1'}).json()['items'][0]['formula'] == 'CH4'
    compatibility=client.post('/assignment/predict',json={'smiles':'CC(=O)Oc1ccccc1C(=O)O'})
    assert compatibility.status_code==200
    assert abs(compatibility.json()['logS_prediction']+1.9999770998376616)<1e-9
    job=client.post('/analyses',json={'smiles':'CCO','mode':'deep'}).json()
    deadline=time.monotonic()+35
    while time.monotonic()<deadline:
        job=client.get('/analyses/'+job['id']).json()
        if job['status'] in {'completed','failed','cancelled'}:break
        time.sleep(.1)
    assert job['status']=='completed',job
    assert job['result']['conformers']['status']=='completed',job['result']['conformers']
    assert client.get(job['result']['conformers']['sdf_url']).status_code==200
    assert client.get('/').status_code == 200
print(json.dumps({'status':'passed','health':health,'source':str(Path(app.__file__).resolve())}))
'''
    environment = dict(os.environ)
    environment.update(OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',PYTHONUTF8='1')
    runtime_bin = Path(sys.base_prefix)/'Library'/'bin'
    if runtime_bin.is_dir():
        environment['PATH'] = str(runtime_bin)+os.pathsep+environment.get('PATH','')
    result = subprocess.run([sys.executable,'-c',smoke_code],cwd=extracted,env=environment,
                            text=True,capture_output=True,encoding='utf-8',timeout=120)
    report = {'checked_utc':datetime.now(timezone.utc).isoformat(),'archive':output.name,
              'bytes':output.stat().st_size,'sha256':digest(output),'files':len(files),
              'crc_check':'passed','extracted_file_hashes':'all matched',
              'smoke_exit_code':result.returncode,'smoke_stdout':result.stdout,'smoke_stderr':result.stderr,
              'scope':'Extracted source/data/models with existing installed dependencies; not a fresh OS/dependency installation.'}
    report_path = output.with_name(output.stem+'_release.json')
    report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if result.returncode:
        raise SystemExit(result.stderr)
    print(json.dumps({k:report[k] for k in ('archive','bytes','sha256','files','crc_check','extracted_file_hashes','smoke_exit_code')},ensure_ascii=False))


if __name__=='__main__':
    main()
