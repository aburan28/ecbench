#!/usr/bin/env python3
"""Rebuild golden payloads from bundled evidence in an isolated directory."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gp',help='PARI/GP path; GLV4 source proof uses this executable.')
    parser.add_argument('--output',type=Path,default=ROOT/'build'/'reproducibility.json')
    args=parser.parse_args()
    environment=dict(os.environ)
    if args.gp:environment['CURVE_REPORT_GP']=str(Path(args.gp).resolve())
    stages=[]
    started=time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='ecbench-rebuild-') as directory:
        isolated=Path(directory)
        for name in ('ecbench','tools','provenance','schemas'):
            shutil.copytree(ROOT/name,isolated/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        (isolated/'build').mkdir()
        steps=[('prime_vectors.py',[]),('binary_vectors.py',[]),('prime_maps.py',[]),('binary_maps.py',[]),
            ('build_dataset.py',[]),('build_algebra_cases.py',[]),('build_algebra_coverage_cases.py',[]),
            ('build_endomorphism_cases.py',[]),('glv4_cases.py',['--pari-crosscheck']),('build_research_manifest.py',[])]
        for name,options in steps:
            stage_started=time.perf_counter()
            process=subprocess.run([sys.executable,str(isolated/'tools'/name),*options],cwd=isolated,
                env=environment,text=True,capture_output=True)
            if process.returncode:
                raise RuntimeError(f'{name} failed: {process.stderr[-4000:]}')
            stages.append({'tool':name,'elapsed_seconds':time.perf_counter()-stage_started,'status':'pass'})
            print(f'{name}: pass',file=sys.stderr,flush=True)
        files=[]
        for expected in sorted((ROOT/'ecbench'/'data').glob('*.json')):
            actual=isolated/'ecbench'/'data'/expected.name
            if not actual.exists() or actual.read_bytes()!=expected.read_bytes():
                raise RuntimeError(f'Rebuilt payload differs: {expected.name}')
            files.append({'path':'ecbench/data/'+expected.name,
                'sha256':hashlib.sha256(actual.read_bytes()).hexdigest(),'byte_identical':True})
        for name in ('golden_ecbench.json','golden_ecbench.jsonl'):
            if (isolated/name).read_bytes()!=(ROOT/name).read_bytes():raise RuntimeError(f'Rebuilt export differs: {name}')
            files.append({'path':name,'sha256':hashlib.sha256((isolated/name).read_bytes()).hexdigest(),'byte_identical':True})
    result={'status':'pass','scope':'Isolated full payload regeneration from bundled source evidence; no original build intermediates copied',
        'elapsed_seconds':time.perf_counter()-started,'python':sys.version.split()[0],
        'pari_metadata_note':'Byte equality of GLV4 source-proof version metadata requires the same PARI version as the release.',
        'stages':stages,'payloads':files}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'pass','payloads_reproduced':len(files),'output':str(args.output)},sort_keys=True))


if __name__=='__main__':main()
