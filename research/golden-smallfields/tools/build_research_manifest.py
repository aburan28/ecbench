#!/usr/bin/env python3
"""Stage deterministic research catalogs and hash their final release payloads."""
import hashlib
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'ecbench'/'data'


def build():
    for name in ('algebra_cases.json','endomorphism_cases.json','glv4.json'):
        shutil.copyfile(ROOT/'build'/name,DATA/name)
    coverage_source=ROOT/'build'/'algebra_coverage_cases.json'
    if coverage_source.exists():shutil.copyfile(coverage_source,DATA/coverage_source.name)
    for name in ('algebra_cases.schema.json','algebra_results.schema.json'):
        source=ROOT/'schemas'/name
        if source.exists():shutil.copyfile(source,DATA/name)
    algebra=json.loads((DATA/'algebra_cases.json').read_text())
    coverage=json.loads((DATA/'algebra_coverage_cases.json').read_text()) if (DATA/'algebra_coverage_cases.json').exists() else {'cases':[]}
    polynomial_cases=algebra['cases']+coverage['cases']
    endo=json.loads((DATA/'endomorphism_cases.json').read_text())
    glv4=json.loads((DATA/'glv4.json').read_text())
    files=sorted(path for path in DATA.glob('*.json') if path.name!='research_manifest.json')
    manifest={
        'schema':'ecbench-research-manifest-v1','release_version':'1.0.0',
        'core_dataset_id':json.loads((DATA/'manifest.json').read_text())['dataset_id'],
        'counts':{'primary_fields':7,'primary_curves':56,'primary_vectors':6931,
            'polynomial_cases':len(polynomial_cases),'main_polynomial_cases':len(algebra['cases']),
            'coverage_polynomial_cases':len(coverage['cases']),
            'fields_with_polynomial_cases':len({row['field']['id'] for row in polynomial_cases}),
            'glv2_frobenius_supported':int(endo['counts']['supported']),
            'glv2_frobenius_unsupported':int(endo['counts']['unsupported']),
            'glv2_frobenius_scalar_vectors':sum(len(row.get('scalar_vectors',[])) for row in endo['cases']),
            'glv_gls4_supplements':1,'glv_gls4_decomposition_vectors':len(glv4['decomposition_vectors'])},
        'sha256':{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
        'boundaries':[
            'Primary 56-curve corpus and its original canonical manifest are unchanged.',
            'Main catalog covers 49 binary and 6 prime two-summand relation-search systems; a separate coverage catalog adds other fields. Not every primary curve has polynomial fixtures.',
            'The GLV/GLS4 extension-field curve is supplemental; no index-calculus polynomial fixtures over that extension field are claimed.',
            'Recorded solver measurements live in build/algebra_buchberger.json and build/algebra_f5b.json; timings are machine-specific.',
            'All-curve endomorphism classification includes unsupported and subgroup-unstable cases explicitly.',
        ],
    }
    destination=DATA/'research_manifest.json'
    destination.write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
    return manifest


if __name__=='__main__':
    result=build()
    print(json.dumps({'status':'pass','counts':result['counts'],'payloads_hashed':len(result['sha256'])},sort_keys=True))
