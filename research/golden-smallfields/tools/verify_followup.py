#!/usr/bin/env python3
"""Verify follow-up proofs, source immutability, hashes, and reported completion."""
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ecbench.cli import catalog
from ecbench.fixtures import verify_cases
from ecbench.preprocess import digest, verify_transformation

OUT = ROOT / 'followups' / '2026-10-15'


def load(name): return json.loads((OUT / name).read_text())


def main():
    cases = {case['id']: case for case in catalog('algebra_cases')['cases']}
    derived = {case['id']: case for case in load('derived_cases.json')['cases']}
    proofs = load('row_equivalence_proofs.json')['transformations']
    for row in proofs:
        original = cases[row['proof']['source_case_id']]
        verify_transformation(original, derived[original['id']+'-gf2-rref'], row['proof'])
    exhaustive = verify_cases(list(derived.values()))
    historical_files = 0
    for line in (ROOT / 'SHA256SUMS').read_text().splitlines():
        expected, name = line.split('  ', 1)
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError('Historical release source changed: ' + name)
        historical_files += 1
    performance = load('preprocessing_benchmark.json')
    metrics = load('complete_mathematical_metrics.json')
    metric_by_id = {row['case_id']: row for row in metrics['rows']}
    completed, timed_out = 0, 0
    for row in performance['rows']:
        case = cases[row['source_case_id']] if row['variant'] == 'original' else derived[row['source_case_id']+'-gf2-rref']
        if digest(case) != row['result']['case_sha256']:
            raise ValueError('Benchmark result refers to a different fixture')
        if row['result']['status'] == 'ok':
            if not row['result']['relations']['truth_matches']:
                raise ValueError('Computed relation truth mismatch')
            if row['result']['groebner_basis'] != metric_by_id[row['source_case_id']]['groebner_basis_metrics']:
                raise ValueError('Equivalent systems returned different reduced-basis metrics')
            completed += 1
        elif row['result']['status'] == 'timeout': timed_out += 1
        else: raise ValueError('Benchmark backend error')
    if {row['case_id'] for row in metrics['rows']} != set(cases):
        raise ValueError('Incomplete original-fixture metric catalog')
    for name, expected in metrics['source_files_sha256'].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError('Changed metric source: ' + name)
    for row in metrics['rows']:
        if row['case_sha256'] != digest(cases[row['case_id']]):
            raise ValueError('Metric case hash mismatch')
        if row['degree_of_regularity'].get('value') is None:
            raise ValueError('Unfinished regularity in complete metric catalog')
    unsafe = []
    def visit(value, location):
        if type(value) is int and abs(value) > 2**53-1: unsafe.append(location)
        elif isinstance(value, dict):
            for key, child in value.items(): visit(child, location+'.'+key)
        elif isinstance(value, list):
            for i, child in enumerate(value): visit(child, location+'['+str(i)+']')
    for path in OUT.glob('*.json'): visit(json.loads(path.read_text()), path.name)
    if unsafe: raise ValueError('Unsafe large JSON integers: '+str(unsafe))
    result = {'status': 'pass', 'row_equivalence_proofs': len(proofs),
        'derived_original_system_checks': exhaustive,
        'historical_release_files_unchanged': historical_files,
        'paired_benchmark_runs': len(performance['rows']), 'completed_benchmark_runs': completed,
        'censored_benchmark_runs': timed_out, 'complete_mathematical_metric_cases': len(metrics['rows']),
        'large_json_integers': 'all mathematical values above 2^53-1 use strings'}
    (OUT / 'verification.json').write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key!='derived_original_system_checks'}, sort_keys=True))


if __name__ == '__main__': main()
