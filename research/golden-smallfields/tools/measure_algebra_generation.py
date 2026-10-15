#!/usr/bin/env python3
"""Measure actual golden-fixture generation without modifying case payloads.

The whole measurement includes factor-base construction/subgroup tests,
polynomial construction, independent golden-truth enumeration, source hashes,
and JSON serialization. It is not a solver benchmark. Timing wrappers alter
no returned values or canonical records. Generation writes into a temporary
directory and must reproduce the canonical case artifact byte for byte.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager, ExitStack
from datetime import datetime, timezone
import hashlib
import json
import platform
from pathlib import Path
import sys
import tempfile
from time import perf_counter

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import build_algebra_cases as builder
import prime_algebra_cases as prime
from ecbench import polynomials as binary

DEFAULT_CASES=ROOT/'build'/'algebra_cases.json'
DEFAULT_OUTPUT=ROOT/'build'/'algebra_generation_benchmark.json'
_ACTIVE_STAGES=set()


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@contextmanager
def _measure(module, name: str, label: str, records: list, describe=None):
    original=getattr(module,name)
    def timed(*args,**kwargs):
        # Golden truth independently reconstructs f3 for point-sum identity
        # checks. Those nested calls belong to truth time, not construction.
        if (label=='prime_semaev_polynomial_construction' and
                'prime_golden_truth_enumeration' in _ACTIVE_STAGES):
            return original(*args,**kwargs)
        started=perf_counter()
        _ACTIVE_STAGES.add(label)
        try:
            result=original(*args,**kwargs)
        finally:
            _ACTIVE_STAGES.discard(label)
        elapsed=perf_counter()-started
        record={'stage':label,'elapsed_seconds':elapsed}
        if describe is not None:
            record.update(describe(args,kwargs,result))
        records.append(record)
        return result
    setattr(module,name,timed)
    try:
        yield
    finally:
        setattr(module,name,original)


def _binary_base(args,kwargs,result):
    base=result[0]
    return {'curve_id':args[1]['id'],'dimension':len(args[2]),'basis_kind':args[3],
            'accepted_x':base['counters']['accepted_x'],
            'scope':'Factor-base span enumeration, field point finding, selected-subgroup membership, sign completion, and metadata'}


def _binary_case(args,kwargs,result):
    return {'case_id':result['id'],'target_kind':result['target']['kind'],
            'ordered_relations':result['expected']['ordered_count'],
            'scope':'Polynomial construction, Boolean/membership equations, exact x and signed-point golden-truth enumeration, and case assembly'}


def _prime_base(args,kwargs,result):
    return {'curve_id':prime.CURVE_ID,'factor_base_size':args[3],
            'scope':'Small subgroup multiples, exact membership/order checks, sign completion, and metadata'}


def _prime_truth(args,kwargs,result):
    expected,checks=result
    return {'curve_id':prime.CURVE_ID,'factor_base_size':len(args[2]),
            'target_x':str(args[1][0]),'ordered_relations':expected['ordered_count'],
            'scope':'Exhaustive polynomial-domain checks and independent signed-point relation enumeration'}


def measure(canonical: Path=DEFAULT_CASES,output: Path=DEFAULT_OUTPUT):
    expected=json.loads(canonical.read_text())
    seed=expected['seed']
    stage_records=[]
    with ExitStack() as stack:
        stack.enter_context(_measure(binary,'_factor_base','binary_factor_base',stage_records,_binary_base))
        stack.enter_context(_measure(binary,'_binary_case','binary_case_construction_and_truth',stage_records,_binary_case))
        stack.enter_context(_measure(prime,'_factor_base','prime_factor_base',stage_records,_prime_base))
        stack.enter_context(_measure(prime,'_expected','prime_golden_truth_enumeration',stage_records,_prime_truth))
        stack.enter_context(_measure(prime,'_semaev_f3','prime_semaev_polynomial_construction',stage_records))
        stack.enter_context(_measure(prime,'_membership_polynomial','prime_membership_polynomial_construction',stage_records))
        with tempfile.TemporaryDirectory(prefix='ecbench-algebra-generation-') as tmp:
            regenerated=Path(tmp)/'algebra_cases.json'
            started=perf_counter()
            generated=builder.build(seed,regenerated)
            elapsed=perf_counter()-started
            regenerated_sha=_sha(regenerated)
            canonical_sha=_sha(canonical)
            if regenerated_sha!=canonical_sha:
                raise RuntimeError('Actual generation did not reproduce the canonical case payload hash')
            if regenerated.read_bytes()!=canonical.read_bytes():
                raise RuntimeError('Hash match did not reproduce identical canonical bytes')
    # Join prime truth measurements to case identifiers after timing, without
    # changing any payload or including this bookkeeping in generation time.
    prime_index={(len(case['factor_base']['x_values']),str(int(case['target']['point']['x'],16))):case['id']
                 for case in generated['cases'] if case['field']['kind']=='prime'}
    for record in stage_records:
        if record['stage']=='prime_golden_truth_enumeration':
            record['case_id']=prime_index[(record['factor_base_size'],record.pop('target_x'))]
    stages={}
    for record in stage_records:
        stage=stages.setdefault(record['stage'],{'calls':0,'elapsed_seconds':0.0})
        stage['calls']+=1
        stage['elapsed_seconds']+=record['elapsed_seconds']
    measured_sum=sum(stage['elapsed_seconds'] for stage in stages.values())
    result={
        'schema':'ecbench-algebra-generation-benchmark-v1',
        'created_at':datetime.now(timezone.utc).isoformat(),
        'status':'passed','seed':seed,'counts':generated['counts'],
        'measurement':{'elapsed_seconds':elapsed,'scope':'Actual full deterministic golden-fixture build, including golden truth',
            'includes':['factor-base generation and selected-subgroup membership checks','polynomial-system construction',
                        'golden finite-domain and signed-point truth enumeration','canonical source hashing and case JSON serialization'],
            'excludes':['Python interpreter startup and module imports','post-build hash comparison and timing-result assembly',
                        'Groebner solving, solver-root enumeration, and relation verification on solver output'],
            'golden_truth_enumeration':'included in whole elapsed time; separately measured for prime cases; combined with construction for binary cases',
            'stage_timing_overhead':'included; lightweight in-process perf_counter wrappers',
            'stage_records_nonoverlapping':True,
            'instrumented_stage_total_seconds':measured_sum,
            'unattributed_pipeline_seconds':elapsed-measured_sum,
            'interpretation':'Fixture preparation cost, not solve-pipeline or index-calculus attack time'},
        'integrity':{'canonical_artifact':str(canonical.relative_to(ROOT)) if canonical.is_relative_to(ROOT) else str(canonical),
                     'canonical_sha256':canonical_sha,'regenerated_sha256':regenerated_sha,
                     'byte_identical':True,'canonical_artifact_modified':False},
        'environment':{'python_version':platform.python_version(),'python_implementation':platform.python_implementation(),
                       'platform':platform.platform(),'processor':platform.processor(),'clock':'time.perf_counter'},
        'source_hashes':{str(path.relative_to(ROOT)):_sha(path) for path in [
            ROOT/'ecbench'/'polynomials.py',ROOT/'tools'/'build_algebra_cases.py',
            ROOT/'tools'/'prime_algebra_cases.py',ROOT/'tools'/'prime_vectors.py']},
        'stage_totals':stages,'stage_records':stage_records}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',type=Path,default=DEFAULT_CASES)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    result=measure(args.cases,args.output)
    print(json.dumps({'status':result['status'],'cases':result['counts']['cases'],
                      'elapsed_seconds':result['measurement']['elapsed_seconds'],
                      'byte_identical':result['integrity']['byte_identical']},sort_keys=True))
    print(f'Wrote {args.output}')


if __name__=='__main__':
    main()
