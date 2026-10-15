#!/usr/bin/env python3
"""Build a deterministic, independently truth-checked two-point algebra battery.

Primary ecbench field/curve data is read-only. Output is a separate experimental
fixture artifact; planted and seeded-target cohorts remain explicitly distinct.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ecbench.polynomials import build_binary_cases
from prime_algebra_cases import build_prime_cases

DEFAULT_OUTPUT=ROOT/'build'/'algebra_cases.json'
DEFAULT_SEED='ecbench-indexcalculus-v1'


def build(seed: str=DEFAULT_SEED,output: Path=DEFAULT_OUTPUT):
    data_dir=ROOT/'ecbench'/'data'
    fields=json.loads((data_dir/'fields.json').read_text())
    curves=json.loads((data_dir/'curves.json').read_text())
    binary=build_binary_cases(fields,curves,seed)
    prime=build_prime_cases(seed)
    cases=binary+prime
    assert len(cases)==55 and len({case['id'] for case in cases})==55
    for case in cases:
        case.setdefault('variable_order',case['variables'])
        for polynomial in case['polynomials']:
            polynomial.setdefault('total_degree',max((sum(term['exponents']) for term in polynomial['terms']),default=0))
            assert all(len(term['exponents'])==len(case['variables']) for term in polynomial['terms'])
            assert all(type(term['coefficient']) is str for term in polynomial['terms'])
        assert len(case['expected']['algebraic_assignments'])==case['expected']['ordered_x_count']
    metadata={'schema':'ecbench-semaev-two-point-cases-v1','seed':seed,
        'coefficient_encoding':'canonical decimal strings modulo coefficient_modulus',
        'exponent_encoding':'nonnegative bounded JSON integers, in variables order',
        'primary_dataset_modified':False,'relation_arity':2,'summation_polynomial_arity':3,
        'counts':{'cases':len(cases),'binary_cases':len(binary),'prime_cases':len(prime),
                  'planted_cases':sum(case['target']['kind']=='planted' for case in cases),
                  'seeded_target_cases':sum(case['target']['kind']=='seeded_subgroup' for case in cases),
                  'empty_factor_base_diagnostics':sum(case['generation'].get('diagnostic_empty_factor_base',False) for case in cases)},
        'source_hashes':{name:hashlib.sha256((data_dir/name).read_bytes()).hexdigest() for name in ('fields.json','curves.json')},
        'boundaries':['All cases are small two-point f3 fixtures, not a complete index-calculus attack.',
            'Planted and subgroup-seeded-basis cases are deliberate performance fixtures, not unbiased relation-yield samples.',
            'Low-basis B83 cases honestly have empty selected-subgroup factor bases and unit membership equations.',
            'Binary strict subgroup membership adds separate indicator equations of degree up to dimension; do not attribute their degree to f3.',
            'Truth is independently enumerated with point arithmetic and supplied only for post-solve conformance checking.',
            'No higher-arity summation systems, generic-hardness, asymptotic degree-of-regularity, or production-security estimate is claimed.'],
        'cases':cases}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(metadata,indent=2)+'\n')
    return metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',default=DEFAULT_SEED)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    result=build(args.seed,args.output)
    print(json.dumps(result['counts'],sort_keys=True))
    print(f'Wrote {args.output}')


if __name__=='__main__':
    main()
