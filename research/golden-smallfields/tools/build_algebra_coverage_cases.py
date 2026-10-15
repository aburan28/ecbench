#!/usr/bin/env python3
"""Build 14 supplementary two-point fixtures covering the other five fields.

The original 55-case catalog and its measurements are read-only. Existing
binary and prime construction helpers are reused without changing their
defaults, parameters or original fixtures.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import hashlib
import itertools
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ecbench import polynomials as binary
from ecbench.algebra import _eval_sparse, _normalize_case, _record_input, _verify_relations
from ecbench.validate import load_dataset
import prime_algebra_cases as prime

DEFAULT_SEED = 'ecbench-indexcalculus-v1'
DEFAULT_OUTPUT = ROOT / 'build' / 'algebra_coverage_cases.json'
BASE_CATALOG = ROOT / 'build' / 'algebra_cases.json'


def _binary_cases(fields, records, seed):
    cases = []
    for curve_id in ('K53-a1', 'K79-a1'):
        record = records[curve_id]
        field_record = fields[record['field_id']]
        field = binary.BinaryField(field_record['degree'], field_record['representation']['polynomial_exponents'])
        curve = binary.BinaryCurve(field, int(record['coefficients']['a2'], 16), int(record['coefficients']['a6'], 16))
        generator = tuple(int(record['subgroup']['generator'][key], 16) for key in ('x', 'y'))
        order = int(record['subgroup']['order'])
        assert curve.on_curve(generator) and curve.mul(order, generator) is None
        for dimension in (2, 3):
            basis_kind = 'subgroup_seeded_basis'
            basis = binary._basis(curve, generator, dimension, basis_kind)
            factor_base, points, masks = binary._factor_base(curve, record, basis, basis_kind)
            planted, pair = binary._planted_target(curve, points, generator, basis_kind)
            assert planted is not None and pair is not None
            cases.append(binary._binary_case(record, field_record, curve, basis, factor_base, points, masks,
                                             planted, 'planted', seed, 3, pair))
            label = f'{curve_id}/{basis_kind}/d{dimension}/seeded_subgroup'
            scalar = binary._seeded_scalar(seed, label, order)
            target = curve.mul(scalar, generator)
            assert target is not None
            cases.append(binary._binary_case(record, field_record, curve, basis, factor_base, points, masks,
                                             target, 'seeded_subgroup', seed, scalar, None))
    return cases


def _prime_target(curve, generator, order, curve_id, size, kind, seed):
    if kind == 'planted':
        # The existing helper's planted branch does not use its default curve ID.
        return prime._target(curve, generator, order, size, kind, seed)
    payload = f'ecbench-prime-f3-v1\0{seed}\0{curve_id}\0{size}\0seeded_subgroup'
    digest = hashlib.sha256(payload.encode('utf-8')).hexdigest()
    scalar = 1 + int(digest, 16) % (order - 1)
    target = curve.mul(scalar, generator)
    assert target is not None and curve.on_curve(target) and curve.mul(order, target) is None
    return target, scalar, {
        'target_method': 'deterministic SHA-256 nonzero subgroup scalar', 'hash': 'SHA-256',
        'hash_input_encoding': 'UTF-8', 'hash_input': payload, 'hash_digest_hex': digest,
        'scalar_formula': '1 + OS2IP(SHA-256(hash_input)) mod (r - 1)',
        'randomness_claim': 'none; deterministic seeded test target',
    }


def _prime_cases(fields, records, seed):
    cases = []
    size = 4
    for curve_id in ('p61-G-0', 'p64-G-0', 'p83s-G-0'):
        record = records[curve_id]
        field = fields[record['field_id']]
        assert field['kind'] == 'prime' and record['model'] == 'short_weierstrass'
        assert all(int(record['coefficients'][key], 16) == 0 for key in ('a1', 'a2', 'a3'))
        curve = prime.PrimeCurve(int(field['characteristic']), int(record['coefficients']['a4'], 16),
                                 int(record['coefficients']['a6'], 16))
        generator = prime._parse_point(record['subgroup']['generator'])
        order, width = int(record['subgroup']['order']), int(field['coordinate_bytes'])
        assert order > 2 * size and record['subgroup']['unique_order_r_subgroup']
        assert record['subgroup']['r_torsion_rank'] == 1
        assert generator is not None and curve.on_curve(generator) and curve.mul(order, generator) is None
        factor_base, xs, points = prime._factor_base(curve, generator, order, size, width)
        for kind in ('planted', 'seeded_subgroup'):
            target, scalar, derivation = _prime_target(curve, generator, order, curve_id, size, kind, seed)
            polynomials = [prime._semaev_f3(curve, target[0]),
                           prime._membership_polynomial(xs, 0, curve.p),
                           prime._membership_polynomial(xs, 1, curve.p)]
            expected, checks = prime._expected(curve, target, xs, points, polynomials, width)
            cases.append({
                'id': f'{curve_id}-f3-fb{size}-{kind}', 'curve_id': curve_id,
                'curve': deepcopy(record), 'field': deepcopy(field), 'coefficient_modulus': str(curve.p),
                'variables': ['x1', 'x2'],
                'polynomials': [prime._serialize_polynomial(poly, name, role) for poly, name, role in zip(
                    polynomials, ('f3', 'membership_x1', 'membership_x2'),
                    ('semaev_f3', 'factor_base_membership', 'factor_base_membership'))],
                'coordinate_groups': [{'name': 'x1', 'variables': ['x1']}, {'name': 'x2', 'variables': ['x2']}],
                'domains': {'x1': list(map(str, xs)), 'x2': list(map(str, xs))},
                'factor_base': deepcopy(factor_base),
                'target': {'point': prime._point_record(target, width), 'kind': kind, 'source': kind, 'scalar': str(scalar)},
                'model_variant': 'strict_selected_subgroup_factor_base', 'expected': expected,
                'arity': 2, 'summation_polynomial_arity': 3,
                'generation': {'method': 'deterministic kG factor base and finite exhaustive verification',
                    'seed': seed, 'factor_base_rule': 'x(kG) for k=1..size; both point signs permitted',
                    'factor_base_scalar_start': '1', 'factor_base_scalar_stop': str(size),
                    'formula': '(x1-x2)^2*z^2 - 2*((x1+x2)*(x1*x2+A)+2*B)*z + (x1*x2-A)^2 - 4*B*(x1+x2)',
                    'checks': checks, **derivation},
            })
    return cases


def _verify(case):
    """Re-evaluate every declared domain assignment, then use the runtime oracle."""
    normalized = _normalize_case(case)
    _, modulus, polynomials, domains, metrics = _record_input(normalized)
    roots = [list(assignment) for assignment in itertools.product(*domains)
             if all(_eval_sparse(list(poly.items()), assignment, modulus) == 0 for poly in polynomials)]
    roots.sort()
    relations = _verify_relations(normalized, roots)
    assert relations['truth_matches'] is True
    assert roots == [[int(value) for value in assignment] for assignment in case['expected']['algebraic_assignments']]
    if case['target']['kind'] == 'planted':
        assert relations['target_success']
    return {'case_id': case['id'], 'status': 'pass', 'domain_assignments_checked': metrics['domain_grid_size'],
            'polynomial_roots': len(roots), 'ordered_relations': relations['ordered_relations'],
            'unordered_relations': relations['unordered_relations'], 'truth_matches': True,
            'selected_subgroup_membership_verified': relations['selected_subgroup_membership_verified']}


def build(seed=DEFAULT_SEED, output=DEFAULT_OUTPUT):
    if not isinstance(seed, str):
        raise TypeError('seed must be a string')
    output = Path(output)
    if output.resolve() == BASE_CATALOG.resolve():
        raise ValueError('The original 55-case catalog must remain unchanged')
    original_bytes = BASE_CATALOG.read_bytes()
    base = json.loads(original_bytes)
    assert len(base['cases']) == 55
    data = load_dataset()
    fields = {record['id']: record for record in data['fields']}
    records = {record['id']: record for record in data['curves']}
    cases = _binary_cases(fields, records, seed) + _prime_cases(fields, records, seed)
    assert len(cases) == 14 and len({case['id'] for case in cases}) == 14
    assert not {case['id'] for case in cases} & {case['id'] for case in base['cases']}
    for case in cases:
        case['field_id'] = case['field']['id']
        case.setdefault('variable_order', case['variables'])
        case['generation']['coverage_cohort'] = 'supplementary all-seven-field coverage'
        for polynomial in case['polynomials']:
            polynomial.setdefault('total_degree', max((sum(term['exponents']) for term in polynomial['terms']), default=0))
            assert all(type(term['coefficient']) is str and str(int(term['coefficient'])) == term['coefficient']
                       and 0 <= int(term['coefficient']) < int(case['coefficient_modulus']) for term in polynomial['terms'])
    verification = [_verify(case) for case in cases]
    original_fields = {case['field']['id'] for case in base['cases']}
    added_fields = {case['field']['id'] for case in cases}
    assert original_fields | added_fields == set(fields)
    counts = Counter(case['field']['id'] for case in cases)
    artifact = {
        'schema': 'ecbench-semaev-two-point-cases-v1', 'catalog_role': 'supplementary all-field coverage',
        'seed': seed, 'coefficient_encoding': 'canonical decimal strings modulo coefficient_modulus',
        'exponent_encoding': 'nonnegative bounded JSON integers, in variables order',
        'primary_dataset_modified': False, 'original_55_catalog_modified': False,
        'relation_arity': 2, 'summation_polynomial_arity': 3,
        'counts': {'cases': 14, 'binary_cases': 8, 'prime_cases': 6, 'planted_cases': 7,
                   'seeded_target_cases': 7, 'empty_factor_base_diagnostics': 0,
                   'cases_by_field': dict(sorted(counts.items())), 'combined_cases': 69, 'combined_fields': 7},
        'coverage': {'original_catalog_fields': sorted(original_fields), 'supplementary_fields': sorted(added_fields),
                     'combined_fields': sorted(original_fields | added_fields)},
        'source_hashes': {name: hashlib.sha256((ROOT / 'ecbench' / 'data' / name).read_bytes()).hexdigest()
                          for name in ('fields.json', 'curves.json')},
        'original_55_catalog_sha256': hashlib.sha256(original_bytes).hexdigest(),
        'validation': {'status': 'pass',
            'method': 'independent helper formula/point enumeration, full sparse-equation domain re-evaluation, and runtime reference signed-point/subgroup checks',
            'cases': verification},
        'boundaries': ['This separate extension adds field coverage to the original 55-case complexity battery.',
                       'The 14 cases are deliberate small two-point fixtures, not an unbiased relation-yield sample or complete ECDLP attack.',
                       'Golden truths are supplied for post-solve conformance; they are not solver roots.',
                       'Generation and validation do not measure Groebner solving or degree of regularity.',
                       'The original 55-case catalog and existing benchmark results are unchanged.'],
        'cases': cases,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2) + '\n')
    assert BASE_CATALOG.read_bytes() == original_bytes
    return artifact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', default=DEFAULT_SEED)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    artifact = build(args.seed, args.output)
    print(json.dumps({'counts': artifact['counts'], 'coverage': artifact['coverage'],
                      'validation': artifact['validation']['status']}, sort_keys=True))
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
