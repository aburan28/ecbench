#!/usr/bin/env python3
"""Reproduce bounded follow-up experiments without replacing the golden release."""
import argparse
from collections import Counter, defaultdict
import datetime
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ecbench.algebra import run_case
from ecbench.cli import catalog
from ecbench.fixtures import verify_cases as verify_fixtures
from ecbench.preprocess import derive_case, digest, verify_transformation
from ecbench.reference import Curve, Field
from ecbench.validate import load_dataset

OUT = ROOT / 'followups' / '2026-10-15'
BENCH_IDS = ['K83-a0-f3-low_monomials-d6-planted',
             'S83-a1-f3-low_monomials-d6-planted',
             'S83-a1-f3-low_monomials-d6-seeded_subgroup',
             'B83-a0-f3-subgroup_seeded_basis-d5-planted']


def write(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')


def stamp():
    return {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'environment': {'python': platform.python_version(), 'platform': platform.platform()},
        'primary_dataset_modified': False}


def derive(cases):
    rows, derived = [], []
    for case in cases:
        if int(case['coefficient_modulus']) != 2:
            continue
        changed, proof = derive_case(case)
        checked = verify_transformation(case, changed, proof)
        rows.append({'proof': proof, 'verification': checked})
        derived.append(changed)
    truth = verify_fixtures(derived)
    write('derived_cases.json', {'schema': 'ecbench-row-reduced-cases-v1', **stamp(),
        'boundaries': 'Exact GF(2) constant-row operations; new case hashes; original fixture catalogs unchanged',
        'counts': {'cases': len(derived)}, 'cases': derived})
    write('row_equivalence_proofs.json', {'schema': 'ecbench-row-equivalence-proofs-v1',
        **stamp(), 'transformations': rows, 'original_system_truth_rechecked': truth})
    print(f'derive: {len(derived)} exact row-span proofs and exhaustive fixture checks pass', flush=True)


def recover(cases, timeout):
    by_id = {case['id']: case for case in cases}
    rows = []
    result = {'schema': 'ecbench-timeout-recovery-v1', **stamp(), 'timeout_seconds': timeout,
        'scope': 'One new run only for each original incomplete solve or timed-out regularity diagnostic; historical runs unchanged',
        'results': rows}
    for method in ('buchberger', 'f5b'):
        original = json.loads((ROOT / 'build' / f'algebra_{method}.json').read_text())
        for old in original['results']:
            if old['status'] == 'ok' and old.get('diagnostic_status') != 'timeout':
                continue
            case = by_id[old['case_id']]
            if old['case_sha256'] != digest(case):
                raise ValueError('Historical case hash changed')
            new = run_case(case, timeout=timeout, method=method, diagnostics=True)
            row = {'original_status': old['status'], 'original_diagnostic_status': old.get('diagnostic_status'),
                'original_wall_budget_seconds': old['timeout_seconds'], 'result': new}
            rows.append(row)
            write('timeout_recovery.json', result)
            print(f'recover {method} {case["id"]}: {new["status"]}; dreg={new["degree_of_regularity"].get("value")}', flush=True)
    result['summary'] = {'attempts': len(rows),
        'completed_solves': sum(row['result']['status'] == 'ok' for row in rows),
        'completed_regularity_diagnostics': sum(row['result']['degree_of_regularity'].get('value') is not None for row in rows),
        'previous_solve_timeouts': sum(row['original_status'] == 'timeout' for row in rows),
        'recovered_solve_timeouts': sum(row['original_status'] == 'timeout' and row['result']['status'] == 'ok' for row in rows),
        'all_completed_truth_matches': all(row['result'].get('relations', {}).get('truth_matches') for row in rows if row['result']['status'] == 'ok')}
    write('timeout_recovery.json', result)


def compare(cases, timeout, repetitions):
    by_id = {case['id']: case for case in cases}
    rows = []
    result = {'schema': 'ecbench-gf2-rref-comparison-v1', **stamp(), 'repetitions': repetitions,
        'timeout_seconds': timeout, 'diagnostics': False,
        'scope': 'Fresh isolated paired runs; exact preprocessing and witness verification included; prebuilt fixtures, no factor-base generation',
        'timing_definitions': {
            'pipeline_wall_seconds': 'One enclosing parent timer from before preprocessing/proof verification through child solver completion; includes child startup and serialization',
            'core_with_preprocessing_seconds': 'Worker core solve time plus parent preprocessing and exact witness verification; excludes process startup and serialization',
            'diagnostic_work': 'Disabled in these performance pairs; separately measured in timeout recovery and regularity follow-up'},
        'rows': rows}
    for repetition in range(repetitions):
        for ci, cid in enumerate(BENCH_IDS):
            for mi, method in enumerate(('buchberger', 'f5b')):
                variants = ['original', 'row_reduced'] if (repetition + ci + mi) % 2 == 0 else ['row_reduced', 'original']
                for variant in variants:
                    original = by_id[cid]
                    started = time.perf_counter()
                    preprocessing = proof_seconds = 0.0
                    case = original
                    if variant == 'row_reduced':
                        stage = time.perf_counter()
                        case, proof = derive_case(original)
                        preprocessing = time.perf_counter() - stage
                        stage = time.perf_counter()
                        verify_transformation(original, case, proof)
                        proof_seconds = time.perf_counter() - stage
                    actual = run_case(case, timeout=timeout, method=method, diagnostics=False)
                    pipeline = time.perf_counter() - started
                    core = actual['stage_times'].get('total_solve_seconds')
                    rows.append({'source_case_id': cid, 'method': method, 'variant': variant,
                        'repetition': repetition, 'execution_index': len(rows),
                        'preprocessing_seconds': preprocessing, 'proof_verification_seconds': proof_seconds,
                        'pipeline_wall_seconds': pipeline,
                        'core_with_preprocessing_seconds': None if core is None else core + preprocessing + proof_seconds,
                        'result': actual})
                    write('preprocessing_benchmark.json', result)
                    print(f'compare r{repetition + 1} {method} {variant} {cid}: {actual["status"]} ({pipeline:.3f}s wall)', flush=True)
    summaries = []
    for cid in BENCH_IDS:
        for method in ('buchberger', 'f5b'):
            selected = [row for row in rows if row['source_case_id'] == cid and row['method'] == method]
            variants = {}
            for variant in ('original', 'row_reduced'):
                chosen = [row for row in selected if row['variant'] == variant]
                completed = [row for row in chosen if row['result']['status'] == 'ok']
                variants[variant] = {'attempts': len(chosen), 'completed': len(completed),
                    'median_pipeline_wall_seconds_completed': statistics.median(row['pipeline_wall_seconds'] for row in completed) if completed else None,
                    'median_core_with_preprocessing_seconds_completed': statistics.median(row['core_with_preprocessing_seconds'] for row in completed) if completed else None}
            matched = [rep for rep in range(repetitions) if all(any(row['repetition'] == rep and row['variant'] == variant and row['result']['status'] == 'ok' for row in selected) for variant in variants)]
            pair_wall_ratios = []
            pair_core_ratios = []
            for rep in matched:
                pair = {row['variant']: row for row in selected if row['repetition'] == rep}
                pair_wall_ratios.append(pair['original']['pipeline_wall_seconds'] / pair['row_reduced']['pipeline_wall_seconds'])
                pair_core_ratios.append(pair['original']['core_with_preprocessing_seconds'] / pair['row_reduced']['core_with_preprocessing_seconds'])
            summaries.append({'source_case_id': cid, 'method': method, 'variants': variants,
                'matched_completed_pairs': len(matched),
                'paired_wall_speedup_ratios': pair_wall_ratios,
                'paired_core_speedup_ratios': pair_core_ratios,
                'median_paired_wall_speedup': statistics.median(pair_wall_ratios) if pair_wall_ratios else None,
                'median_paired_core_speedup': statistics.median(pair_core_ratios) if pair_core_ratios else None})
    result['summary'] = summaries
    result['all_completed_truth_matches'] = all(row['result'].get('relations', {}).get('truth_matches') for row in rows if row['result']['status'] == 'ok')
    write('preprocessing_benchmark.json', result)


def rational(value):
    return {'numerator': str(value.numerator), 'denominator': str(value.denominator), 'approximation': float(value)}


def regularity(cases, timeout):
    by_id = {case['id']: case for case in cases}
    rows = []
    worker = '''import json,sys,time
import sympy
from ecbench.algebra import _normalize_case,_record_input,_degree_regularity
case=json.load(sys.stdin)
names,modulus,polynomials,domains,metrics=_record_input(_normalize_case(case))
symbols=sympy.symbols(names)
started=time.perf_counter()
value=_degree_regularity(sympy,polynomials,symbols,modulus,'buchberger','grevlex')
print(json.dumps({'status':'ok','degree_of_regularity':value,'diagnostic_seconds':time.perf_counter()-started}))
'''
    result = {'schema': 'ecbench-regularity-generator-dependence-v1', **stamp(),
        'method': 'buchberger', 'order': 'grevlex', 'timeout_seconds': timeout,
        'scope': 'Separate homogeneous-top-part diagnostics. Original and row-reduced affine ideals are exactly equal, but their highest-part generator ideals need not be.', 'rows': rows}
    for cid in BENCH_IDS:
        original = by_id[cid]
        changed, proof = derive_case(original)
        verify_transformation(original, changed, proof)
        for variant, case in [('original', original), ('row_reduced', changed)]:
            started = time.perf_counter()
            try:
                process = subprocess.run([sys.executable, '-c', worker], input=json.dumps(case),
                    capture_output=True, text=True, cwd=ROOT, timeout=timeout)
                actual = json.loads(process.stdout) if process.returncode == 0 else {'status': 'error', 'error': process.stderr[-2000:]}
            except subprocess.TimeoutExpired:
                actual = {'status': 'timeout', 'degree_of_regularity': {'value': None, 'reason': 'diagnostic_timeout'}}
            rows.append({'source_case_id': cid, 'case_sha256': digest(case), 'variant': variant,
                'observed_process_wall_seconds': time.perf_counter() - started, 'result': actual})
            write('regularity_preprocessing_comparison.json', result)
            print(f'regularity {variant} {cid}: {actual["status"]} {actual.get("degree_of_regularity",{}).get("value")}', flush=True)


def density(cases):
    groups = defaultdict(list)
    for case in cases:
        encoded = json.dumps({'curve_id': case['curve_id'], 'points': case['factor_base']['points_by_x']}, sort_keys=True, separators=(',', ':'))
        groups[encoded].append(case)
    rows = []
    model_total = Fraction(0)
    seeded_count = 0
    for encoded, family in groups.items():
        case = family[0]
        curve = Curve(case['curve'], Field(case['field']))
        pool = sorted({curve.parse_point(point) for points in case['factor_base']['points_by_x'].values() for point in points})
        if any(point is None or not curve.on_curve(point) or not curve.subgroup_contains(point) for point in pool):
            raise ValueError('Invalid factor base in density enumeration')
        ordered, unordered = Counter(), Counter()
        for i, p in enumerate(pool):
            for q in pool[i:]:
                target = curve.add(p, q)
                unordered[target] += 1
                ordered[target] += 1 if p == q else 2
        nonidentity = [point for point in unordered if point is not None]
        success = Fraction(len(nonidentity), curve.r - 1)
        expected_unordered = Fraction(sum(unordered[point] for point in nonidentity), curve.r - 1)
        checks = []
        for member in family:
            target = curve.parse_point(member['target']['point'])
            if ordered[target] != int(member['expected']['ordered_count']) or unordered[target] != int(member['expected']['swap_orbit_count']):
                raise ValueError('Pair-sum census disagrees with golden relation counts')
            checks.append({'case_id': member['id'], 'ordered_relations': ordered[target], 'unordered_relations': unordered[target]})
            if member['target']['kind'] == 'seeded_subgroup':
                model_total += success
                seeded_count += 1
        support = [{'target': curve.point(point), 'ordered_relations': ordered[point], 'unordered_relations': unordered[point]}
                   for point in sorted(nonidentity)]
        rows.append({'curve_id': case['curve_id'], 'field_id': case['field']['id'],
            'factor_base_id': hashlib.sha256(encoded.encode()).hexdigest(),
            'case_ids': [member['id'] for member in family], 'factor_base_points': len(pool),
            'subgroup_order': str(curve.r), 'nonzero_target_population': str(curve.r - 1),
            'distinct_nonidentity_pair_sums': len(nonidentity),
            'identity_unordered_pairs_excluded_from_nonzero_model': unordered[None],
            'uniform_nonzero_target_success_probability': rational(success),
            'uniform_nonzero_target_expected_unordered_relations': rational(expected_unordered),
            'support': support, 'golden_target_checks': checks})
    result = {'schema': 'ecbench-exact-pair-sum-density-v1', **stamp(), 'distinct_factor_bases': len(rows),
        'golden_cases_checked': sum(len(row['case_ids']) for row in rows),
        'cohort': {'seeded_targets': seeded_count,
            'uniform_model_expected_successes': rational(model_total),
            'uniform_model_any_success_union_bound': rational(min(Fraction(1), model_total))},
        'model_boundary': 'Exact support fractions for hypothetical uniform NONZERO subgroup targets; the existing deterministic hash/modulo targets are not random trials, so this is not a p-value or a measured attack success probability',
        'factor_bases': rows}
    write('relation_density.json', result)
    print(f'density: {len(rows)} factor bases; all {len(cases)} golden relation counts match; uniform-model expected successes={float(model_total):.6g}', flush=True)


def isogeny():
    data = load_dataset()
    field = Field(next(row for row in data['fields'] if row['id'] == 'p61'))
    source_record = next(row for row in data['curves'] if row['id'] == 'p61-E4-0')
    target_record = next(row for row in data['curves'] if row['id'] == 'p61-E4-1')
    source, target = Curve(source_record, field), Curve(target_record, field)
    p = field.q
    ratio = field.parse(target_record['coefficients']['a4']) * pow(4, -1, p) % p
    inner = field.sqrt(ratio)
    for value in (inner, -inner % p):
        try:
            unit = field.sqrt(value)
            break
        except ValueError:
            pass
    if pow(unit, 4, p) != ratio:
        raise ValueError('Invalid quotient-to-representative isomorphism')
    u2, u3 = pow(unit, 2, p), pow(unit, 3, p)
    ui2, ui3 = pow(u2, -1, p), pow(u3, -1, p)
    inv4, inv8 = pow(4, -1, p), pow(8, -1, p)

    def forward(point):
        if point is None or point[0] == 0:
            return None
        x, y = point
        inv = pow(x, -1, p)
        return (u2 * (x - inv) % p, u3 * y * (1 + inv * inv) % p)

    def dual(point):
        if point is None or point[0] == 0:
            return None
        x, y = point[0] * ui2 % p, point[1] * ui3 % p
        inv = pow(x, -1, p)
        return ((x + 4 * inv) * inv4 % p, y * (1 - 4 * inv * inv) * inv8 % p)

    import sympy as sp
    x = sp.Symbol('x')
    X = x - 1 / x
    curve_identity = sp.cancel((x**3 - x) * (1 + 1 / x**2)**2 - X**3 - 4 * X)
    Xback = (x + 4 / x) / 4
    dual_identity = sp.cancel((x**3 + 4*x) * (1 - 4 / x**2)**2 / 64 - Xback**3 + Xback)
    if curve_identity != 0 or dual_identity != 0:
        raise ValueError('Symbolic isogeny equation failed')
    y = sp.Symbol('y')
    xf, yf = x - 1/x, y*(1 + 1/x**2)
    xd, yd = (x + 4/x)/4, y*(1 - 4/x**2)/8
    compositions = [((xf + 4/xf)/4, yf*(1 - 4/xf**2)/8, -1),
                    (xd - 1/xd, yd*(1 + 1/xd**2), 4)]
    for xx, yy, aa in compositions:
        slope = (3*x*x + aa)/(2*y)
        double_x = slope*slope - 2*x
        double_y = slope*(x - double_x) - y
        relation = sp.Poly(y*y - x*x*x - aa*x, y)
        for difference in (xx - double_x, yy - double_y):
            numerator = sp.together(difference).as_numer_denom()[0]
            if not sp.Poly(numerator, y).rem(relation).is_zero:
                raise ValueError('Symbolic dual composition is not doubling')
    vectors, pair_checks = [], 0
    for curve, other, mapping, reverse, direction in [(source, target, forward, dual, 'forward'), (target, source, dual, forward, 'dual')]:
        points = {None, (0, 0)}
        points.update(curve.parse_point(row['point']) for row in curve.record['named_points'].values())
        points.update(curve.mul(k, curve.parse_point(curve.record['named_points']['R']['point'])) for k in [-7, -2, -1, 1, 2, 3, 17, 2**59, 2**60-1])
        ordered_points = [None] + sorted(point for point in points if point is not None)
        for point in ordered_points:
            image = mapping(point)
            if not other.on_curve(image) or reverse(image) != curve.double(point):
                raise ValueError('Isogeny target or dual-composition check failed')
            vectors.append({'direction': direction, 'point': curve.point(point),
                'image': other.point(image), 'dual_composition': curve.point(reverse(image)), 'expected_double': curve.point(curve.double(point))})
        for left in ordered_points:
            for right in ordered_points:
                if mapping(curve.add(left, right)) != other.add(mapping(left), mapping(right)):
                    raise ValueError('Isogeny homomorphism check failed')
                pair_checks += 1
    result = {'schema': 'ecbench-represented-dual-isogeny-v1', **stamp(), 'status': 'pass',
        'source_curve_id': source_record['id'], 'target_curve_id': target_record['id'], 'degree': 2,
        'field_prime': str(p), 'representative_scale_u': str(unit),
        'source_equation': 'y^2=x^3-x', 'target_equation': 'y^2=x^3-37*x',
        'auxiliary_quotient_equation': 'Y^2=X^3+4*X',
        'forward_formula': '(u^2*(x-1/x), u^3*y*(1+1/x^2)); infinity and (0,0) -> infinity',
        'dual_formula': 'Undo (u^2,u^3) scaling, then ((X+4/X)/4,Y*(1-4/X^2)/8); infinity and (0,0) -> infinity',
        'source_group': source_record['invariants']['group_invariant_factors'],
        'target_group': target_record['invariants']['group_invariant_factors'],
        'source_order': str(math.prod(map(int, source_record['invariants']['group_invariant_factors']))),
        'target_order': str(math.prod(map(int, target_record['invariants']['group_invariant_factors']))),
        'symbolic_curve_equations_verified': True, 'symbolic_dual_compositions_verified': True,
        'homomorphism_pair_checks': pair_checks,
        'vectors': vectors, 'boundary': 'Constructive degree-two path between existing representatives; equal point counts do not imply isomorphic rational point groups; this does not compute a general isogeny volcano'}
    write('represented_isogeny.json', result)
    print(f'isogeny: {len(vectors)} point/dual vectors, {pair_checks} homomorphism checks', flush=True)


def consolidate(cases):
    """Collect complete mathematical diagnostics without fabricating timings."""
    providers = []
    for filename in ['algebra_buchberger.json', 'algebra_f5b.json', 'algebra_coverage_buchberger.json', 'algebra_coverage_f5b.json']:
        path = ROOT / 'build' / filename
        for index, row in enumerate(json.loads(path.read_text())['results']):
            providers.append((str(path.relative_to(ROOT)), f'/results/{index}', row))
    path = OUT / 'timeout_recovery.json'
    for index, row in enumerate(json.loads(path.read_text())['results']):
        providers.append((str(path.relative_to(ROOT)), f'/results/{index}/result', row['result']))
    by_id = defaultdict(list)
    for path, pointer, row in providers:
        by_id[row['case_id']].append((path, pointer, row))
    rows = []
    for case in cases:
        matches = by_id[case['id']]
        if any(row['case_sha256'] != digest(case) for path, pointer, row in matches):
            raise ValueError('Metric provider has the wrong case hash')
        completed = [(path, pointer, row) for path, pointer, row in matches
                     if row['status'] == 'ok' and row.get('relations', {}).get('truth_matches')]
        regularities = [(path, pointer, row) for path, pointer, row in matches
                        if row['degree_of_regularity'].get('value') is not None]
        if not completed or not regularities:
            raise ValueError('Incomplete metric coverage: ' + case['id'])
        if len({(row['degree_of_regularity']['value'], tuple(row['degree_of_regularity'].get('hilbert_function', []))) for path, pointer, row in regularities}) != 1:
            raise ValueError('Conflicting homogeneous diagnostics: ' + case['id'])
        if len({json.dumps(row['groebner_basis'], sort_keys=True) for path, pointer, row in completed}) != 1:
            raise ValueError('Conflicting reduced-basis metrics: ' + case['id'])
        bp, bi, basis = completed[0]
        rp, ri, regular = regularities[0]
        rows.append({'case_id': case['id'], 'case_sha256': digest(case),
            'groebner_basis_metrics': basis['groebner_basis'],
            'degree_of_regularity': regular['degree_of_regularity'],
            'ordered_relations': basis['relations']['ordered_relations'],
            'unordered_relations': basis['relations']['unordered_relations'],
            'basis_source': {'file': bp, 'json_pointer': bi},
            'regularity_source': {'file': rp, 'json_pointer': ri},
            'algorithm_maximum_working_degree': None})
    source_names = {name for row in rows for name in [row['basis_source']['file'], row['regularity_source']['file']]}
    write('complete_mathematical_metrics.json', {'schema': 'ecbench-complete-mathematical-metrics-v1',
        **stamp(), 'cases': len(rows), 'monomial_order': 'grevlex',
        'scope': 'Complete final-basis metrics and exact highest-input-part Hilbert diagnostics for original fixtures. Mathematical diagnostics are algorithm-independent. No method-specific timing is filled from another method.',
        'source_files_sha256': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(source_names)},
        'rows': rows})
    print(f'consolidate: all {len(rows)} original systems have verified final-basis metrics and exact regularity diagnostics', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase', choices=['derive', 'recover', 'compare', 'regularity', 'density', 'isogeny', 'consolidate', 'all'], default='all')
    parser.add_argument('--timeout', type=float, default=20)
    parser.add_argument('--repetitions', type=int, default=3)
    args = parser.parse_args()
    if args.repetitions < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        raise ValueError('Positive repetitions and finite positive timeout required')
    cases = catalog('algebra_cases')['cases']
    phases = ['derive', 'recover', 'compare', 'regularity', 'density', 'isogeny', 'consolidate'] if args.phase == 'all' else [args.phase]
    for phase in phases:
        if phase == 'derive': derive(cases)
        elif phase == 'recover': recover(cases, args.timeout)
        elif phase == 'compare': compare(cases, args.timeout, args.repetitions)
        elif phase == 'regularity': regularity(cases, args.timeout)
        elif phase == 'density': density(cases)
        elif phase == 'consolidate': consolidate(cases)
        else: isogeny()


if __name__ == '__main__': main()
