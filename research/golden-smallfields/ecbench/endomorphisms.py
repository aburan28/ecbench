"""Verified experimental GLV2 and binary Frobenius scalar decompositions.

The primary curve corpus is never changed. Maps act on the selected cyclic
subgroup only after their scalar action is checked on its generator. Rational
r-torsion membership alone does not establish that this subgroup is stable.
All arithmetic is variable-time reference code, unsuitable for secrets.
"""
from __future__ import annotations

from collections import Counter
import datetime
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import time

from .extension import babai_decompose
from .reference import Curve, Field

SCHEMA = 'ecbench-endomorphism-glv2-cases-v1'
DEFAULT_CASES = Path(__file__).resolve().parent / 'data' / 'endomorphism_cases.json'


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _point(point):
    return None if point is None else {'x': str(point[0]), 'y': str(point[1])}


def _parse_point(point, field):
    if point is None:
        return None
    _require(isinstance(point, dict) and set(point) == {'x', 'y'}, 'Invalid decimal point')
    result = tuple(_integer(point[key]) for key in ('x', 'y'))
    field.check_elements(*result)
    return result


def _integer(value):
    _require(type(value) is str, 'Mathematical integers must be decimal strings')
    result = int(value)
    _require(str(result) == value, 'Noncanonical decimal integer')
    return result


def _source_hashes(data):
    return {name + '.json': hashlib.sha256(
        (json.dumps(data[name], sort_keys=True, indent=2, ensure_ascii=True) + '\n').encode()).hexdigest()
        for name in ('fields', 'curves')}


def _sqrt_mod(value, prime):
    """Tonelli-Shanks, returning the canonical smaller root or None."""
    value %= prime
    if prime == 2 or value == 0:
        return value
    if pow(value, (prime - 1) // 2, prime) != 1:
        return None
    if prime % 4 == 3:
        root = pow(value, (prime + 1) // 4, prime)
        return min(root, prime - root)
    odd, exponent = prime - 1, 0
    while odd % 2 == 0:
        odd //= 2
        exponent += 1
    nonresidue = 2
    while pow(nonresidue, (prime - 1) // 2, prime) != prime - 1:
        nonresidue += 1
    c = pow(nonresidue, odd, prime)
    root = pow(value, (odd + 1) // 2, prime)
    residue, order = pow(value, odd, prime), exponent
    while residue != 1:
        power, trial = 0, residue
        while trial != 1:
            trial = trial * trial % prime
            power += 1
            _require(power < order, 'Tonelli-Shanks invariant failed')
        unit = pow(c, 1 << (order - power - 1), prime)
        root = root * unit % prime
        c = unit * unit % prime
        residue = residue * c % prime
        order = power
    _require(root * root % prime == value, 'Invalid square root')
    return min(root, prime - root)


def _quadratic_roots(linear, constant, prime):
    if prime == 2:
        return [root for root in range(2) if (root * root + linear * root + constant) % 2 == 0]
    discriminant_root = _sqrt_mod(linear * linear - 4 * constant, prime)
    if discriminant_root is None:
        return []
    inverse_two = pow(2, -1, prime)
    return sorted({(-linear + sign * discriminant_root) * inverse_two % prime for sign in (-1, 1)})


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _nearest(value):
    return (2 * value.numerator + value.denominator) // (2 * value.denominator)


def _gauss_reduce(order, eigenvalue):
    """Exact Lagrange-Gauss reduction of the two-dimensional kernel lattice."""
    first, second = [order, 0], [-eigenvalue, 1]
    while True:
        if _dot(second, second) < _dot(first, first):
            first, second = second, first
        coefficient = _nearest(Fraction(_dot(first, second), _dot(first, first)))
        if coefficient == 0:
            return [first, second]
        second = [b - coefficient * a for a, b in zip(first, second)]


def _component_bounds(basis):
    """Exact nearest-plane residual bounds, derived from Gram-Schmidt."""
    first, second = basis
    star1 = list(map(Fraction, first))
    mu = Fraction(_dot(first, second), _dot(first, first))
    star2 = [Fraction(y) - mu * x for x, y in zip(first, second)]
    bounds = [(abs(star1[j]) + abs(star2[j])) / 2 for j in range(2)]
    ceilings = [(value.numerator + value.denominator - 1) // value.denominator for value in bounds]
    return {'coordinate_absolute_bounds': [
        {'numerator': str(value.numerator), 'denominator': str(value.denominator)} for value in bounds],
        'coordinate_absolute_bound_ceilings': list(map(str, ceilings)),
        'guaranteed_component_bits': str(max(value.bit_length() for value in ceilings))}


def apply_endomorphism(curve, description, point):
    parameter = description.get('parameter')
    return curve.map_point(description['operation'], point, None if parameter is None else _integer(parameter))


def joint_mul(curve, scalars, points):
    """Signed interleaved binary multiplication, with no implicit order reduction."""
    _require(len(scalars) == len(points), 'Scalar/point count mismatch')
    prepared = []
    for scalar, point in zip(scalars, points):
        _require(type(scalar) is int, 'Invalid signed scalar')
        curve.check(point)
        prepared.append((abs(scalar), curve._neg(point) if scalar < 0 else point))
    result = None
    for bit in reversed(range(max((scalar.bit_length() for scalar, _ in prepared), default=0))):
        result = curve._add(result, result)
        for scalar, point in prepared:
            if scalar >> bit & 1:
                result = curve._add(result, point)
    return result


def scalar_mul(curve, description, basis, scalar):
    components = babai_decompose(scalar, curve.r, basis)
    image = apply_endomorphism(curve, description, curve.G)
    return joint_mul(curve, components, [curve.G, image])


def _map_spec(curve):
    f = curve.field
    if f.kind == 'binary':
        if curve.a not in (0, 1) or curve.b != 1:
            return None, None, 'relative_frobenius_targets_distinct_curve'
        return {'operation': 'relative_frobenius', 'parameter': None,
                'formula': 'tau(x,y)=(x^2,y^2)', 'base_field': 'F2'}, [1, -(2 * curve.a - 1), 2], None
    if curve.record['family'] == 'E4':
        if f.q % 4 != 1:
            return None, None, 'gaussian_unit_not_in_base_field'
        unit = f.sqrt(f.q - 1)
        return {'operation': 'gaussian_automorphism', 'parameter': str(unit),
                'formula': 'phi(x,y)=(-x,i*y); i^2=-1',
                'unit_selection': 'canonical smaller square root of -1'}, [1, 0, 1], None
    if curve.record['family'] == 'E3':
        if f.q % 3 != 1:
            return None, None, 'eisenstein_unit_not_in_base_field'
        primitive = int(f.record['representation']['multiplicative_primitive_root'])
        unit = pow(primitive, (f.q - 1) // 3, f.q)
        return {'operation': 'eisenstein_automorphism', 'parameter': str(unit),
                'formula': 'phi(x,y)=(zeta*x,y); zeta^2+zeta+1=0',
                'unit_selection': 'recorded primitive root raised to (p-1)/3'}, [1, 1, 1], None
    return None, None, 'no_explicit_nontrivial_endomorphism_in_this_battery'


def _map_relation(curve, description, polynomial, point):
    image = apply_endomorphism(curve, description, point)
    second = apply_endomorphism(curve, description, image)
    return curve.add(curve.add(second, curve.mul(polynomial[1], image)),
                     curve.mul(polynomial[2], point)) is None


def _scalars(order, curve_id):
    values = [('zero', 0), ('one', 1), ('negative_one', -1), ('r_minus_one', order - 1),
              ('r', order), ('r_plus_one', order + 1), ('fixed_80_bit', (1 << 79) + 0x12345),
              ('fixed_128_bit', (1 << 127) + 0x123456789abcdef)]
    values += [('seeded_subgroup_' + str(index), int.from_bytes(hashlib.sha256(
        f'ecbench-glv2-v1:{curve_id}:{index}'.encode()).digest(), 'big') % order) for index in range(8)]
    return values


def build_cases(data, binary_provenance):
    """Generate deterministic cases using only the independent reference oracle."""
    fields = {record['id']: Field(record) for record in data['fields']}
    provenance = {record['label']: record for record in binary_provenance['curves']}
    cases = []
    for record in data['curves']:
        curve = Curve(record, fields[record['field_id']])
        description, polynomial, reason = _map_spec(curve)
        row = {'id': record['id'] + '-endomorphism', 'curve_id': record['id'], 'field_id': record['field_id'],
               'subgroup_order': str(curve.r), 'r_torsion_rank': str(curve.rank), 'generator': _point(curve.G)}
        if reason is not None:
            row.update(status='unsupported', reason=reason)
            cases.append(row)
            continue
        image = apply_endomorphism(curve, description, curve.G)
        roots = _quadratic_roots(polynomial[1], polynomial[2], curve.r)
        eigenvalues = [root for root in roots if curve.mul(root, curve.G) == image]
        row.update(map=description, characteristic_polynomial=list(map(str, polynomial)),
                   characteristic_polynomial_coefficient_order='descending',
                   characteristic_roots_mod_r=list(map(str, roots)), image_of_generator=_point(image),
                   map_image_on_curve=curve.on_curve(image), map_image_is_r_torsion=curve.r_torsion(image),
                   selected_subgroup_is_stable=curve.subgroup_contains(image))
        if not eigenvalues:
            _require(not row['selected_subgroup_is_stable'], 'Stable cyclic subgroup has no verified eigenvalue')
            row.update(status='unsupported', reason='selected_cyclic_subgroup_not_stable')
            cases.append(row)
            continue
        _require(len(eigenvalues) == 1, 'Nonunique scalar action on the selected generator')
        eigenvalue = eigenvalues[0]
        if curve.field.kind == 'binary':
            source = provenance[record['id']]
            _require(int(source['largest_prime_subgroup']) == curve.r, 'Binary provenance subgroup mismatch')
            _require(int(source['F2_frobenius_scalar_on_largest_subgroup']) == eigenvalue, 'Binary provenance eigenvalue mismatch')
            row['eigenvalue_source'] = 'provenance/binary_data.json:F2_frobenius_scalar_on_largest_subgroup'
            row['frobenius_iteration_degree'] = str(curve.field.bits)
        else:
            row['eigenvalue_source'] = 'root of characteristic polynomial verified by [lambda]G=phi(G)'
        basis = _gauss_reduce(curve.r, eigenvalue)
        determinant = basis[0][0] * basis[1][1] - basis[0][1] * basis[1][0]
        _require(abs(determinant) == curve.r, 'Kernel lattice determinant mismatch')
        vectors = []
        for name, scalar in _scalars(curve.r, record['id']):
            components = babai_decompose(scalar, curve.r, basis)
            _require((components[0] + eigenvalue * components[1] - scalar) % curve.r == 0, 'Scalar congruence failed')
            expected = curve.mul(scalar, curve.G)
            signed = curve.add(curve.mul(components[0], curve.G), curve.mul(components[1], image))
            _require(signed == expected == joint_mul(curve, components, [curve.G, image]), 'Signed reconstruction failed')
            vectors.append({'name': name, 'scalar': str(scalar), 'scalar_mod_r': str(scalar % curve.r),
                            'components': list(map(str, components)),
                            'max_component_bits': str(max(abs(component).bit_length() for component in components)),
                            'expected_point': _point(expected)})
        row.update(status='supported', decomposition='GLV2' if curve.field.kind == 'prime' else 'binary-Frobenius-two-term',
                   eigenvalue=str(eigenvalue), dimension='2', lattice_seed=[[str(curve.r), '0'], [str(-eigenvalue), '1']],
                   lattice_basis=[[str(value) for value in vector] for vector in basis], lattice_determinant=str(determinant),
                   lattice_reduction='exact Lagrange-Gauss', rounding='nearest integer; ties toward +infinity',
                   component_bounds=_component_bounds(basis), scalar_vectors=vectors,
                   observed_max_component_bits=str(max(int(vector['max_component_bits']) for vector in vectors)))
        cases.append(row)
    counts = Counter(row['status'] for row in cases)
    return {'schema': SCHEMA, 'source_dataset_id': data['manifest']['dataset_id'],
            'source_hashes': _source_hashes(data), 'integer_encoding': 'canonical decimal strings',
            'primary_dataset_modified': False, 'counts': {'curves': str(len(cases)),
            'supported': str(counts['supported']), 'unsupported': str(counts['unsupported']),
            'binary_frobenius': str(sum(row.get('decomposition') == 'binary-Frobenius-two-term' for row in cases))},
            'boundaries': ['All maps and decompositions are checked on the selected cyclic subgroup.',
                           'Rational r-torsion membership alone does not prove membership in the selected cyclic subgroup.',
                           'Relative Frobenius on binary curves with coefficients outside F2 is not an endomorphism of the same curve.',
                           'This battery has two-term decompositions only; powers of Frobenius do not establish a GLV4 construction.',
                           'The genuine supplementary quadratic-extension GLV+GLS case is a separate artifact.',
                           'Timings use variable-time Python reference arithmetic and provide no production-security claim.'],
            'cases': cases}


def _load_cases(cases):
    if cases is None:
        cases = DEFAULT_CASES
    return json.loads(Path(cases).read_text()) if isinstance(cases, (str, Path)) else cases


def verify_cases(data, cases=None):
    """Independently verify maps, stable-subgroup evidence, lattices and vectors."""
    artifact = _load_cases(cases)
    _require(artifact['schema'] == SCHEMA, 'Unsupported endomorphism artifact schema')
    _require(artifact['source_hashes'] == _source_hashes(data), 'Endomorphism artifact source hash mismatch')
    _require(artifact['source_dataset_id'] == data['manifest']['dataset_id'], 'Dataset identity mismatch')
    fields = {record['id']: Field(record) for record in data['fields']}
    records = {record['id']: record for record in data['curves']}
    _require(len(artifact['cases']) == len(records), 'Endomorphism classification must cover every curve')
    _require({row['curve_id'] for row in artifact['cases']} == set(records), 'Missing or duplicate curve classification')
    counts = Counter()
    vectors_checked = 0
    for row in artifact['cases']:
        record = records[row['curve_id']]
        curve = Curve(record, fields[record['field_id']])
        _require(row['field_id'] == record['field_id'], 'Field identity mismatch')
        _require(_integer(row['subgroup_order']) == curve.r and _integer(row['r_torsion_rank']) == curve.rank, 'Subgroup parameters mismatch')
        _require(_parse_point(row['generator'], curve.field) == curve.G, 'Generator mismatch')
        description, polynomial, reason = _map_spec(curve)
        counts[row['status']] += 1
        if reason is not None:
            _require(row['status'] == 'unsupported' and row['reason'] == reason, 'Incorrect unavailable-map classification')
            continue
        _require(row['map'] == description, 'Endomorphism map parameter mismatch')
        _require([_integer(value) for value in row['characteristic_polynomial']] == polynomial, 'Characteristic polynomial mismatch')
        roots = _quadratic_roots(polynomial[1], polynomial[2], curve.r)
        _require([_integer(value) for value in row['characteristic_roots_mod_r']] == roots, 'Incorrect characteristic roots')
        image = apply_endomorphism(curve, description, curve.G)
        _require(_parse_point(row['image_of_generator'], curve.field) == image, 'Generator image mismatch')
        _require(row['map_image_on_curve'] == curve.on_curve(image), 'Image curve membership mismatch')
        _require(row['map_image_is_r_torsion'] == curve.r_torsion(image), 'Image torsion membership mismatch')
        stable = curve.subgroup_contains(image)
        _require(row['selected_subgroup_is_stable'] == stable, 'Selected-subgroup stability mismatch')
        matching = [root for root in roots if curve.mul(root, curve.G) == image]
        for point in [None, curve.G, curve.parse_point(record['named_points']['R']['point'])]:
            _require(_map_relation(curve, description, polynomial, point), 'Endomorphism polynomial identity failed')
        if not matching:
            _require(row['status'] == 'unsupported' and row['reason'] == 'selected_cyclic_subgroup_not_stable' and not stable,
                     'Incorrect unstable-subgroup classification')
            continue
        _require(row['status'] == 'supported' and len(matching) == 1, 'Incorrect supported classification')
        eigenvalue = _integer(row['eigenvalue'])
        _require(eigenvalue == matching[0] and 0 <= eigenvalue < curve.r, 'Incorrect selected eigenvalue')
        if curve.field.kind == 'binary':
            _require(_integer(row['frobenius_iteration_degree']) == curve.field.bits, 'Incorrect Frobenius iteration degree')
            repeated = curve.G
            for _ in range(curve.field.bits):
                repeated = apply_endomorphism(curve, description, repeated)
            _require(repeated == curve.G, 'tau^m(G) != G')
            counts['binary_frobenius'] += 1
        _require(_integer(row['dimension']) == 2, 'Incorrect decomposition dimension')
        seed = [[_integer(value) for value in vector] for vector in row['lattice_seed']]
        _require(seed == [[curve.r, 0], [-eigenvalue, 1]], 'Incorrect kernel lattice seed')
        basis = [[_integer(value) for value in vector] for vector in row['lattice_basis']]
        _require(len(basis) == 2 and all(len(vector) == 2 for vector in basis), 'Incorrect lattice shape')
        determinant = basis[0][0] * basis[1][1] - basis[0][1] * basis[1][0]
        _require(abs(determinant) == curve.r and _integer(row['lattice_determinant']) == determinant, 'Incorrect lattice determinant')
        _require(all((vector[0] + eigenvalue * vector[1]) % curve.r == 0 for vector in basis), 'Lattice row outside kernel')
        first_norm, second_norm = _dot(basis[0], basis[0]), _dot(basis[1], basis[1])
        _require(first_norm <= second_norm and 2 * abs(_dot(basis[0], basis[1])) <= first_norm,
                 'Kernel lattice basis is not Gauss reduced')
        _require(row['component_bounds'] == _component_bounds(basis), 'Incorrect nearest-plane component bounds')
        _require([(vector['name'], _integer(vector['scalar'])) for vector in row['scalar_vectors']] == _scalars(curve.r, row['curve_id']), 'Scalar-vector battery mismatch')
        observed = 0
        bound_ceilings = list(map(_integer, row['component_bounds']['coordinate_absolute_bound_ceilings']))
        for vector in row['scalar_vectors']:
            scalar = _integer(vector['scalar'])
            components = list(map(_integer, vector['components']))
            _require(components == babai_decompose(scalar, curve.r, basis), 'Incorrect scalar decomposition')
            _require(_integer(vector['scalar_mod_r']) == scalar % curve.r, 'Incorrect scalar residue')
            _require((components[0] + eigenvalue * components[1] - scalar) % curve.r == 0, 'Scalar congruence failed')
            _require(all(abs(value) <= bound for value, bound in zip(components, bound_ceilings)), 'Component exceeds proven bound')
            expected = curve.mul(scalar, curve.G)
            signed = curve.add(curve.mul(components[0], curve.G), curve.mul(components[1], image))
            _require(signed == expected == joint_mul(curve, components, [curve.G, image]), 'Signed point reconstruction failed')
            _require(_parse_point(vector['expected_point'], curve.field) == expected, 'Expected reconstruction point mismatch')
            bits = max(abs(value).bit_length() for value in components)
            _require(_integer(vector['max_component_bits']) == bits, 'Incorrect observed component bits')
            observed = max(observed, bits)
            vectors_checked += 1
        _require(_integer(row['observed_max_component_bits']) == observed, 'Incorrect observed maximum component bits')
    _require(artifact['counts'] == {'curves': str(len(records)), 'supported': str(counts['supported']),
             'unsupported': str(counts['unsupported']), 'binary_frobenius': str(counts['binary_frobenius'])}, 'Case counts mismatch')
    return {'status': 'pass', 'schema': SCHEMA, 'curves_classified': len(records),
            'supported_cases': counts['supported'], 'unsupported_cases': counts['unsupported'],
            'binary_frobenius_cases': counts['binary_frobenius'], 'scalar_vectors_checked': vectors_checked,
            'checks': ['canonical map units', 'characteristic roots', 'selected cyclic subgroup stability',
                       'endomorphism polynomial identities', 'binary tau^m', 'kernel lattice determinant',
                       'Gauss-reduced basis', 'exact nearest-plane bounds', 'signed scalar congruence and point reconstruction']}


def benchmark_cases(data, cases=None, ids=None, samples=3, iterations=1):
    """Timing helper; expensive fixture verification is outside measured calls."""
    _require(samples > 0 and iterations > 0, 'Samples and iterations must be positive')
    artifact = _load_cases(cases)
    preflight = verify_cases(data, artifact)
    fields = {record['id']: Field(record) for record in data['fields']}
    records = {record['id']: record for record in data['curves']}
    selected = set(ids or [])
    _require(not selected - set(records), 'Unknown curve ID in endomorphism benchmark')
    results = []
    for row in artifact['cases']:
        if row['status'] != 'supported' or (selected and row['curve_id'] not in selected):
            continue
        record = records[row['curve_id']]
        curve = Curve(record, fields[record['field_id']])
        basis = [[_integer(value) for value in vector] for vector in row['lattice_basis']]
        scalar = ((1 << 127) + 0x123456789abcdef) % curve.r or 1
        components = babai_decompose(scalar, curve.r, basis)
        image = apply_endomorphism(curve, row['map'], curve.G)
        expected = curve.mul(scalar, curve.G)
        _require(scalar_mul(curve, row['map'], basis, scalar) == expected, 'Benchmark reconstruction mismatch')
        operations = [
            ('endomorphism', lambda: apply_endomorphism(curve, row['map'], curve.G)),
            ('scalar_decompose', lambda: babai_decompose(scalar, curve.r, basis)),
            ('scalar_mul_subgroup_baseline', lambda: curve.mul(scalar, curve.G)),
            ('glv2_joint_mul_precomputed', lambda: joint_mul(curve, components, [curve.G, image])),
            ('glv2_scalar_mul_full', lambda: scalar_mul(curve, row['map'], basis, scalar)),
        ]
        for name, operation in operations:
            operation()
            elapsed = []
            for _ in range(samples):
                started = time.perf_counter_ns()
                for _ in range(iterations):
                    operation()
                elapsed.append((time.perf_counter_ns() - started) / iterations)
            ordered = sorted(elapsed)
            results.append({'curve_id': row['curve_id'], 'field_id': row['field_id'], 'operation': name,
                'scalar_mod_r': str(scalar), 'components': list(map(str, components)),
                'samples': samples, 'iterations_per_sample': iterations, 'unit': 'nanoseconds_per_operation',
                'min': min(elapsed), 'median': statistics.median(elapsed), 'mean': statistics.mean(elapsed),
                'p95': ordered[max(0, math.ceil(.95 * len(ordered)) - 1)], 'sample_values': elapsed})
    return {'schema': SCHEMA, 'preflight': preflight, 'implementation': 'variable-time Python reference',
            'dataset_id': data['manifest']['dataset_id'], 'fixture_source_hashes': artifact['source_hashes'],
            'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'environment': {'python': platform.python_version(), 'python_implementation': platform.python_implementation(),
                            'platform': platform.platform(), 'machine': platform.machine(), 'processor': platform.processor()},
            'timing_scope': 'arithmetic call only; dataset loading and fixture verification excluded',
            'precomputed_scope': 'glv2_joint_mul_precomputed excludes decomposition and endomorphism; glv2_scalar_mul_full includes both',
            'results': results}

