"""Integrity, proof, parameter, and known-answer validation."""
from pathlib import Path
import hashlib
import json
import math
import re
import time
from .reference import Field, Curve, ECError
DATA = Path(__file__).resolve().parent / 'data'
INTEGER = re.compile('(?:0|[1-9][0-9]*|-[1-9][0-9]*)\\Z')

def integer(value):
    if not isinstance(value, str) or not INTEGER.fullmatch(value):
        raise ValueError(f'Noncanonical decimal integer: {value!r}')
    return int(value)

def canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True) + '\n').encode()

def load_dataset(data_dir=None, bundle=None):
    if bundle:
        data = json.loads(Path(bundle).read_text())
    else:
        directory = Path(data_dir) if data_dir else DATA
        data = {name: json.loads((directory / (name + '.json')).read_text()) for name in ['manifest', 'fields', 'curves', 'vectors', 'relations', 'certificates']}
    for name, digest in data['manifest']['files'].items():
        if not hashlib.sha256(canonical_bytes(data[name[:-5]])).hexdigest() == digest:
            raise ValueError(f'Integrity mismatch: {name}')
    return data

def factorization(factors, primes):
    result = 1
    seen = set()
    for row in factors:
        p = integer(row['prime'])
        e = row['exponent']
        if not (p in primes and p not in seen and (type(e) is int) and (e > 0)):
            raise ValueError('Verification condition failed')
        seen.add(p)
        result *= p ** e
    return result

def certify(records):
    proofs = {integer(c['n']): c for c in records}
    if not len(proofs) == len(records):
        raise ValueError('Duplicate primality certificate')
    valid = set()
    active = set()

    def visit(n):
        if n in valid:
            return
        if not (n not in active and n in proofs):
            raise ValueError('Missing or cyclic certificate')
        c = proofs[n]
        if n == 2:
            if not (c['method'] == 'base' and (not c['n_minus_1_factorization'])):
                raise ValueError('Verification condition failed')
            valid.add(n)
            return
        if not (n > 2 and n % 2 == 1 and (c['method'] == 'lucas_complete_n_minus_1')):
            raise ValueError('Verification condition failed')
        active.add(n)
        for row in c['n_minus_1_factorization']:
            visit(integer(row['prime']))
        if not factorization(c['n_minus_1_factorization'], valid) == n - 1:
            raise ValueError('Verification condition failed')
        if not set(c['witnesses']) == {row['prime'] for row in c['n_minus_1_factorization']}:
            raise ValueError('Verification condition failed')
        for row in c['n_minus_1_factorization']:
            q = integer(row['prime'])
            a = integer(c['witnesses'][str(q)])
            if not (1 < a < n and pow(a, n - 1, n) == 1):
                raise ValueError('Verification condition failed')
            if not math.gcd(pow(a, (n - 1) // q, n) - 1, n) == 1:
                raise ValueError('Verification condition failed')
        active.remove(n)
        valid.add(n)
    for n in proofs:
        visit(n)
    return (valid, proofs)

def poly_mod(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << a.bit_length() - b.bit_length()
    return a

def poly_gcd(a, b):
    while b:
        a, b = (b, poly_mod(a, b))
    return a

def small_prime_divisors(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out

def validate_field(record, primes, proofs):
    f = Field(record)
    p = integer(record['characteristic'])
    degree = record['degree']
    if not (p in primes and integer(record['cardinality']) == p ** degree):
        raise ValueError('Verification condition failed')
    if not f.width == (f.bits + 7) // 8:
        raise ValueError('Verification condition failed')
    if not integer(record['representation']['modulus_integer']) == f.modulus:
        raise ValueError('Verification condition failed')
    if f.kind == 'prime':
        if not (degree == 1 and p > 3 and (f.q == f.modulus) and (f.bits == p.bit_length())):
            raise ValueError('Verification condition failed')
        if not record['proof']['certificate_n'] == str(p):
            raise ValueError('Verification condition failed')
        g = integer(record['representation']['multiplicative_primitive_root'])
        if not (1 < g < p and pow(g, p - 1, p) == 1):
            raise ValueError('Verification condition failed')
        for row in proofs[p]['n_minus_1_factorization']:
            if not pow(g, (p - 1) // integer(row['prime']), p) != 1:
                raise ValueError('Verification condition failed')
    else:
        if not (f.kind == 'binary' and p == 2 and (degree == f.bits) and (degree % 2 == 1)):
            raise ValueError('Verification condition failed')
        exps = record['representation']['polynomial_exponents']
        if not (exps == sorted(set(exps), reverse=True) and exps[0] == degree and (exps[-1] == 0)):
            raise ValueError('Verification condition failed')
        if not f.modulus == sum((1 << e for e in exps)):
            raise ValueError('Verification condition failed')
        if not f.parse(record['representation']['basis_generator']) == 2:
            raise ValueError('Verification condition failed')
        checks = {degree // ell for ell in small_prime_divisors(degree)}
        x = 2
        for i in range(1, degree + 1):
            x = f.square(x)
            if i in checks:
                if not poly_gcd(x ^ 2, f.modulus) == 1:
                    raise ValueError('Verification condition failed')
        if not x == 2:
            raise ValueError('Rabin irreducibility check failed')
    return f

def fundamental_discriminant(ring, primes):
    dk = integer(ring['field_discriminant'])
    fs = {integer(v['prime']): v['exponent'] for v in ring['field_discriminant_factorization']}
    if not (dk < 0 and factorization(ring['field_discriminant_factorization'], primes) == -dk):
        raise ValueError('Verification condition failed')
    if dk % 4 == 1:
        if not all((e == 1 for e in fs.values())):
            raise ValueError('Verification condition failed')
    else:
        if not (dk % 4 == 0 and dk // 4 % 4 in (2, 3)):
            raise ValueError('Verification condition failed')
        fs[2] = fs.get(2, 0) - 2
        if not fs[2] >= 0:
            raise ValueError('Verification condition failed')
        if fs[2] == 0:
            del fs[2]
        if not all((e == 1 for e in fs.values())):
            raise ValueError('Verification condition failed')
    return dk

def validate_curve(record, f, primes, proofs):
    c = Curve(record, f)
    co = c.coefficients
    inv = record['invariants']
    if f.kind == 'prime':
        if not (record['model'] == 'short_weierstrass' and co['a1'] == co['a2'] == co['a3'] == 0):
            raise ValueError('Verification condition failed')
        denominator = (4 * c.a ** 3 + 27 * c.b ** 2) % f.q
        if not denominator != 0:
            raise ValueError('Verification condition failed')
        delta = -16 * denominator % f.q
        j = 1728 * 4 * c.a ** 3 * pow(denominator, -1, f.q) % f.q
    else:
        if not (record['model'] == 'ordinary_binary' and co['a1'] == 1 and (co['a3'] == co['a4'] == 0) and (c.b != 0)):
            raise ValueError('Verification condition failed')
        delta = c.b
        j = f.inverse(c.b)
    if not (f.parse(inv['delta']) == delta and f.parse(inv['j']) == j):
        raise ValueError('Verification condition failed')
    n, t = (integer(inv['order']), integer(inv['trace']))
    if not (n == f.q + 1 - t and t * t <= 4 * f.q):
        raise ValueError('Verification condition failed')
    if not inv['frobenius_polynomial'] == ['1', str(-t), str(f.q)]:
        raise ValueError('Verification condition failed')
    if not factorization(inv['order_factorization'], primes) == n:
        raise ValueError('Verification condition failed')
    d1, d2 = map(integer, inv['group_invariant_factors'])
    if not (d1 * d2 == n and d2 % d1 == 0 and ((f.q - 1) % d1 == 0)):
        raise ValueError('Verification condition failed')
    if not (inv['group_cyclic'] == (d1 == 1) and integer(inv['group_exponent']) == d2):
        raise ValueError('Verification condition failed')
    if not inv['ordinary'] == (t % integer(f.record['characteristic']) != 0):
        raise ValueError('Verification condition failed')
    dpi = integer(inv['frobenius_discriminant'])
    if not dpi == t * t - 4 * f.q:
        raise ValueError('Verification condition failed')
    end = inv['endomorphism']
    ring = end['base_field']
    dk = fundamental_discriminant(ring, primes)
    fpi = integer(end['frobenius_order']['conductor'])
    fe = integer(ring['order_conductor'])
    if not (fpi > 0 and fe > 0 and (fpi % fe == 0) and (dpi == dk * fpi * fpi)):
        raise ValueError('Verification condition failed')
    if not integer(end['frobenius_order']['discriminant']) == dpi:
        raise ValueError('Verification condition failed')
    if not (integer(ring['order_discriminant']) == dk * fe * fe and ring['maximal'] == (fe == 1)):
        raise ValueError('Verification condition failed')
    if inv['ordinary']:
        if not (end['geometric'] == ring and ring['type'] == 'imaginary_quadratic_order'):
            raise ValueError('Verification condition failed')
    else:
        if not end['geometric']['type'] == 'maximal_quaternion_order':
            raise ValueError('Verification condition failed')
        if not integer(end['geometric']['ramified_prime']) == integer(f.record['characteristic']):
            raise ValueError('Verification condition failed')
    if ring['class_number'] is not None:
        if not (integer(ring['class_number']) == 1 and inv['ordinary'] and (dk * fe * fe in (-3, -4, -7))):
            raise ValueError('Verification condition failed')
    subgroup = record['subgroup']
    r = integer(subgroup['order'])
    h = integer(subgroup['cofactor'])
    if not (r in primes and r == max((integer(v['prime']) for v in inv['order_factorization'])) and (h * r == n)):
        raise ValueError('Verification condition failed')
    if not c.rank == (2 if d1 % r == 0 else 1):
        raise ValueError('Verification condition failed')
    if not subgroup['unique_order_r_subgroup'] == (c.rank == 1):
        raise ValueError('Verification condition failed')
    if not subgroup['cofactor_coprime_to_r'] == (math.gcd(h, r) == 1):
        raise ValueError('Verification condition failed')
    if not integer(subgroup['projection_scalar']) == d2 // r:
        raise ValueError('Verification condition failed')
    source = c.parse_point(subgroup['source_point'])
    if not (c.G is not None and c.mul(r, c.G) is None):
        raise ValueError('Verification condition failed')
    if not c.mul(d2 // r, source) == c.G:
        raise ValueError('Verification condition failed')
    emb = inv['subgroup_embedding_degrees']
    if not {v['prime'] for v in emb} == {v['prime'] for v in inv['order_factorization']}:
        raise ValueError('Verification condition failed')
    for row in emb:
        rr = integer(row['prime'])
        if rr == integer(f.record['characteristic']):
            if not row['degree'] is None:
                raise ValueError('Verification condition failed')
            continue
        k = integer(row['degree'])
        if not (k > 0 and (rr - 1) % k == 0 and (pow(f.q, k, rr) == 1)):
            raise ValueError('Verification condition failed')
        divisor_primes = {v['prime'] for v in proofs[rr]['n_minus_1_factorization'] if k % integer(v['prime']) == 0}
        if not set(row['minimality_prime_divisors']) == divisor_primes:
            raise ValueError('Verification condition failed')
        for ell in divisor_primes:
            if not pow(f.q, k // integer(ell), rr) != 1:
                raise ValueError('Verification condition failed')
        if rr == r:
            if not k == integer(subgroup['embedding_degree']):
                raise ValueError('Verification condition failed')
    for name, row in record['named_points'].items():
        p = c.parse_point(row['point'])
        order = integer(row['order'])
        if not (c.on_curve(p) and n % order == 0 and (factorization(row['order_factorization'], primes) == order)):
            raise ValueError('Verification condition failed')
        if not c.mul(order, p) is None:
            raise ValueError('Verification condition failed')
        for ell in row['order_factorization']:
            if not c.mul(order // integer(ell['prime']), p) is not None:
                raise ValueError('Verification condition failed')
    if not integer(record['named_points']['R']['order']) == d2:
        raise ValueError('Verification condition failed')
    return c

def evaluate(vector, fields, curves):
    inputs = vector['inputs']
    op = vector['operation']
    domain = vector['domain']
    if domain == 'field':
        f = fields[vector['field_id']]

        def operand(name):
            return integer(inputs['raw_' + name]) if 'raw_' + name in inputs else f.parse(inputs[name])
        if op == 'canonical':
            return {'boolean': f.canonical(integer(inputs['value']))}
        if op == 'reduce':
            out = f.reduce(integer(inputs['value']))
        elif op in ('add', 'sub', 'mul', 'div'):
            out = getattr(f, op)(operand('a'), operand('b'))
        elif op in ('pow', 'frobenius'):
            out = (f.power if op == 'pow' else f.frobenius)(operand('a'), integer(inputs['exponent']))
        else:
            out = getattr(f, op)(operand('a'))
        return {'integer': str(out)} if op == 'trace' else {'element': f.hex(out)}
    c = curves[vector['curve_id']]

    def point():
        if 'raw_point' in inputs:
            p = inputs['raw_point']
            return None if p is None else (integer(p['x']), integer(p['y']))
        return c.parse_point(inputs['P'])
    if domain == 'curve':
        p = point()
        if op == 'add':
            out = c.add(p, c.parse_point(inputs['Q']))
        elif op == 'mul':
            out = c.mul(integer(inputs['scalar']), p)
        else:
            out = getattr(c, op)(p)
        return {'point': c.point(out)}
    if domain == 'map':
        parameter = c.field.parse(inputs['parameter']) if 'parameter' in inputs else None
        out = c.map_point(op, point(), parameter)
        target = {k: c.field.parse(v) for k, v in inputs['target_coefficients'].items()}
        if op == 'degree2_isogeny':
            expected = {'a1': 0, 'a2': 0, 'a3': 0, 'a4': 4 % c.field.q, 'a6': 0}
        elif op == 'relative_frobenius':
            expected = {k: c.field.square(v) for k, v in c.coefficients.items()}
        else:
            expected = c.coefficients
        if not target == expected:
            raise ValueError('Map target coefficient mismatch')
        auxiliary = Curve({**c.record, 'coefficients': inputs['target_coefficients'], 'subgroup': {**c.record['subgroup'], 'generator': None}}, c.field)
        if not auxiliary.on_curve(out):
            raise ValueError('Map image is not on its target curve')
        return {'point': c.point(out)}
    if domain == 'membership':
        function = {'validate_point': c.on_curve, 'is_r_torsion': c.r_torsion, 'subgroup_contains': c.subgroup_contains, 'public_key_valid': c.public_key_valid}[op]
        return {'boolean': function(point())}
    if domain == 'encoding':
        if op == 'encode':
            return {'octets': c.encode(point(), inputs['format'] == 'compressed')}
        if op == 'decode':
            return {'point': c.point(c.decode(inputs['octets']))}
    raise ValueError(f'Unsupported vector operation {domain}.{op}')

def verify(data, with_pari=False, gp_path=None):
    started = time.perf_counter()
    m = data['manifest']
    counts = m['counts']
    if not (m['schema_version'] == '1.0.0' and m['dataset_id'] == 'ecbench-golden-smallfields-v1'):
        raise ValueError('Verification condition failed')
    primes, proofs = certify(data['certificates'])
    fields = {f['id']: validate_field(f, primes, proofs) for f in data['fields']}
    if not len(fields) == len(data['fields']) == counts['fields'] == 7:
        raise ValueError('Verification condition failed')
    curves = {c['id']: validate_curve(c, fields[c['field_id']], primes, proofs) for c in data['curves']}
    if not len(curves) == len(data['curves']) == counts['curves'] == 56:
        raise ValueError('Verification condition failed')
    for c in data['curves']:
        digest = hashlib.sha256(json.dumps(c, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        if not digest == m['record_hashes'][c['id']]:
            raise ValueError('Verification condition failed')
        cc = curves[c['id']]
        partner = curves[c['twist']['partner_id']]
        f = cc.field
        if not partner.field.record['id'] == f.record['id']:
            raise ValueError('Verification condition failed')
        if not c['invariants']['j'] == partner.record['invariants']['j']:
            raise ValueError('Verification condition failed')
        if not integer(c['invariants']['trace']) == -integer(partner.record['invariants']['trace']):
            raise ValueError('Verification condition failed')
        d = f.parse(c['twist']['parameter'])
        if f.kind == 'binary':
            if not (f.trace(d) == 1 and partner.a == cc.a ^ d and (partner.b == cc.b)):
                raise ValueError('Verification condition failed')
        else:
            if not pow(d, (f.q - 1) // 2, f.q) == f.q - 1:
                raise ValueError('Verification condition failed')
            scale = f.parse(c['twist']['representative_isomorphism_scale'])
            if not f.mul(f.power(scale, 4), partner.a) == f.mul(f.power(d, 2), cc.a):
                raise ValueError('Verification condition failed')
            if not f.mul(f.power(scale, 6), partner.b) == f.mul(f.power(d, 3), cc.b):
                raise ValueError('Verification condition failed')
    ids = set()
    for v in data['vectors']:
        if not v['id'] not in ids:
            raise ValueError('Verification condition failed')
        ids.add(v['id'])
        try:
            result = evaluate(v, fields, curves)
        except ECError as error:
            result = {'error': error.code}
        if not result == v['expected']:
            raise ValueError(f"Vector {v['id']} failed: got {result!r}, expected {v['expected']!r}")
    if not (len(ids) == counts['vectors'] and len(primes) == counts['primality_certificates']):
        raise ValueError('Verification condition failed')
    if not sum((len(c['named_points']) for c in data['curves'])) == counts['named_points']:
        raise ValueError('Verification condition failed')
    classes = data['relations']['isogeny_classes']
    assigned = set()
    for row in classes:
        for cid in row['curve_ids']:
            if not cid not in assigned:
                raise ValueError('Verification condition failed')
            assigned.add(cid)
            c = curves[cid]
            if not (c.field.record['id'] == row['field_id'] and c.record['invariants']['trace'] == row['trace']):
                raise ValueError('Verification condition failed')
    if not (assigned == set(curves) and len(classes) == counts['represented_isogeny_classes']):
        raise ValueError('Verification condition failed')
    assigned = set()
    for row in data['relations']['geometric_classes_over_base_field']:
        for cid in row['curve_ids']:
            if not cid not in assigned:
                raise ValueError('Verification condition failed')
            assigned.add(cid)
            c = curves[cid]
            if not (c.field.record['id'] == row['field_id'] and c.record['invariants']['j'] == row['j']):
                raise ValueError('Verification condition failed')
    if not assigned == set(curves):
        raise ValueError('Verification condition failed')
    assigned = set()
    for row in data['relations']['ordinary_endomorphism_order_classes']:
        for cid in row['curve_ids']:
            if not cid not in assigned:
                raise ValueError('Verification condition failed')
            assigned.add(cid)
            c = curves[cid]
            ring = c.record['invariants']['endomorphism']['geometric']
            if not (c.record['invariants']['ordinary'] and ring['field_discriminant'] == row['field_discriminant'] and (ring['order_conductor'] == row['conductor'])):
                raise ValueError('Verification condition failed')
    if not assigned == {cid for cid, c in curves.items() if c.record['invariants']['ordinary']}:
        raise ValueError('Verification condition failed')
    expected_pairs = {tuple(sorted([c['id'], c['twist']['partner_id']])) for c in data['curves']}
    actual_pairs = {tuple(row['curve_ids']) for row in data['relations']['twist_pairs']}
    if not (expected_pairs == actual_pairs and len(actual_pairs) == len(data['relations']['twist_pairs'])):
        raise ValueError('Verification condition failed')
    result = {'status': 'pass', 'dataset_id': m['dataset_id'], 'fields': len(fields), 'curves': len(curves), 'vectors': len(ids), 'primality_certificates': len(primes), 'named_points': counts['named_points'], 'point_counts': 'stored_verified_values; not independently recomputed in core mode',
              'endomorphism_orders':'source-proven metadata; discriminant consistency checked, model-specific ring proofs not recomputed',
              'class_numbers':'known small endomorphism-order class numbers checked; remaining evidence retained as source metadata',
              'elapsed_seconds': round(time.perf_counter() - started, 4)}
    if with_pari:
        from .pari import verify_with_pari
        result['pari'] = verify_with_pari(data['fields'], data['curves'], data['vectors'], gp_path)
        result['point_counts'] = 'independently recomputed with PARI/GP'
    return result
