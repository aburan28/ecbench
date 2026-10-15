"""Proof-carrying GF(2) row reduction of exact polynomial generators.

Only constant-coefficient linear combinations are used. No monomial
multiples, Boolean power reductions, root guesses, or point relations are
added. Equality of row spaces proves equality of the polynomial ideals.
"""
import copy
import hashlib
import json

from .algebra import _normalize_case, _record_input


def digest(case):
    return hashlib.sha256(json.dumps(case, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _input(case):
    names, modulus, polynomials, domains, metrics = _record_input(_normalize_case(case))
    if modulus != 2:
        raise ValueError('GF(2) row reduction requires coefficient characteristic two')
    return names, polynomials


def _column_key(exponents):
    return sum(exponents), tuple(-value for value in reversed(exponents))


def _encode(polynomials, columns):
    positions = {exponents: i for i, exponents in enumerate(columns)}
    return [sum(1 << positions[exponents] for exponents in polynomial) for polynomial in polynomials]


def _echelon(rows):
    pivots = {}
    for index, row in enumerate(rows):
        witness = 1 << index
        while row:
            pivot = (row & -row).bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = [row, witness]
                break
            row ^= pivots[pivot][0]
            witness ^= pivots[pivot][1]
    for pivot in sorted(pivots, reverse=True):
        for earlier in pivots:
            if earlier < pivot and (pivots[earlier][0] >> pivot) & 1:
                pivots[earlier][0] ^= pivots[pivot][0]
                pivots[earlier][1] ^= pivots[pivot][1]
    return pivots


def derive_case(case):
    """Return an equivalent fixture and exact XOR-combination witnesses."""
    names, polynomials = _input(case)
    columns = sorted(set().union(*(set(poly) for poly in polynomials)), key=_column_key, reverse=True)
    encoded = _encode(polynomials, columns)
    pivots = _echelon(encoded)
    # A discovered unit goes first so solvers can immediately detect it.
    rows = sorted(pivots.values(), key=lambda row: (0 if row[0] and all(
        not any(columns[i]) for i in range(len(columns)) if row[0] >> i & 1) else 1,
        (row[0] & -row[0]).bit_length()))
    output = []
    for index, (row, witness) in enumerate(rows):
        output.append({'id': f'gf2-rref-{index}', 'role': 'gf2_row_reduced_generator',
            'terms': [{'coefficient': '1', 'exponents': list(exponents)}
                      for i, exponents in enumerate(columns) if row >> i & 1]})
    derived = copy.deepcopy(case)
    derived['id'] = case['id'] + '-gf2-rref'
    derived['polynomials'] = output
    proof = {'schema': 'ecbench-gf2-row-equivalence-v1', 'source_case_id': case['id'],
        'source_case_sha256': digest(case), 'derived_case_sha256': digest(derived),
        'column_order': 'descending grevlex on the exact recorded monomials',
        'input_generator_count': len(polynomials), 'output_generator_count': len(output),
        'monomial_columns': len(columns),
        'output_xor_input_masks': [hex(witness) for row, witness in rows],
        'new_linear_generators': sum(max((sum(term['exponents']) for term in poly['terms']), default=0) == 1 for poly in output),
        'unit_generator_found': any(poly['terms'] == [{'coefficient': '1', 'exponents': [0] * len(names)}] for poly in output),
        'proof_scope': 'Both inclusions of constant-coefficient polynomial row spans, hence equality of generated affine polynomial ideals; highest-part ideals can change'}
    return derived, proof


def verify_transformation(original, derived, proof):
    """Reconstruct each witness and check reverse row-span inclusion."""
    if proof.get('schema') != 'ecbench-gf2-row-equivalence-v1':
        raise ValueError('Unsupported row-equivalence proof')
    if proof['source_case_sha256'] != digest(original) or proof['derived_case_sha256'] != digest(derived):
        raise ValueError('Preprocessing case hash mismatch')
    for key, value in original.items():
        if key not in {'id', 'polynomials'} and derived.get(key) != value:
            raise ValueError('Preprocessing changed non-polynomial fixture data')
    if derived['id'] != original['id'] + '-gf2-rref':
        raise ValueError('Unexpected derived case identity')
    names, inputs = _input(original)
    derived_names, outputs = _input(derived)
    if names != derived_names:
        raise ValueError('Preprocessing changed variables')
    columns = sorted(set().union(*(set(poly) for poly in inputs + outputs)), key=_column_key, reverse=True)
    source_rows = _encode(inputs, columns)
    output_rows = _encode(outputs, columns)
    actual_linear = sum(max(map(sum, poly), default=0) == 1 for poly in outputs)
    actual_unit = any(poly == {(0,) * len(names): 1} for poly in outputs)
    if (proof.get('source_case_id') != original['id'] or
            proof.get('input_generator_count') != len(inputs) or
            proof.get('output_generator_count') != len(outputs) or
            proof.get('monomial_columns') != len(columns) or
            proof.get('new_linear_generators') != actual_linear or
            proof.get('unit_generator_found') != actual_unit):
        raise ValueError('Preprocessing metadata does not match actual generators')
    witnesses = proof['output_xor_input_masks']
    if len(witnesses) != len(outputs):
        raise ValueError('Missing row-equivalence witness')
    for expected, encoded_witness in zip(output_rows, witnesses):
        witness = int(encoded_witness, 16)
        if witness < 0 or witness.bit_length() > len(inputs):
            raise ValueError('Invalid row-equivalence witness mask')
        actual = 0
        for i, row in enumerate(source_rows):
            if witness >> i & 1:
                actual ^= row
        if actual != expected:
            raise ValueError('Derived generator does not match its source combination')
    echelon = _echelon(output_rows)
    for row in source_rows:
        while row:
            pivot = (row & -row).bit_length() - 1
            if pivot not in echelon:
                raise ValueError('Original generator is outside derived row span')
            row ^= echelon[pivot][0]
    return {'status': 'pass', 'source_case_id': original['id'],
        'input_generators': len(inputs), 'output_generators': len(outputs),
        'proof': 'Exact row-span equality checked in both directions; affine ideals equal'}
