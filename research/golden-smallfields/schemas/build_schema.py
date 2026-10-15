"""Build strict Draft 2020-12 release schema, retaining semantic check boundaries."""
from pathlib import Path
import json
import re
import jsonschema

ROOT = Path(__file__).resolve().parent.parent
DATA = json.loads((ROOT / 'golden_ecbench.json').read_text())
DEFS = {}

def ref(name): return {'$ref': '#/$defs/' + name}
def obj(properties, required=None):
    return {'type': 'object', 'properties': properties,
            'required': list(properties) if required is None else required,
            'additionalProperties': False}
def arr(items, minimum=0, maximum=None):
    out = {'type': 'array', 'items': items, 'minItems': minimum}
    if maximum is not None: out['maxItems'] = maximum
    return out
def choices(*names): return {'oneOf': [ref(name) for name in names]}
def unique(schemas):
    result = {}
    for schema in schemas: result[json.dumps(schema, sort_keys=True)] = schema
    values = list(result.values())
    return values[0] if len(values) == 1 else {'anyOf': values}
def constant(value): return {'const': value}
def condition(key, value, target):
    return {'if': {'properties': {key: constant(value)}, 'required': [key]}, 'then': target}

DEFS.update({
    'integer': {'type': 'string', 'pattern': r'^(?:0|-?[1-9][0-9]*)$'},
    'unsigned': {'type': 'string', 'pattern': r'^(?:0|[1-9][0-9]*)$'},
    'positive': {'type': 'string', 'pattern': r'^[1-9][0-9]*$'},
    'negative': {'type': 'string', 'pattern': r'^-[1-9][0-9]*$'},
    'element': {'type': 'string', 'pattern': r'^0x(?:[0-9a-f]{2})+$',
                '$comment': 'Field-specific byte width and range require semantic verification.'},
    'mask': {'type': 'string', 'pattern': r'^0x[1-9a-f][0-9a-f]*$'},
    'octets': {'type': 'string', 'pattern': r'^(?:[0-9a-f]{2})*$',
               '$comment': 'Empty octets are permitted for negative decode vectors.'},
    'hash': {'type': 'string', 'pattern': r'^[0-9a-f]{64}$'},
    'fieldId': {'enum': [f['id'] for f in DATA['fields']]},
    'curveId': {'enum': [c['id'] for c in DATA['curves']]},
    'finitePoint': obj({'x': ref('element'), 'y': ref('element')}),
    'rawFinitePoint': obj({'x': ref('integer'), 'y': ref('integer')}),
})
DEFS['point'] = {'oneOf': [{'type': 'null'}, ref('finitePoint')]}
DEFS['rawPoint'] = {'oneOf': [{'type': 'null'}, ref('rawFinitePoint')]}
DEFS['factor'] = obj({'prime': ref('positive'), 'exponent': {'type': 'integer', 'minimum': 1, 'maximum': 4096}})
DEFS['factorization'] = arr(ref('factor'))
DEFS['coefficients'] = obj({name: ref('element') for name in ('a1', 'a2', 'a3', 'a4', 'a6')})

ENUM_KEYS = {'kind', 'model', 'family', 'basis', 'type', 'method', 'status',
             'class_number_status', 'identification_status', 'membership_rule',
             'partner_relation', 'projection_formula', 'projector_formula',
             'legacy_source_file'}
POSITIVE_KEYS = {'cardinality', 'characteristic', 'modulus_integer',
                 'multiplicative_primitive_root', 'n', 'prime', 'order', 'cofactor',
                 'group_exponent', 'conductor', 'order_conductor', 'class_number',
                 'projection_scalar', 'projector', 'embedding_degree',
                 'scan_points_examined', 'ramified_prime', 'certificate_n'}
PAIRS = {'group_invariant_factors', 'G_source', 'H_source'}

def infer(values, path):
    """Infer only shapes; semantic mathematical consistency remains separate."""
    key = path.rsplit('.', 1)[-1]
    nonnull = [v for v in values if v is not None]
    if len(nonnull) != len(values):
        if not nonnull: return {'type': 'null'}
        return {'anyOf': [{'type': 'null'}, infer(nonnull, path)]}
    if key in {'coefficients', 'target_coefficients'}: return ref('coefficients')
    if key in {'point', 'generator', 'source_point'} and isinstance(values[0], dict): return ref('finitePoint')
    if key == 'named_points' and isinstance(values[0], dict):
        names = sorted(set().union(*(set(v) for v in values)))
        DEFS['namedPoint'] = obj({'point': ref('point'), 'order': ref('positive'), 'order_factorization': ref('factorization')})
        return obj({name: ref('namedPoint') for name in names}, ['G', 'O', 'P', 'Q', 'R'])
    if all(isinstance(v, dict) for v in values):
        # Preserve disjoint authoritative variants where their discriminators differ.
        for discriminator in ('kind', 'type', 'method'):
            if all(discriminator in v for v in values) and len({str(v[discriminator]) for v in values}) > 1:
                groups = {}
                for value in values: groups.setdefault(value[discriminator], []).append(value)
                return {'oneOf': [infer(group, path) for group in groups.values()]}
        if key == 'witnesses':
            return {'type': 'object', 'propertyNames': ref('positive'), 'additionalProperties': ref('positive')}
        if key == 'record_hashes':
            return obj({cid: ref('hash') for cid in DEFS['curveId']['enum']})
        if key == 'files': return obj({name: ref('hash') for name in values[0]})
        names = sorted(set().union(*(set(v) for v in values)))
        common = set(values[0]).intersection(*(set(v) for v in values[1:]))
        properties = {name: infer([v[name] for v in values if name in v], path + '.' + name) for name in names}
        if key == 'vectors_by_domain': properties['map'] = {'type': 'integer', 'minimum': 0, 'maximum': 10000000}
        output = obj(properties, sorted(common))
        if key in {'base_field', 'geometric'} and 'class_number_status' in properties:
            output['allOf'] = [condition('class_number_status', 'known_exact', {'properties': {'class_number': ref('positive')}}),
                               condition('class_number_status', 'not_computed', {'properties': {'class_number': {'type': 'null'}}})]
        return output
    if all(isinstance(v, list) for v in values):
        if key.endswith('factorization'): return ref('factorization')
        children = [child for value in values for child in value]
        item = infer(children, path + '[]') if children else False
        output = arr(item)
        if key in PAIRS or all(len(v) == 2 and all(isinstance(x, str) and re.fullmatch(r'\d+', x) for x in v) for v in values):
            output.update(minItems=2, maxItems=2)
        if key in PAIRS or key in {'generator', 'source_point'} and all(isinstance(v, list) for v in values):
            output['items'] = ref('unsigned')
        if key == 'frobenius_polynomial':
            output = {'type': 'array', 'prefixItems': [constant('1'), ref('integer'), ref('positive')], 'items': False, 'minItems': 3, 'maxItems': 3}
        if key in {'polynomial_exponents', 'minimality_prime_divisors', 'curve_ids'} and key != 'curve_ids': output['uniqueItems'] = True
        return output
    if all(type(v) is bool for v in values):
        if key in {'ramified_at_infinity', 'target_is_auxiliary_model', 'explicit_maps_computed'}: return constant(values[0])
        return {'type': 'boolean'}
    if all(type(v) is int for v in values):
        return {'type': 'integer', 'minimum': 1 if key in {'exponent', 'degree', 'element_bits', 'coordinate_bytes'} else 0, 'maximum': 512 if key == 'coordinate_bytes' else 4096 if key in {'exponent', 'degree', 'element_bits'} else 10000000}
    if all(isinstance(v, str) for v in values):
        if key in {'field_id'}: return ref('fieldId')
        if key in {'partner_id', 'source_label'}: return ref('curveId')
        if key == 'id':
            return ref('fieldId' if path.startswith('fields') else 'curveId')
        if key == 'modulus': return ref('mask')
        if all(v.startswith('0x') for v in values): return ref('element')
        if all(re.fullmatch(r'(?:0|-?[1-9][0-9]*)', v) for v in values):
            if 'discriminant' in key: return ref('negative')
            if key.endswith('scan_x') or key == 'scan_final_x': return ref('unsigned')
            return ref('positive' if key in POSITIVE_KEYS else 'integer')
        if key in ENUM_KEYS: return {'type': 'string', 'enum': sorted(set(values))}
        return {'type': 'string', 'minLength': 1}
    raise TypeError((path, {type(v).__name__ for v in values}))

for name, values in [('field', DATA['fields']), ('curve', DATA['curves']), ('certificate', DATA['certificates'])]:
    DEFS[name] = infer(values, name + 's[]')
DEFS['relations'] = infer([DATA['relations']], 'relations')
DEFS['manifest'] = infer([DATA['manifest']], 'manifest')
for key in ('release_version', 'schema_version', 'dataset_id', 'integer_encoding', 'field_element_encoding', 'sec1_octet_encoding'):
    DEFS['manifest']['properties'][key] = constant(DATA['manifest'][key])
DEFS['manifest']['properties']['created_date'] = {'type': 'string', 'format': 'date'}
# Current fields have exact object variants, including null binary proof subjects.
for variant in DEFS['field']['oneOf']:
    kind = variant['properties']['kind']['enum'][0]
    variant['properties']['kind'] = constant(kind)
    if kind == 'prime': variant['properties']['degree'] = constant(1)
    else: variant['properties']['characteristic'] = constant('2')
for variant in DEFS['certificate']['oneOf']:
    method = variant['properties']['method']['enum'][0]
    variant['properties']['method'] = constant(method)
    if method == 'base':
        variant['properties']['n'] = constant('2')
        variant['properties']['n_minus_1_factorization'] = arr(ref('factor'), 0, 0)
        variant['properties']['witnesses'] = {'type': 'object', 'maxProperties': 0, 'additionalProperties': False}
    else:
        variant['properties']['n_minus_1_factorization'] = arr(ref('factor'), 1)
        variant['properties']['witnesses']['minProperties'] = 1

# Expected values are exactly one tagged result, never an ambiguous combination.
for tag, definition in [('point', ref('point')), ('element', ref('element')),
                        ('integer', ref('integer')), ('boolean', {'type': 'boolean'}),
                        ('octets', ref('octets')), ('error', {'type': 'string', 'pattern': r'^[a-z][a-z0-9_]*$'})]:
    DEFS['expected_' + tag] = obj({tag: definition})
DEFS['expected_octets']['properties']['octets'] = {'allOf': [ref('octets'), {'minLength': 2}]}
DEFS['expected'] = choices(*('expected_' + tag for tag in ('point', 'element', 'integer', 'boolean', 'octets', 'error')))
for tag in ('point', 'element', 'integer', 'boolean', 'octets'):
    DEFS[tag + 'Result'] = choices('expected_' + tag, 'expected_error')

DEFS['unaryElementInput'] = {'oneOf': [obj({'a': ref('element')}), obj({'raw_a': ref('integer')})]}
DEFS['binaryElementInput'] = {'oneOf': [obj({a: ref('element' if a == 'a' else 'integer'), b: ref('element' if b == 'b' else 'integer')}) for a in ('a', 'raw_a') for b in ('b', 'raw_b')]}
DEFS['integerInput'] = obj({'value': ref('integer')})
DEFS['powerInput'] = obj({'a': ref('element'), 'exponent': ref('integer')})
DEFS['pointInput'] = obj({'P': ref('point')})
DEFS['pointsInput'] = obj({'P': ref('point'), 'Q': ref('point')})
DEFS['scalarInput'] = obj({'P': ref('point'), 'scalar': ref('integer')})
DEFS['encodeInput'] = obj({'P': ref('point'), 'format': {'enum': ['compressed', 'uncompressed']}})
DEFS['decodeInput'] = obj({'octets': ref('octets')})
DEFS['membershipInput'] = {'oneOf': [obj({'P': ref('point')}), obj({'raw_point': ref('rawPoint')})]}
DEFS['mapInput'] = obj({'P': ref('point'), 'target_coefficients': ref('coefficients'), 'parameter': ref('element')}, ['P', 'target_coefficients'])

def vector(domain, cases):
    properties = {'id': {'type': 'string', 'pattern': r'^(?:field|curve|encoding|membership|map):[A-Za-z0-9][A-Za-z0-9_-]*:[0-9]{4,}$'},
                  'domain': constant(domain), 'operation': {'enum': list(cases)},
                  'purpose': {'enum': ['known_answer', 'input_rejection', 'subgroup_conformance', 'map_conformance', 'explicit_map_known_answer']},
                  'inputs': {}, 'expected': ref('expected')}
    properties['field_id' if domain == 'field' else 'curve_id'] = ref('fieldId' if domain == 'field' else 'curveId')
    output = obj(properties)
    output['allOf'] = [condition('operation', op, {'properties': {'inputs': ref(inputs), 'expected': ref(expected)}}) for op, (inputs, expected) in cases.items()]
    return output
field_cases = {}
for operation in ('add', 'sub', 'mul', 'div', 'square', 'inverse', 'sqrt', 'pow', 'frobenius', 'trace', 'halftrace', 'reduce', 'canonical'):
    inputs = 'binaryElementInput' if operation in ('add', 'sub', 'mul', 'div') else 'powerInput' if operation in ('pow', 'frobenius') else 'integerInput' if operation in ('reduce', 'canonical') else 'unaryElementInput'
    result = 'integerResult' if operation == 'trace' else 'booleanResult' if operation == 'canonical' else 'elementResult'
    field_cases[operation] = inputs, result
DEFS['fieldVector'] = vector('field', field_cases)
DEFS['curveVector'] = vector('curve', {'add': ('pointsInput', 'pointResult'), 'double': ('pointInput', 'pointResult'), 'neg': ('pointInput', 'pointResult'), 'mul': ('scalarInput', 'pointResult')})
DEFS['encodingVector'] = vector('encoding', {'encode': ('encodeInput', 'octetsResult'), 'decode': ('decodeInput', 'pointResult')})
DEFS['membershipVector'] = vector('membership', {op: ('membershipInput', 'booleanResult') for op in ('validate_point', 'is_r_torsion', 'subgroup_contains', 'public_key_valid')})
DEFS['mapVector'] = vector('map', {op: ('mapInput', 'pointResult') for op in ('degree2_isogeny', 'gaussian_automorphism', 'eisenstein_automorphism', 'relative_frobenius')})
DEFS['mapVector']['allOf'] += [condition('operation', op, {'properties': {'inputs': {'required': ['parameter']}}}) for op in ('gaussian_automorphism', 'eisenstein_automorphism')]
DEFS['vector'] = choices('fieldVector', 'curveVector', 'encodingVector', 'membershipVector', 'mapVector')

SCHEMA = {'$schema': 'https://json-schema.org/draft/2020-12/schema', '$id': 'urn:ecbench:golden:dataset:1.0.0',
          'title': 'ECBench golden small-field dataset 1.0.0',
          'description': 'Strict structural schema. Fixed widths, field ranges, references, counts, hashes, proofs, and arithmetic require the semantic verifier.',
          '$defs': DEFS, **obj({'manifest': ref('manifest'), 'fields': arr(ref('field'), 7, 7),
                              'curves': arr(ref('curve'), 56, 56), 'vectors': arr(ref('vector'), 1),
                              'relations': ref('relations'), 'certificates': arr(ref('certificate'), 1)})}
jsonschema.Draft202012Validator.check_schema(SCHEMA)
validator = jsonschema.Draft202012Validator(SCHEMA, format_checker=jsonschema.FormatChecker())
errors = list(validator.iter_errors(DATA))
if errors:
    for error in errors[:5]: print(list(error.absolute_path), error.message[:300])
    raise SystemExit(f'{len(errors)} schema errors')
(ROOT / 'schemas/dataset.schema.json').write_text(json.dumps(SCHEMA, indent=2, sort_keys=True) + '\n')
print('Wrote strict Draft 2020-12 schema; current complete golden bundle validates.')
