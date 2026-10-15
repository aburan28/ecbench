#!/usr/bin/env python3
"""Independent relative-Frobenius known answers for the fourteen binary curves.

Output rows retain raw integer field elements. The target is the auxiliary
model with b squared; this adds no new canonical curve records. Koblitz b=1
makes that target the source curve, giving a same-curve endomorphism.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from binary_vectors import BinaryCurve, BinaryField, DEFAULT_SOURCE, PROJECT_ROOT

DEFAULT_VECTORS = PROJECT_ROOT / 'build' / 'binary_vectors.json'
DEFAULT_OUTPUT = PROJECT_ROOT / 'build' / 'binary_maps.json'


def mapped_point(field: BinaryField, point):
    if point is None:
        return None
    return tuple(field.square(coordinate) for coordinate in point)


def build(source: Path = DEFAULT_SOURCE, vectors: Path = DEFAULT_VECTORS,
          output: Path = DEFAULT_OUTPUT):
    normalized = json.loads(source.read_text())
    raw = json.loads(vectors.read_text())
    fields = {field['id']: BinaryField(field['bits'], field['exponents'])
              for field in normalized['fields'] if field['kind'] == 'binary'}
    rows = []
    for info in normalized['curves']:
        if info['kind'] != 'binary':
            continue
        field = fields[info['field']]
        data = raw['curves'][info['id']]
        coefficients = data['coefficients']
        assert coefficients['a1'] == 1 and coefficients['a3'] == coefficients['a4'] == 0
        assert coefficients['a2'] in (0, 1)
        source_curve = BinaryCurve(field, coefficients['a2'], coefficients['a6'])
        target_b = field.square(coefficients['a6'])
        target_curve = BinaryCurve(field, coefficients['a2'], target_b)
        target_coefficients = {'a1': 1, 'a2': coefficients['a2'], 'a3': 0,
                               'a4': 0, 'a6': target_b}
        named_points = {'O': None, 'G': data['subgroup']['generator'],
                        'R': data['R'], 'P': data['P'], 'T': data['two_torsion_point']}
        for name, point in named_points.items():
            point = None if point is None else tuple(point)
            assert source_curve.on_curve(point)
            result = mapped_point(field, point)
            assert target_curve.on_curve(result)
            rows.append({'curve_id': info['id'], 'operation': 'relative_frobenius',
                         'point': point, 'parameter': None,
                         'target_coefficients': target_coefficients, 'result': result})
            # Independently check compatibility with negation and doubling.
            assert mapped_point(field, source_curve.neg(point)) == target_curve.neg(result)
            assert mapped_point(field, source_curve.double(point)) == target_curve.double(result)
        p, q = tuple(data['P']), tuple(data['Q'])
        assert mapped_point(field, source_curve.add(p, q)) == target_curve.add(
            mapped_point(field, p), mapped_point(field, q))
        if info['family'] == 'K':
            assert target_coefficients == coefficients
            point = tuple(data['subgroup']['generator'])
            for _ in range(field.degree):
                point = mapped_point(field, point)
            assert point == tuple(data['subgroup']['generator'])
        else:
            assert target_b != coefficients['a6']
    assert len(rows) == 70
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(rows, indent=2) + '\n')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--vectors', type=Path, default=DEFAULT_VECTORS)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    rows = build(args.source, args.vectors, args.output)
    print(f'Wrote {len(rows)} independently checked relative-Frobenius vectors to {args.output}')


if __name__ == '__main__':
    main()
