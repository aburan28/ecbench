#!/usr/bin/env python3
"""Independent two-point Semaev cases over the canonical 83-bit prime field.

The only point-arithmetic dependency is ``tools/prime_vectors.py``.  Polynomial
construction, evaluation, and exhaustive answer generation use local standard
library code and do not import ecbench's runtime arithmetic or validator.

For each factor-base size 4, 8, and 16, the x-coordinates are those of kG for
1 <= k <= size.  Both signs of every point are permitted.  The planted target
is G + 2G.  The other target is a deterministic SHA-256-derived nonzero scalar
multiple of G; it is a reproducible test target, not a claim of true randomness.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Optional

try:
    from .prime_vectors import PrimeCurve
except ImportError:  # Direct execution: python tools/prime_algebra_cases.py.
    from prime_vectors import PrimeCurve


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CURVE_ID = "p83-S-0"
DEFAULT_SEED = "ecbench-indexcalculus-v1"
FACTOR_BASE_SIZES = (4, 8, 16)
Point = Optional[tuple[int, int]]
SparsePolynomial = dict[tuple[int, int], int]


def _point_record(point: Point, width: int) -> Optional[dict[str, str]]:
    if point is None:
        return None
    return {"x": f"0x{point[0]:0{2 * width}x}",
            "y": f"0x{point[1]:0{2 * width}x}"}


def _parse_point(record: Optional[dict[str, str]]) -> Point:
    return None if record is None else (int(record["x"], 16),
                                        int(record["y"], 16))


def _load_curve() -> tuple[dict, dict, PrimeCurve, Point, int]:
    data_dir = PROJECT_ROOT / "ecbench" / "data"
    curves = json.loads((data_dir / "curves.json").read_text())
    fields = json.loads((data_dir / "fields.json").read_text())
    record = next(curve for curve in curves if curve["id"] == CURVE_ID)
    field = next(field for field in fields if field["id"] == record["field_id"])
    assert field["kind"] == "prime" and record["model"] == "short_weierstrass"
    p = int(field["characteristic"])
    coefficients = record["coefficients"]
    assert all(int(coefficients[name], 16) == 0 for name in ("a1", "a2", "a3"))
    curve = PrimeCurve(p, int(coefficients["a4"], 16),
                       int(coefficients["a6"], 16))
    generator = _parse_point(record["subgroup"]["generator"])
    subgroup_order = int(record["subgroup"]["order"])
    assert generator is not None and curve.on_curve(generator)
    assert subgroup_order > 2 * max(FACTOR_BASE_SIZES)
    assert curve.mul(subgroup_order, generator) is None
    assert record["subgroup"]["unique_order_r_subgroup"]
    assert record["subgroup"]["r_torsion_rank"] == 1
    return record, field, curve, generator, subgroup_order


def _semaev_f3(curve: PrimeCurve, z: int) -> SparsePolynomial:
    """Expand f3(x1,x2,z), with canonical nonzero coefficients modulo p."""
    p, a, b = curve.p, curve.A, curve.B
    raw = {
        (2, 2): 1,
        (2, 1): -2 * z,
        (1, 2): -2 * z,
        (2, 0): z * z,
        (0, 2): z * z,
        (1, 1): -2 * z * z - 2 * a,
        (1, 0): -2 * a * z - 4 * b,
        (0, 1): -2 * a * z - 4 * b,
        (0, 0): a * a - 4 * b * z,
    }
    return {exponents: coefficient % p for exponents, coefficient in raw.items()
            if coefficient % p}


def _f3_direct(curve: PrimeCurve, x: int, y: int, z: int) -> int:
    """Evaluate the unexpanded Semaev formula as a construction cross-check."""
    a, b = curve.A, curve.B
    return ((x - y)**2 * z**2
            - 2 * ((x + y) * (x * y + a) + 2 * b) * z
            + (x * y - a)**2 - 4 * b * (x + y)) % curve.p


def _membership_polynomial(values: list[int], variable: int,
                           p: int) -> SparsePolynomial:
    """Expand the monic polynomial product(u - x) for factor-base x values."""
    coefficients = [1]  # Ascending powers of u.
    for value in values:
        expanded = [0] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            expanded[degree] = (expanded[degree] - value * coefficient) % p
            expanded[degree + 1] = (expanded[degree + 1] + coefficient) % p
        coefficients = expanded
    return {((degree, 0) if variable == 0 else (0, degree)): coefficient
            for degree, coefficient in enumerate(coefficients) if coefficient}


def _evaluate(polynomial: SparsePolynomial, x: int, y: int, p: int) -> int:
    return sum(coefficient * pow(x, exponents[0], p)
               * pow(y, exponents[1], p)
               for exponents, coefficient in polynomial.items()) % p


def _serialize_polynomial(polynomial: SparsePolynomial, name: str,
                          role: str) -> dict:
    terms = sorted(polynomial.items(),
                   key=lambda item: (sum(item[0]), item[0]), reverse=True)
    return {"id": name, "role": role,
            "terms": [{"coefficient": str(coefficient),
                       "exponents": list(exponents)}
                      for exponents, coefficient in terms]}


def _factor_base(curve: PrimeCurve, generator: Point, order: int,
                 size: int, width: int) -> tuple[dict, list[int], list[Point]]:
    positive_points = [curve.mul(k, generator) for k in range(1, size + 1)]
    assert all(point is not None and curve.on_curve(point)
               and curve.mul(order, point) is None for point in positive_points)
    xs = [point[0] for point in positive_points]
    assert len(xs) == len(set(xs)) == size
    signed_points: list[Point] = []
    points_by_x = {}
    scalar_by_x = {}
    for scalar, point in enumerate(positive_points, 1):
        negative = curve.neg(point)
        assert negative is not None and negative != point
        pair = sorted([point, negative])
        signed_points.extend(pair)
        key = f"0x{point[0]:0{2 * width}x}"
        points_by_x[key] = [_point_record(signed, width) for signed in pair]
        scalar_by_x[key] = str(scalar)
    assert len(signed_points) == 2 * size == len(set(signed_points))
    return ({
        "size_requested": size,
        "x_values": [f"0x{x:0{2 * width}x}" for x in xs],
        "points_by_x": points_by_x,
        "positive_scalar_by_x": scalar_by_x,
        "generator": _point_record(generator, width),
        "subgroup_order": str(order),
        "counters": {
            "scalar_candidates_considered": size,
            "accepted_distinct_x_values": size,
            "signed_points": len(signed_points),
            "infinity_rejections": 0,
            "duplicate_x_rejections": 0,
            "subgroup_annihilation_checks": size,
        },
        "membership_policy": "selected cyclic subgroup",
        "basis_kind": "small_generator_multiples",
    }, xs, sorted(signed_points))


def _target(curve: PrimeCurve, generator: Point, order: int, size: int,
            kind: str, seed: str) -> tuple[Point, int, dict]:
    if kind == "planted":
        first = curve.mul(1, generator)
        second = curve.mul(2, generator)
        point = curve.add(first, second)
        scalar = 3
        derivation = {"target_method": "sum of first two positive kG points",
                      "planted_scalars": ["1", "2"]}
    else:
        payload = f"ecbench-prime-f3-v1\0{seed}\0{CURVE_ID}\0{size}\0seeded_subgroup"
        digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
        scalar = 1 + int(digest, 16) % (order - 1)
        point = curve.mul(scalar, generator)
        derivation = {
            "target_method": "deterministic SHA-256 nonzero subgroup scalar",
            "hash": "SHA-256",
            "hash_input_encoding": "UTF-8",
            "hash_input": payload,
            "hash_digest_hex": digest,
            "scalar_formula": "1 + OS2IP(SHA-256(hash_input)) mod (r - 1)",
            "randomness_claim": "none; deterministic seeded test target",
        }
    assert point is not None and curve.on_curve(point)
    assert point == curve.mul(scalar, generator)
    assert curve.mul(order, point) is None
    return point, scalar, derivation


def _expected(curve: PrimeCurve, target: Point, xs: list[int],
              points: list[Point], polynomials: list[SparsePolynomial],
              width: int) -> tuple[dict, dict]:
    assert target is not None
    z = target[0]
    algebraic = []
    for x in xs:
        for y in xs:
            expanded = _evaluate(polynomials[0], x, y, curve.p)
            assert expanded == _f3_direct(curve, x, y, z)
            if all(_evaluate(poly, x, y, curve.p) == 0 for poly in polynomials):
                algebraic.append((x, y))
    algebraic.sort()

    decompositions = []
    finite_sums = 0
    infinity_sums = 0
    for first in points:
        for second in points:
            result = curve.add(first, second)
            if result is None:
                infinity_sums += 1
            else:
                finite_sums += 1
                # Check the polynomial identity for every supplied point pair,
                # including pairs whose sum is not this case's fixed target.
                assert _f3_direct(curve, first[0], second[0], result[0]) == 0
                assert _evaluate(_semaev_f3(curve, result[0]),
                                 first[0], second[0], curve.p) == 0
            if result == target:
                decompositions.append((first, second))
    decompositions.sort()
    x_relations = sorted({(first[0], second[0])
                          for first, second in decompositions})
    # All x values lift to signed rational subgroup points.  Exhaustive point
    # enumeration therefore detects every f3 solution and every decomposition.
    assert algebraic == x_relations
    point_orbits = {tuple(sorted((first, second)))
                    for first, second in decompositions}
    x_orbits = {tuple(sorted(pair)) for pair in algebraic}
    return ({
        "algebraic_assignments": [[str(x), str(y)] for x, y in algebraic],
        "ordered_point_decompositions": [
            {"P": _point_record(first, width), "Q": _point_record(second, width)}
            for first, second in decompositions],
        "ordered_count": len(decompositions),
        "swap_orbit_count": len(point_orbits),
        "ordered_x_count": len(algebraic),
        "swap_orbit_x_count": len(x_orbits),
    }, {
        "x_assignments_enumerated": len(xs)**2,
        "ordered_signed_point_pairs_enumerated": len(points)**2,
        "finite_point_sums_identity_checked": finite_sums,
        "infinity_point_sums": infinity_sums,
        "expanded_formula_matches_direct_formula": True,
        "polynomial_zeros_match_point_x_relations": True,
    })


def build_prime_cases(seed: str = DEFAULT_SEED) -> list[dict]:
    """Return all six verified prime-field finite factor-base cases."""
    if not isinstance(seed, str):
        raise TypeError("seed must be a string")
    record, field, curve, generator, order = _load_curve()
    width = int(field["coordinate_bytes"])
    cases = []
    for size in FACTOR_BASE_SIZES:
        factor_base, xs, points = _factor_base(curve, generator, order, size, width)
        for kind in ("planted", "seeded_subgroup"):
            target, scalar, derivation = _target(curve, generator, order,
                                                 size, kind, seed)
            polynomials = [_semaev_f3(curve, target[0]),
                           _membership_polynomial(xs, 0, curve.p),
                           _membership_polynomial(xs, 1, curve.p)]
            expected, checks = _expected(curve, target, xs, points,
                                         polynomials, width)
            cases.append({
                "id": f"{CURVE_ID}-f3-fb{size}-{kind}",
                "curve_id": record["id"],
                "curve": deepcopy(record),
                "field": deepcopy(field),
                "coefficient_modulus": str(curve.p),
                "variables": ["x1", "x2"],
                "polynomials": [
                    _serialize_polynomial(polynomials[0], "f3", "semaev_f3"),
                    _serialize_polynomial(polynomials[1], "membership_x1",
                                          "factor_base_membership"),
                    _serialize_polynomial(polynomials[2], "membership_x2",
                                          "factor_base_membership")],
                "coordinate_groups": [
                    {"name": "x1", "variables": ["x1"]},
                    {"name": "x2", "variables": ["x2"]}],
                "domains": {"x1": list(map(str, xs)), "x2": list(map(str, xs))},
                "factor_base": deepcopy(factor_base),
                "target": {"point": _point_record(target, width),
                           "kind": kind, "source": kind, "scalar": str(scalar)},
                "model_variant": "strict_selected_subgroup_factor_base",
                "expected": expected,
                "arity": 2,
                "summation_polynomial_arity": 3,
                "generation": {
                    "method": "deterministic kG factor base and finite exhaustive verification",
                    "seed": seed,
                    "factor_base_rule": "x(kG) for k=1..size; both point signs permitted",
                    "factor_base_scalar_start": "1",
                    "factor_base_scalar_stop": str(size),
                    "formula": "(x1-x2)^2*z^2 - 2*((x1+x2)*(x1*x2+A)+2*B)*z + (x1*x2-A)^2 - 4*B*(x1+x2)",
                    "checks": checks,
                    **derivation,
                },
            })
    return cases


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", default=DEFAULT_SEED)
    parser.add_argument("--output", type=Path,
                        default=PROJECT_ROOT / "build" / "prime_algebra_cases.json")
    args = parser.parse_args()
    cases = build_prime_cases(args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(cases, indent=2) + "\n")
    print(f"Wrote {len(cases)} cases to {args.output}")
    for case in cases:
        expected = case["expected"]
        print(f"{case['id']}: {expected['ordered_count']} ordered point decompositions, "
              f"{expected['ordered_x_count']} ordered x assignments")


if __name__ == "__main__":
    main()
