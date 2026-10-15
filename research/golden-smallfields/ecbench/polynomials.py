"""Independent, deterministic small Semaev systems for algebra experiments.

Only the standard library is needed for these binary generators. Polynomial
coefficients are canonical decimal strings; exponents and bounded structural
metadata are JSON integers. These deliberately small systems model two-point
relations, not a complete index-calculus attack or a security estimate.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Optional

Point = Optional[tuple[int, int]]
class BinaryField:
    def __init__(self, degree: int, exponents: list[int]):
        if degree % 2 != 1:
            raise ValueError('Half trace needs an odd field degree')
        if max(exponents) != degree or min(exponents) != 0:
            raise ValueError('Expected a monic polynomial with constant one')
        self.degree = degree
        self.exponents = list(exponents)
        self.modulus = sum(1 << e for e in exponents)
        self.q = 1 << degree
        self.mask = self.q - 1

    def reduce(self, polynomial: int) -> int:
        if polynomial < 0:
            raise ValueError('A polynomial bitmask must be nonnegative')
        while polynomial.bit_length() > self.degree:
            polynomial ^= self.modulus << (polynomial.bit_length() - self.degree - 1)
        return polynomial

    def add(self, a: int, b: int) -> int:
        return a ^ b

    sub = add

    def mul(self, a: int, b: int) -> int:
        result = 0
        a = self.reduce(a)
        b = self.reduce(b)
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= self.modulus
        return result

    def square(self, a: int) -> int:
        a = self.reduce(a)
        result = 0
        while a:
            bit = a & -a
            result ^= 1 << (2 * (bit.bit_length() - 1))
            a ^= bit
        return self.reduce(result)

    def inv(self, a: int) -> int:
        a = self.reduce(a)
        if a == 0:
            raise ZeroDivisionError('zero has no field inverse')
        # Polynomial extended Euclid; subtraction is XOR.
        u, v = a, self.modulus
        s, t = 1, 0
        while u != 1:
            if u == 0:
                raise ValueError('Modulus is reducible or input is not a unit')
            shift = u.bit_length() - v.bit_length()
            if shift < 0:
                u, v = v, u
                s, t = t, s
                shift = -shift
            u ^= v << shift
            s ^= t << shift
        return self.reduce(s)

    def div(self, a: int, b: int) -> int:
        return self.mul(a, self.inv(b))

    def pow(self, a: int, exponent: int) -> int:
        if exponent < 0:
            a, exponent = self.inv(a), -exponent
        result = 1
        a = self.reduce(a)
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            exponent >>= 1
            if exponent:
                a = self.square(a)
        return result

    def frobenius(self, a: int, iterations: int = 1) -> int:
        if iterations < 0:
            raise ValueError('Frobenius iteration must be nonnegative')
        a = self.reduce(a)
        for _ in range(iterations % self.degree):
            a = self.square(a)
        return a

    def sqrt(self, a: int) -> int:
        return self.frobenius(a, self.degree - 1)

    def trace(self, a: int) -> int:
        a = self.reduce(a)
        value = 0
        for _ in range(self.degree):
            value ^= a
            a = self.square(a)
        if value not in (0, 1):
            raise AssertionError('Absolute trace must lie in F_2')
        return value

    def halftrace(self, a: int) -> int:
        a = self.reduce(a)
        value = 0
        for _ in range((self.degree + 1) // 2):
            value ^= a
            a = self.square(self.square(a))
        return value


class BinaryCurve:
    def __init__(self, field: BinaryField, a: int, b: int):
        self.field = field
        self.a = field.reduce(a)
        self.b = field.reduce(b)
        if not self.b:
            raise ValueError('Singular ordinary binary model: b=0')

    def on_curve(self, point: Point) -> bool:
        if point is None:
            return True
        x, y = point
        if not (0 <= x < self.field.q and 0 <= y < self.field.q):
            return False
        f = self.field
        lhs = f.square(y) ^ f.mul(x, y)
        x2 = f.square(x)
        rhs = f.mul(x2, x) ^ f.mul(self.a, x2) ^ self.b
        return lhs == rhs

    def point_at_x(self, x: int) -> Point:
        f = self.field
        if not 0 <= x < f.q:
            raise ValueError('Noncanonical x encoding')
        if x == 0:
            point = (0, f.sqrt(self.b))
        else:
            # z=y/x solves z^2+z=x+a+b/x^2.
            beta = x ^ self.a ^ f.div(self.b, f.square(x))
            if f.trace(beta):
                return None
            z = f.halftrace(beta)
            if f.square(z) ^ z != beta:
                raise AssertionError('Half trace failed to solve the equation')
            y = f.mul(x, z)
            point = (x, min(y, y ^ x))
        if not self.on_curve(point):
            raise AssertionError('Point finder produced an invalid point')
        return point

    def points(self, start_x: int = 0):
        for x in range(start_x, self.field.q):
            point = self.point_at_x(x)
            if point is not None:
                yield point

    def neg(self, point: Point) -> Point:
        if point is None:
            return None
        x, y = point
        return x, x ^ y

    def double(self, point: Point) -> Point:
        if point is None or point[0] == 0:
            return None
        x, y = point
        f = self.field
        slope = x ^ f.div(y, x)
        xx = f.square(slope) ^ slope ^ self.a
        yy = f.square(x) ^ f.mul(slope ^ 1, xx)
        return xx, yy

    def add(self, p: Point, q: Point) -> Point:
        if p is None:
            return q
        if q is None:
            return p
        x1, y1 = p
        x2, y2 = q
        if x1 == x2:
            if y1 == y2:
                return self.double(p)
            if y1 ^ y2 != x1:
                raise ValueError('Same-x inputs are not valid curve points')
            return None
        f = self.field
        slope = f.div(y1 ^ y2, x1 ^ x2)
        xx = f.square(slope) ^ slope ^ x1 ^ x2 ^ self.a
        yy = f.mul(slope, x1 ^ xx) ^ xx ^ y1
        return xx, yy

    def mul(self, scalar: int, point: Point) -> Point:
        if scalar < 0:
            return self.mul(-scalar, self.neg(point))
        result = None
        while scalar:
            if scalar & 1:
                result = self.add(result, point)
            scalar >>= 1
            if scalar:
                point = self.double(point)
        return result




def binary_f3_value(field: BinaryField, x: int, y: int, z: int, b: int) -> int:
    """Semaev f3 for y^2+xy=x^3+a*x^2+b, with a in F_2."""
    pair_sum=field.mul(x,y)^field.mul(x,z)^field.mul(y,z)
    return field.square(pair_sum)^field.mul(field.mul(x,y),z)^b


def _hex(field: BinaryField, value: int) -> str:
    return '0x'+format(value,f'0{2*((field.degree+7)//8)}x')


def _point(field: BinaryField, point: Point):
    if point is None:
        return None
    return {'x':_hex(field,point[0]),'y':_hex(field,point[1])}


def _bit_span(basis: list[int], mask: int) -> int:
    value=0
    for i,element in enumerate(basis):
        if mask>>i&1:
            value ^= element
    return value


def _independent(elements: list[int]) -> bool:
    pivots={}
    for element in elements:
        value=element
        while value:
            pivot=value.bit_length()-1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot]=value
                break
        if value==0:
            return False
    return True


def _basis(curve: BinaryCurve, generator: Point, dimension: int, kind: str) -> list[int]:
    if kind=='low_monomials':
        return [1<<i for i in range(dimension)]
    if kind!='subgroup_seeded_basis':
        raise ValueError('Unknown factor-base basis style')
    assert generator is not None
    doubled=curve.double(generator)
    assert doubled is not None
    basis=[generator[0],doubled[0]]
    assert _independent(basis)
    for i in range(curve.field.degree):
        if len(basis)==dimension:
            break
        candidate=1<<i
        if _independent(basis+[candidate]):
            basis.append(candidate)
    assert len(basis)==dimension
    return basis


def _xor_poly(left: dict[int,int],right: dict[int,int]) -> dict[int,int]:
    result=dict(left)
    for monomial,coefficient in right.items():
        result[monomial]=result.get(monomial,0)^coefficient
        if not result[monomial]:
            del result[monomial]
    return result


def _field_poly_product(field: BinaryField, left: dict[int,int],right: dict[int,int]) -> dict[int,int]:
    """Product in F_(2^m)[v]/(v_i^2-v_i), using bitmask monomials."""
    result={}
    for lm,lc in left.items():
        for rm,rc in right.items():
            monomial=lm|rm
            result[monomial]=result.get(monomial,0)^field.mul(lc,rc)
            if not result[monomial]:
                del result[monomial]
    return result


def _field_poly_scale(field: BinaryField, polynomial: dict[int,int],constant: int) -> dict[int,int]:
    return {m:c for m,v in polynomial.items() if (c:=field.mul(v,constant))}


def _f3_field_polynomial(field: BinaryField,basis: list[int],z: int,b: int) -> dict[int,int]:
    d=len(basis)
    x={1<<i:value for i,value in enumerate(basis)}
    y={1<<(d+i):value for i,value in enumerate(basis)}
    xy=_field_poly_product(field,x,y)
    summed=_xor_poly(xy,_xor_poly(_field_poly_scale(field,x,z),_field_poly_scale(field,y,z)))
    squared={monomial:field.square(value) for monomial,value in summed.items()}
    return _xor_poly(squared,_xor_poly(_field_poly_scale(field,xy,z),{0:b}))


def _field_poly_evaluate(polynomial: dict[int,int],assignment_mask: int) -> int:
    value=0
    for monomial,coefficient in polynomial.items():
        if assignment_mask&monomial==monomial:
            value ^= coefficient
    return value


def _bool_poly_product(left: dict[int,int],right: dict[int,int]) -> dict[int,int]:
    result={}
    for lm,lc in left.items():
        for rm,rc in right.items():
            mask=lm|rm
            result[mask]=result.get(mask,0)^(lc&rc)
            if result[mask]==0:
                del result[mask]
    return result


def _membership_indicator(dimension: int,allowed_masks: list[int],offset: int) -> dict[int,int]:
    """Polynomial 1+sum delta_allowed, zero exactly on permitted bit tuples."""
    result={0:1}
    for value in allowed_masks:
        delta={0:1}
        for i in range(dimension):
            bit=1<<(offset+i)
            factor={bit:1} if value>>i&1 else {0:1,bit:1}
            delta=_bool_poly_product(delta,factor)
        result=_xor_poly(result,delta)
    for assignment in range(1<<dimension):
        actual=_field_poly_evaluate(result,assignment<<offset)
        assert actual==int(assignment not in allowed_masks)
    return result


def _sparse_boolean_terms(polynomial: dict[int,int],variables: int):
    return [{'coefficient':str(coefficient),
             'exponents':[(mask>>i)&1 for i in range(variables)]}
            for mask,coefficient in sorted(polynomial.items()) if coefficient]


def _poly_degree(terms: list[dict]) -> int:
    return max((sum(term['exponents']) for term in terms),default=0)


def _factor_base(curve: BinaryCurve,curve_record: dict,basis: list[int],kind: str):
    field=curve.field
    r=int(curve_record['subgroup']['order'])
    assert curve_record['subgroup']['r_torsion_rank']==1
    allowed_masks=[]
    raw_points={}
    counters={'candidate_x':1<<len(basis),'x_without_curve_points':0,
              'x_with_non_r_points':0,'accepted_x':0,'accepted_points':0}
    for mask in range(1<<len(basis)):
        x=_bit_span(basis,mask)
        p=curve.point_at_x(x)
        if p is None:
            counters['x_without_curve_points']+=1
            continue
        if curve.mul(r,p) is not None:
            counters['x_with_non_r_points']+=1
            continue
        assert p is not None and curve.on_curve(p)
        points=[p]
        neg=curve.neg(p)
        if neg!=p:
            points.append(neg)
        assert all(curve.mul(r,point) is None for point in points)
        allowed_masks.append(mask)
        raw_points[x]=points
        counters['accepted_x']+=1
        counters['accepted_points']+=len(points)
    assert sum(counters[key] for key in ('x_without_curve_points','x_with_non_r_points','accepted_x'))==counters['candidate_x']
    serialized={'dimension':len(basis),'basis_kind':kind,
        'basis_elements':[_hex(field,value) for value in basis],
        'span_x_values':[_hex(field,_bit_span(basis,mask)) for mask in range(1<<len(basis))],
        'x_values':[_hex(field,x) for x in sorted(raw_points)],
        'x_encoding_masks':{_hex(field,_bit_span(basis,mask)):str(mask) for mask in allowed_masks},
        'points_by_x':{_hex(field,x):[_point(field,p) for p in points] for x,points in sorted(raw_points.items())},
        'membership_policy':'on curve and in the selected rank-one rational prime-order subgroup; both signs allowed',
        'counters':counters,
        'sampling_note':('Deliberately includes x(G),x([2]G); not an unbiased factor-base or relation-yield sample'
                         if kind=='subgroup_seeded_basis' else 'Fixed low polynomial basis; no selected-subgroup point existence assumed')}
    return serialized,raw_points,allowed_masks


def _seeded_scalar(seed: str,label: str,r: int) -> int:
    # Deterministic pseudorandom scalar; reduction bias is documented, not hidden.
    digest=hashlib.sha256((seed+'\0'+label).encode()).digest()
    return 1+int.from_bytes(digest,'big')%(r-1)


def _point_key(point: Point):
    return (-1,-1) if point is None else point


def _planted_target(curve: BinaryCurve,raw_points: dict[int,list[Point]],
                    generator: Point,basis_kind: str):
    if not raw_points:
        return None,None
    if basis_kind=='subgroup_seeded_basis':
        p=generator
        q=curve.double(generator)
        assert p is not None and q is not None
        assert p in raw_points[p[0]] and q in raw_points[q[0]]
        target=curve.add(p,q)
        assert target is not None
        return target,(p,q)
    xs=sorted(raw_points)
    p=raw_points[xs[0]][0]
    q=raw_points[xs[1]][0] if len(xs)>1 else p
    target=curve.add(p,q)
    assert target is not None
    return target,(p,q)


def _binary_case(curve_record: dict,field_record: dict,curve: BinaryCurve,basis: list[int],
                 factor_base: dict,raw_points: dict,allowed_masks: list[int],
                 target: Point,target_kind: str,seed: str,scalar: int|None,
                 planted_pair: tuple[Point,Point]|None):
    field=curve.field
    d=len(basis)
    variables=[f'x1_{i}' for i in range(d)]+[f'x2_{i}' for i in range(d)]
    assert target is not None and curve.on_curve(target)
    z=target[0]
    expanded=_f3_field_polynomial(field,basis,z,curve.b)
    polynomials=[]
    for coordinate in range(field.degree):
        terms=_sparse_boolean_terms({mask:1 for mask,coefficient in expanded.items() if coefficient>>coordinate&1},2*d)
        polynomials.append({'id':f'f3_coordinate_{coordinate}','role':'semaev_f3_coordinate','coordinate':coordinate,
                            'terms':terms,'total_degree':_poly_degree(terms)})
    for i,name in enumerate(variables):
        square=[0]*(2*d);square[i]=2
        linear=[0]*(2*d);linear[i]=1
        polynomials.append({'id':f'boolean_{name}','role':'boolean_field_equation',
            'terms':[{'coefficient':'1','exponents':square},{'coefficient':'1','exponents':linear}],
            'total_degree':2})
    for number,offset in [(1,0),(2,d)]:
        terms=_sparse_boolean_terms(_membership_indicator(d,allowed_masks,offset),2*d)
        polynomials.append({'id':f'membership_x{number}','role':'factor_base_membership',
                            'terms':terms,'total_degree':_poly_degree(terms)})
    algebraic=[]
    x_pairs=[]
    for first in allowed_masks:
        x=_bit_span(basis,first)
        for second in allowed_masks:
            y=_bit_span(basis,second)
            mask=first|(second<<d)
            value=binary_f3_value(field,x,y,z,curve.b)
            assert value==_field_poly_evaluate(expanded,mask)
            if value==0:
                algebraic.append([str((mask>>i)&1) for i in range(2*d)])
                x_pairs.append((x,y))
    algebraic.sort(key=lambda row:tuple(map(int,row)))
    ordered=[]
    point_pairs=[]
    points=[p for x in sorted(raw_points) for p in raw_points[x]]
    for p in points:
        for q in points:
            summed=curve.add(p,q)
            if summed is not None:
                # This verifies the symbolic f3 formula using independent ECC operations.
                assert binary_f3_value(field,p[0],q[0],summed[0],curve.b)==0
            if summed==target:
                ordered.append({'P':_point(field,p),'Q':_point(field,q)})
                point_pairs.append((p,q))
    relation_x=set((p[0],q[0]) for p,q in point_pairs)
    assert set(x_pairs)==relation_x
    # No rational odd-prime-order factor-base point is its own negative.
    assert len(point_pairs)==len(x_pairs)
    orbit_points={tuple(sorted((p,q))) for p,q in point_pairs}
    orbit_x={tuple(sorted((x,y))) for x,y in x_pairs}
    target_record={'point':_point(field,target),'kind':target_kind,'source':target_kind,
                   'scalar':str(scalar) if scalar is not None else None,
                   'sampling_method':('SHA256(seed NUL case-label), reduced modulo r-1, plus one; deterministic pseudorandom, not exact uniform sampling'
                                      if target_kind=='seeded_subgroup' else 'Sum of two explicitly retained factor-base points')}
    if planted_pair is not None:
        target_record['planted_decomposition']={'P':_point(field,planted_pair[0]),'Q':_point(field,planted_pair[1])}
    kind=factor_base['basis_kind']
    case_id=f"{curve_record['id']}-f3-{kind}-d{d}-{target_kind}"
    degree_roles={role:max((p['total_degree'] for p in polynomials if p['role']==role),default=0)
                  for role in ['semaev_f3_coordinate','boolean_field_equation','factor_base_membership']}
    return {'id':case_id,'curve_id':curve_record['id'],'curve':curve_record,'field':field_record,
        'coefficient_modulus':'2','variables':variables,'variable_order':variables,
        'polynomials':polynomials,'coordinate_groups':[
            {'name':'x1','variables':variables[:d],'basis':factor_base['basis_elements'],'encoding':'binary_linear'},
            {'name':'x2','variables':variables[d:],'basis':factor_base['basis_elements'],'encoding':'binary_linear'}],
        'domains':{name:['0','1'] for name in variables},'factor_base':factor_base,
        'target':target_record,'arity':2,'summation_polynomial_arity':3,
        'model_variant':'strict_selected_subgroup_factor_base',
        'expected':{'algebraic_assignments':algebraic,'ordered_point_decompositions':ordered,
                    'ordered_count':len(point_pairs),'swap_orbit_count':len(orbit_points),
                    'ordered_x_count':len(x_pairs),'swap_orbit_x_count':len(orbit_x),
                    'algebraic_assignment_count':len(algebraic)},
        'generation':{'seed':seed,'method':'Independent binary polynomial arithmetic and affine ECC; exact finite factor-base enumeration',
            'f3_formula':'(x1*x2 + x1*z + x2*z)^2 + x1*x2*z + b',
            'original_f3_total_degree':4,'boolean_reduced_f3_max_degree':degree_roles['semaev_f3_coordinate'],
            'polynomial_degrees_by_role':degree_roles,
            'membership_filter_note':'Additional Boolean truth-table equations restrict x to actual selected-subgroup points; their degree is not the degree of Semaev f3',
            'diagnostic_empty_factor_base':not bool(raw_points),
            'truth_method':'Exhaustive selected factor-base x-pair f3 evaluation, independently matched to signed point-pair sums',
            'claims':'Small two-point fixture only; no higher-arity relation or cryptographic hardness estimate'}}


def build_binary_cases(fields: list[dict],curves: list[dict],seed: str='ecbench-indexcalculus-v1') -> list[dict]:
    field_by_id={record['id']:record for record in fields}
    curve_by_id={record['id']:record for record in curves}
    cases=[]
    for cid in ['K83-a0','S83-a1','B83-a0']:
        record=curve_by_id[cid]
        field_record=field_by_id[record['field_id']]
        field=BinaryField(field_record['degree'],field_record['representation']['polynomial_exponents'])
        coefficients=record['coefficients']
        curve=BinaryCurve(field,int(coefficients['a2'],16),int(coefficients['a6'],16))
        generator=tuple(int(record['subgroup']['generator'][name],16) for name in ('x','y'))
        r=int(record['subgroup']['order'])
        assert curve.on_curve(generator) and curve.mul(r,generator) is None
        for basis_kind,dimensions in [('low_monomials',[2,3,4,5,6]),('subgroup_seeded_basis',[2,3,4,5])]:
            for dimension in dimensions:
                basis=_basis(curve,generator,dimension,basis_kind)
                factor_base,raw_points,allowed_masks=_factor_base(curve,record,basis,basis_kind)
                planted,pair=_planted_target(curve,raw_points,generator,basis_kind)
                if planted is not None:
                    cases.append(_binary_case(record,field_record,curve,basis,factor_base,raw_points,allowed_masks,
                                              planted,'planted',seed,3 if basis_kind=='subgroup_seeded_basis' else None,pair))
                label=f'{cid}/{basis_kind}/d{dimension}/seeded_subgroup'
                scalar=_seeded_scalar(seed,label,r)
                target=curve.mul(scalar,generator)
                assert target is not None
                cases.append(_binary_case(record,field_record,curve,basis,factor_base,raw_points,allowed_masks,
                                          target,'seeded_subgroup',seed,scalar,None))
    assert len(cases)==49
    return cases
