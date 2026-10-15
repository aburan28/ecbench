#!/usr/bin/env python3
"""Independent, integer-only polynomial-basis binary-field and ECC vectors.

The internal format uses raw Python integers and affine [x,y] pairs, with None
for the point at infinity. This module uses neither PARI nor any other ECC
implementation. Point finding uses the odd-degree half trace, and scans the
integer encoding of x in increasing order. It handles the general ordinary
model y^2 + x*y = x^3 + a*x^2 + b.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

Point = Optional[tuple[int, int]]
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = PROJECT_ROOT / 'provenance' / 'normalized_data.json'
if not DEFAULT_SOURCE.is_file():
    session_source = Path('/workspace/curve_report/normalized_data.json')
    if session_source.is_file():
        DEFAULT_SOURCE = session_source
DEFAULT_OUTPUT = PROJECT_ROOT / 'build' / 'binary_vectors.json'


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



class SEC1DecodeError(ValueError):
    def __init__(self, code: str):
        super().__init__(code)
        self.code = code


def sec1_encode(curve: BinaryCurve, point: Point, compressed: bool) -> str:
    if not curve.on_curve(point):
        raise ValueError('Cannot encode an off-curve point')
    if point is None:
        return '00'
    f = curve.field
    width = (f.degree + 7) // 8
    x,y = point
    xb = x.to_bytes(width,'big')
    if compressed:
        ybit = 0 if x == 0 else f.div(y,x)&1
        return (bytes([2+ybit])+xb).hex()
    return (b'\x04'+xb+y.to_bytes(width,'big')).hex()


def sec1_decode(curve: BinaryCurve, encoded: str) -> Point:
    try:
        payload=bytes.fromhex(encoded)
    except ValueError:
        raise SEC1DecodeError('invalid_hex') from None
    if not payload:
        raise SEC1DecodeError('invalid_length')
    prefix=payload[0]
    f=curve.field
    width=(f.degree+7)//8
    if prefix==0:
        if len(payload)!=1:
            raise SEC1DecodeError('invalid_length')
        return None
    if prefix not in (2,3,4):
        raise SEC1DecodeError('invalid_prefix')
    if len(payload)!=(1+width if prefix in (2,3) else 1+2*width):
        raise SEC1DecodeError('invalid_length')
    x=int.from_bytes(payload[1:1+width],'big')
    if x>=f.q:
        raise SEC1DecodeError('noncanonical_coordinate')
    if prefix==4:
        y=int.from_bytes(payload[1+width:],'big')
        if y>=f.q:
            raise SEC1DecodeError('noncanonical_coordinate')
        point=(x,y)
        if not curve.on_curve(point):
            raise SEC1DecodeError('point_not_on_curve')
        return point
    sign=prefix&1
    if x==0:
        if sign:
            raise SEC1DecodeError('noncanonical_x0_sign')
        return (0,f.sqrt(curve.b))
    beta=x ^ curve.a ^ f.div(curve.b,f.square(x))
    if f.trace(beta):
        raise SEC1DecodeError('no_quadratic_solution')
    z=f.halftrace(beta)
    if z&1 != sign:
        z ^= 1
    point=(x,f.mul(x,z))
    if not curve.on_curve(point):
        raise AssertionError('Compressed decoding produced an invalid point')
    return point


def sec1_vectors(curve: BinaryCurve, p: Point, q: Point, g: Point, t: Point, full_generator: Point) -> dict:
    valid=[]
    for label,point in [('P',p),('Q',q),('G',g),('O',None),('T',t),('R',full_generator)]:
        compressed=sec1_encode(curve,point,True)
        uncompressed=sec1_encode(curve,point,False)
        assert sec1_decode(curve,compressed)==point
        assert sec1_decode(curve,uncompressed)==point
        valid.append({'point_label':label,'point':point,
                      'compressed':compressed,'uncompressed':uncompressed})
    f=curve.field
    width=(f.degree+7)//8
    zeros=bytes(width)
    high=f.q.to_bytes(width,'big')
    comp=sec1_encode(curve,g,True)
    uncomp=sec1_encode(curve,g,False)
    px,py=p
    off_y=py ^ (1 if px!=1 else 2)
    assert not curve.on_curve((px,off_y))
    missing_x=next(x for x in range(1,10000) if curve.point_at_x(x) is None)
    negatives=[
        ('empty','', 'invalid_length'),
        ('infinity_with_payload','0000','invalid_length'),
        ('unsupported_prefix',(b'\x05'+zeros).hex(),'invalid_prefix'),
        ('hybrid_prefix',(b'\x06'+zeros+zeros).hex(),'invalid_prefix'),
        ('compressed_short',comp[:-2],'invalid_length'),
        ('compressed_long',comp+'00','invalid_length'),
        ('uncompressed_short',uncomp[:-2],'invalid_length'),
        ('uncompressed_long',uncomp+'00','invalid_length'),
        ('compressed_high_x',(b'\x02'+high).hex(),'noncanonical_coordinate'),
        ('uncompressed_high_x',(b'\x04'+high+zeros).hex(),'noncanonical_coordinate'),
        ('uncompressed_high_y',(b'\x04'+zeros+high).hex(),'noncanonical_coordinate'),
        ('x0_prefix03',(b'\x03'+zeros).hex(),'noncanonical_x0_sign'),
        ('uncompressed_offcurve',(b'\x04'+px.to_bytes(width,'big')+off_y.to_bytes(width,'big')).hex(),'point_not_on_curve'),
        ('compressed_trace_one',(b'\x02'+missing_x.to_bytes(width,'big')).hex(),'no_quadratic_solution'),
    ]
    invalid=[]
    for label,encoded,expected in negatives:
        try:
            sec1_decode(curve,encoded)
        except SEC1DecodeError as err:
            assert err.code==expected
        else:
            raise AssertionError(f'Invalid SEC1 payload accepted: {label}')
        invalid.append({'id':label,'encoding':encoded,'error':expected})
    return {'coordinate_bytes':width,'compressed_sign_rule':'LSB(y/x) for x != 0; 0 for x=0',
            'valid':valid,'invalid':invalid}

def field_vectors(field: BinaryField) -> list[dict]:
    f = field
    a = f.reduce(f.mask ^ 0x123456789abcdef)
    b = f.reduce(0xfedcba9876543210 ^ (1 << (f.degree-1)) ^ (1 << (f.degree//2)))
    assert a and b
    wide = (1 << (2*f.degree+7)) ^ (1 << (f.degree+13)) ^ 0xdeadbeef12345
    vectors = [
        {'operation':'add', 'inputs':[a,b], 'output':f.add(a,b)},
        {'operation':'subtract', 'inputs':[a,b], 'output':f.sub(a,b)},
        {'operation':'multiply', 'inputs':[a,b], 'output':f.mul(a,b)},
        {'operation':'square', 'inputs':[a], 'output':f.square(a)},
        {'operation':'inverse', 'inputs':[a], 'output':f.inv(a)},
        {'operation':'divide', 'inputs':[a,b], 'output':f.div(a,b)},
        {'operation':'reduce', 'inputs':[wide], 'output':f.reduce(wide)},
        {'operation':'reduce_modulus', 'inputs':[f.modulus], 'output':0},
        {'operation':'frobenius', 'inputs':[a,3], 'output':f.frobenius(a,3)},
        {'operation':'frobenius_period', 'inputs':[a,f.degree], 'output':a},
        {'operation':'sqrt', 'inputs':[b], 'output':f.sqrt(b)},
        {'operation':'pow', 'inputs':[a,17], 'output':f.pow(a,17)},
        {'operation':'trace', 'inputs':[a], 'output':f.trace(a)},
        {'operation':'halftrace', 'inputs':[a], 'output':f.halftrace(a)},
    ]
    assert f.mul(a,f.inv(a))==1
    assert f.square(f.sqrt(b))==b
    assert f.square(a)==f.mul(a,a)
    assert f.square(f.halftrace(a)) ^ f.halftrace(a)==a ^ f.trace(a)
    assert f.frobenius(a,f.degree)==a
    return vectors



def point_order(curve: BinaryCurve, point: Point, group_order: int, factors: dict) -> int:
    """Exact order using the supplied fully factored rational group order."""
    if not curve.on_curve(point):
        raise ValueError('Point is not on the curve')
    if curve.mul(group_order,point) is not None:
        raise ValueError('Supplied group order does not annihilate point')
    order=group_order
    primes=sorted(int(p) for p in factors)
    for prime in primes:
        while order%prime==0 and curve.mul(order//prime,point) is None:
            order//=prime
    assert curve.mul(order,point) is None
    assert all(curve.mul(order//prime,point) is not None for prime in primes if order%prime==0)
    return order

def curve_record(curve_data: dict, field: BinaryField) -> dict:
    exponent = curve_data['b_exponent']
    b = 1 if exponent == 0 else field.pow(2,exponent)
    e = BinaryCurve(field,curve_data['a'],b)
    n, r = curve_data['N'], curve_data['r']
    d1,d2 = curve_data['group']
    assert d1*d2==n and d2%r==0
    projection = d2//r
    attempts = 0
    # Critically, projection uses the group exponent, not automatically N/r.
    for source in e.points():
        attempts += 1
        g = e.mul(projection,source)
        if g is not None and e.mul(r,g) is None:
            break
        if attempts > 10000:
            raise RuntimeError(f"Subgroup scan did not finish for {curve_data['id']}")
    else:
        raise RuntimeError('No nonzero subgroup generator found')
    assert g is not None and e.on_curve(g)
    p = source
    q = next(e.points(p[0]+1))
    t = e.point_at_x(0)
    assert t is not None
    r_attempts=0
    for full_generator in e.points():
        r_attempts+=1
        if e.mul(d2,full_generator) is not None:
            raise AssertionError('Group exponent did not annihilate point')
        if all(e.mul(d2//int(prime),full_generator) is not None
               for prime in curve_data['factors'] if d2%int(prime)==0):
            break
        if r_attempts>10000:
            raise RuntimeError(f"Maximal-order scan did not finish for {curve_data['id']}")
    else:
        raise RuntimeError('No maximal-order point found')
    assert d1==1 and d2==n  # All fourteen supplied binary rational groups are cyclic.
    named_points={'P':p,'Q':q,'G':g,'T':t,'R':full_generator,'O':None}
    orders={label:point_order(e,point,n,curve_data['factors'])
            for label,point in named_points.items()}
    assert orders['R']==d2 and orders['G']==r and orders['T']==2 and orders['O']==1
    vectors = []
    def add_vector(operation, inputs, output):
        if not e.on_curve(output):
            raise AssertionError(f'Invalid {operation} output')
        vectors.append({'operation':operation,'inputs':inputs,'output':output})
    add_vector('add',[p,q],e.add(p,q))
    add_vector('double',[p],e.double(p))
    add_vector('neg',[p],e.neg(p))
    add_vector('add',[p,e.neg(q)],e.add(p,e.neg(q)))
    width_scalar=(1 << (field.degree-1)) | (0x123456789abcdef & ((1 << (field.degree-1))-1))
    scalars=[-1,0,1,2,3,r-1,r,r+1,width_scalar,(1 << 128)+1,n]
    for scalar in scalars:
        add_vector('mul',[scalar,p],e.mul(scalar,p))
    add_vector('add',[p,None],p)
    add_vector('add',[None,p],p)
    add_vector('add',[p,e.neg(p)],None)
    add_vector('add',[p,p],e.double(p))
    add_vector('double',[t],None)
    add_vector('add',[t,t],None)
    add_vector('neg',[t],t)
    add_vector('mul',[r,t],t)
    add_vector('add',[p,t],e.add(p,t))
    for scalar in [r,n,width_scalar]:
        add_vector('mul',[scalar,full_generator],e.mul(scalar,full_generator))
    subgroup_vectors=[]
    for scalar in scalars[:-1]:
        result=e.mul(scalar,g)
        assert e.on_curve(result)
        subgroup_vectors.append({'operation':'mul','inputs':[scalar,g],'output':result})
    assert e.mul(n,p) is None and e.mul(n,q) is None
    assert e.add(p,q)==e.add(q,p)
    assert e.add(e.add(p,q),g)==e.add(p,e.add(q,g))
    assert e.double(p)==e.add(p,p)
    assert e.add(p,e.neg(p)) is None
    assert e.mul(r,g) is None
    return {
        'id':curve_data['id'], 'field':curve_data['field'],
        'equation':'y^2 + x*y = x^3 + a2*x^2 + a6',
        'coefficients':{'a1':1,'a2':e.a,'a3':0,'a4':0,'a6':e.b},
        'delta':e.b,'j':field.inv(e.b),
        'P':p,'Q':q,'R':full_generator,'two_torsion_point':t,'point_orders':orders,
        'R_scan':{'points_examined':r_attempts,'final_x':full_generator[0],
                  'rule':'ascending canonical x integers; min y; first point of exact group exponent d2',
                  'order':d2,'exact_order_prime_checks':True},
        'subgroup':{'order':r,'projection_scalar':projection,'source_point':source,
                    'generator':g,'scan_points_examined':attempts,'scan_final_x':source[0],
                    'generator_rule':'ascending canonical x integers; min y; first nonzero [d2/r]P of order r'},
        'curve_vectors':vectors,'subgroup_vectors':subgroup_vectors,
        'sec1_vectors':sec1_vectors(e,p,q,g,t,full_generator),
        'verification':{'P_and_Q_on_curve':True,'N_P_and_N_Q_are_infinity':True,
                        'r_G_is_infinity':True,'G_nonzero':True,'arithmetic_identity_checks':True,
                        'all_named_point_orders_verified':True,'R_exact_order_is_group_exponent':True},
    }


def build(source: Path=DEFAULT_SOURCE, output: Path=DEFAULT_OUTPUT) -> dict:
    normalized=json.loads(source.read_text())
    fields={}
    objects={}
    for info in normalized['fields']:
        if info['kind']!='binary':
            continue
        f=BinaryField(info['bits'],info['exponents'])
        objects[info['id']]=f
        fields[info['id']]={'id':info['id'],'degree':f.degree,'order':f.q,
            'modulus':f.modulus,'modulus_exponents':f.exponents,
            'encoding':'polynomial-basis bitmask; bit i is coefficient of u^i',
            'field_vectors':field_vectors(f)}
    curves={}
    for info in normalized['curves']:
        if info['kind']=='binary':
            curves[info['id']]=curve_record(info,objects[info['field']])
    assert len(fields)==3 and len(curves)==14
    data={'schema':'binary raw vectors v1','source':str(source),
          'point_encoding':'[x,y] raw integers; null is infinity',
          'fields':fields,'curves':curves}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(data,indent=2)+'\n')
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    data=build(args.source,args.output)
    for c in data['curves'].values():
        print(c['id'],'b=',c['coefficients']['a6'],'G=',c['subgroup']['generator'],
              'source_x=',c['subgroup']['scan_final_x'])
    print(f"Wrote {len(data['fields'])} fields and {len(data['curves'])} curves to {args.output}")

if __name__=='__main__':
    main()
