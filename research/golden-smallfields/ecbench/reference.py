"""Independent, variable-time, standard-library reference arithmetic.

This code is a correctness oracle and benchmark baseline, not a production
cryptographic implementation. Scalar multiplication never reduces a scalar by
the selected subgroup order when operating on an arbitrary curve point.
"""
from math import isqrt
import re


class ECError(ValueError):
    def __init__(self, code):
        self.code = code
        super().__init__(code)


class Field:
    def __init__(self, record):
        self.record = record
        self.kind = record['kind']
        self.q = int(record['cardinality'])
        self.bits = record['element_bits']
        self.width = record['coordinate_bytes']
        self.modulus = int(record['representation']['modulus'], 16)

    def canonical(self, a):
        return type(a) is int and 0 <= a < self.q

    def check_elements(self, *values):
        if not all(self.canonical(value) for value in values):
            raise ECError('noncanonical_element')

    def reduce(self, a):
        if type(a) is not int:
            raise ECError('invalid_integer')
        if self.kind == 'prime':
            return a % self.q
        if a < 0:
            raise ECError('negative_polynomial')
        while a.bit_length() > self.bits:
            a ^= self.modulus << (a.bit_length() - self.bits - 1)
        return a

    def add(self, a, b):
        self.check_elements(a,b)
        return (a+b) % self.q if self.kind == 'prime' else a ^ b

    def sub(self, a, b):
        self.check_elements(a,b)
        return (a-b) % self.q if self.kind == 'prime' else a ^ b

    def mul(self, a, b):
        self.check_elements(a,b)
        if self.kind == 'prime':
            return a*b % self.q
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= self.modulus
        return out

    def square(self, a):
        self.check_elements(a)
        if self.kind == 'prime':
            return a*a % self.q
        expanded = 0
        for i in range(a.bit_length()):
            if a >> i & 1:
                expanded |= 1 << (2*i)
        return self.reduce(expanded)

    def inverse(self, a):
        self.check_elements(a)
        if not a:
            raise ECError('inverse_zero')
        if self.kind == 'prime':
            return pow(a, -1, self.q)
        # Polynomial extended Euclid, independent of the vector generators.
        u, v, s, t = a, self.modulus, 1, 0
        while u != 1:
            if not u:
                raise ECError('noninvertible')
            shift = u.bit_length()-v.bit_length()
            if shift < 0:
                u, v, s, t = v, u, t, s
                shift = -shift
            u ^= v << shift
            s ^= t << shift
        return self.reduce(s)

    def div(self, a, b):
        return self.mul(a, self.inverse(b))

    def power(self, a, exponent):
        self.check_elements(a)
        if type(exponent) is not int:
            raise ECError('invalid_integer')
        if exponent < 0:
            a, exponent = self.inverse(a), -exponent
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            a = self.square(a)
            exponent >>= 1
        return result

    def frobenius(self, a, exponent):
        self.check_elements(a)
        if type(exponent) is not int:
            raise ECError('invalid_integer')
        if self.kind != 'binary' or exponent < 0:
            raise ECError('unsupported_operation')
        for _ in range(exponent):
            a = self.square(a)
        return a

    def trace(self, a):
        self.check_elements(a)
        if self.kind != 'binary':
            raise ECError('unsupported_operation')
        result = 0
        for _ in range(self.bits):
            result ^= a
            a = self.square(a)
        return result

    def halftrace(self, a):
        self.check_elements(a)
        if self.kind != 'binary' or self.bits % 2 != 1:
            raise ECError('unsupported_operation')
        result, term = 0, a
        for _ in range((self.bits+1)//2):
            result ^= term
            term = self.square(self.square(term))
        return result

    def sqrt(self, a):
        self.check_elements(a)
        if self.kind == 'binary':
            return self.frobenius(a, self.bits-1)
        if a == 0:
            return 0
        p = self.q
        if pow(a, (p-1)//2, p) != 1:
            raise ECError('point_not_on_curve')
        if p % 4 == 3:
            root = pow(a, (p+1)//4, p)
        else:
            odd, s = p-1, 0
            while odd % 2 == 0:
                odd //= 2
                s += 1
            z = 2
            while pow(z, (p-1)//2, p) != p-1:
                z += 1
            c, root, residue, m = pow(z, odd, p), pow(a, (odd+1)//2, p), pow(a, odd, p), s
            while residue != 1:
                i, trial = 0, residue
                while trial != 1:
                    trial = trial*trial % p
                    i += 1
                    if i == m:
                        raise ECError('point_not_on_curve')
                b = pow(c, 1 << (m-i-1), p)
                root = root*b % p
                c = b*b % p
                residue = residue*c % p
                m = i
        return min(root, p-root)

    def hex(self, a):
        if not self.canonical(a):
            raise ECError('noncanonical_coordinate')
        return f'0x{a:0{2*self.width}x}'

    def parse(self, a):
        if not isinstance(a, str) or not re.fullmatch('0x[0-9a-f]{'+str(2*self.width)+'}', a):
            raise ECError('noncanonical_storage')
        value = int(a, 16)
        if not self.canonical(value):
            raise ECError('noncanonical_coordinate')
        return value


class Curve:
    def __init__(self, record, field):
        self.record, self.field = record, field
        self.coefficients = {k:field.parse(v) for k,v in record['coefficients'].items()}
        self.a = self.coefficients['a4'] if field.kind == 'prime' else self.coefficients['a2']
        self.b = self.coefficients['a6']
        self.r = int(record['subgroup']['order'])
        self.G = self.parse_point(record['subgroup']['generator'])
        self.rank = record['subgroup']['r_torsion_rank']
        self._bsgs = None

    def parse_point(self, point):
        if point is None:
            return None
        if not isinstance(point, dict) or set(point) != {'x','y'}:
            raise ECError('invalid_point_shape')
        return self.field.parse(point['x']), self.field.parse(point['y'])

    def point(self, p):
        if p is None:
            return None
        return {'x':self.field.hex(p[0]),'y':self.field.hex(p[1])}

    def on_curve(self, p):
        if p is None:
            return True
        if not isinstance(p, (tuple,list)) or len(p) != 2 or not all(self.field.canonical(a) for a in p):
            return False
        f = self.field
        x,y = p
        x2,y2 = f.square(x),f.square(y)
        if f.kind == 'prime':
            return y2 == f.add(f.add(f.mul(x2,x),f.mul(self.a,x)),self.b)
        return f.add(y2,f.mul(x,y)) == f.add(f.add(f.mul(x2,x),f.mul(self.a,x2)),self.b)

    def check(self, p):
        if not self.on_curve(p):
            raise ECError('point_not_on_curve')

    def _neg(self, p):
        if p is None:
            return None
        return (p[0],(-p[1]) % self.field.q) if self.field.kind == 'prime' else (p[0],p[0]^p[1])

    def neg(self, p):
        self.check(p)
        return self._neg(p)

    def _add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        f = self.field
        x,y = p;u,v = q
        if f.kind == 'prime':
            if x == u:
                if (y+v) % f.q == 0:
                    return None
                slope = (3*x*x+self.a)*pow(2*y,-1,f.q) % f.q
            else:
                slope = (v-y)*pow(u-x,-1,f.q) % f.q
            outx = (slope*slope-x-u) % f.q
            return outx,(slope*(x-outx)-y) % f.q
        if x == u:
            if y ^ v == x:
                return None
            slope = x ^ f.div(y,x)
            outx = f.square(slope) ^ slope ^ self.a
            return outx,f.square(x) ^ f.mul(slope^1,outx)
        slope = f.div(y^v,x^u)
        outx = f.square(slope) ^ slope ^ x ^ u ^ self.a
        return outx,f.mul(slope,x^outx) ^ outx ^ y

    def add(self, p, q):
        self.check(p);self.check(q)
        p = tuple(p) if p is not None else None
        q = tuple(q) if q is not None else None
        return self._add(p,q)

    def double(self, p):
        return self.add(p,p)

    def mul(self, scalar, p):
        if type(scalar) is not int:
            raise ECError('invalid_scalar')
        self.check(p)
        p = tuple(p) if p is not None else None
        if scalar < 0:
            scalar,p = -scalar,self._neg(p)
        # No reduction modulo r: p can be outside the selected subgroup.
        out = None
        while scalar:
            if scalar & 1:
                out = self._add(out,p)
            p = self._add(p,p)
            scalar >>= 1
        return out

    def r_torsion(self, p):
        return self.on_curve(p) and self.mul(self.r,p) is None

    def subgroup_contains(self, p):
        if not self.r_torsion(p):
            return False
        if self.rank == 1:
            return True
        p = tuple(p) if p is not None else None
        # The dataset's three rank-two cases have r<=65537. Exact BSGS
        # distinguishes the chosen cyclic subgroup from all rational r-torsion.
        if self.r > 1000000:
            raise ECError('membership_search_limit')
        if self._bsgs is None:
            step = isqrt(self.r)+1
            baby, value = {},None
            for j in range(step):
                baby.setdefault(value,j)
                value = self._add(value,self.G)
            self._bsgs = step,baby,self._neg(self.mul(step,self.G))
        step,baby,giant = self._bsgs
        value = p
        for i in range(step+1):
            if value in baby and i*step+baby[value] < self.r:
                return True
            value = self._add(value,giant)
        return False

    def public_key_valid(self, p):
        return p is not None and self.subgroup_contains(p)

    def map_point(self, operation, p, parameter=None):
        """Apply one explicitly supplied isogeny/automorphism to a source point."""
        self.check(p)
        f=self.field
        if operation=='degree2_isogeny':
            if f.kind!='prime' or self.a!=f.q-1 or self.b!=0:
                raise ECError('map_source_mismatch')
            if p is None or p[0]==0:return None
            x,y=p;inverse=f.inverse(x)
            return f.sub(x,inverse),f.mul(y,f.add(1,f.square(inverse)))
        if operation=='relative_frobenius':
            if f.kind!='binary':raise ECError('map_source_mismatch')
            return None if p is None else (f.square(p[0]),f.square(p[1]))
        if operation=='gaussian_automorphism':
            if f.kind!='prime' or self.b!=0 or not f.canonical(parameter) or f.square(parameter)!=f.q-1:
                raise ECError('map_source_mismatch')
            return None if p is None else ((-p[0])%f.q,f.mul(parameter,p[1]))
        if operation=='eisenstein_automorphism':
            if f.kind!='prime' or self.a!=0 or not f.canonical(parameter) or parameter==1 or f.power(parameter,3)!=1:
                raise ECError('map_source_mismatch')
            return None if p is None else (f.mul(parameter,p[0]),p[1])
        raise ECError('unsupported_operation')

    def encode(self, p, compressed=True):
        self.check(p)
        if p is None:
            return '00'
        f = self.field
        x,y = p
        xb,yb = x.to_bytes(f.width,'big'),y.to_bytes(f.width,'big')
        if not compressed:
            return (b'\x04'+xb+yb).hex()
        if f.kind == 'prime':
            bit = y & 1
        else:
            bit = f.div(y,x) & 1 if x else 0
        return (bytes([2+bit])+xb).hex()

    def decode(self, octets):
        if not isinstance(octets,str) or not re.fullmatch(r'(?:[0-9a-f]{2})*',octets):
            raise ECError('invalid_hex')
        raw = bytes.fromhex(octets)
        if not raw:
            raise ECError('invalid_length')
        f = self.field
        prefix = raw[0]
        if prefix == 0:
            if len(raw) != 1:
                raise ECError('invalid_length')
            return None
        if prefix not in (2,3,4):
            raise ECError('invalid_prefix')
        length = 1+(2 if prefix==4 else 1)*f.width
        if len(raw) != length:
            raise ECError('invalid_length')
        x = int.from_bytes(raw[1:1+f.width],'big')
        if not f.canonical(x):
            raise ECError('noncanonical_coordinate')
        if prefix == 4:
            y = int.from_bytes(raw[1+f.width:],'big')
            if not f.canonical(y):
                raise ECError('noncanonical_coordinate')
        elif f.kind == 'prime':
            rhs = (x*x*x+self.a*x+self.b) % f.q
            y = f.sqrt(rhs)
            if y & 1 != prefix & 1:
                y = (-y) % f.q
            if y & 1 != prefix & 1:
                raise ECError('invalid_y_selector')
        elif not x:
            if prefix & 1:
                raise ECError('invalid_y_selector')
            y = f.sqrt(self.b)
        else:
            rhs = x ^ self.a ^ f.div(self.b,f.square(x))
            if f.trace(rhs):
                raise ECError('point_not_on_curve')
            z = f.halftrace(rhs)
            if z & 1 != prefix & 1:
                z ^= 1
            y = f.mul(x,z)
        p = x,y
        self.check(p)
        return p
