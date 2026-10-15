"""Independent experimental arithmetic over a quadratic prime-field extension.

Elements of F_p[u]/(u^2-d) are packed as c0+p*c1, with 0<=c0,c1<p.
This is a numerical representation, not a SEC1 or byte-encoding convention.
All routines are variable-time correctness oracles, not production code.
"""
from fractions import Fraction


def babai_decompose(scalar,subgroup_order,basis):
    """Exact nearest-plane reduction in a supplied subgroup kernel lattice."""
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    basis=[list(map(int,row)) for row in basis]
    stars=[];norms=[]
    for row in basis:
        v=list(map(Fraction,row))
        for star,norm in zip(stars,norms):
            mu=dot(row,star)/norm
            v=[x-mu*y for x,y in zip(v,star)]
        stars.append(v);norms.append(dot(v,v))
    residual=[int(scalar)%int(subgroup_order)]+[0]*(len(basis)-1)
    for row,star,norm in reversed(list(zip(basis,stars,norms))):
        coefficient=dot(residual,star)/norm
        rounded=(2*coefficient.numerator+coefficient.denominator)//(2*coefficient.denominator)
        residual=[x-rounded*y for x,y in zip(residual,row)]
    return residual


class Fp2:
    def __init__(self,p,d):
        self.p,self.d=int(p),int(d)
        self.q=self.p*self.p
        if self.p<=3 or not 0<self.d<self.p:
            raise ValueError("invalid quadratic extension parameters")
        if pow(self.d,(self.p-1)//2,self.p)!=self.p-1:
            raise ValueError("d must be a nonsquare in Fp")

    def pack(self,c0,c1=0):
        return c0%self.p+self.p*(c1%self.p)

    def unpack(self,a):
        if type(a) is not int or not 0<=a<self.q:
            raise ValueError("noncanonical packed Fp2 element")
        return a%self.p,a//self.p

    def add(self,a,b):
        a0,a1=self.unpack(a);b0,b1=self.unpack(b)
        return self.pack(a0+b0,a1+b1)

    def neg(self,a):
        a0,a1=self.unpack(a)
        return self.pack(-a0,-a1)

    def sub(self,a,b):
        return self.add(a,self.neg(b))

    def mul(self,a,b):
        a0,a1=self.unpack(a);b0,b1=self.unpack(b)
        return self.pack(a0*b0+self.d*a1*b1,a0*b1+a1*b0)

    def square(self,a):
        a0,a1=self.unpack(a)
        return self.pack(a0*a0+self.d*a1*a1,2*a0*a1)

    def inverse(self,a):
        a0,a1=self.unpack(a)
        norm=(a0*a0-self.d*a1*a1)%self.p
        if not norm:
            raise ValueError("inverse_zero")
        inverse=pow(norm,-1,self.p)
        return self.pack(a0*inverse,-a1*inverse)

    def div(self,a,b):
        return self.mul(a,self.inverse(b))

    def power(self,a,n):
        if n<0:
            return self.power(self.inverse(a),-n)
        result=1
        while n:
            if n&1:
                result=self.mul(result,a)
            a=self.square(a)
            n>>=1
        return result

    def frobenius(self,a):
        a0,a1=self.unpack(a)
        return self.pack(a0,-a1)

    def sqrt(self,a):
        """Tonelli-Shanks in Fp2; choose the smaller packed square root."""
        self.unpack(a)
        if not a:
            return 0
        minus_one=self.p-1
        if self.power(a,(self.q-1)//2)!=1:
            return None
        odd=self.q-1;s=0
        while not odd&1:
            odd>>=1;s+=1
        # Base-field elements are all squares in Fp2; start with 1+u.
        z=self.pack(1,1)
        while self.power(z,(self.q-1)//2)!=minus_one:
            z+=1
            if z>=self.q:
                raise ArithmeticError("nonsquare search exhausted")
        c=self.power(z,odd);x=self.power(a,(odd+1)//2)
        t=self.power(a,odd);m=s
        while t!=1:
            i=1;tt=self.square(t)
            while tt!=1:
                tt=self.square(tt);i+=1
                if i>=m:
                    raise ArithmeticError("Tonelli-Shanks invariant failed")
            b=self.power(c,1<<(m-i-1))
            x=self.mul(x,b);c=self.square(b);t=self.mul(t,c);m=i
        assert self.square(x)==a
        return min(x,self.neg(x))


class ExtensionCurve:
    """Short-Weierstrass affine curve over Fp2; None denotes infinity."""
    def __init__(self,field,A,B):
        self.field=field;self.A=A;self.B=B
        field.unpack(A);field.unpack(B)
        if not field.add(field.mul(field.pack(4),field.power(A,3)),field.mul(field.pack(27),field.square(B))):
            raise ValueError("singular curve")

    def on_curve(self,P):
        if P is None:return True
        try:
            x,y=P;f=self.field
            f.unpack(x);f.unpack(y)
            return f.square(y)==f.add(f.add(f.power(x,3),f.mul(self.A,x)),self.B)
        except (ValueError,TypeError):
            return False

    def neg(self,P):
        return None if P is None else (P[0],self.field.neg(P[1]))

    def add(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;z,w=Q;f=self.field
        if x==z:
            if f.add(y,w)==0:return None
            slope=f.div(f.add(f.mul(3,f.square(x)),self.A),f.mul(2,y))
        else:
            slope=f.div(f.sub(w,y),f.sub(z,x))
        xx=f.sub(f.sub(f.square(slope),x),z)
        result=(xx,f.sub(f.mul(slope,f.sub(x,xx)),y))
        assert self.on_curve(result)
        return result

    def mul(self,n,P):
        if n<0:return self.mul(-n,self.neg(P))
        result=None
        while n:
            if n&1:result=self.add(result,P)
            P=self.add(P,P);n>>=1
        return result

    def scan(self,limit=100000):
        f=self.field
        for x in range(min(limit,f.q)):
            rhs=f.add(f.add(f.power(x,3),f.mul(self.A,x)),self.B)
            y=f.sqrt(rhs)
            if y is not None:
                yield x,y

    def joint_mul(self,scalars,points):
        """Signed interleaved binary multi-scalar multiplication (Straus)."""
        if len(scalars)!=len(points):
            raise ValueError("scalar/point count mismatch")
        prepared=[(abs(int(n)),self.neg(P) if int(n)<0 else P) for n,P in zip(scalars,points)]
        result=None
        bits=max((n.bit_length() for n,_ in prepared),default=0)
        for bit in reversed(range(bits)):
            result=self.add(result,result)
            for n,P in prepared:
                if n>>bit&1:
                    result=self.add(result,P)
        return result
