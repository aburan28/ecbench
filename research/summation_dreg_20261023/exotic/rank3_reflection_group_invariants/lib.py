import sympy as sp, random, itertools, hashlib, re, subprocess, os, tempfile
from sympy.polys.polyfuncs import symmetrize
import semaev_Y as SY

# ---------- F_{p^n} = F_p[t]/(g), elements are tuples of length n ----------
class Fq:
    def __init__(s,p,n,g):
        s.p,s.n,s.g=p,n,g   # g: list of n lower coeffs of monic g (t^n = -sum g_i t^i)
        s.zero=(0,)*n; s.one=(1,)+(0,)*(n-1); s.q=p**n
    def add(s,x,y): return tuple((u+v)%s.p for u,v in zip(x,y))
    def sub(s,x,y): return tuple((u-v)%s.p for u,v in zip(x,y))
    def neg(s,x): return tuple((-u)%s.p for u in x)
    def sc(s,k,x): return tuple((k*u)%s.p for u in x)
    def mul(s,x,y):
        n,p=s.n,s.p; r=[0]*(2*n-1)
        for i,u in enumerate(x):
            if u:
                for j,v in enumerate(y): r[i+j]+=u*v
        for k in range(2*n-2,n-1,-1):
            co=r[k]%p
            if co:
                for i in range(n): r[k-n+i]-=co*s.g[i]
            r[k]=0
        return tuple(v%p for v in r[:n])
    def pw(s,x,e):
        r=s.one
        while e:
            if e&1: r=s.mul(r,x)
            x=s.mul(x,x); e>>=1
        return r
    def inv(s,x): return s.pw(x,s.q-2)
    def const(s,k): return (k%s.p,)+(0,)*(s.n-1)
    def is_sq(s,x): return x==s.zero or s.pw(x,(s.q-1)//2)==s.one
    def sqrt(s,x,rng):
        if x==s.zero: return x
        # Cipolla-free: Tonelli-Shanks
        q=s.q; Q=q-1; S=0
        while Q%2==0: Q//=2; S+=1
        while True:
            z=tuple(rng.randrange(s.p) for _ in range(s.n))
            if z!=s.zero and not s.is_sq(z): break
        M=S; cc=s.pw(z,Q); t=s.pw(x,Q); R=s.pw(x,(Q+1)//2)
        while t!=s.one:
            i=0; tt=t
            while tt!=s.one: tt=s.mul(tt,tt); i+=1
            b=cc
            for _ in range(M-i-1): b=s.mul(b,b)
            M=i; cc=s.mul(b,b); t=s.mul(t,cc); R=s.mul(R,b)
        assert s.mul(R,R)==x
        return R

def find_field(p,n,rng):
    while True:
        g=[rng.randrange(p) for _ in range(n)]
        # irreducible iff no root (n<=3)
        if all(sum(g[i]*pow(r,i,p) for i in range(n))+pow(r,n,p)!=0 for r in range(0,p) for _ in [0]) if False else \
           all((sum(g[i]*pow(r,i,p) for i in range(n))+pow(r,n,p))%p for r in range(p)):
            return Fq(p,n,g)

class Curve:
    def __init__(s,F,a,c): s.F=F; s.a=a; s.c=c; s.b=F.const(c*c)
    def f(s,x):
        F=s.F; x2=F.mul(x,x)
        return F.add(F.add(F.mul(x2,x),F.mul(s.a,x2)),F.mul(s.b,x))
    def add(s,P,Q):
        F=s.F
        if P is None: return Q
        if Q is None: return P
        (x1,y1),(x2,y2)=P,Q
        if x1==x2:
            if F.add(y1,y2)==F.zero: return None
            x12=F.mul(x1,x1)
            num=F.add(F.add(F.sc(3,x12),F.sc(2,F.mul(s.a,x1))),s.b); den=F.sc(2,y1)
        else: num=F.sub(y2,y1); den=F.sub(x2,x1)
        l=F.mul(num,F.inv(den))
        x3=F.sub(F.sub(F.sub(F.mul(l,l),s.a),x1),x2)
        y3=F.sub(F.mul(l,F.sub(x1,x3)),y1)
        return (x3,y3)
    def neg(s,P): return None if P is None else (P[0],s.F.neg(P[1]))
    def Y(s,x):
        F=s.F; cc=F.const(s.c)
        return F.mul(F.sub(x,cc),F.inv(F.add(x,cc)))

def factor_base(C,rng):
    """points with x in F_p (x != 0, +-c) and y in F_q: list of (x,y) (one of +-)."""
    F=C.F; p=F.p; out=[]
    for xv in range(p):
        if xv in (0,C.c%p,(-C.c)%p): continue
        x=F.const(xv); fx=C.f(x)
        if fx!=F.zero and F.is_sq(fx): out.append((x,F.sqrt(fx,rng)))
    return out

# ---------- polynomials in Y as dict {exp tuple: Fq element} ----------
_S={}
def semaev_dict(m):
    key='S3' if m==2 else 'S4'
    if key in _S: return _S[key]
    d=SY.load()[key]
    vars_=(SY.Y1,SY.Y2,SY.Y3) if m==2 else (SY.Y1,SY.Y2,SY.Y3,SY.YR)  # m=2: Y3 plays YR
    P=sp.Poly(d,*vars_,SY.a,SY.c)
    _S[key]=P.terms(); return _S[key]

def eval_S(C,m,YR):
    """S_{m+1}(Y_1..Y_m; YR) as dict {exp tuple(m): Fq elt}"""
    F=C.F; terms=semaev_dict(m); out={}
    apow={0:F.one}; ypow={0:F.one}
    def pw(base,tab,k):
        while k not in tab:
            j=max(tab); tab[j+1]=F.mul(tab[j],base)
        return tab[k]
    apow={0:F.one}; ypow={0:F.one}
    for mon,co in terms:
        ex=mon[:m]; ey=mon[m]; ea=mon[m+1]; ec=mon[m+2]
        v=F.sc(int(co)*pow(C.c,ec,F.p)%F.p,F.mul(pw(C.a,apow,ea),pw(YR,ypow,ey)))
        out[ex]=F.add(out.get(ex,F.zero),v)
    return {k:v for k,v in out.items() if v!=F.zero}

def pmul(F,A,B):
    out={}
    for ea,va in A.items():
        for eb,vb in B.items():
            e=tuple(x+y for x,y in zip(ea,eb)); out[e]=F.add(out.get(e,F.zero),F.mul(va,vb))
    return {k:v for k,v in out.items() if v!=F.zero}

def components(F,D):
    """descend: dict exp-> Fq elt  ->  list of n dicts exp->int mod p"""
    comps=[{} for _ in range(F.n)]
    for e,v in D.items():
        for i,u in enumerate(v):
            if u: comps[i][e]=u
    return comps

# ---------- invariant rewriting ----------
def sym_rewrite(poly_int, m):
    """poly_int: dict exp->int in m vars, symmetric. returns dict {exp over e1..em: int}"""
    ys=sp.symbols('y1:%d'%(m+1)); expr=sum(c*sp.prod([y**k for y,k in zip(ys,e)]) for e,c in poly_int.items())
    if expr==0: return {}
    sym,rem,defs=symmetrize(sp.expand(expr),*ys,formal=True)
    assert rem==0, 'not symmetric'
    ss=[d[0] for d in defs]; assert len(ss)==m
    P=sp.Poly(sym,*ss)
    return {mon:int(c) for mon,c in P.terms()}

def lift(d,p):  # centered lift not required
    return {k:int(v) for k,v in d.items()}

def to_Zvars(d, parity):
    """d: dict exp->int, all exps have given parity; return dict in Z=Y^2 exps"""
    out={}
    for e,c in d.items():
        assert all(x%2==parity for x in e), e
        out[tuple((x-parity)//2 for x in e)]=c
    return out

def modp(d,p): 
    r={}
    for k,v in d.items():
        v%=p
        if v: r[k]=v
    return r

def poly_text(d,vars_,p):
    if not d: return '0'
    terms=[]
    for e,c in sorted(d.items(),reverse=True):
        mon='*'.join(f'{v}^{k}' if k>1 else v for v,k in zip(vars_,e) if k)
        terms.append(f'{c}'+('*'+mon if mon else ''))
    return '+'.join(terms)

def tdeg(d): return max(sum(e) for e in d) if d else -1
