"""Weil restriction of binary Semaev systems to boolean (F2) polynomial systems,
with exact monomial algebra (c^2=c built in). Independent of the crypto repo.

Field element over GF(2^n) whose coords are F2-polynomials in boolean vars c_*:
represented as a list of n 'bpoly's. A bpoly is a set (frozenset) of monomials;
a monomial is a frozenset of integer var-ids. XOR add = symmetric_difference;
mult = XOR-accumulate of unions (c_i^2=c_i => exponent capped at 1)."""
from itertools import product, combinations
from gf2 import IRRED

def irred_coords(n):
    m=IRRED[n]; return [ (m>>i)&1 for i in range(n+1) ]  # bit i

# bpoly ops (set of frozenset monomials, over F2)
def b_add(A,B): return A ^ B
def b_mul(A,B):
    R=set()
    for ma in A:
        for mb in B:
            mon=ma|mb
            if mon in R: R.discard(mon)
            else: R.add(mon)
    return R
ZERO=frozenset  # not used
def const1(): return {frozenset()}   # the constant 1
def b_is_zero(A): return len(A)==0

# field element = list of n bpolys (coord of t^i)
def f_zero(n): return [set() for _ in range(n)]
def f_const(c,n):  # c in {0,1}
    e=f_zero(n)
    if c&1: e[0]=const1()
    return e
def f_add(X,Y,n): return [ X[i]^Y[i] for i in range(n) ]
def f_mul(X,Y,n):
    # schoolbook into degree 2n-2 then reduce mod irreducible
    acc=[set() for _ in range(2*n-1)]
    for i in range(n):
        if not X[i]: continue
        for j in range(n):
            if not Y[j]: continue
            p=b_mul(X[i],Y[j])
            acc[i+j]^=p
    ic=irred_coords(n)
    for k in range(2*n-2, n-1, -1):
        if acc[k]:
            # t^k = t^{k-n} * (t^n), t^n = sum_{i<n} ic[i] t^i
            carry=acc[k]; acc[k]=set()
            for i in range(n):
                if ic[i]:
                    acc[k-n+i]^=carry
    return acc[:n]
def f_sq(X,n): return f_mul(X,X,n)

def var_elt(prefix,i,basis,n):
    """field element x_i = sum_j c_{i,j} * basis[j], basis[j] a concrete GF(2^n) int,
    c_{i,j} boolean vars. Returns (elt, list_of_varids)."""
    l=len(basis); elt=f_zero(n); vids=[]
    for j in range(l):
        vid=(prefix,i,j); vids.append(vid)
        bj=basis[j]
        for k in range(n):
            if (bj>>k)&1:
                elt[k]=elt[k]^{frozenset([vid])}
    return elt,vids

def S3_bin_sym(X1,X2,X3,a6,n):
    e2=f_add(f_add(f_mul(X1,X2,n),f_mul(X1,X3,n),n),f_mul(X2,X3,n),n)
    e3=f_mul(f_mul(X1,X2,n),X3,n)
    return f_add(f_add(f_sq(e2,n),e3,n),f_const(a6,n),n)

def weil_equations(elt,n):
    """elt=0 in GF(2^n) <=> each of n coords = 0. Return list of bpolys."""
    return [elt[i] for i in range(n) if elt[i]]

if __name__=='__main__':
    # smoke test: S3_bin_sym with all-constant inputs must match numeric S3_bin
    from gf2 import S3_bin
    import random
    n=5; a6=1
    for _ in range(200):
        x1,x2,x3=[random.randrange(1<<n) for _ in range(3)]
        X=[]
        for x in (x1,x2,x3):
            e=f_zero(n)
            for k in range(n):
                if (x>>k)&1: e[k]=const1()
            X.append(e)
        s=S3_bin_sym(X[0],X[1],X[2],a6,n)
        val=0
        for k in range(n):
            if s[k]:  # constant coord (set containing empty monomial => 1)
                if frozenset() in s[k] and len(s[k])==1: val|=(1<<k)
                else: raise SystemExit("non-constant where constant expected")
        if val!=S3_bin(x1,x2,x3,a6,n): raise SystemExit("MISMATCH")
    print("SMOKE_OK symbolic S3 matches numeric over 200 random triples, n=5")

# ---- unsplit S4 via resultant of two quadratics (char 2), in boolean-Weil ----
def S3_coeffs_in_X(X1,X2,a6,n):
    """S3_bin(X1,X2,X) = A*X^2 + B*X + C with A=(X1+X2)^2, B=X1X2, C=(X1X2)^2+a6."""
    s=f_add(X1,X2,n); A=f_sq(s,n)
    B=f_mul(X1,X2,n)
    C=f_add(f_sq(B,n),f_const(a6,n),n)
    return A,B,C
def S4_unsplit(X1,X2,X3,X4,a6,n):
    A1,B1,C1=S3_coeffs_in_X(X1,X2,a6,n)
    A2,B2,C2=S3_coeffs_in_X(X3,X4,a6,n)
    # Res = (A1 C2 + A2 C1)^2 + (A1 B2 + A2 B1)(B1 C2 + B2 C1)  (char 2)
    t1=f_sq(f_add(f_mul(A1,C2,n),f_mul(A2,C1,n),n),n)
    t2=f_mul(f_add(f_mul(A1,B2,n),f_mul(A2,B1,n),n),
             f_add(f_mul(B1,C2,n),f_mul(B2,C1,n),n),n)
    return f_add(t1,t2,n)
