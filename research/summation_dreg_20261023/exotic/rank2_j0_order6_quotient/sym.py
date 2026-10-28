"""Symbolic gate: semi-invariance of S3 (a=0) and u-polynomiality of the orbit product."""
import sympy as sp, json, sys
x1,x2,X,b,w,u1,u2=sp.symbols('x1 x2 X b w u1 u2')
def S3(a1,a2,a3):
    return (a1-a2)**2*a3**2-2*((a1+a2)*a1*a2+2*b)*a3+(a1*a2)**2-4*b*(a1+a2)
cyc=w**2+w+1
def red(e): return sp.rem(sp.expand(e),cyc,w)
# 1. semi-invariance
si=red(S3(w*x1,w*x2,w*X)-w*S3(x1,x2,X))
print("semi-invariance residue:",sp.simplify(si))
assert si==0
# 2. orbit product
Q=1
for k1 in range(3):
    for k2 in range(3):
        Q=red(sp.expand(Q*S3(w**k1*x1,w**k2*x2,X)))
        Q=sp.expand(Q)
P=sp.Poly(Q,x1,x2)
assert not Q.has(w), "w remains"
for (i,j),c in P.terms():
    assert i%3==0 and j%3==0,(i,j)
Qu=sp.Poly(sp.expand(sum(c*u1**(i//3)*u2**(j//3) for (i,j),c in P.terms())),u1,u2)
print("Q degs in u1,u2:",Qu.degree(u1),Qu.degree(u2),"total",Qu.total_degree(),"terms",len(Qu.terms()))
# 3. target-orbit variant, checked numerically over F_p (w in F_p): Q27 == unit*Q^3
import random
random.seed(1)
ok=True
for p in (61,97,1009,4093,65521):
    wv=[pow(g,(p-1)//3,p) for g in range(2,p) if pow(g,(p-1)//3,p)!=1][0]
    for _ in range(30):
        vx1,vx2,vX,vb=[random.randrange(1,p) for _ in range(4)]
        def s3(a1,a2,a3): return ((a1-a2)**2*a3**2-2*((a1+a2)*a1*a2+2*vb)*a3+(a1*a2)**2-4*vb*(a1+a2))%p
        Qv=1
        for k1 in range(3):
            for k2 in range(3): Qv=Qv*s3(pow(wv,k1,p)*vx1%p,pow(wv,k2,p)*vx2%p,vX)%p
        Q27=1
        for k in range(3):
            for k1 in range(3):
                for k2 in range(3): Q27=Q27*s3(pow(wv,k1,p)*vx1%p,pow(wv,k2,p)*vx2%p,pow(wv,k,p)*vX%p)%p
        ok&= any(Q27==pow(wv,e,p)*pow(Qv,3,p)%p for e in range(3))
diff=0 if ok else 1; diff2=diff3=1
# export: dict (i,j)-> poly in X,b as string
out={}
for (i,j),c in Qu.terms():
    out[f"{i},{j}"]=str(sp.expand(c))
json.dump(out,open('raw/Q_u_poly.json','w'),indent=1,sort_keys=True)
json.dump(dict(semi_invariant=True,u_polynomial=True,deg_u1=int(Qu.degree(u1)),deg_u2=int(Qu.degree(u2)),total=int(Qu.total_degree()),
  target_orbit_is_Q_cubed=bool(diff==0 or diff2==0 or diff3==0)),open('raw/sym_gate.json','w'))
