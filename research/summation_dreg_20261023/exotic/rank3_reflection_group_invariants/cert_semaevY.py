"""Brute-force certificate: S_{m+1} in Y coordinate vanishes <=> exist signs with sum eps_i P_i = +-R,
over factor-base points (x in F_p) of curves over F_{p^n}. 0 false pos / 0 false neg required."""
import random, itertools, json, sys
from lib import *
def run(p,n,m,seed):
    rng=random.Random(seed); F=find_field(p,n,rng)
    while True:
        a=tuple(rng.randrange(p) for _ in range(n)); c=rng.randrange(1,p)
        C=Curve(F,a,c)
        disc=F.sub(F.mul(a,a),F.sc(4,C.b))
        if disc!=F.zero: break
    fb=factor_base(C,rng); ok=0;fp=0;fn=0;tot=0
    pts=[pt for pt in fb]
    for tup in itertools.product(range(len(pts)),repeat=m+1):
        Ps=[pts[i] for i in tup[:m]]; R=pts[tup[m]]
        # brute: exists signs
        found=False
        for eps in itertools.product([1,-1],repeat=m):
            S=None
            for e,P in zip(eps,Ps): S=C.add(S,P if e==1 else C.neg(P))
            if S is not None and S[0]==R[0]: found=True;break
        D=eval_S(C,m,C.Y(R[0]))
        val=F.zero
        ys=[C.Y(P[0]) for P in Ps]
        for e,v in D.items():
            t=v
            for y,k in zip(ys,e): t=F.mul(t,F.pw(y,k))
            val=F.add(val,t)
        zero=(val==F.zero)
        tot+=1
        if zero and not found: fp+=1
        if found and not zero: fn+=1
    return dict(p=p,n=n,m=m,seed=seed,fb=len(pts),tuples=tot,false_pos=fp,false_neg=fn)
if __name__=='__main__':
    res=[]
    for (p,n,m,seed) in [(11,2,2,1),(13,2,2,2),(11,3,2,3),(13,3,2,4),(11,2,3,5),(13,3,3,6),(11,3,3,7),(17,2,3,8),(31,2,2,9),(23,2,3,10),(19,3,3,11)]:
        r=run(p,n,m,seed); print(r,flush=True); res.append(r)
    json.dump(res,open('raw/cert_semaevY.json','w'),indent=1)
    assert all(r['false_pos']==0 and r['false_neg']==0 for r in res); print("CERT OK")
