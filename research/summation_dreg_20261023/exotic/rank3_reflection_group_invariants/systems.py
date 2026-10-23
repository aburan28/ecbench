import random, itertools, hashlib, json
from lib import *

def make_instance(p,n,m,seed):
    rng=random.Random(seed); F=find_field(p,n,rng)
    while True:
        a=tuple(rng.randrange(p) for _ in range(n))
        if n>1 and all(v==0 for v in a[1:]): continue   # a not in F_p
        c=rng.randrange(1,p); C=Curve(F,a,c)
        if F.sub(F.mul(a,a),F.sc(4,C.b))==F.zero: continue
        fb=factor_base(C,rng)
        if len(fb)<m+1: continue
        for _ in range(100):
            pick=rng.sample(fb,m)
            pick=[(x,y if rng.random()<.5 else F.neg(y)) for x,y in pick]
            R=None
            for P in pick: R=C.add(R,P)
            if R is None or R[0]==F.const(-c): continue
            if n>1 and all(v==0 for v in C.Y(R[0])[1:]): continue   # Y_R in F_p degenerate (descent loses equations)
            return F,C,pick,R
def evalpoly(F,D,ys):
    val=F.zero
    for e,v in D.items():
        t=v
        for y,k in zip(ys,e): t=F.mul(t,F.pw(y,k))
        val=F.add(val,t)
    return val

def build_systems(p,n,m,seed):
    F,C,pick,R=make_instance(p,n,m,seed)
    YR=C.Y(R[0]); yint=[C.Y(P[0])[0] for P in pick]
    ys=[F.const(v) for v in yint]; assert all(C.Y(P[0])==F.const(C.Y(P[0])[0]) for P in pick)
    D=eval_S(C,m,YR); assert evalpoly(F,D,ys)==F.zero, "planted not a root"
    comps=[modp(d,p) for d in components(F,D)]
    out={}
    # planted values
    Ys=yint
    ek=lambda vals:[sum(__import__('math').prod(c) for c in itertools.combinations(vals,k))%p for k in range(1,m+1)]
    prod=1
    for v in Ys: prod=prod*v%p
    pl_a=Ys; pl_b=ek(Ys); pl_u=ek([v*v%p for v in Ys])
    # (a)
    out['a']=dict(vars=[f'y{i+1}' for i in range(m)],eqs=comps,planted=pl_a)
    # (b)
    out['b']=dict(vars=[f'e{i+1}' for i in range(m)],eqs=[modp(sym_rewrite(c,m),p) for c in comps],planted=pl_b)
    # (c1),(c1')
    c1=[];c1p=[]
    for comp in comps:
        ev={e:v for e,v in comp.items() if all(x%2==0 for x in e)}
        od={e:v for e,v in comp.items() if all(x%2==1 for x in e)}
        assert len(ev)+len(od)==len(comp), 'not D_m invariant'
        g=sym_rewrite(to_Zvars(ev,0),m) if ev else {}
        h=sym_rewrite(to_Zvars(od,1),m) if od else {}
        full={}   # exps (u1..um,pi)
        for e,v in g.items(): full[e+(0,)]=(full.get(e+(0,),0)+v)%p
        for e,v in h.items(): full[e+(1,)]=(full.get(e+(1,),0)+v)%p
        full=modp(full,p); c1p.append(full)
        red={}
        for e,v in full.items():
            ee=e[:m-1]+(e[m]+2*e[m-1],)
            red[ee]=(red.get(ee,0)+v)%p
        c1.append(modp(red,p))
    out['c1']=dict(vars=[f'u{i+1}' for i in range(m-1)]+['pi'],eqs=c1,planted=pl_u[:m-1]+[prod])
    rel={tuple([0]*(m-1)+[1,0]):p-1, (0,)*m+(2,):1}   # pi^2 - u_m
    out['c1p']=dict(vars=[f'u{i+1}' for i in range(m)]+['pi'],eqs=c1p+[modp(rel,p)],planted=pl_u+[prod])
    # (c2) B_m product
    D2=eval_S(C,m,F.neg(YR)); Pm=pmul(F,D,D2)
    c2=[]
    for comp in components(F,Pm):
        comp=modp(comp,p); assert all(all(x%2==0 for x in e) for e in comp), 'product not even'
        c2.append(modp(sym_rewrite(to_Zvars(comp,0),m),p))
    out['c2']=dict(vars=[f'u{i+1}' for i in range(m)],eqs=c2,planted=pl_u)
    # check planted satisfies all
    for k,s in out.items():
        for eq in s['eqs']:
            val=0
            for e,v in eq.items():
                t=v
                for x,kx in zip(s['planted'],e): t=t*pow(x,kx,p)%p
                val=(val+t)%p
            assert val==0,(k,'planted not root')
    meta=dict(p=p,n=n,m=m,seed=seed,field_g=F.g,a=C.a,c=C.c,planted_Y=yint,YR=YR)
    return out,meta
