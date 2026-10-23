import sys, os, re, json, random, hashlib, subprocess, tempfile, itertools, time
import numpy as np
import sympy as sp
sys.path.insert(0,'/home/user/ecbench/research/summation_dreg_20261023/harness')
from semaev import ec_add, points

def sqrt_set(p): return {(y*y)%p for y in range(p)}

def S3n(x1,x2,x3,a,b,p):
    return ((x1-x2)**2*x3**2-2*((x1+x2)*(x1*x2+a)+2*b)*x3+(x1*x2-a)**2-4*b*(x1+x2))%p

def pick_instance(p,d,seed):
    rng=random.Random(seed)
    g=next(g for g in range(2,p) if all(pow(g,(p-1)//q,p)!=1 for q in sp.factorint(p-1)))
    while True:
        a=rng.randrange(1,p); b=rng.randrange(1,p)
        if (4*a**3+27*b*b)%p==0: continue
        c=pow(g,d*rng.randrange(1,(p-1)//d),p)
        V=sorted(x for x in range(p) if pow(x,d,p)==c)
        if len(V)!=d: continue
        pts=points(a,b,p)
        Pv=[P for P in pts if P and P[0] in V and P[1]!=0]
        xs={P[0] for P in Pv}
        if len(xs)<4: continue
        for _ in range(50):
            P1,P2,P3=rng.sample(Pv,3)
            if len({P1[0],P2[0],P3[0]})<3: continue
            R=ec_add(ec_add(P1,P2,a,p),P3,a,p)
            if R is None or R[0] in V: continue
            return dict(p=p,d=d,a=a,b=b,c=c,V=V,planted=[P1[0],P2[0],P3[0]],xR=R[0],pts=pts)

# ---------- symbolic systems ----------
X1,X2,X3,U,U2,V2,XR=sp.symbols('x1 x2 x3 u u2 v xR')
def S3s(A,B,C,a,b): return sp.expand((A-B)**2*C**2-2*((A+B)*(A*B+a)+2*b)*C+(A*B-a)**2-4*b*(A+B))
def vieta(A,B,s,t,a,b):
    D=(A-B)**2
    return [sp.expand(D*(s+t)-2*((A*B+a)*(A+B)+2*b)), sp.expand(D*s*t-((A*B-a)**2-4*b*(A+B)))]
def flat(A,B,s,t,a,b,tag):
    e,w,dl,sg,pi=[sp.Symbol(f'{n}{tag}') for n in ('e','w','dl','sg','pi')]
    return [e-(A+B), w-A*B, dl-(e**2-4*w), sg-(s+t), pi-s*t,
            sp.expand(dl*sg-2*((w+a)*e+2*b)), sp.expand(dl*pi-((w-a)**2-4*b*e))]
def build(inst,m):
    p,d,a,b,c,xR=inst['p'],inst['d'],inst['a'],inst['b'],inst['c'],inst['xR']
    xs=[X1,X2,X3][:m]
    field=[x**d-c for x in xs]
    xr=sp.Integer(xR)
    out={}
    if m==3:
        S4=sp.resultant(S3s(X1,X2,U,a,b),S3s(U,X3,xr,a,b),U)
        out['a']=[sp.expand(S4)]+field
        out['b']=[S3s(X1,X2,U,a,b),S3s(U,X3,xr,a,b)]+field
        out['c1']=vieta(X1,X2,U,U2,a,b)+vieta(U,X3,xr,V2,a,b)+field
        out['c2']=flat(X1,X2,U,U2,a,b,'A')+flat(U,X3,xr,V2,a,b,'B')+field
    else:
        out['a']=[S3s(X1,X2,xr,a,b)]+field
        out['b']=out['a']
        out['c1']=vieta(X1,X2,xr,V2,a,b)+field
        out['c2']=flat(X1,X2,xr,V2,a,b,'A')+field
    res={}
    for k,eqs in out.items():
        eqs=[sp.Poly(e,*sorted(set().union(*[q.free_symbols for q in eqs]),key=str),modulus=p).as_expr() if True else e for e in eqs]
        res[k]=eqs
    return res

def varsof(eqs): return sorted(set().union(*[e.free_symbols for e in eqs]),key=str)
def tdeg(e,vs,p): return sp.Poly(e,*vs,modulus=p).total_degree()

def fmt(e,p): 
    s=str(sp.expand(e)).replace('**','^')
    return s
def to_ms(eqs,vs,p): return ",".join(str(v) for v in vs)+f"\n{p}\n"+",\n".join(fmt(e,p) for e in eqs)+"\n"
def to_sing(eqs,vs,p,extra=None,tail=""):
    gens=[fmt(e,p) for e in eqs]+(extra or [])
    return f"ring r={p},({','.join(str(v) for v in vs)}),dp;\noption(prot);\nideal I={','.join(gens)};\n"+tail

def run(cmd,inp=None):
    c=['nice','-n','10','timeout','300','prlimit','--as=2500000000']+cmd
    try:
        r=subprocess.run(c,input=inp,capture_output=True,text=True,timeout=320)
        return r.returncode,r.stdout+r.stderr
    except subprocess.TimeoutExpired: return 124,''

def msolve_deg(eqs,vs,p):
    with tempfile.NamedTemporaryFile('w',suffix='.ms',delete=False) as fh:
        fh.write(to_ms(eqs,vs,p)); path=fh.name
    rc,out=run(['msolve','-t','1','-v','2','-f',path,'-g','2','-o','/dev/null'])
    os.unlink(path)
    degs=[int(m.group(1)) for l in out.splitlines() for m in [re.match(r'\s*(\d+)\s+\d+\s+\d+\s+\d+ x \d+',l)] if m]
    return (max(degs) if degs else None), rc

def singular(eqs,vs,p,planted):
    # d_solv(Singular) := least degBound D such that std(I) with degBound=D already has vdim == vdim of full std (rank-3 metric)
    ex=",".join(f"{v}-{planted[str(v)]}" for v in vs if str(v) in planted)
    tail=('degBound=0;ideal G=std(I);int V=vdim(G);print("VDIM "+string(V));print("DIM "+string(dim(G)));\n'
          +(f'ideal Mp={ex};ideal H=std(G+Mp);print("PLANTED "+string(H[1]!=1));\n' if ex else '')
          +'int d=1;int found=0;ideal Gd;\nwhile(found==0 and d<=60){degBound=d;Gd=std(I);if(vdim(Gd)==V){found=d;}d=d+1;}\nprint("DSOLV "+string(found));\n')
    rc,out=run(['Singular','-q'],to_sing(eqs,vs,p,tail=tail).replace('option(prot);\n','')+'$\n')
    vd=re.search(r'VDIM (-?\d+)',out); dm=re.search(r'DIM (-?\d+)',out); ds=re.search(r'DSOLV (\d+)',out); pl=re.search(r'PLANTED (\d)',out)
    return dict(sing_deg=(int(ds.group(1)) or None) if ds else None, vdim=int(vd.group(1)) if vd else None,
                dim=int(dm.group(1)) if dm else None, sing_rc=rc,
                planted_ok=bool(int(pl.group(1))) if pl else None)

def bezout(eqs,vs,p):
    ds=sorted(tdeg(e,vs,p) for e in eqs)
    n=len(vs); sm=ds[:n] if len(ds)>=n else ds
    return int(np.prod(sm,dtype=object)), int(np.prod(ds,dtype=object)), ds

# ---------- brute-force correctness ----------
def truth(inst):
    p,a,V=inst['p'],inst['a'],inst['V']; pts=inst['pts']; xR=inst['xR']
    Pv=[P for P in pts if P and P[0] in V]
    neg=lambda P:(P[0],(-P[1])%p)
    T=set()
    for P1 in Pv:
      for P2 in Pv:
        for P3 in Pv:
          for s2 in (0,1):
            for s3 in (0,1):
              Q2=neg(P2) if s2 else P2; Q3=neg(P3) if s3 else P3
              S=ec_add(ec_add(P1,Q2,a,p),Q3,a,p)
              if S is not None and S[0]==xR: T.add((P1[0],P2[0],P3[0]))
    return T
def fsols(inst,m):
    """F_p-solvability of a,b,c1 over V^m, for on-curve x's; also for c2 (equiv to c1 by algebra, checked numerically)."""
    p,a,b,V,xR=inst['p'],inst['a'],inst['b'],inst['V'],inst['xR']
    fset={x for x in range(p) if pow((x**3+a*x+b)%p,(p-1)//2,p)==1 or (x**3+a*x+b)%p==0}
    on=[x for x in V if x in fset]
    allp=np.arange(p)
    def linkpairs(A,B):  # all (s,t) in F_p^2 with Vieta eqs
        D=(A-B)**2%p; N=2*((A*B+a)*(A+B)+2*b)%p; Bq=((A*B-a)**2-4*b*(A+B))%p
        sols=[]
        if D:
            Di=pow(D,-1,p); sg=N*Di%p; pr=Bq*Di%p
            for s in range(p):
                t=(sg-s)%p
                if s*t%p==pr: sols.append((s,t))
        else:
            if N==0 and Bq==0: sols=[(s,t) for s in range(p) for t in range(p)]
        return sols
    def link_ok_const(A,B,xr):  # exists t with Vieta for (s=xr)
        return any(s==xr for s,t in linkpairs(A,B))
    cache={}
    def lp(A,B):
        k=(min(A,B),max(A,B))
        if k not in cache: cache[k]=linkpairs(*k)
        return cache[k]
    A_,C1=set(),set()
    for x1 in on:
      for x2 in on:
        if m==2:
            if S3n(x1,x2,xR,a,b,p)==0: A_.add((x1,x2))
            if any(s==xR for s,t in lp(x1,x2)): C1.add((x1,x2))
            continue
        for x3 in on:
            # b: exists u in F_p
            okb=any(S3n(x1,x2,u,a,b,p)==0 and S3n(u,x3,xR,a,b,p)==0 for u in range(p))
            if okb: A_.add((x1,x2,x3))
            okc=any(any(s==xR for s,t in lp(u,x3)) for (u,u2) in lp(x1,x2))
            if okc: C1.add((x1,x2,x3))
    return A_,C1
def s4_zero(inst,m):
    p,a,b,V,xR=inst['p'],inst['a'],inst['b'],inst['V'],inst['xR']
    if m==2: return None
    S4=sp.resultant(S3s(X1,X2,U,a,b),S3s(U,X3,sp.Integer(xR),a,b),U)
    f=sp.lambdify((X1,X2,X3),sp.expand(S4),'math')
    fset={x for x in range(p) if pow((x**3+a*x+b)%p,(p-1)//2,p)==1}
    on=[x for x in V if x in fset]
    return {(x1,x2,x3) for x1 in on for x2 in on for x3 in on if int(f(x1,x2,x3))%p==0}

def correctness(inst,m):
    T=truth(inst)
    if m==2:
        # truth for m=2: x(+-P1+-P2)=xR
        p,a,V=inst['p'],inst['a'],inst['V']; Pv=[P for P in inst['pts'] if P and P[0] in V]
        T=set()
        for P1 in Pv:
            for P2 in Pv:
                for s in (0,1):
                    Q=(P2[0],(-P2[1])%p) if s else P2
                    S=ec_add(P1,Q,a,p)
                    if S and S[0]==inst['xR']: T.add((P1[0],P2[0]))
    A_,C1=fsols(inst,m)
    S4z=s4_zero(inst,m)
    ontrue=T
    def cmp(S,name):
        return dict(name=name,n=len(S),false_pos=len(S-ontrue),false_neg=len(ontrue-S))
    rep=[cmp(A_,'a/b F_p'),cmp(C1,'c1 F_p')]
    if S4z is not None: rep.append(cmp(S4z,'a (S4 vanishing)'))
    return dict(truth_n=len(T),checks=rep)
