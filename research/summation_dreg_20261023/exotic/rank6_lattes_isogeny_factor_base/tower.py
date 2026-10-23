import sys, random, subprocess, re, os, json, hashlib, tempfile, math
sys.path.insert(0,'/home/user/ecbench/research/summation_dreg_20261023/harness')
import sympy as sp
from semaev import points, ec_add, S3_closed

# ---------- rational maps over F_p: (num, den) as lists of coeffs low->high mod p ----------
def trim(a):
    while a and a[-1]==0: a.pop()
    return a
def padd(a,b,p):
    n=max(len(a),len(b)); return trim([((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0))%p for i in range(n)])
def pmul(a,b,p):
    if not a or not b: return []
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%p
    return trim(r)
def pscal(a,c,p): return trim([x*c%p for x in a])
def ppow(a,e,p):
    r=[1]
    for _ in range(e): r=pmul(r,a,p)
    return r
def peval(a,x,p):
    r=0
    for c in reversed(a): r=(r*x+c)%p
    return r
def pgcd_deg(a,b,p):
    F=sp.GF(p)
    A=sp.Poly(list(reversed(a)),sp.Symbol('x'),modulus=p); B=sp.Poly(list(reversed(b)),sp.Symbol('x'),modulus=p)
    return sp.gcd(A,B).degree()
def rdeg(f,p):
    return max(len(f[0]),len(f[1]))-1
def compose(f,g,p):
    """f o g, f=(N1,D1), g=(N2,D2): returns (N,D) = f(g)."""
    N1,D1=f; N2,D2=g; d=max(len(N1),len(D1))-1
    def hom(c):
        r=[]
        for i,co in enumerate(c):
            t=pscal(pmul(ppow(N2,i,p),ppow(D2,d-i,p),p),co,p); r=padd(r,t,p)
        return r
    N=hom(N1); D=hom(D1)
    g_=pgcd_deg(N,D,p)
    if g_>0:
        A=sp.Poly(list(reversed(N)),sp.Symbol('x'),modulus=p); B=sp.Poly(list(reversed(D)),sp.Symbol('x'),modulus=p)
        G=sp.gcd(A,B); A=sp.quo(A,G); B=sp.quo(B,G)
        N=[int(c)%p for c in reversed(A.all_coeffs())]; D=[int(c)%p for c in reversed(B.all_coeffs())]
    return (trim(N),trim(D))
def rapply(f,x,p):
    d=peval(f[1],x,p)
    return None if d==0 else peval(f[0],x,p)*pow(d,-1,p)%p   # None = infinity

# ---------- layer families ----------
def kummer(d,p): return ([0]*d+[1],[1])
def generic(d,p,rng):
    while True:
        N=[rng.randrange(p) for _ in range(d+1)]; D=[rng.randrange(p) for _ in range(d+1)]
        if N[d]==0 and D[d]==0: continue
        N=trim(N); D=trim(D)
        if not N or not D: continue
        if max(len(N),len(D))-1==d and pgcd_deg(N,D,p)==0: return (N,D)
def lattes(a,b,p):
    # x([2]) on y^2=x^3+ax+b: (x^4-2a x^2-8b x+a^2)/(4(x^3+ax+b))
    return ([a*a%p,(-8*b)%p,(-2*a)%p,0,1],[4*b%p,4*a%p,0,4])
def rand_curve(p,rng):
    while True:
        a,b=rng.randrange(p),rng.randrange(p)
        if (4*a**3+27*b*b)%p: return a,b
def velu2(a,b,x0,p):
    t=(3*x0*x0+a)%p; w=x0*t%p
    f=([(x0*1 - 0)*0+ (-x0*x0+t*0)%p*0,0],[1])  # placeholder (overwritten)
    # x' = x + t/(x-x0) = (x^2 - x0 x + t)/(x - x0)
    return (( (t)%p, (-x0)%p, 1), ([(-x0)%p,1])), ((a-5*t)%p,(b-7*w)%p)
def velu3(a,b,x0,p):
    t=(6*x0*x0+2*a)%p; u=4*(x0**3+a*x0+b)%p
    # x' = x + t/(x-x0) + u/(x-x0)^2 = ((x)(x-x0)^2 + t(x-x0) + u)/(x-x0)^2
    xm=[(-x0)%p,1]; xm2=pmul(xm,xm,p)
    N=padd(padd(pmul([0,1],xm2,p),pscal(xm,t,p),p),[u],p)
    return (N,xm2),((a-5*t)%p,(b-7*(u+x0*t))%p)
def roots_2(a,b,p): return [x for x in range(p) if (x**3+a*x+b)%p==0]
def roots_3(a,b,p): return [x for x in range(p) if (3*x**4+6*a*x*x+12*b*x-a*a)%p==0]
def isogeny_chain(ell,n,p,rng,tries=20000):
    """chain of n rational ell-isogenies via Velu from a random curve; returns (list of maps, curves)."""
    for _ in range(tries):
        a,b=rand_curve(p,rng); maps=[]; ok=True; cur=(a,b); curves=[cur]
        for _i in range(n):
            rts=roots_2(*cur,p) if ell==2 else roots_3(*cur,p)
            if not rts: ok=False;break
            x0=rng.choice(rts)
            f,cur=(velu2 if ell==2 else velu3)(cur[0],cur[1],x0,p)
            if (4*cur[0]**3+27*cur[1]**2)%p==0: ok=False;break
            maps.append((trim(list(f[0])),trim(list(f[1])))); curves.append(cur)
        if not ok: continue
        # check composite has full degree ell^n (no backtracking)
        F=maps[0]
        for m in maps[1:]: F=compose(m,F,p)
        if rdeg(F,p)==ell**n: return maps,curves
    raise RuntimeError('no isogeny chain')
def check_isog_maps(maps,curves,p):
    """x-map sends x of points on curve i to x of points on curve i+1."""
    for f,(c0,c1) in zip(maps,zip(curves,curves[1:])):
        X0={x for x in range(p) if pow((x**3+c0[0]*x+c0[1])%p,(p-1)//2,p)==1 or (x**3+c0[0]*x+c0[1])%p==0}
        X1={x for x in range(p) if pow((x**3+c1[0]*x+c1[1])%p,(p-1)//2,p)==1 or (x**3+c1[0]*x+c1[1])%p==0}
        for x in X0:
            y=rapply(f,x,p)
            if y is not None and y not in X1: return False
    return True

def build_layers(family,k,d,p,rng):
    if family=='kummer': return [kummer(d,p)]*k
    if family=='generic': return [generic(d,p,rng) for _ in range(k)]
    if family=='lattes':
        assert d==4; a,b=rand_curve(p,rng); L=lattes(a,b,p)
        assert rdeg(L,p)==4; return [L]*k
    if family=='isogeny':
        ell,per=(2,1) if d==2 else (3,1) if d==3 else (2,2)
        maps,curves=isogeny_chain(ell,per*k,p,rng); assert check_isog_maps(maps,curves,p)
        layers=[]
        for j in range(k):
            F=maps[per*j]
            for m in maps[per*j+1:per*(j+1)]: F=compose(m,F,p)
            assert rdeg(F,p)==d; layers.append(F)
        return layers
    raise ValueError(family)

def applyL(layers,x,p):
    for f in layers:
        if x is None: return None
        x=rapply(f,x,p)
    return x

# ---------- instance ----------
def fibre_info(layers,p):
    img={}
    for x in range(p):
        v=applyL(layers,x,p)
        if v is not None: img.setdefault(v,[]).append(x)
    return img

def decomp_count(F,ptsx,a,b,p,R_list):
    """avg # ordered (x1,x2) in F x F with S3(x1,x2,xR)=0 (algebraic), restricted to rational-point x's"""
    Fs=[x for x in F if x in ptsx]; tot=0
    for xR in R_list:
        for i,x1 in enumerate(Fs):
            for x2 in Fs:
                if S3_closed(x1,x2,xR,a,b)%p==0: tot+=1
    return tot/len(R_list)

def make_instance(p,family,k,d,s,seed,a=None,b=None):
    rng=random.Random(f"{p}-{family}-{k}-{d}-{s}-{seed}")
    # target curve (shared per p,seed)
    trng=random.Random(f"E-{p}-{seed}"); a,b=rand_curve(p,trng)
    pts=[P for P in points(a,b,p) if P is not None]; ptsx=sorted({P[0] for P in pts})
    layers=build_layers(family,k,d,p,rng)
    img=fibre_info(layers,p); keys=sorted(img)
    s_=min(s,len(keys)); C=sorted(rng.sample(keys,s_))
    F=sorted(x for c in C for x in img[c])
    # planted decomposition: P1,P2 with x in F, P1 != +-P2
    Fp=[P for P in pts if P[0] in set(F)]
    plant=None
    for _ in range(1000):
        if len(Fp)<2: break
        P1,P2=rng.choice(Fp),rng.choice(Fp)
        if P1[0]!=P2[0]:
            R=ec_add(P1,P2,a,p)
            if R is not None: plant=(P1,P2,R); break
    return dict(p=p,a=a,b=b,family=family,k=k,d=d,s=s_,seed=seed,layers=layers,C=C,F=F,ptsx=ptsx,plant=plant,rng=rng)

def poly_str(c,var,p):
    terms=[]
    for e,co in enumerate(c):
        co%=p
        if co==0: continue
        terms.append(f"{co}*{var}^{e}" if e>1 else (f"{co}*{var}" if e==1 else f"{co}"))
    return "+".join(terms) if terms else "0"

def build_system(inst):
    p,a,b,k=inst['p'],inst['a'],inst['b'],inst['k']
    xR=inst['plant'][2][0]
    vars_=[f"y{i}_{j}" for j in range(k+1) for i in (1,2)]
    eqs=[]
    for i in (1,2):
        for j,(N,D) in enumerate(inst['layers']):
            yj=f"y{i}_{j}"; yn=f"y{i}_{j+1}"
            # yn*D(yj)-N(yj)
            Dn=[(-c)%p for c in []]
            eqs.append(f"{yn}*({poly_str(D,yj,p)})-({poly_str(N,yj,p)})")
        top=[1]
        for c in inst['C']: top=pmul(top,[(-c)%p,1],p)
        eqs.append(poly_str(top,f"y{i}_{k}",p))
    X,Y=sp.symbols('y1_0 y2_0')
    S=sp.Poly(sp.expand(S3_closed(X,Y,xR,a,b)),X,Y,modulus=p)
    eqs.append(str(S.as_expr()).replace('**','^'))
    degs=[]
    for e in eqs: pass
    return vars_,eqs,xR

def planted_values(inst):
    P1,P2,R=inst['plant']; p=inst['p']; k=inst['k']; vals={}
    for i,P in ((1,P1),(2,P2)):
        x=P[0]
        for j in range(k+1):
            vals[f"y{i}_{j}"]=x
            if j<k: x=rapply(inst['layers'][j],x,p)
            if x is None: return None
    return vals

def check_roots(inst,vars_,eqs):
    vals=planted_values(inst); p=inst['p']
    if vals is None: return False
    for e in eqs:
        v=eval(e.replace('^','**'),{},dict(vals))%p
        if v: return False
    return True

def bezout(inst):
    k,d,s=inst['k'],inst['d'],inst['s']
    return ((d+1)**k*s)**2*4

def sys_hash(vars_,eqs,p):
    return hashlib.sha256((",".join(vars_)+str(p)+"|".join(eqs)).encode()).hexdigest()[:16]

# ---------- engines ----------
LIM=['nice','-n','10','timeout','300','prlimit','--as=2500000000']
def run_msolve(vars_,eqs,p):
    inp=",".join(vars_)+"\n%d\n"%p+",\n".join(eqs)+"\n"
    with tempfile.NamedTemporaryFile('w',suffix='.ms',delete=False) as fh: fh.write(inp); path=fh.name
    try:
        r=subprocess.run(LIM+['msolve','-t','1','-v','2','-f',path,'-g','2','-o','/dev/null'],capture_output=True,text=True)
        out=r.stdout+r.stderr
        degs=[int(m.group(1)) for line in out.splitlines() for m in [re.match(r'\s*(\d+)\s+\d+\s+\d+\s+\d+ x \d+',line)] if m]
        if r.returncode==124 or (not degs and r.returncode!=0):
            return dict(status='timeout_or_fail',rc=r.returncode,deg=max(degs) if degs else None)
        return dict(status='ok',deg=max(degs) if degs else None,rc=r.returncode)
    finally: os.unlink(path)
def run_singular(vars_,eqs,p,plantvals):
    v=",".join(vars_)
    lin=",".join(f"{n}-{x}" for n,x in plantvals.items())
    script=(f"ring r={p},({v}),dp;\noption(prot);\nideal I={','.join(eqs)};\nideal G=std(I);\n"
            f"print(\"DIM \"+string(dim(G)));\nprint(\"VDIM \"+string(vdim(G)));\n"
            f"ideal J=I,{lin};\nideal H=std(J);\nprint(\"PLANT \"+string(size(H)>1 or H[1]!=1));\n$\n")
    with tempfile.NamedTemporaryFile('w',suffix='.sing',delete=False) as fh: fh.write(script); path=fh.name
    try:
        r=subprocess.run(LIM+['Singular','-q',path],capture_output=True,text=True,stdin=subprocess.DEVNULL)
        out=r.stdout+r.stderr
        if 'VDIM' not in out: return dict(status='timeout_or_fail',rc=r.returncode,deg=None)
        # protocol is printed for the first std and the second; take max over the first std only
        part=out.split('DIM ')[0]
        part=re.sub(r'\[\d+:\d+\]','',part)
        degs=[int(x) for x in re.findall(r'(?<![\d(:])(\d+)',part)]
        vd=int(re.search(r'VDIM (-?\d+)',out).group(1)); dm=int(re.search(r'\bDIM (-?\d+)',out).group(1))
        pl=re.search(r'PLANT (\d)',out)
        return dict(status='ok',deg=max(degs) if degs else None,vdim=vd,dim=dm,plant_not_eliminated=(pl.group(1)=='1') if pl else None)
    finally: os.unlink(path)
