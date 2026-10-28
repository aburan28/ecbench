"""Main experiment. usage: python3 run.py <outdir> """
import sys, os, json, re, random, hashlib, subprocess, tempfile, itertools
OUT=sys.argv[1]; os.makedirs(OUT,exist_ok=True)
Qc={tuple(map(int,k.split(','))):v for k,v in json.load(open('raw/Q_u_poly.json')).items()}
import sympy as sp
Xs,bs=sp.symbols('X b')
Qpolys={k:sp.Poly(sp.sympify(v),Xs,bs) for k,v in Qc.items()}

def S3n(a1,a2,a3,b,p):
    return ((a1-a2)**2*a3**2-2*((a1+a2)*a1*a2+2*b)*a3+(a1*a2)**2-4*b*(a1+a2))%p
def ec_add(P,Q,b,p):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P;x2,y2=Q
    if x1==x2 and (y1+y2)%p==0: return None
    lam=(3*x1*x1)*pow(2*y1,-1,p)%p if P==Q else (y2-y1)*pow(x2-x1,-1,p)%p
    x3=(lam*lam-x1-x2)%p; y3=(lam*(x1-x3)-y1)%p
    return (x3,y3)
def sqrt_mod(n,p):
    n%=p
    if n==0: return 0
    if pow(n,(p-1)//2,p)!=1: return None
    return sp.sqrt_mod(n,p)
def omega(p):
    for g in range(2,p):
        w=pow(g,(p-1)//3,p)
        if w!=1: return min(w,pow(w,2,p))
def Qnum_coeffs(Xv,b,p):
    return {k:int(P.eval({Xs:Xv,bs:b}))%p if False else int(sum(c*pow(Xv,m[0],p)*pow(b,m[1],p) for m,c in P.terms()))%p for k,P in Qpolys.items()}
def Qval(coef,u1,u2,p):
    return sum(c*pow(u1,i,p)*pow(u2,j,p) for (i,j),c in coef.items())%p

def poly_str_S3(p,b,Xv):
    # S3(x1,x2,Xv) as Singular/msolve string
    return f"(x1-x2)^2*{Xv*Xv%p}-2*((x1+x2)*x1*x2+{2*b%p})*{Xv}+(x1*x2)^2-{4*b%p}*(x1+x2)"
def poly_str_Q(coef,p):
    return "+".join(f"{c}*u1^{i}*u2^{j}" for (i,j),c in sorted(coef.items()) if c)

def run_singular(vars_,p,eqs,extra=""):
    code=f'ring r={p},({",".join(vars_)}),dp;\nideal I={",".join(eqs)};\noption(prot);\nideal G=std(I);\noption(noprot);\nprint("VDIM="+string(vdim(G)));\n{extra}quit;\n'
    with tempfile.NamedTemporaryFile('w',suffix='.sing',delete=False) as f: f.write(code); path=f.name
    try:
        r=subprocess.run(f'nice -n 10 timeout 300 prlimit --as=2500000000 Singular -q {path}',shell=True,capture_output=True,text=True,timeout=320)
        out=r.stdout
    finally: os.unlink(path)
    if 'VDIM=' not in out: return dict(status='timeout_or_fail',deg=None,vdim=None,rad=None)
    brk=[int(x) for x in re.findall(r'\[\d+:(\d+)\]',out)]
    plines=[l for l in out.splitlines() if l.startswith('[')]
    prot=re.sub(r'\(\d+\)',' ',re.sub(r'\[\d+:\d+\]',' '," ".join(plines)))
    lv=[]; prev=0
    for t in re.findall(r'\d+',prot):
        vals=[int(t)]
        if len(t)>=4 and int(t)>prev+8:   # consecutive pair-levels printed with no separator, e.g. 147148
            for i in range(2,len(t)-1):
                a,b2=int(t[:i]),int(t[i:])
                if a<=prev+8 and a<=b2<=a+8: vals=[a,b2]; break
        lv+=vals; prev=max(prev,max(vals))
    degs=lv  # degree = max pair-level number in the protocol (bracket value is uninformative in this build; recorded separately)
    vd=int(re.search(r'VDIM=(\d+)',out).group(1))
    rad=re.search(r'RADVDIM=(\d+)',out)
    return dict(status='ok',deg=max(degs) if degs else 0,bracket_deg=max(brk) if brk else None,vdim=vd,rad=int(rad.group(1)) if rad else None)

def run_msolve(vars_,p,eqs):
    txt=",".join(vars_)+"\n"+str(p)+"\n"+",\n".join(eqs)+"\n"
    with tempfile.NamedTemporaryFile('w',suffix='.ms',delete=False) as f: f.write(txt); path=f.name
    try:
        r=subprocess.run(f'nice -n 10 timeout 300 prlimit --as=2500000000 msolve -t 1 -v 2 -g 2 -f {path} -o /dev/null',shell=True,capture_output=True,text=True,timeout=320)
        out=r.stdout+r.stderr
    finally: os.unlink(path)
    degs=[]
    for line in out.splitlines():
        m=re.match(r'\s*(\d+)\s+\S+\s+\d+\s+\d+\s+x\s+\d+',line)
        if m: degs.append(int(m.group(1)))
    if r.returncode!=0 or not degs: return dict(status='timeout_or_fail' if r.returncode in (124,137) else 'nodeg',deg=max(degs) if degs else None,rc=r.returncode)
    return dict(status='ok',deg=max(degs))

def divisors3(n): return [d for d in range(3,n+1,3) if n%d==0]

def instance(p,d,seed):
    rnd=random.Random(f"{p}-{d}-{seed}")
    w=omega(p)
    while True:
        b=rnd.randrange(1,p)
        # V coset: pick x0 with rational y
        x0=rnd.randrange(1,p); y0=sqrt_mod(x0**3+b,p)
        if y0 in (None,0): continue
        c=pow(x0,d,p)
        V=[x for x in range(1,p) if pow(x,d,p)==c]
        assert len(V)==d
        cand=[x for x in V if sqrt_mod(x**3+b,p) not in (None,0)]
        x1=x0; x2=rnd.choice(cand)
        P1=(x1,y0); P2=(x2,sqrt_mod(x2**3+b,p))
        if rnd.random()<.5: P2=(P2[0],(-P2[1])%p)
        R=ec_add(P1,P2,b,p)
        if R is None or R[1]==0: continue
        if any(S3n(a,bb,R[0],b,p) for a,bb in [(x1,x2)]): continue  # sanity
        return dict(p=p,d=d,b=b,c=c,w=w,V=V,planted=(x1,x2),xR=R[0])

def measure(inst,tag):
    p,d,b,c,V,xR=inst['p'],inst['d'],inst['b'],inst['c'],inst['V'],inst['xR']
    w=inst['w']
    # brute force
    plain=[(a,bb) for a in V for bb in V if S3n(a,bb,xR,b,p)==0]
    coef=Qnum_coeffs(xR,b,p)
    U=sorted({pow(x,3,p) for x in V}); assert len(U)==d//3
    quot=[(a,bb) for a in U for bb in U if Qval(coef,a,bb,p)==0]
    # expansion check
    lifts=set()
    for (a,bb) in quot:
        for x in V:
            if pow(x,3,p)!=a: continue
            for y in V:
                if pow(y,3,p)==bb and S3n(x,y,xR,b,p)==0: lifts.add((x,y))
    exp_ok=(lifts==set(plain))
    planted_ok=(inst['planted'] in set(plain)) and ((pow(inst['planted'][0],3,p),pow(inst['planted'][1],3,p)) in set(quot))
    # also: Q zero iff some of the 9 rotated S3 zero (all U^2 pairs, checked via V lifts)
    nine_ok=all((Qval(coef,pow(x,3,p),pow(y,3,p),p)==0)==any(S3n(pow(w,k1,p)*x%p,pow(w,k2,p)*y%p,xR,b,p)==0 for k1 in range(3) for k2 in range(3)) for x in V[:min(d,12)] for y in V[:min(d,12)])
    s3=poly_str_S3(p,b,xR)
    P_eqs=[f"x1^{d}-{c}",f"x2^{d}-{c}",s3]
    Q_eqs=[f"u1^{d//3}-{c}",f"u2^{d//3}-{c}",poly_str_Q(coef,p)]
    degP=[d,d,4]; degQ=[d//3,d//3,max(i+j for (i,j),cc in coef.items() if cc)]
    def bez(ds): return dict(B_all=ds[0]*ds[1]*ds[2],B_ci=min(ds[0]*ds[1],ds[0]*ds[2],ds[1]*ds[2]))
    rad="LIB \"primdec.lib\"; ideal Rd=radical(I); ideal GR=std(Rd); print(\"RADVDIM=\"+string(vdim(GR)));\n"
    res=dict(p=p,d=d,b=b,c=c,xR=xR,N_plain=len(plain),N_quot=len(quot),expansion_ok=exp_ok,planted_ok=planted_ok,nine_factor_ok=nine_ok)
    for name,vars_,eqs,ds in (('plain',['x1','x2'],P_eqs,degP),('quot',['u1','u2'],Q_eqs,degQ)):
        h=hashlib.sha256(("|".join(eqs)+f"|{p}").encode()).hexdigest()[:16]
        sg=run_singular(vars_,p,eqs,rad); ms=run_msolve(vars_,p,eqs)
        res[name]=dict(hash=h,nvars=2,maxdeg=max(ds),eqdegs=ds,**bez(ds),singular=sg,msolve=ms)
    return res

if __name__=='__main__':
    grid={61:[6,12,30,60],97:[6,12,24,48,96],193:[6,12,24,48,96],1009:[9,18,36,63,126],4093:[6,12,33,66],65521:[6,12,18,36,45,90]}
    results=[]
    for p,ds in grid.items():
        for d in ds:
            assert (p-1)%d==0 and d%3==0
            for seed in (0,1):
                inst=instance(p,d,seed)
                r=measure(inst,'')
                r['seed']=seed
                results.append(r)
                print(p,d,seed,r['N_plain'],r['N_quot'],r['expansion_ok'],r['plain']['singular'].get('deg'),r['plain']['msolve'].get('deg'),r['quot']['singular'].get('deg'),r['quot']['msolve'].get('deg'),flush=True)
    json.dump(results,open(os.path.join(OUT,'results.json'),'w'),indent=1,sort_keys=True)
