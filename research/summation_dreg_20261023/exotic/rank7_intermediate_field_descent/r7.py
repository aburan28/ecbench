import sys, os, re, json, random, hashlib, subprocess, tempfile, itertools
HD='/home/user/ecbench/research/summation_dreg_20261023/harness'
sys.path.insert(0,HD)
import gf2, boolean_weil as bw, pdp, exp_chain as EC
from gf2 import gmul, gsq, bin_add, IRRED
bw.f_const=lambda c,n: EC.const_elt(c,n)   # allow general GF(2^n) constants (harness f_const only takes 0/1); local monkeypatch
LIM=['nice','-n','10','timeout','300','prlimit','--as=2500000000']

def gpow(x,e,n):
    r=1
    while e:
        if e&1: r=gmul(r,x,n)
        x=gsq(x,n); e>>=1
    return r

def setup(n,k,m,a2,a6,seed):
    q=1<<k; assert n==k*m
    # alpha: smallest root of IRRED[k] in GF(2^n)
    ir=IRRED[k]
    def ev(x):
        r=0
        for i in range(k+1):
            if (ir>>i)&1: r^=gpow(x,i,n)
        return r
    alpha=min(x for x in range(2,1<<n) if ev(x)==0)
    basis=[gpow(alpha,j,n) for j in range(k)]
    V={}
    for bits in range(q):
        v=0
        for j in range(k):
            if (bits>>j)&1: v^=basis[j]
        V[v]=bits           # harness elt -> a-poly bits (a^j <-> alpha^j)
    assert len(V)==q and all(gpow(v,q,n)==v for v in V)
    Vinv={b:v for v,b in V.items()}
    z=2
    # min poly of z over GF(q): prod (t - z^{q^i}); coefficients in V
    poly=[1]  # poly[i] coeff of t^i, harness elts
    for i in range(m):
        c=gpow(z,q**i,n)
        new=[0]*(len(poly)+1)
        for d,co in enumerate(poly):
            new[d+1]^=co; new[d]^=gmul(co,c,n)
        poly=new
    assert all(c in V for c in poly) and poly[m]==1
    g=[V[c] for c in poly]   # GF(q) a-bits
    # coordinates of harness element in basis z^j over GF(q)
    zp=[gpow(z,j,n) for j in range(m)]
    coord={}
    for tup in itertools.product(range(q),repeat=m):
        v=0
        for j,c in enumerate(tup): v^=gmul(Vinv[c],zp[j],n)
        coord[v]=tup
    assert len(coord)==1<<n
    return dict(n=n,k=k,m=m,q=q,alpha=alpha,basis=basis,V=V,Vinv=Vinv,g=g,coord=coord,a2=a2,a6=a6)

def fb_points(S):
    n,a2,a6=S['n'],S['a2'],S['a6']; pts=[]
    for x in S['V']:
        rhs=gmul(gsq(x,n),x,n)^gmul(a2,gsq(x,n),n)^a6
        for y in range(1<<n):
            if gsq(y,n)^gmul(x,y,n)==rhs: pts.append((x,y))
    return pts

def plant(S,seed):
    rng=random.Random(seed); fb=fb_points(S); n,a2,a6,m=S['n'],S['a2'],S['a6'],S['m']
    for _try in range(5000):
        ch=[rng.choice(fb) for _ in range(m)]
        R=None
        for P in ch: R=bin_add(R,P,a2,a6,n)
        if R is not None and R[0] not in S['V'] and len(set(P[0] for P in ch))==m: break   # distinct x (generic decomposition)
    else: return None
    S['fb']=fb; S['xs']=[P[0] for P in ch]; S['xR']=R[0]; S['R']=R
    return S

def bruteforce(S):
    """ordered x-tuples in V^m with some points P_i (x(P_i)=x_i) and sum = +-R."""
    n,a2,a6,m=S['n'],S['a2'],S['a6'],S['m']
    byx={}
    for P in S['fb']: byx.setdefault(P[0],[]).append(P)
    R=S['R']; negR=(R[0],R[0]^R[1])
    sols=set()
    for tup in itertools.product(list(byx),repeat=m):
        ok=False
        for ps in itertools.product(*[byx[x] for x in tup]):
            s=None
            for P in ps: s=bin_add(s,P,a2,a6,n)
            if s==R or s==negR: ok=True;break
        if ok: sols.add(tup)
    return sols

def numeric_S(S,xs):
    """harness-formula S_{m+1}(xs,xR) evaluated numerically."""
    n,a6,m=S['n'],S['a6'],S['m']; xR=S['xR']
    if m==2: return gf2.S3_bin(xs[0],xs[1],xR,a6,n)
    X=[EC.const_elt(x,n) for x in list(xs)+[xR]]
    s=bw.S4_unsplit(X[0],X[1],X[2],X[3],a6,n)
    return sum((1<<i) for i in range(n) if s[i])  # constants only; frozenset() monomial
# ---------------- A: GF(q) descent in Singular ----------------
def aelt(bits,k):
    if bits==0: return "0"
    return "("+"+".join("1" if j==0 else f"a^{j}" for j in range(k) if (bits>>j)&1)+")"
def tower_poly(S,val):   # harness elt -> poly in t over GF(q)
    k,m=S['k'],S['m']; c=S['coord'][val]
    return "("+"+".join(f"{aelt(ci,k)}*t^{j}" for j,ci in enumerate(c))+")"

def singular_script_A(S,fe,mode):
    k,m,q=S['k'],S['m'],S['q']
    xv=[f"x{i+1}" for i in range(m)]
    g="+".join(f"{aelt(S['g'][i],k)}*t^{i}" for i in range(m+1))
    s=f"ring r=(2^{k},a),({','.join(xv)},t),dp;\n"
    s+=f"poly g={g};\nideal gi=std(g);\n"
    s+=f"poly a6={tower_poly(S,S['a6'])};\npoly XR={tower_poly(S,S['xR'])};\n"
    if m==2:
        s+="poly e2=x1*x2+x1*XR+x2*XR;\npoly Sx=e2^2+x1*x2*XR+a6;\n"
    else:
        s+="poly A1=(x1+x2)^2; poly B1=x1*x2; poly C1=B1^2+a6;\npoly A2=(x3+XR)^2; poly B2=x3*XR; poly C2=B2^2+a6;\n"
        s+="poly Sx=(A1*C2+A2*C1)^2+(A1*B2+A2*B1)*(B1*C2+B2*C1);\n"
    s+="Sx=reduce(Sx,gi);\nmatrix M=coeffs(Sx,t);\nideal E;\nint i;\nfor(i=1;i<=nrows(M);i++){ if(M[i,1]!=0){E=E,M[i,1];} }\nE=simplify(E,2);\n"
    s+=f'ring r2=(2^{k},a),({",".join(xv)}),dp;\nsetring r2;\nideal E=imap(r,E);\n'
    s+='print("NEQ "+string(size(E)));\nfor(i=1;i<=size(E);i++){print("EQDEG "+string(deg(E[i])));}\n'
    s+=f'ideal F={",".join(f"x{i+1}^{q}+x{i+1}" for i in range(m))};\n'
    # planted check
    pl="".join(f"p=subst(p,x{i+1},{aelt(S['V'][S['xs'][i]],k)});" for i in range(m))
    s+='int ok=1; poly p; for(i=1;i<=size(E);i++){p=E[i];'+pl+'if(p!=0){ok=0;}} print("PLANT "+string(ok));\n'
    if fe: s+="E=E+F;\n"
    s+='option(prot);\nideal G=std(E);\noption(noprot);\nprint("DIM "+string(dim(G)));\nif(dim(G)==0){print("VDIM "+string(vdim(G)));}\nprint("GDEG "+string(deg(G[size(G)])));\n$\n'
    return s

def run_cmd(cmd,inp=None,timeout=320):
    try:
        p=subprocess.run(LIM+cmd,input=inp,capture_output=True,text=True,timeout=timeout)
        return p.returncode,p.stdout+p.stderr
    except subprocess.TimeoutExpired:
        return 124,"TIMEOUT"

def prot_labels(out):
    """max degree label in the Singular protocol (bare integers outside (..) and [..]); the harness
    convention '[k:d]' is also returned separately."""
    body=out.split('product criterion')[0]
    body=re.sub(r'\[\d+:\d+\]','',body)
    body=re.sub(r'\(\d+\)','',body)
    body=re.sub(r'\b(NEQ|EQDEG|PLANT|DIM|VDIM|GDEG)\b[^\n]*','',body)
    lab=[int(x) for x in re.findall(r'\d+',body)]
    return max(lab) if lab else None

def parse_A(rc,out):
    if rc!=0: return dict(status="timeout/fail rc=%s"%rc)
    degs=[int(x) for x in re.findall(r'\[\d+:(\d+)\]',out)]
    d=lambda key:(lambda mm:int(mm.group(1)) if mm else None)(re.search(key+r' (-?\d+)',out))
    eqd=[int(x) for x in re.findall(r'EQDEG (\d+)',out)]
    return dict(status="ok",solving_degree=max(degs) if degs else None,max_label_degree=prot_labels(out),dim=d('DIM'),vdim=d('VDIM'),
                eq_degrees=eqd,bezout=int(__import__('math').prod(eqd)) if eqd else None,
                planted_ok=d('PLANT')==1,neq=d('NEQ'))

# ---------------- B: F2 descent ----------------
def build_B(S):
    n,k,m=S['n'],S['k'],S['m']; a6=S['a6']
    Xs=[];vars_=[]
    for i in range(m):
        X,v=bw.var_elt('x',i+1,S['basis'],n); Xs.append(X); vars_+=v
    XR=EC.const_elt(S['xR'],n)
    if m==2: s=bw.S3_bin_sym(Xs[0],Xs[1],XR,a6,n)
    else: s=bw.S4_unsplit(Xs[0],Xs[1],Xs[2],XR,a6,n)
    eqs=bw.weil_equations(s,n)
    return eqs,sorted(vars_)

def planted_assign_B(S):
    asg={}
    for i,x in enumerate(S['xs']):
        bits=S['V'][x]
        for j in range(S['k']): asg[('x',i+1,j)]=(bits>>j)&1
    return asg
def eval_b(f,asg):
    r=0
    for mon in f:
        v=1
        for u in mon: v&=asg[u]
        r^=v
    return r
def msolve_deg(eqs,vars_):
    vs=sorted(vars_); vmap={v:f"c{i}" for i,v in enumerate(vs)}
    gens=[to_poly(f,vmap) for f in eqs]+[f"{vmap[v]}^2+{vmap[v]}" for v in vs]
    inp=",".join(vmap[v] for v in vs)+"\n2\n"+",\n".join(gens)+"\n"
    with tempfile.NamedTemporaryFile('w',suffix='.ms',delete=False) as fh:
        fh.write(inp); path=fh.name
    rc,out=run_cmd(['msolve','-t','1','-v','2','-f',path,'-g','2','-o','/dev/null'])
    os.unlink(path)
    if rc!=0: return dict(status="timeout/fail rc=%s"%rc)
    degs=[int(mm.group(1)) for l in out.splitlines() for mm in [re.match(r'\s*(\d+)\s+\d+\s+\d+\s+\d+ x \d+',l)] if mm]
    return dict(status="ok",solving_degree=max(degs) if degs else None)
def singular_B(eqs,vars_):
    vs=sorted(vars_); vmap={v:f"c{i}" for i,v in enumerate(vs)}
    gens=[to_poly(f,vmap) for f in eqs]+[f"{vmap[v]}^2+{vmap[v]}" for v in vs]
    sc=f"ring r=2,({','.join(vmap[v] for v in vs)}),dp;\noption(prot);\nideal I={','.join(gens)};\nideal G=std(I);\noption(noprot);\nprint(\"DIM \"+string(dim(G)));\nprint(\"VDIM \"+string(vdim(G)));\n$\n"
    rc,out=run_cmd(['Singular','-q'],sc)
    if rc!=0: return dict(status="timeout/fail rc=%s"%rc)
    degs=[int(x) for x in re.findall(r'\[\d+:(\d+)\]',out)]
    vd=re.search(r'VDIM (-?\d+)',out)
    return dict(status="ok",solving_degree=max(degs) if degs else None,max_label_degree=prot_labels(out),vdim=int(vd.group(1)) if vd else None)

def to_poly(f,vmap):
    if not f: return "0"
    terms=sorted("1" if not mon else "*".join(vmap[v] for v in sorted(mon)) for mon in f)
    return "+".join(terms)
def h(s): return hashlib.sha256(s.encode()).hexdigest()[:16]

def run_config(n,k,m,a2,a6name,seed=1):
    S=setup(n,k,m,a2,1,seed)
    a6=1 if a6name=='koblitz' else 7  # generic: z^2+z+1, z the harness generator (not in the subfield)
    S['a6']=a6
    if plant(S,seed) is None: return dict(n=n,k=k,m=m,a6=a6name,status='no generic planted decomposition (xR in V always)')
    bf=bruteforce(S)
    pl=tuple(S['xs'])
    # numeric check of the S formula against brute force
    num=set(t for t in itertools.product(list(S['V']),repeat=m) if numeric_S(S,t)==0)
    rec=dict(n=n,k=k,m=m,a2=a2,a6=a6name,a6_val=a6,xs=S['xs'],xR=S['xR'],nfb=len(S['fb']),
             brute_count=len(bf),numeric_S_count=len(num),numeric_eq_brute=(num==bf),planted_in_brute=pl in bf)
    # A
    for tag,fe in (('A0',False),('A1',True)):
        sc=singular_script_A(S,fe,tag)
        rec[tag+'_hash']=h(sc)
        rc,out=run_cmd(['Singular','-q'],sc)
        rec[tag]=parse_A(rc,out)
    # B
    eqs,vars_=build_B(S)
    asg=planted_assign_B(S)
    vmap={v:f"c{i}" for i,v in enumerate(sorted(vars_))}
    canon="|".join(sorted(to_poly(f,vmap) for f in eqs))
    rec['B_hash']=h(canon); rec['B_nvars']=len(vars_); rec['B_neq']=len(eqs)
    rec['B_eq_degrees']=sorted(max(len(mn) for mn in f) for f in eqs)
    rec['B_planted_ok']=all(eval_b(f,asg)==0 for f in eqs)
    rec['B_msolve']=msolve_deg(eqs,vars_)
    rec['B_singular']=singular_B(eqs,vars_)
    return rec

CONFIGS=[(4,2,2),(6,3,2),(8,4,2),(6,2,3),(9,3,3),(12,4,3)]
if __name__=='__main__':
    outp=sys.argv[1]; res=[]
    sel=sys.argv[2:] 
    for (n,k,m) in CONFIGS:
        for cv in (('koblitz',1),('generic',1)):
            r=run_config(n,k,m,1,cv[0]); res.append(r); print(json.dumps(r),flush=True)
    json.dump(res,open(outp,'w'),indent=1)
