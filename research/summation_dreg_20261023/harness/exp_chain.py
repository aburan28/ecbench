"""Semaev S3-chain PDP systems over GF(2^n), Weil-restricted, with levers:
 - field_eqs: include c^2=c as explicit generators (Singular/msolve always quotient)
 - aux_in_V: constrain auxiliary chain vars to subspace V (vs free in GF(2^n))
Measures: independent Macaulay first-fall + rank profile; Singular solving degree; msolve correctness.
All deterministic outputs are hashed for reproducibility. Independent of crypto repo."""
import sys, random, json, hashlib, time
H='/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness'
sys.path.insert(0,H)
import boolean_weil as bw, pdp
from gf2 import bin_points, bin_add

def const_elt(val,n):
    e=bw.f_zero(n)
    for k in range(n):
        if (val>>k)&1: e[k]=bw.const1()
    return e

def free_elt(prefix,i,n):
    e=bw.f_zero(n); vids=[]
    for k in range(n):
        vid=(prefix,i,k); vids.append(vid); e[k]={frozenset([vid])}
    return e,vids

def build_chain(a2,a6,n,l,m,seed, aux_in_V=False):
    rng=random.Random(seed)
    B=pdp.random_subspace_basis(n,l,rng)
    dec=pdp.find_decomposition(a2,a6,n,B,m,rng)
    if dec is None: return None
    xR=dec['xR']; xs=dec['xs']
    # factor-base point x-vars x1..xm in V
    X=[]; allv=[]
    for i in range(m):
        xi,vi=bw.var_elt('x',i+1,B,n); X.append(xi); allv+=vi
    # auxiliary chain vars A_1..A_{m-2}
    # chain: S3(x1,x2,A1); S3(A1,x3,A2); ...; S3(A_{m-2}, x_m, xR)
    eqs=[]
    aux=[]
    if m==2:
        S=bw.S3_bin_sym(X[0],X[1],const_elt(xR,n),a6,n)
        eqs+=bw.weil_equations(S,n)
    else:
        for j in range(m-2):
            if aux_in_V:
                Aj,av=bw.var_elt('a',j+1,B,n)
            else:
                Aj,av=free_elt('a',j+1,n)
            aux.append(Aj); allv+=av
        # first link
        S=bw.S3_bin_sym(X[0],X[1],aux[0],a6,n); eqs+=bw.weil_equations(S,n)
        for j in range(1,m-2):
            S=bw.S3_bin_sym(aux[j-1],X[j+1],aux[j],a6,n); eqs+=bw.weil_equations(S,n)
        # last link to target
        S=bw.S3_bin_sym(aux[-1],X[m-1],const_elt(xR,n),a6,n); eqs+=bw.weil_equations(S,n)
    allv=sorted(set(allv))
    return dict(B=B,dec=dec,eqs=eqs,vars=allv,xR=xR,xs=xs,m=m,n=n,l=l)

def sys_hash(eqs,vars_):
    vidx={v:i for i,v in enumerate(sorted(vars_))}
    norm=sorted(tuple(sorted(vidx[v] for v in mon)) for f in eqs for mon in f)
    return hashlib.sha256(repr(norm).encode()).hexdigest()[:16]

def run(a2,a6,n,l,m,seed,Dmax=7,aux_in_V=False,do_msolve=True):
    S=build_chain(a2,a6,n,l,m,seed,aux_in_V=aux_in_V)
    if S is None: return dict(n=n,l=l,m=m,seed=seed,status='no_decomp')
    t0=time.time()
    mac=pdp.macaulay_rank_seq(S['eqs'],S['vars'],Dmax)
    t1=time.time()
    sg=pdp.parse_singular_degree(pdp.run_singular(S['eqs'],S['vars']))
    t2=time.time()
    res=dict(a2=a2,a6=a6,n=n,l=l,m=m,seed=seed,aux_in_V=aux_in_V,
             nvars=len(S['vars']),neqs=len(S['eqs']),xR=S['xR'],
             sys_sha=sys_hash(S['eqs'],S['vars']),
             first_fall_degree=mac['first_fall_degree'],
             macaulay=mac['profile'], singular=sg,
             t_macaulay=round(t1-t0,2),t_singular=round(t2-t1,2))
    if do_msolve:
        mv=pdp.run_msolve_check(S['eqs'],S['vars'],{})
        res['msolve']=mv; res['t_msolve']=round(time.time()-t2,2)
    return res

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--grid',default='m2')
    ap.add_argument('--out',default=H+'/out/exp_chain.json')
    a=ap.parse_args()
    jobs=[]
    if a.grid=='m2':
        for (n,l) in [(6,3),(8,4),(10,5)]: jobs.append((n,l,2,False))
    elif a.grid=='m3free':
        for (n,l) in [(6,2),(8,2),(8,3),(10,3)]: jobs.append((n,l,3,False))
    elif a.grid=='m3auxV':
        for (n,l) in [(6,2),(8,2),(8,3),(10,3)]: jobs.append((n,l,3,True))
    out=[]
    for (n,l,m,auxV) in jobs:
        r=run(1,1,n,l,m,seed=1,aux_in_V=auxV)
        pr=r.get('singular'); ff=r.get('first_fall_degree')
        print(json.dumps({k:r[k] for k in ('n','l','m','aux_in_V','nvars','neqs','first_fall_degree','singular') if k in r}))
        out.append(r)
    json.dump(out,open(a.out,'w'),indent=1)
    print("WROTE",a.out)

def _span(B,n):
    V={0}
    for b in B: V|={v^b for v in V}
    return V
