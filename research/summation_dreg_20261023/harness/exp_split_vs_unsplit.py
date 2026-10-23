"""Direct test of Semaev's splitting claim: decompose R into 3 factor-base points.
 UNSPLIT: single equation WeilRestrict(S4(x1,x2,x3,xR)) in 3l vars, degree up to 4.
 SPLIT  : WeilRestrict(S3(x1,x2,A)) + WeilRestrict(S3(A,x3,xR)) in 3l+n vars, degree 2.
Same curve/subspace/target/seed so the comparison is apples-to-apples.
Deterministic outputs (sys_sha, first_fall, solving_degree, vdim) are the reproducible invariants."""
import sys,random,json,hashlib,time
H='/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness'
sys.path.insert(0,H)
import boolean_weil as bw, pdp, exp_chain as EC
from gf2 import bin_points, bin_add

def sysh(eqs,vars_):
    vidx={v:i for i,v in enumerate(sorted(vars_))}
    norm=sorted(tuple(sorted(vidx[v] for v in mon)) for f in eqs for mon in f)
    return hashlib.sha256(repr(norm).encode()).hexdigest()[:16]

def measure(eqs,vars_,Dmax):
    mac=pdp.macaulay_rank_seq(eqs,vars_,Dmax)
    sg=pdp.parse_singular_degree(pdp.run_singular(eqs,vars_))
    return dict(nvars=len(vars_),neqs=len(eqs),sys_sha=sysh(eqs,vars_),
                first_fall_degree=mac['first_fall_degree'],
                macaulay=mac['profile'],singular=sg)

def run(a2,a6,n,l,seed,Dmax_unsplit=6,Dmax_split=5):
    rng=random.Random(seed)
    B=pdp.random_subspace_basis(n,l,rng)
    dec=pdp.find_decomposition(a2,a6,n,B,3,rng)
    if dec is None: return dict(status='no_decomp',n=n,l=l)
    xR=dec['xR']
    X1,v1=bw.var_elt('x',1,B,n); X2,v2=bw.var_elt('x',2,B,n); X3,v3=bw.var_elt('x',3,B,n)
    XR=EC.const_elt(xR,n)
    # UNSPLIT
    S4=bw.S4_unsplit(X1,X2,X3,XR,a6,n)
    eqU=bw.weil_equations(S4,n); vU=sorted(set(v1)|set(v2)|set(v3))
    tu=time.time(); U=measure(eqU,vU,Dmax_unsplit); tu=round(time.time()-tu,1)
    # SPLIT
    A,va=EC.free_elt('a',1,n)
    S_a=bw.S3_bin_sym(X1,X2,A,a6,n); S_b=bw.S3_bin_sym(A,X3,XR,a6,n)
    eqS=bw.weil_equations(S_a,n)+bw.weil_equations(S_b,n); vS=sorted(set(v1)|set(v2)|set(v3)|set(va))
    ts=time.time(); S=measure(eqS,vS,Dmax_split); ts=round(time.time()-ts,1)
    return dict(n=n,l=l,a2=a2,a6=a6,seed=seed,xR=xR,
                unsplit=U,split=S,t_unsplit=tu,t_split=ts)

if __name__=='__main__':
    out=[]
    for (n,l) in [(6,2),(8,2),(9,3),(10,2),(8,3)]:
        r=run(1,1,n,l,seed=1)
        if r.get('status'): print(json.dumps(r)); out.append(r); continue
        U,S=r['unsplit'],r['split']
        print(json.dumps(dict(n=n,l=l,
            unsplit=dict(nvars=U['nvars'],neqs=U['neqs'],ff=U['first_fall_degree'],sg=U['singular']),
            split=dict(nvars=S['nvars'],neqs=S['neqs'],ff=S['first_fall_degree'],sg=S['singular']))),flush=True)
        out.append(r)
        json.dump(out,open(H+'/out/split_vs_unsplit.json','w'),indent=1)
    print("DONE")
