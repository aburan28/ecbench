import sys, random, json, hashlib, time
H='/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness'
sys.path.insert(0,H)
import boolean_weil as bw, pdp
from gf2 import bin_points, bin_add

def build_m2(a2,a6,n,l,seed):
    rng=random.Random(seed)
    B=pdp.random_subspace_basis(n,l,rng)
    dec=pdp.find_decomposition(a2,a6,n,B,2,rng)
    if dec is None: return None
    xR=dec['xR']
    # symbolic x1,x2 in V
    X1,v1=bw.var_elt('x',1,B,n)
    X2,v2=bw.var_elt('x',2,B,n)
    XR=bw.f_zero(n)
    for k in range(n):
        if (xR>>k)&1: XR[k]=bw.const1()
    S=bw.S3_bin_sym(X1,X2,XR,a6,n)
    eqs=bw.weil_equations(S,n)
    allvars=sorted(set(v1)|set(v2))
    return dict(B=B,dec=dec,eqs=eqs,vars=allvars,xR=xR)

def sys_hash(eqs,vars_):
    vidx={v:i for i,v in enumerate(sorted(vars_))}
    norm=sorted(tuple(sorted(vidx[v] for v in mon)) for f in eqs for mon in f)
    return hashlib.sha256(repr(norm).encode()).hexdigest()[:16]

def run(a2,a6,n,l,seed,Dmax=6):
    S=build_m2(a2,a6,n,l,seed)
    if S is None: return dict(n=n,l=l,seed=seed,status='no_decomp')
    t0=time.time()
    mac=pdp.macaulay_rank_seq(S['eqs'],S['vars'],Dmax)
    t1=time.time()
    sgout=pdp.run_singular(S['eqs'],S['vars'])
    sg=pdp.parse_singular_degree(sgout)
    t2=time.time()
    # first fall: smallest D where rank < cols-1 but a degree<D elt appears => proxy:
    # smallest D where (rows>0 and rank deficient vs trivial) -> we report full seq
    solved_D=next((r['D'] for r in mac if r['solved']),None)
    return dict(a2=a2,a6=a6,n=n,l=l,seed=seed,nvars=len(S['vars']),neqs=len(S['eqs']),
                xR=S['xR'], sys_sha=sys_hash(S['eqs'],S['vars']),
                macaulay=mac, macaulay_solved_D=solved_D,
                singular=sg, t_macaulay=round(t1-t0,2), t_singular=round(t2-t1,2))

if __name__=='__main__':
    a2,a6=1,1  # Koblitz-type
    out=[]
    for (n,l) in [(6,3),(8,4)]:
        r=run(a2,a6,n,l,seed=1)
        out.append(r)
        print(json.dumps({k:v for k,v in r.items() if k!='macaulay'}))
        for m in r.get('macaulay',[]): print("   ",m)
    json.dump(out,open(H+'/out/exp_m2_validate.json','w'),indent=1)
