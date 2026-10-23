"""Cross-check survey_A's central claim with independent tooling:
SAT solving degree (decomposable target) vs UNSAT refutation degree
(non-decomposable target), same (n,ell,m,V), measured by msolve/Singular.
If UNSAT refutes at a lower/flatter degree than SAT solves, that corroborates
survey_A: the 'degree' depends on the engine's job, and fall-reuse keeps
refutation low."""
import sys,random,json
import os; H=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,H)
import boolean_weil as bw, pdp, exp_chain as EC
from gf2 import bin_points, bin_add

def decomposable_xr(a2,a6,n,B,m,rng):
    d=pdp.find_decomposition(a2,a6,n,B,m,rng,tries=40000)
    return d['xR'] if d else None

def all_decomposable(a2,a6,n,B,m):
    """set of x_R that ARE decomposable into m FB points."""
    pts=[P for P in bin_points(a2,a6,n) if P is not None]
    V=pdp.subspace_elements(B,n)
    fb=[P for P in pts if P[0] in V]
    import itertools
    good=set()
    for combo in itertools.product(fb,repeat=m):
        S=None
        for P in combo: S=bin_add(S,P,a2,a6,n)
        if S is not None: good.add(S[0])
    return good, [P[0] for P in pts]

def build_split(a2,a6,n,B,m,xR):
    X=[];allv=[]
    for i in range(m):
        xi,vi=bw.var_elt('x',i+1,B,n);X.append(xi);allv+=vi
    aux=[]
    if m==2:
        eqs=bw.weil_equations(bw.S3_bin_sym(X[0],X[1],EC.const_elt(xR,n),a6,n),n)
    else:
        for j in range(m-2):
            Aj,av=EC.free_elt('a',j+1,n);aux.append(Aj);allv+=av
        eqs=bw.weil_equations(bw.S3_bin_sym(X[0],X[1],aux[0],a6,n),n)
        for j in range(1,m-2):
            eqs+=bw.weil_equations(bw.S3_bin_sym(aux[j-1],X[j+1],aux[j],a6,n),n)
        eqs+=bw.weil_equations(bw.S3_bin_sym(aux[-1],X[m-1],EC.const_elt(xR,n),a6,n),n)
    return eqs,sorted(set(allv))

def run(a2,a6,n,ell,m,nsat=6,nunsat=6,seed=1):
    rng=random.Random(seed)
    B=pdp.random_subspace_basis(n,ell,rng)
    good,allx=all_decomposable(a2,a6,n,B,m)
    bad=[x for x in allx if x not in good]
    sat=[]; unsat=[]
    for xR in list(good)[:nsat]:
        eqs,vars_=build_split(a2,a6,n,B,m,xR)
        sg=pdp.parse_singular_degree(pdp.run_singular(eqs,vars_))
        sat.append((xR,sg['solving_degree'],sg['vdim']))
    for xR in bad[:nunsat]:
        eqs,vars_=build_split(a2,a6,n,B,m,xR)
        sg=pdp.parse_singular_degree(pdp.run_singular(eqs,vars_))
        unsat.append((xR,sg['solving_degree'],sg['vdim']))
    return dict(n=n,ell=ell,m=m,nB=len(good),
                sat_degs=sorted(set(d for _,d,_ in sat)),sat=sat,
                unsat_degs=sorted(set(d for _,d,_ in unsat)),unsat=unsat)

if __name__=='__main__':
    out=[]
    for (n,ell,m) in [(8,2,3),(8,3,3),(10,2,3)]:
        r=run(1,1,n,ell,m)
        print(json.dumps({k:r[k] for k in ('n','ell','m','nB','sat_degs','unsat_degs')}),flush=True)
        out.append(r)
    json.dump(out,open(H+'/out/sat_vs_unsat.json','w'),indent=1)
    print("DONE")
