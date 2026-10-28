"""Canonical reproducibility runner. Emits ONLY deterministic invariants
(system hash, first-fall degree, Singular solving degree, vdim, dim, Macaulay
rank profile) for a fixed grid. Run twice; the two JSON files must be identical
(timings excluded by construction). This is the reproducibility certificate."""
import sys,json,hashlib
import os; H=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,H)
import boolean_weil as bw, pdp, exp_chain as EC
from gf2 import bin_points

def sysh(eqs,vars_):
    vidx={v:i for i,v in enumerate(sorted(vars_))}
    norm=sorted(tuple(sorted(vidx[v] for v in mon)) for f in eqs for mon in f)
    return hashlib.sha256(repr(norm).encode()).hexdigest()[:16]

def inv_measure(eqs,vars_,Dmax):
    mac=pdp.macaulay_rank_seq(eqs,vars_,Dmax)
    sg=pdp.parse_singular_degree(pdp.run_singular(eqs,vars_))
    prof=[(r['D'],r['rows'],r['cols'],r['rank'],r['falls']) for r in mac['profile']]
    return dict(nvars=len(vars_),neqs=len(eqs),sys_sha=sysh(eqs,vars_),
                first_fall=mac['first_fall_degree'],
                solving_degree=sg['solving_degree'],vdim=sg['vdim'],dim=sg['dim'],
                macaulay_profile=prof)

def chain_case(n,l,m,auxV,Dmax):
    S=EC.build_chain(1,1,n,l,m,seed=1,aux_in_V=auxV)
    if S is None: return dict(kind='chain',n=n,l=l,m=m,auxV=auxV,status='no_decomp')
    r=inv_measure(S['eqs'],S['vars'],Dmax); r.update(kind='chain',n=n,l=l,m=m,auxV=auxV,xR=S['xR'])
    return r

def svu_case(n,l,Dmax_u,Dmax_s):
    import random
    rng=random.Random(1); B=pdp.random_subspace_basis(n,l,rng)
    dec=pdp.find_decomposition(1,1,n,B,3,rng)
    if dec is None: return [dict(kind='unsplit',n=n,l=l,status='no_decomp'),dict(kind='split',n=n,l=l,status='no_decomp')]
    xR=dec['xR']
    X1,v1=bw.var_elt('x',1,B,n); X2,v2=bw.var_elt('x',2,B,n); X3,v3=bw.var_elt('x',3,B,n); XR=EC.const_elt(xR,n)
    S4=bw.S4_unsplit(X1,X2,X3,XR,1,n); eqU=bw.weil_equations(S4,n); vU=sorted(set(v1)|set(v2)|set(v3))
    U=inv_measure(eqU,vU,Dmax_u); U.update(kind='unsplit',n=n,l=l,xR=xR)
    A,va=EC.free_elt('a',1,n)
    eqS=bw.weil_equations(bw.S3_bin_sym(X1,X2,A,1,n),n)+bw.weil_equations(bw.S3_bin_sym(A,X3,XR,1,n),n)
    vS=sorted(set(v1)|set(v2)|set(v3)|set(va)); Sp=inv_measure(eqS,vS,Dmax_s); Sp.update(kind='split',n=n,l=l,xR=xR)
    return [U,Sp]

GRID_CHAIN=[(6,3,2,False,6),(8,4,2,False,6),(10,5,2,False,6),(12,6,2,False,5),
            (6,2,3,False,6),(8,2,3,False,5),(10,2,3,False,5),(8,3,3,False,5),(9,3,3,False,5),
            (6,2,3,True,6),(8,2,3,True,5),(8,3,3,True,5)]
GRID_SVU=[(6,2,6,5),(8,2,6,5),(10,2,6,5),(9,3,6,5)]

if __name__=='__main__':
    tag=sys.argv[1] if len(sys.argv)>1 else 'run1'
    rows=[]
    for (n,l,m,auxV,Dm) in GRID_CHAIN: rows.append(chain_case(n,l,m,auxV,Dm))
    for (n,l,Du,Ds) in GRID_SVU: rows+=svu_case(n,l,Du,Ds)
    path=f"{H}/out/canonical_{tag}.json"
    json.dump(rows,open(path,'w'),indent=1,sort_keys=True)
    # also a flat digest of just the headline invariants
    digest=[{k:row.get(k) for k in ('kind','n','l','m','auxV','nvars','neqs','sys_sha','first_fall','solving_degree','vdim','dim')} for row in rows]
    json.dump(digest,open(f"{H}/out/digest_{tag}.json",'w'),indent=1,sort_keys=True)
    print(f"WROTE {path}  ({len(rows)} rows)")
