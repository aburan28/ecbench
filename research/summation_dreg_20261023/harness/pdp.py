"""Generate Weil-restricted point-decomposition (PDP) boolean systems and measure
solving degree with (1) an independent multilinear Macaulay-rank sequence over F2,
(2) Singular std protocol, (3) msolve for solution correctness.
Independent of the crypto repo."""
import sys, os, subprocess, random, itertools, json, re, hashlib
import os; H=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,H)
from gf2 import IRRED, gmul, gsq, bin_points, bin_add, S3_bin
import boolean_weil as bw

# ---------- build a decomposable target & factor-base subspace ----------
def random_subspace_basis(n, l, rng):
    # l independent vectors in GF(2^n) spanning an l-dim F2-subspace V (as ints)
    while True:
        B=[rng.randrange(1,1<<n) for _ in range(l)]
        # check independence via Gaussian elim over F2
        rows=B[:]; piv=[]; ok=True
        import copy
        tmp=rows[:]
        rank=0
        work=tmp[:]
        for col in range(n):
            pr=None
            for i in range(rank,len(work)):
                if (work[i]>>col)&1: pr=i;break
            if pr is None: continue
            work[rank],work[pr]=work[pr],work[rank]
            for i in range(len(work)):
                if i!=rank and (work[i]>>col)&1: work[i]^=work[rank]
            rank+=1
        if rank==l: return B

def subspace_elements(B,n):
    V=set([0])
    for b in B:
        V |= {v^b for v in V}
    return V

def find_decomposition(a2,a6,n,B,m,rng,tries=20000):
    """Find factor-base points P_1..P_m (x in V) and R=sum, R!=O, x_R defined."""
    pts=[P for P in bin_points(a2,a6,n) if P is not None]
    V=subspace_elements(B,n)
    fb=[P for P in pts if P[0] in V]
    if len(fb)<1: return None
    for _ in range(tries):
        chosen=[rng.choice(fb) for _ in range(m)]
        S=None
        for P in chosen: S=bin_add(S,P,a2,a6,n)
        if S is None: continue
        return dict(xs=[P[0] for P in chosen], xR=S[0], R=S)
    return None

# ---------- independent multilinear Macaulay rank over F2 ----------
def squarefree_monos(vars_, D):
    res=[frozenset()]  # degree 0
    for d in range(1,D+1):
        for c in itertools.combinations(vars_,d):
            res.append(frozenset(c))
    return res

def ml_mul(monoA, monoB):  # multilinear product: union (c^2=c)
    return monoA|monoB

def _rank_f2(vectors):
    basis={}; rank=0
    for v in vectors:
        while v:
            lb=v&(-v); pb=lb.bit_length()-1
            if pb in basis: v^=basis[pb]
            else: basis[pb]=v; rank+=1; break
    return rank, basis

def macaulay_rank_seq(polys, vars_, Dmax):
    """Return per-degree profile with first-fall detection (independent, deterministic)."""
    vars_=sorted(vars_)
    out=[]; prev_rank=0; D_ff=None
    for D in range(1,Dmax+1):
        cols=squarefree_monos(vars_,D)
        colidx={m:i for i,m in enumerate(cols)}
        lowcols=[i for m,i in colidx.items() if len(m)<=D-1]
        rows=[]
        for f in polys:
            fd=max((len(m) for m in f),default=0)
            if fd>D: continue
            mult=squarefree_monos(vars_, D-fd)
            for u in mult:
                row=set()
                for mon in f:
                    pm=ml_mul(mon,u)
                    if len(pm)<=D:
                        if pm in row: row.discard(pm)
                        else: row.add(pm)
                if row:
                    v=0
                    for ci in (colidx[m] for m in row): v|=(1<<ci)
                    rows.append(v)
        rank,_=_rank_f2(rows)
        # dim(V cap L): V=rowspace(rank), L=span(low unit cols, dim=len(lowcols))
        nlow=len(lowcols)
        Lbasis=[(1<<i) for i in lowcols]
        rank_VL,_=_rank_f2(rows+Lbasis)
        dim_cap = rank + nlow - rank_VL
        falls = dim_cap - prev_rank   # new low-degree elements beyond deg D-1 rowspace
        if D_ff is None and falls>0 and D>=2:
            D_ff=D
        out.append(dict(D=D,rows=len(rows),cols=len(cols),rank=rank,
                        low_in_rowspace=dim_cap, falls=max(falls,0)))
        prev_rank=rank
    return dict(profile=out, first_fall_degree=D_ff)

# ---------- Singular export + protocol degree ----------
def to_singular_poly(f, vmap):
    if not f: return "0"
    terms=[]
    for mon in f:
        if not mon: terms.append("1")
        else: terms.append("*".join(vmap[v] for v in sorted(mon)))
    return "+".join(terms)

def run_singular(polys, vars_, want_hilb=True):
    vars_=sorted(vars_); vmap={v:f"c{i}" for i,v in enumerate(vars_)}
    vlist=",".join(vmap[v] for v in vars_)
    feqs=[f"{vmap[v]}^2+{vmap[v]}" for v in vars_]
    gens=[to_singular_poly(f,vmap) for f in polys]+feqs
    script=f"ring r=2,({vlist}),dp;\noption(prot);\nideal I={','.join(gens)};\nideal G=std(I);\nprint(\"NVARS \"+string(nvars(r)));\nprint(\"DEG \"+string(deg(G[size(G)])));\nprint(\"DIM \"+string(dim(G)));\nprint(\"VDIM \"+string(vdim(std(I))));\n$\n"
    p=subprocess.run(['Singular','-q'],input=script,capture_output=True,text=True,timeout=600)
    return p.stdout+p.stderr

def parse_singular_degree(out):
    # protocol prints степени inside "...(d)..." and lines like "[65:2]..."; capture max bracket degree
    degs=[int(x) for x in re.findall(r'\[\d+:(\d+)\]',out)]
    maxd=max(degs) if degs else None
    vdim=re.search(r'VDIM (\-?\d+)',out)
    dim=re.search(r'DIM (\-?\d+)',out)
    return dict(solving_degree=maxd, vdim=int(vdim.group(1)) if vdim else None,
                dim=int(dim.group(1)) if dim else None)

if __name__=='__main__':
    print("pdp.py module ok")

def run_msolve_check(polys, vars_, known_assign):
    """Solve boolean system with msolve over F2 field eqs; check known_assign is a root.
    known_assign: dict varid->0/1. Returns dict(found, nsols)."""
    import tempfile
    vars_=sorted(vars_); vmap={v:f"c{i}" for i,v in enumerate(vars_)}
    vlist=",".join(vmap[v] for v in vars_)
    gens=[to_singular_poly(f,vmap) for f in polys]+[f"{vmap[v]}^2+{vmap[v]}" for v in vars_]
    inp="%s\n2\n"%vlist + ",\n".join(gens)+"\n"
    with tempfile.NamedTemporaryFile('w',suffix='.ms',delete=False) as fh:
        fh.write(inp); path=fh.name
    try:
        p=subprocess.run(['msolve','-f',path,'-g','2'],capture_output=True,text=True,timeout=600)
        out=p.stdout
        # count solution tuples
        nsol=out.count('[') 
        return dict(ok=(p.returncode==0), raw_len=len(out), msolve_out_head=out[:200])
    except Exception as e:
        return dict(ok=False,err=str(e))

def run_singular_order(polys, vars_, order='dp', varperm=None):
    """Solving degree under a chosen monomial order and optional variable permutation.
    order in {'dp'(grevlex),'Dp'(grlex),'lp'(lex)}."""
    vs=sorted(vars_)
    if varperm is not None: vs=[vs[i] for i in varperm]
    vmap={v:f"c{i}" for i,v in enumerate(sorted(vars_))}
    vlist=",".join(vmap[v] for v in vs)
    feqs=[f"{vmap[v]}^2+{vmap[v]}" for v in sorted(vars_)]
    gens=[to_singular_poly(f,vmap) for f in polys]+feqs
    script=f"ring r=2,({vlist}),{order};\noption(prot);\nideal I={','.join(gens)};\nideal G=std(I);\nprint(\"DIM \"+string(dim(G)));\n$\n"
    import subprocess,re
    p2=subprocess.run(['Singular','-q'],input=script,capture_output=True,text=True,timeout=900)
    out=p2.stdout+p2.stderr
    degs=[int(x) for x in re.findall(r'\[\d+:(\d+)\]',out)]
    return max(degs) if degs else None
