"""Same ideal, two engines: does msolve-F4 (fall reuse) reach a lower max degree
than Singular std on MY descended PDP systems? Tests survey_A's claim that the
measured 'degree' is engine-dependent (fall reuse), independent of the repo."""
import sys,random,json,subprocess,tempfile,re,os
import os; H=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,H)
import boolean_weil as bw, pdp, exp_chain as EC

def to_msolve(eqs,vars_):
    vs=sorted(vars_); vmap={v:f"c{i}" for i,v in enumerate(vs)}
    gens=[pdp.to_singular_poly(f,vmap) for f in eqs]+[f"{vmap[v]}^2+{vmap[v]}" for v in vs]
    return ",".join(vmap[v] for v in vs)+"\n2\n"+",\n".join(gens)+"\n"

def msolve_maxdeg(eqs,vars_):
    inp=to_msolve(eqs,vars_)
    with tempfile.NamedTemporaryFile('w',suffix='.ms',delete=False) as fh:
        fh.write(inp); path=fh.name
    try:
        p=subprocess.run(['msolve','-t','1','-v','2','-f',path,'-g','2','-o','/dev/null'],
                         capture_output=True,text=True,timeout=300)
        out=p.stdout+p.stderr
        degs=[]
        for line in out.splitlines():
            m=re.match(r'\s*(\d+)\s+\d+\s+\d+\s+\d+ x \d+',line)
            if m: degs.append(int(m.group(1)))
        return max(degs) if degs else None
    finally:
        os.unlink(path)

def build(n,ell,m,model,seed=1):
    rng=random.Random(seed); B=pdp.random_subspace_basis(n,ell,rng)
    dec=pdp.find_decomposition(1,1,n,B,m,rng); xR=dec['xR']
    if model=='split':
        S=EC.build_chain(1,1,n,ell,m,seed=seed,aux_in_V=False); return S['eqs'],S['vars']
    X1,v1=bw.var_elt('x',1,B,n);X2,v2=bw.var_elt('x',2,B,n);X3,v3=bw.var_elt('x',3,B,n)
    S4=bw.S4_unsplit(X1,X2,X3,EC.const_elt(xR,n),1,n)
    return bw.weil_equations(S4,n),sorted(set(v1)|set(v2)|set(v3))

if __name__=='__main__':
    out=[]
    for (n,ell,m,model) in [(8,2,3,'split'),(10,2,3,'split'),(8,3,3,'split'),
                            (8,2,3,'unsplit'),(10,2,3,'unsplit'),(9,3,3,'unsplit')]:
        eqs,vars_=build(n,ell,m,model)
        sing=pdp.parse_singular_degree(pdp.run_singular(eqs,vars_))['solving_degree']
        mso=msolve_maxdeg(eqs,vars_)
        r=dict(n=n,ell=ell,m=m,model=model,nvars=len(vars_),singular_std=sing,msolve_f4=mso)
        print(json.dumps(r),flush=True); out.append(r)
    json.dump(out,open(H+'/out/engine_compare.json','w'),indent=1)
    print("DONE")
