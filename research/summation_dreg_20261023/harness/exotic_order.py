"""EXOTIC LEVER (Rank 8, survey_E): monomial/variable-order independence of the
solving degree. Theory (Huang-Kosters-Yeo last fall degree) says the floor is
order-independent. We run grevlex/grlex/lex + random variable permutations on a
fixed descended PDP system and check the max degree reached is invariant.
Falsifier of 'ordering is a lever': all orders reach the SAME max degree."""
import sys,random,json
H='/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness'
sys.path.insert(0,H)
import boolean_weil as bw, pdp, exp_chain as EC

def build(n,ell,m,model,seed=1):
    S=EC.build_chain(1,1,n,ell,m,seed=seed,aux_in_V=False) if model=='split' else None
    if model=='split': return S['eqs'],S['vars']
    rng=random.Random(seed); B=pdp.random_subspace_basis(n,ell,rng)
    dec=pdp.find_decomposition(1,1,n,B,m,rng); xR=dec['xR']
    X1,v1=bw.var_elt('x',1,B,n);X2,v2=bw.var_elt('x',2,B,n);X3,v3=bw.var_elt('x',3,B,n)
    S4=bw.S4_unsplit(X1,X2,X3,EC.const_elt(xR,n),1,n)
    return bw.weil_equations(S4,n),sorted(set(v1)|set(v2)|set(v3))

if __name__=='__main__':
    out=[]
    for (n,ell,m,model) in [(8,2,3,'unsplit'),(10,2,3,'unsplit'),(8,2,3,'split')]:
        eqs,vars_=build(n,ell,m,model)
        res={'n':n,'ell':ell,'m':m,'model':model,'nvars':len(vars_)}
        for order in ['dp','Dp','lp']:
            res[f'order_{order}']=pdp.run_singular_order(eqs,vars_,order=order)
        rng=random.Random(7); perms=[]
        for _ in range(4):
            p=list(range(len(vars_))); rng.shuffle(p); perms.append(tuple(p))
        res['grevlex_randperm']=[pdp.run_singular_order(eqs,vars_,order='dp',varperm=list(p)) for p in perms]
        degs=[res['order_dp'],res['order_Dp'],res['order_lp']]+res['grevlex_randperm']
        res['all_equal']=(len(set(d for d in degs if d is not None))==1)
        res['max_over_orders']=max(d for d in degs if d is not None)
        res['min_over_orders']=min(d for d in degs if d is not None)
        print(json.dumps(res),flush=True)
        out.append(res)
        json.dump(out,open(H+'/out/exotic_order.json','w'),indent=1)
    print("DONE")
