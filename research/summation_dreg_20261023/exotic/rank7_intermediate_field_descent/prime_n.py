import sys, json, random
sys.path.insert(0,'.')
import r7, pdp
from gf2 import gmul, gsq
def subfields(n):
    divs=[d for d in range(1,n+1) if n%d==0]
    out={}
    for d in divs:
        cnt=sum(1 for x in range(1<<n) if r7.gpow(x,1<<d,n)==x) if n<=13 else None
        out[d]=cnt
    return divs,out
def prime_run(n,l,m,seed):
    rng=random.Random(seed); B=pdp.random_subspace_basis(n,l,rng)
    V={v:0 for v in pdp.subspace_elements(B,n)}
    S=dict(n=n,k=l,m=m,a2=1,a6=1,V=V,basis=B)
    # random non-subfield subspace: build V as a dict(x->bits)
    elts=sorted(V); idx={}
    for bits in range(1<<l):
        v=0
        for j in range(l):
            if (bits>>j)&1: v^=B[j]
        idx[v]=bits
    S['V']=idx
    fb=r7.fb_points(S); R=None
    for _ in range(5000):
        ch=[rng.choice(fb) for _ in range(m)]
        R=None
        for P in ch: R=r7.bin_add(R,P,1,1,n)
        if R is not None: break
    S['xR']=R[0]; S['xs']=[P[0] for P in ch]
    eqs,vs=r7.build_B(S)
    asg={}
    for i,x in enumerate(S['xs']):
        for j in range(l): asg[('x',i+1,j)]=(idx[x]>>j)&1
    return dict(n=n,l=l,m=m,nvars=len(vs),planted_ok=all(r7.eval_b(f,asg)==0 for f in eqs),
                msolve=r7.msolve_deg(eqs,vs),singular=r7.singular_B(eqs,vs))
if __name__=='__main__':
    res={}
    for n in (11,13):
        divs,cnt=subfields(n); res[f'n{n}']=dict(divisors=divs,fixed_pts_by_d=cnt)
        res[f'n{n}']['prime_run']=prime_run(n,3,2,1)
    for n in (17,131,233,283):
        res[f'n{n}']=dict(divisors=[d for d in range(1,n+1) if n%d==0])
    json.dump(res,open(sys.argv[1],'w'),indent=1); print(json.dumps(res,indent=1))
