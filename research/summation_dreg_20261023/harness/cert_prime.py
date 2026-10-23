import sys, itertools, json
sys.path.insert(0,'/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness')
from semaev import ec_add, points

def S3(x1,x2,x3,a,b,p):
    return ((x1-x2)**2*x3**2
            - 2*((x1+x2)*(x1*x2+a)+2*b)*x3
            + (x1*x2-a)**2 - 4*b*(x1+x2))%p

def cert(a,b,p):
    pts=[P for P in points(a,b,p) if P is not None]
    tp=fp=fn=tn=0
    for P,Q,R in itertools.product(pts,repeat=3):
        summ = ec_add(ec_add(P,Q,a,p),R,a,p)
        is_zero_sum = (summ is None)
        s = (S3(P[0],Q[0],R[0],a,b,p)==0)
        # S3 vanishing <=> exists sign choice with P±Q±R=O; group by x only.
        # Correct statement: S3(xP,xQ,xR)=0 iff (P+Q+R=O or P+Q-R=O or P-Q+R=O or P-Q-R=O)
        any_sum = any(ec_add(ec_add(sP,sQ,a,p),sR,a,p) is None
                      for sP in (P,(P[0],(-P[1])%p))
                      for sQ in (Q,(Q[0],(-Q[1])%p))
                      for sR in (R,(R[0],(-R[1])%p)))
        if s and any_sum: tp+=1
        elif s and not any_sum: fp+=1
        elif not s and any_sum: fn+=1
        else: tn+=1
    return dict(p=p,a=a,b=b,n_aff=len(pts),tp=tp,fp=fp,fn=fn,tn=tn,ok=(fp==0 and fn==0))

res=[]
for (a,b,p) in [(2,3,11),(1,6,13),(0,7,17),(3,5,19),(1,1,23),(2,1,29)]:
    # skip singular curves (4a^3+27b^2=0)
    if (4*a**3+27*b**2)%p==0: continue
    res.append(cert(a,b,p))
for r in res: print(r)
print("ALL_OK", all(r['ok'] for r in res))
json.dump(res,open('/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness/out/cert_prime.json','w'),indent=0)
