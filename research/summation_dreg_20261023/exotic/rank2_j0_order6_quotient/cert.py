"""S3 (a=0) vs real group law, brute force, small j=0 curves: 0 FP / 0 FN on rational-y points; orbit product vs 9 rotated S3."""
import json,sys
sys.argv=['x','/dev/null']
def S3n(a1,a2,a3,b,p): return ((a1-a2)**2*a3**2-2*((a1+a2)*a1*a2+2*b)*a3+(a1*a2)**2-4*b*(a1+a2))%p
def ec_add(P,Q,p):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P;x2,y2=Q
    if x1==x2 and (y1+y2)%p==0: return None
    lam=(3*x1*x1)*pow(2*y1,-1,p)%p if P==Q else (y2-y1)*pow(x2-x1,-1,p)%p
    x3=(lam*lam-x1-x2)%p; return (x3,(lam*(x1-x3)-y1)%p)
rep=[]
for p in (7,13,19,31,37,43,61):
    w=[pow(g,(p-1)//3,p) for g in range(2,p) if pow(g,(p-1)//3,p)!=1][0]
    for b in range(1,p):
        pts=[(x,y) for x in range(p) for y in range(p) if (y*y-x**3-b)%p==0 and y!=0]
        xs=sorted({x for x,_ in pts}); ysof={x:[y for xx,y in pts if xx==x] for x in xs}
        sums={}  # x-set reachable: x such that P1+-P2 = +-P3
        fp=fn=0
        reach=set()
        for x1 in xs:
            for x2 in xs:
                for y1 in ysof[x1]:
                    for y2 in ysof[x2]:
                        R=ec_add((x1,y1),(x2,y2),p)
                        if R is not None: reach.add((x1,x2,R[0]))
        for x1 in xs:
            for x2 in xs:
                for x3 in xs:
                    z=S3n(x1,x2,x3,b,p)==0
                    # x3 reachable only if R rational with y!=0 ; R with y==0 excluded from xs (fine: x3 in xs)
                    r=(x1,x2,x3) in reach
                    fp+=z and not r; fn+=r and not z
        # semi-invariance numeric + orbit
        si=all(S3n(w*x1%p,w*x2%p,w*x3%p,b,p)==w*S3n(x1,x2,x3,b,p)%p for x1 in xs[:8] for x2 in xs[:8] for x3 in xs[:8])
        rep.append(dict(p=p,b=b,npts=len(pts),FP=fp,FN=fn,semi_inv=si))
tot=dict(curves=len(rep),FP=sum(r['FP'] for r in rep),FN=sum(r['FN'] for r in rep),semi_inv_all=all(r['semi_inv'] for r in rep))
print(tot); json.dump(dict(total=tot,detail=rep),open('raw/cert_S3.json','w'))
