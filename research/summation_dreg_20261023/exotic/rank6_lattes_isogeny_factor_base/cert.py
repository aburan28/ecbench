# correctness certificate: S3 vs brute force (p=101), isogeny/Lattes maps send points to points
import random
from tower import *
p=101; fp=fn=0; n=0
for seed in range(3):
    a,b=rand_curve(p,random.Random(seed)); pts=[P for P in points(a,b,p) if P is not None]
    xs=sorted({P[0] for P in pts})
    neg=lambda P:(P[0],(-P[1])%p)
    for x1 in xs:
        for x2 in xs:
            # set of x(R) for R=+-P1+-P2 over actual points
            P1=[P for P in pts if P[0]==x1][0]; P2=[P for P in pts if P[0]==x2][0]
            Rx=set()
            for A in (P1,neg(P1)):
                for B in (P2,neg(P2)):
                    R=ec_add(A,B,a,p)
                    if R is not None: Rx.add(R[0])
            for xR in xs:
                n+=1; z=S3_closed(x1,x2,xR,a,b)%p==0
                if z and xR not in Rx: fp+=1
                if (not z) and xR in Rx: fn+=1
print("S3 brute-force checks",n,"false_pos",fp,"false_neg",fn)
rng=random.Random(5)
for q in (211,521):
    for ell,nn in ((2,3),(3,2),(2,4)):
        maps,curves=isogeny_chain(ell,nn,q,rng); print("p",q,"ell",ell,"chain",nn,"maps-points-to-points",check_isog_maps(maps,curves,q))
    a,b=rand_curve(q,rng); L=lattes(a,b,q); pts=[P for P in points(a,b,q) if P is not None]
    ok=all((lambda R:(R is None and peval(L[1],P[0],q)==0) or (R is not None and rapply(L,P[0],q)==R[0]))(ec_add(P,P,a,q)) for P in pts)
    print("p",q,"Lattes = x([2])",ok)
