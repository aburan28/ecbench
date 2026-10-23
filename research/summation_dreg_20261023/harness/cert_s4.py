import sys,random,itertools
H='/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness'
sys.path.insert(0,H)
import boolean_weil as bw
from gf2 import bin_points, bin_add

def const(x,n):
    e=bw.f_zero(n)
    for k in range(n):
        if (x>>k)&1: e[k]=bw.const1()
    return e
def S4num(x1,x2,x3,x4,a6,n):
    r=bw.S4_unsplit(const(x1,n),const(x2,n),const(x3,n),const(x4,n),a6,n)
    val=0
    for k in range(n):
        if r[k]:
            assert frozenset() in r[k] and len(r[k])==1
            val|=1<<k
    return val

# certify: S4(x1,x2,x3,x4)=0 iff exists sign choice P1±P2±P3±P4=O
def cert(a2,a6,n,samples=4000):
    pts=[P for P in bin_points(a2,a6,n) if P is not None]
    negs={P:(P[0],P[0]^P[1]) for P in pts}
    rng=random.Random(7); fp=fn=tp=tn=0
    for _ in range(samples):
        P1,P2,P3,P4=[rng.choice(pts) for _ in range(4)]
        s=(S4num(P1[0],P2[0],P3[0],P4[0],a6,n)==0)
        anysum=False
        for s1 in (P1,negs[P1]):
         for s2 in (P2,negs[P2]):
          for s3 in (P3,negs[P3]):
           for s4 in (P4,negs[P4]):
            acc=None
            for Q in (s1,s2,s3,s4): acc=bin_add(acc,Q,a2,a6,n)
            if acc is None: anysum=True
        if s and anysum:tp+=1
        elif s and not anysum:fp+=1
        elif not s and anysum:fn+=1
        else:tn+=1
    return dict(n=n,a2=a2,a6=a6,tp=tp,fp=fp,fn=fn,tn=tn,ok=(fp==0 and fn==0))
for (a2,a6,n) in [(1,1,4),(1,1,5),(0,1,5),(1,1,6)]:
    print(cert(a2,a6,n))
