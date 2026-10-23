"""GF(2^n) polynomial-basis arithmetic + binary EC group law + binary S3 cert.
Independent of the crypto repo."""
import itertools, json
# irreducible polynomials (as bitmask incl top bit) for small n
IRRED={2:0b111,3:0b1011,4:0b10011,5:0b100101,6:0b1000011,7:0b10001001,8:0b100011011,
       9:0b1000010001,10:0b10000001001,11:0b100000000101,12:0b1000001010011,13:0b10000000011011}
def gmul(x,y,n):
    m=IRRED[n]; r=0
    while y:
        if y&1: r^=x
        y>>=1; x<<=1
        if x>>n & 1: x^=m
    return r
def gsq(x,n): return gmul(x,x,n)
def ginv(x,n):
    # x^(2^n-2)
    r=1; e=(1<<n)-2; b=x
    while e:
        if e&1: r=gmul(r,b,n)
        b=gsq(b,n); e>>=1
    return r
def gdiv(x,y,n): return gmul(x,ginv(y,n),n)

def bin_add(P,Q,a2,a6,n):
    # E: y^2+xy=x^3+a2 x^2+a6 ; -P=(x,x+y)
    if P is None: return Q
    if Q is None: return P
    x1,y1=P; x2,y2=Q
    if x1==x2:
        if (y1^y2)==x1 or x1==0: # Q=-P  (or 2-torsion x=0)
            if (y1^y2)==x1: return None
        if P==Q:
            if x1==0: return None
            lam=x1 ^ gdiv(y1,x1,n)
            x3=gsq(lam,n)^lam^a2
            y3=gsq(x1,n) ^ gmul(lam^1,x3,n)
            return (x3,y3)
        # x1==x2, not equal, not inverse -> inverse actually
        return None
    lam=gdiv(y1^y2,x1^x2,n)
    x3=gsq(lam,n)^lam^x1^x2^a2
    y3=gmul(lam,x1^x3,n)^x3^y1
    return (x3,y3)

def bin_points(a2,a6,n):
    pts=[None]
    for x in range(1<<n):
        rhs=gmul(gsq(x,n),x,n)^gmul(a2,gsq(x,n),n)^a6  # x^3+a2 x^2+a6
        for y in range(1<<n):
            if gsq(y,n)^gmul(x,y,n)==rhs: pts.append((x,y))
    return pts

def S3_bin(x1,x2,x3,a6,n):
    # candidate: (x1x2+x1x3+x2x3)^2 + x1x2x3 + a6
    e2=gmul(x1,x2,n)^gmul(x1,x3,n)^gmul(x2,x3,n)
    return gsq(e2,n) ^ gmul(gmul(x1,x2,n),x3,n) ^ a6

def cert(a2,a6,n):
    pts=[P for P in bin_points(a2,a6,n) if P is not None]
    fp=fn=0; tp=tn=0
    negs={}
    for P in pts: negs[P]=(P[0],P[0]^P[1])
    for P,Q,R in itertools.product(pts,repeat=3):
        s=(S3_bin(P[0],Q[0],R[0],a6,n)==0)
        anysum=False
        for sP in (P,negs[P]):
            for sQ in (Q,negs[Q]):
                for sR in (R,negs[R]):
                    if bin_add(bin_add(sP,sQ,a2,a6,n),sR,a2,a6,n) is None: anysum=True;break
                if anysum:break
            if anysum:break
        if s and anysum:tp+=1
        elif s and not anysum:fp+=1
        elif not s and anysum:fn+=1
        else:tn+=1
    return dict(n=n,a2=a2,a6=a6,npts=len(pts),tp=tp,fp=fp,fn=fn,tn=tn,ok=(fp==0 and fn==0))

if __name__=='__main__':
    res=[]
    for (a2,a6,n) in [(0,1,3),(1,1,3),(0,1,4),(1,1,4),(1,1,5),(0,1,5)]:
        res.append(cert(a2,a6,n))
    for r in res:print(r)
    print("ALL_OK",all(r['ok'] for r in res))
    json.dump(res,open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'out','cert_binary.json'),'w'))
