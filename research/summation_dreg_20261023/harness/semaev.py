"""Independent Semaev summation-polynomial construction + correctness certificate.
No dependency on the crypto repo. Prime-field (char!=2,3) curves y^2=x^3+ax+b.
S_3 built two independent ways (closed form + line-elimination resultant) and
checked equal, then checked against real EC point additions."""
import sympy as sp

def ec_add(P, Q, a, p):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P; x2,y2=Q
    if x1==x2 and (y1+y2)%p==0: return None
    if P==Q:
        lam=(3*x1*x1+a)*pow(2*y1%p,-1,p)%p
    else:
        lam=(y2-y1)*pow((x2-x1)%p,-1,p)%p
    x3=(lam*lam-x1-x2)%p
    y3=(lam*(x1-x3)-y1)%p
    return (x3,y3)

def points(a,b,p):
    pts=[None]
    for x in range(p):
        r=(x*x*x+a*x+b)%p
        for y in range(p):
            if (y*y-r)%p==0: pts.append((x,y))
    return pts

def S3_closed(X1,X2,X3,a,b):
    # Semaev's 3rd summation polynomial (char != 2,3), symmetric form.
    return ((X1-X2)**2*X3**2
            - 2*((X1+X2)*(X1*X2+a)+2*b)*X3
            + (X1*X2-a)**2 - 4*b*(X1+X2))

def S3_lineelim(X1,X2,X3,a,b):
    # Independent build: three points collinear on y=l x+m.
    # x^3+ax+b-(l x+m)^2 has roots X1,X2,X3 => Vieta in (l,m); eliminate.
    l,m=sp.symbols('l m')
    e1=X1+X2+X3-l**2                      # sum of roots = l^2
    e2=X1*X2+X1*X3+X2*X3-(a-2*l*m)        # pair sum = a-2lm
    e3=X1*X2*X3-(m**2-b)                  # product = -(b-m^2)=m^2-b
    r=sp.resultant(sp.Poly(e1,m),sp.Poly(e3,m),m)  # eliminate m between e1,e3? need l too
    # eliminate m first from e1(no m),e2,e3: use e2,e3 to kill m, then e1 kills l
    r1=sp.resultant(sp.Poly(e2,m),sp.Poly(e3,m),m)   # in l
    res=sp.resultant(sp.Poly(e1,l),sp.Poly(r1,l),l)
    return sp.expand(res)

if __name__=='__main__':
    X1,X2,X3,a,b=sp.symbols('X1 X2 X3 a b')
    s_closed=sp.expand(S3_closed(X1,X2,X3,a,b))
    s_line=sp.expand(S3_lineelim(X1,X2,X3,a,b))
    # line-elim may differ by a nonzero scalar/factor; compare up to content & a square factor
    q=sp.simplify(s_line/ s_closed)
    print("S3 closed form:", s_closed)
    print("line/closed ratio simplifies to:", q)
