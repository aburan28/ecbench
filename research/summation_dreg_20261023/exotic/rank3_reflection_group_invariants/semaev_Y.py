"""Semaev S3,S4 for y^2=x^3+a x^2+c^2 x in the coordinate Y=(x-c)/(x+c); symbolic in a,c (sympy). Cached."""
import sympy as sp, pickle, os
H=os.path.dirname(os.path.abspath(__file__))
a,c=sp.symbols('a c')
Y1,Y2,Y3,YR,z=sp.symbols('Y1 Y2 Y3 YR z')

def xof(Y): return c*(1+Y)/(1-Y)

def S3x(x1,x2,x3):
    f=lambda x:x**3+a*x**2+c**2*x
    N=((x3+x1+x2+a)*(x2-x1)**2-f(x1)-f(x2))**2-4*f(x1)*f(x2)
    q,r=sp.div(sp.expand(N),sp.expand((x1-x2)**2),x1,x2,x3,a,c)
    assert r==0
    return sp.expand(q)

def S3Y(A,B,C):
    e=S3x(xof(A),xof(B),xof(C))
    e=sp.together(e)
    num,den=sp.fraction(sp.cancel(e))
    return sp.expand(num)

def build():
    X1,X2,X3=sp.symbols('X1 X2 X3')
    s3=S3x(X1,X2,X3)          # polynomial in x
    # transform
    def toY(A,B,C):
        e=sp.together(s3.subs({X1:xof(A),X2:xof(B),X3:xof(C)},simultaneous=True))
        num,den=sp.fraction(sp.cancel(e)); return sp.expand(num)
    P3=toY(Y1,Y2,Y3)
    P4=sp.expand(sp.resultant(toY(Y1,Y2,z),toY(Y3,YR,z),z))
    return {'S3':P3,'S4':P4,'S3x':s3}

def load():
    f=H+'/raw/semaevY.pkl'
    if os.path.exists(f): return pickle.load(open(f,'rb'))
    d=build(); pickle.dump(d,open(f,'wb')); return d

if __name__=='__main__':
    d=load()
    for k in ('S3','S4'):
        P=sp.Poly(d[k],Y1,Y2,Y3,YR,a,c)
        print(k,'nterms',len(P.terms()),'deg per var',[P.degree(v) for v in (Y1,Y2,Y3,YR)],'total',P.total_degree())
    P=sp.Poly(d['S3x'],*sp.symbols('X1 X2 X3')); print('S3x deg per var',[P.degree(v) for v in sp.symbols('X1 X2 X3')])
