#!/usr/bin/env python3
"""Generate a genuine GLV+GLS four-term supplementary Fp2 case.

The four maps are 1, phi, psi, phi*psi. They are not four powers of one CM
automorphism. The ordinary endomorphism algebra is still quadratic; four-term
decomposition concerns a finite subgroup and a rank-four congruence lattice.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ecbench.extension import Fp2,ExtensionCurve
from prime_vectors import sqrt_mod
from sympy import Matrix,Rational,factorint,n_order,prevprime

P=2199023254699
TRACE=-1456168
R=1208925818676176351078857


def point_json(point):
    return None if point is None else {"x":str(point[0]),"y":str(point[1])}


def strings(value):
    if isinstance(value,bool) or value is None:return value
    if isinstance(value,int):return str(value)
    if isinstance(value,dict):return {str(k):strings(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [strings(v) for v in value]
    return value


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def gram_schmidt(basis):
    stars=[];norms=[]
    for row in basis:
        v=list(map(Fraction,row))
        for star,norm in zip(stars,norms):
            mu=dot(row,star)/norm
            v=[x-mu*y for x,y in zip(v,star)]
        stars.append(v);norms.append(dot(v,v))
    assert all(norm>0 for norm in norms)
    return stars,norms


def babai(k,r,basis):
    residual=[k%r]+[0]*(len(basis)-1)
    stars,norms=gram_schmidt(basis)
    for row,star,norm in reversed(list(zip(basis,stars,norms))):
        coefficient=dot(residual,star)/norm
        # Exact nearest-integer rounding, ties towards +infinity.
        rounded=(coefficient.numerator*2+coefficient.denominator)//(2*coefficient.denominator)
        residual=[x-rounded*y for x,y in zip(residual,row)]
    return residual


def lucas_certificates(roots):
    certificates={}
    def certify(n):
        if str(n) in certificates:return
        if n==2:
            certificates[str(n)]={"method":"base","n":n}
            return
        factors={int(q):int(e) for q,e in factorint(n-1).items()}
        assert math.prod(q**e for q,e in factors.items())==n-1
        for q in factors:certify(q)
        for a in range(2,10000):
            if pow(a,n-1,n)==1 and all(math.gcd(pow(a,(n-1)//q,n)-1,n)==1 for q in factors):
                certificates[str(n)]={"method":"Lucas","n":n,"witness":a,
                                      "n_minus_one_factorization":[{"prime":q,"exponent":e} for q,e in factors.items()]}
                return
        raise RuntimeError("No Lucas primality witness found")
    for root in roots:certify(root)
    return certificates


def build():
    p,r,t=P,R,TRACE;q=p*p;N=(p-1)**2+t*t
    assert N==4*r
    f=Fp2(p,2);eta=f.pack(1,1);B=f.power(eta,3);E=ExtensionCurve(f,0,B)
    assert f.unpack(B)==(7,5)
    assert f.power(eta,(q-1)//2)==p-1
    source=next(E.scan());G=E.mul(2,source)
    assert G is not None and E.mul(r,G) is None
    assert E.mul(2*r,source) is None and E.mul(r,source) is not None
    T2=E.mul(r,source)
    zeta=pow(2,(p-1)//3,p)
    p_minus_factors={int(prime):int(e) for prime,e in factorint(p-1).items()}
    assert all(pow(2,(p-1)//prime,p)!=1 for prime in p_minus_factors)
    cx=f.power(eta,1-p);cy=f.power(eta,3*(1-p)//2)
    def phi(point):
        return None if point is None else (f.mul(zeta,point[0]),point[1])
    def psi(point):
        return None if point is None else (f.mul(cx,f.frobenius(point[0])),f.mul(cy,f.frobenius(point[1])))
    lambda_psi=(p-1)*pow(t,-1,r)%r
    square_root=sqrt_mod(-3,r)
    roots=[(-1+square_root)*pow(2,-1,r)%r,(-1-square_root)*pow(2,-1,r)%r]
    lambda_phi=next(value for value in roots if E.mul(value,G)==phi(G))
    assert E.mul(lambda_psi,G)==psi(G)
    assert (lambda_phi*lambda_phi+lambda_phi+1)%r==0
    assert lambda_psi*lambda_psi%r==r-1
    eigenvalues=[1,lambda_phi,lambda_psi,lambda_phi*lambda_psi%r]
    assert len(set(eigenvalues))==4
    basis_points=[G,phi(G),psi(G),phi(psi(G))]
    assert phi(psi(G))==psi(phi(G))
    # A balanced exact CM seed avoids an assertion bug in SymPy's LLL on the
    # very unbalanced 80-bit canonical rank-four kernel basis.
    cm_c=math.isqrt((4*p-t*t)//3)
    assert t*t-4*p==-3*cm_c*cm_c
    cm_a=(t+cm_c)//2
    assert cm_a*cm_a-cm_a*cm_c+cm_c*cm_c==p
    cm_mu=(lambda_psi-cm_a)*pow(cm_c,-1,r)%r
    actual_cm_a,actual_cm_c=cm_a,cm_c
    h=(cm_a-1)//2;j=(cm_a+1)//2;c=cm_c//2
    seed=[[h,c,-j,-c],[-j,-c,-h,-c],[c-j,-j,c-h,-h],[c-h,-h,-(c-j),j]]
    if lambda_phi!=cm_mu:
        assert lambda_phi==(-1-cm_mu)%r
        seed=[[a-b,-b,c-d,-d] for a,b,c,d in seed]
        actual_cm_a,actual_cm_c=cm_a-cm_c,-cm_c
    assert (actual_cm_a+actual_cm_c*lambda_phi)%r==lambda_psi
    basis4=[[int(x) for x in row] for row in Matrix(seed).lll(delta=Rational(3,4)).tolist()]
    basis2=[[int(x) for x in row] for row in Matrix([[r,0],[-lambda_phi,1]]).lll(delta=Rational(3,4)).tolist()]
    for basis,values in [(basis4,eigenvalues),(basis2,eigenvalues[:2])]:
        assert abs(int(Matrix(basis).det()))==r
        assert all(dot(row,values)%r==0 for row in basis)
    universal_bounds={}
    for name,basis in [("glv4",basis4),("glv2",basis2)]:
        stars,_=gram_schmidt(basis)
        bounds=[sum(abs(row[j]) for row in stars)/2 for j in range(len(basis))]
        ceilings=[(value.numerator+value.denominator-1)//value.denominator for value in bounds]
        universal_bounds[name]={"coordinate_absolute_bound_rationals":
                               [{"numerator":value.numerator,"denominator":value.denominator} for value in bounds],
                               "coordinate_absolute_bound_ceilings":ceilings,
                               "guaranteed_component_bits":max(value.bit_length() for value in ceilings)}
    # Each scalar is reconstructed with signed multiplications of genuinely
    # different efficient maps, and also checked by ordinary scalar multiply.
    scalars=[0,1,-1,r-1,r,r+1,(1<<79)+0x12345,(1<<128)+1]
    scalars += [int.from_bytes(hashlib.sha256(f"ecbench-glv4-{i}".encode()).digest(),"big")%r for i in range(8)]
    decomp=[]
    for k in scalars:
        c4=babai(k,r,basis4);c2=babai(k,r,basis2)
        assert dot(c4,eigenvalues)%r==k%r and dot(c2,eigenvalues[:2])%r==k%r
        Q4=None;Q2=None
        for coefficient,point in zip(c4,basis_points):Q4=E.add(Q4,E.mul(coefficient,point))
        for coefficient,point in zip(c2,basis_points[:2]):Q2=E.add(Q2,E.mul(coefficient,point))
        expected=E.mul(k,G)
        assert Q4==Q2==expected
        assert E.joint_mul(c4,basis_points)==expected
        assert E.joint_mul(c2,basis_points[:2])==expected
        decomp.append({"scalar":k,"glv4_components":c4,"glv2_components":c2,
                       "glv4_max_component_bits":max(abs(x).bit_length() for x in c4),
                       "glv2_max_component_bits":max(abs(x).bit_length() for x in c2),
                       "expected_point":point_json(expected)})
    points={"O":None,"P":source,"R":source,"G":G,"T2":T2,
            "phi_G":phi(G),"psi_G":psi(G),"phi_psi_G":phi(psi(G))}
    maps=[]
    for name in ["O","G","R","T2"]:
        point=points[name]
        assert E.on_curve(phi(point)) and E.on_curve(psi(point))
        assert phi(phi(phi(point)))==point
        assert psi(psi(point))==E.neg(point)
        assert phi(psi(point))==psi(phi(point))
        assert E.add(E.mul(actual_cm_a,point),E.mul(actual_cm_c,phi(point)))==psi(point)
        for op,image in [("phi",phi(point)),("psi",psi(point)),
                         ("phi_psi",phi(psi(point))),("psi_phi",psi(phi(point))),
                         ("phi_cubed",phi(phi(phi(point)))),("psi_squared",psi(psi(point)))]:
            maps.append({"operation":op,"point_name":name,"point":point_json(point),"expected_point":point_json(image)})
    k=int(n_order(q,r))
    assert pow(q,k,r)==1
    for prime in factorint(k):assert pow(q,k//int(prime),r)!=1
    p_list=[];candidate=int(prevprime(1<<41))
    while len(p_list)<13:
        if candidate%3==1:p_list.append(candidate)
        candidate=int(prevprime(candidate))
    assert p_list[-1]==p
    field_values=[("eta",eta),("twist_coefficient",B),("psi_x_multiplier",cx),("psi_y_multiplier",cy)]
    field_vectors=[]
    for name,value in field_values:
        for operation,answer in [("square",f.square(value)),("inverse",f.inverse(value)),("frobenius",f.frobenius(value))]:
            field_vectors.append({"operation":operation,"name":name,"input":value,"expected":answer})
    certs=lucas_certificates([p,r,2]+list(map(int,factorint(r-1))))
    T=2*p-t*t
    result={"schema":"ecbench-glv-gls-supplement-v1","supplement_id":"j0-glv-gls4-p41",
        "experimental":True,"primary_dataset_modified":False,
        "field":{"kind":"quadratic_prime_extension","base_prime":p,"cardinality":q,"element_bits":q.bit_length(),
                 "nonresidue_d":2,"polynomial":"u^2-2",
                 "packing":"c0 + p*c1, 0 <= c0,c1 < p; packed coordinates are decimal strings; this is not SEC1",
                 "eta":eta,"eta_coefficients":[1,1],"eta_norm":p-1},
        "base_curve":{"equation":"y^2=x^3+1 over Fp","b":1,"trace":t,"order":p+1-t,
                      "Frobenius_discriminant":t*t-4*p,"CM_discriminant":-3},
        "curve":{"equation":"y^2=x^3+(7+5u) over Fp2","a4":0,"a6":B,"a6_coefficients":[7,5],
                 "trace":T,"order":N,"order_factorization":[{"prime":2,"exponent":2},{"prime":r,"exponent":1}],
                 "group_invariant_factors":[2,2*r],"Frobenius_discriminant":T*T-4*q,
                 "Frobenius_order_conductor":abs(t)*cm_c,"endomorphism_discriminant":-3,
                 "endomorphism_ring":"Z[zeta3], maximal; j=0 and ordinary",
                 "subgroup_order":r,"subgroup_bits":r.bit_length(),"cofactor":4,
                 "generator":point_json(G),"generator_projector":2,"generator_source":point_json(source),
                 "generator_rule":"ascending packed x from zero; minimum packed square root; first nonzero [d2/r]P",
                 "embedding_degree":k},
        "endomorphisms":{"phi":{"formula":"(zeta3*x,y)","zeta3":zeta,"base_primitive_root":2,"degree":1,"subgroup_eigenvalue":lambda_phi},
                         "psi":{"formula":"(eta^(1-p)*x^p, eta^(3*(1-p)/2)*y^p)","x_multiplier":cx,"y_multiplier":cy,
                                  "degree":p,"subgroup_eigenvalue":lambda_psi,
                                  "global_relation":"psi^2=-pi_q; psi^2-t0*psi+[p]=0",
                                  "actual_CM_linear_relation":{"formula":"psi=[a]+[c]*phi","a":actual_cm_a,"c":actual_cm_c},
                                  "rational_point_relation":"psi^2=[-1] on E(Fp2)",
                                  "eigenvalue_derivation":"lambda_psi=(p-1)/t0 mod r"},
                         "commute":True,"four_term_maps":["1","phi","psi","phi*psi"],"four_term_eigenvalues":eigenvalues,
                         "proof_boundary":"The geometric endomorphism algebra remains quadratic. Four-dimensional refers to a rank-four congruence lattice for distinct efficient GLV and GLS maps on a finite subgroup."},
        "lattices":{"definition":"rows b satisfy sum b_i*lambda_i = 0 mod r; determinant magnitude r",
                    "LLL_delta":"3/4","glv4":{"balanced_seed":seed,"basis":basis4,"determinant":int(Matrix(basis4).det()),"cm_a":cm_a,"cm_c":cm_c},
                    "glv2":{"basis":basis2,"determinant":int(Matrix(basis2).det())},
                    "Babai_universal_bounds":universal_bounds,
                    "decomposition":"exact rational Gram-Schmidt and Babai nearest plane; ties to +infinity; input normalized mod r"},
        "named_points":{name:{"point":point_json(point),"order":1 if point is None else (2*r if name in ["P","R"] else 2 if name=="T2" else r)} for name,point in points.items()},
        "field_vectors":field_vectors,"endomorphism_vectors":maps,"decomposition_vectors":decomp,
        "certificates":certs,
        "search":{"rule":"descending primes below 2^41 with p=1 mod3; six j0 twist representatives b=g^j; accept first N=(p-1)^2+t0^2 with h<=12 and prime r>=2^74",
                  "fields_examined":13,"examined_base_primes":p_list,"selected_b_index":0},
        "verification":{"independent_Fp2_and_affine_ECC":True,"generator_exact_order":True,
                        "eigenvalue_point_checks":True,"phi_psi_commute":True,"unit_order_restriction_checks":True,
                        "kernel_lattice_determinants":True,"Babai_reconstruction_vectors":len(decomp),
                        "max_glv4_component_bits":max(v["glv4_max_component_bits"] for v in decomp),
                        "max_glv2_component_bits":max(v["glv2_max_component_bits"] for v in decomp)},
        "limitations":["About 40 bits of generic subgroup security; experimental, not a production parameter set.",
                       "No SEC1 convention is asserted for the packed Fp2 representation.",
                       "The small CM discriminant -3 is intentional; this case compares scalar decompositions, not generic-curve arithmetic.",
                       "For j=1728, a GLS order-four restriction collapses with the Gaussian unit, and the analogous order factors into roughly p-sized pieces."]}
    return result


def pari_crosscheck(data):
    p=P;r=R;N=4*r
    executable=os.environ.get("CURVE_REPORT_GP") or shutil.which("gp") or "/tmp/pari-local/usr/bin/gp"
    def element(value):
        value=int(value);return f"({value%p}+{value//p}*u)"
    def point(record):
        return "[0]" if record is None else "["+element(record["x"])+","+element(record["y"])+"]"
    lines=[f"p={p};u=ffgen(Mod(z^2-2,p),'u);E=ellinit([0*u,7+5*u]);checks=0;"]
    checks=0
    def check(expression,label):
        nonlocal checks
        lines.append(f"if({expression},checks++,print({json.dumps('FAIL '+label)});quit(1));");checks+=1
    check(f"ellcard(E)=={N}","twist cardinality")
    check(f"ellgroup(E)==[{2*r},2]","twist point group")
    lines.append("E0=ellinit([Mod(0,p),Mod(1,p)]);")
    check(f"ellcard(E0)=={p+1-TRACE}","base cardinality")
    for name,record in data["named_points"].items():
        check(f"ellorder(E,{point(record['point'])},{N})=={record['order']}","order "+name)
    G=point(data["curve"]["generator"])
    phi=data["endomorphisms"]["phi"];psi=data["endomorphisms"]["psi"]
    for vector in data["decomposition_vectors"]:
        expected=point(vector["expected_point"]);k=int(vector["scalar"])
        check(f"ellmul(E,{G},{k})=={expected}","scalar known answer")
    for vector in data["endomorphism_vectors"]:
        Q=vector["point"];expected=point(vector["expected_point"])
        if Q is None:expression="[0]"
        else:
            x=element(Q["x"]);y=element(Q["y"]);op=vector["operation"]
            if op in ["phi","phi_psi","psi_phi"]:
                if op=="phi":expression=f"[{phi['zeta3']}*{x},{y}]"
                else:expression=f"[{phi['zeta3']}*{element(psi['x_multiplier'])}*{x}^p,{element(psi['y_multiplier'])}*{y}^p]"
            elif op=="psi":expression=f"[{element(psi['x_multiplier'])}*{x}^p,{element(psi['y_multiplier'])}*{y}^p]"
            elif op=="phi_cubed":expression=point(Q)
            else:expression=f"ellneg(E,{point(Q)})"
        check(f"{expression}=={expected}","map "+vector["operation"])
    lines.append("print([checks,version()]);")
    result=subprocess.run([executable,"-f","-q","-s","256M"],input="\n".join(lines)+"\n",text=True,capture_output=True,check=True)
    assert not result.stderr,result.stderr
    actual,version=json.loads(result.stdout)
    assert actual==checks
    data["verification"]["PARI_independent_checks"]=checks
    data["verification"]["PARI_version"]=".".join(map(str,version))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=ROOT/"build"/"glv4.json")
    parser.add_argument("--pari-crosscheck",action="store_true")
    args=parser.parse_args()
    data=build()
    if args.pari_crosscheck:pari_crosscheck(data)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(strings(data),indent=2)+"\n")
    print(f"Wrote genuine GLV+GLS4 case, max component bits {data['verification']['max_glv4_component_bits']} versus GLV2 {data['verification']['max_glv2_component_bits']}")


if __name__=="__main__":main()
