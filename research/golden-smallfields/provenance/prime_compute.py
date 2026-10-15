#!/usr/bin/env python3
"""Reproducible, integer-only certificates for CM examples over four prime fields.

Point counts use the CM trace candidate theorem and finite scalar-multiplication
witnesses eliminating all wrong candidates. Primes use full-factorization Lucas
certificates; factor discovery uses SymPy but primality certification does not.
"""
import json
import math
import subprocess
import os
import shutil
from pathlib import Path
from functools import lru_cache
from sympy import factorint, sqrt_mod, primitive_root, n_order

OUT = Path(__file__).resolve().parent
FIELDS = [
    ("p61", "2^61 - 1", 2**61-1),
    ("p64", "2^64 - 2^32 + 1", 2**64-2**32+1),
    ("p83s", "2^83 - 2^59 - 1", 2**83-2**59-1),
    ("p83", "2^83 - 55", 2**83-55),
]
CERTS = {}
GP = os.environ.get("CURVE_REPORT_GP") or shutil.which("gp") or "/tmp/pari-local/usr/bin/gp"

@lru_cache(None)
def factors(n):
    return {int(k): int(v) for k,v in factorint(int(n)).items()}

def certify_prime(n):
    """Lucas criterion: complete certified factorization n-1 plus a witness."""
    key=str(n)
    if key in CERTS:
        return
    if n == 2:
        CERTS[key] = {"method": "base", "n": n}
        return
    if n < 2 or n%2 == 0:
        raise ValueError(f"not prime: {n}")
    fs=factors(n-1)
    assert math.prod(q**e for q,e in fs.items()) == n-1
    for q in fs:
        certify_prime(q)
    for a in range(2, 10000):
        if pow(a,n-1,n)==1 and all(math.gcd(pow(a,(n-1)//q,n)-1,n)==1 for q in fs):
            CERTS[key] = {"method":"Lucas", "n":n, "a":a,
                          "n_minus_1_factors":{str(q):e for q,e in fs.items()}}
            return
    raise RuntimeError(f"Lucas witness search failed for {n}")

def cornacchia(p,d):
    for r in sqrt_mod(-d,p,all_roots=True):
        a,b=p,int(r)
        while b*b > p:
            a,b=b,a%b
        rem=p-b*b
        if rem%d == 0:
            c=math.isqrt(rem//d)
            if b*b+d*c*c == p:
                return b,c
    raise ValueError((p,d))

def add(P,Q,A,p):
    if P is None: return Q
    if Q is None: return P
    x,y=P; z,w=Q
    if x==z:
        if (y+w)%p==0: return None
        slope=((3*x*x+A)*pow(2*y,-1,p))%p
    else:
        slope=((w-y)*pow(z-x,-1,p))%p
    r=(slope*slope-x-z)%p
    return r,(slope*(x-r)-y)%p

def mul(k,P,A,p):
    R=None
    while k:
        if k&1: R=add(R,P,A,p)
        P=add(P,P,A,p)
        k>>=1
    return R

def count_cm(p,A,B,candidates):
    candidates=set(candidates)
    assert all(t*t <= 4*p for t in candidates)
    eliminated=[]
    for x in range(1000):
        rhs=(x*x*x+A*x+B)%p
        if rhs==0:
            y=0
        elif pow(rhs,(p-1)//2,p)!=1:
            continue
        else:
            y=int(sqrt_mod(rhs,p))
        P=(x,y)
        for t in sorted(list(candidates)):
            R=mul(p+1-t,P,A,p)
            if R is not None:
                eliminated.append({"rejected_trace":t,"point":list(P),"multiple":list(R)})
                candidates.remove(t)
        if len(candidates)==1:
            return candidates.pop(),eliminated
    raise RuntimeError(f"Could not distinguish CM traces: {p,A,B,candidates}")

def order_data(N,p):
    fs=factors(N)
    for r in fs: certify_prime(r)
    r=max(fs)
    k=None if r==p else int(n_order(p,r))
    order_factor=factors(r-1)
    for ell in order_factor: certify_prime(ell)
    # Verify order exactly using its complete factorization.
    if k is not None:
        assert pow(p,k,r)==1
        for ell in factors(k):
            assert pow(p,k//ell,r)!=1
    return {"order_factors":{str(r):e for r,e in fs.items()},
            "largest_prime_factor":r,"largest_prime_bits":r.bit_length(),
            "cofactor_for_largest_prime":N//r,"embedding_degree":k,
            "r_minus_1_factors":{str(q):e for q,e in order_factor.items()},
            "generic_security_bits_approx":math.log2(r)/2}

def pari(lines):
    result=subprocess.run([GP,"-f","-q","-s","256M"],
                          input="\n".join(lines)+"\n",text=True,
                          capture_output=True,check=True)
    if result.stderr.strip():
        raise RuntimeError(result.stderr)
    return [json.loads(line) for line in result.stdout.splitlines() if line.strip()]

def fundamental_discriminant(D):
    assert D<0 and D%4 in (0,1)
    fs=factors(-D)
    for q in fs: certify_prime(q)
    d=-math.prod(q for q,e in fs.items() if e%2)
    DK=d if d%4==1 else 4*d
    f=math.isqrt(D//DK)
    assert DK*f*f==D
    return DK,f,fs

def make_generic(p,g,B,N,group,proof=None):
    t=p+1-N
    D=t*t-4*p
    DK,f,fs=fundamental_discriminant(D)
    assert t%p != 0 # ordinary
    A=(-3)%p
    delta=(-16*(4*(-3)**3+27*B**2))%p
    j=(1728*4*(-3)**3*pow(4*(-3)**3+27*B**2,-1,p))%p
    row={"index":0,"A":A,"B":B,"trace":t,"order":N,
         "curve_discriminant_mod_p":delta,"j_mod_p":j,
         "Frobenius_discriminant":D,"Frobenius_discriminant_absolute_factors":
             {str(q):e for q,e in fs.items()},
         "quadratic_field_discriminant":DK,"Frobenius_order_conductor":f,
         "status":"ordinary", "count_method":"PARI/GP ellcard 2.17.2",
         "group_invariant_factors":[group[1] if len(group)>1 else 1,group[0]],
         "endomorphism_ring_over_Fp":"O_K", "endomorphism_discriminant_over_Fp":DK,
         "geometric_endomorphism_ring":"O_K",
         "endomorphism_conductor":1,
         "endomorphism_proof": proof or "Frobenius discriminant is fundamental (f_pi=1)."}
    row.update(order_data(N,p))
    return row

def enrich_pari(result):
    lines=[];cm_rows=[]
    for field in result["fields"]:
        p=field["p"]
        for fam in field["families"]:
            for v in fam["twists"]:
                lines.append(f'E=ellinit([Mod({v["A"]},{p}),Mod({v["B"]},{p})]);print([ellcard(E),ellgroup(E)]);')
                cm_rows.append(v)
    checks=pari(lines)
    assert len(checks)==len(cm_rows)
    for v,(N,group) in zip(cm_rows,checks):
        assert v["order"]==N
        assert v["group_invariant_factors"]==[group[1] if len(group)>1 else 1,group[0]]
        v["independent_crosscheck"]="PARI/GP 2.17.2 ellcard and ellgroup agree"
    # Generic control, together with its central quadratic twist by nonsquare g.
    lines=[]
    for field in result["fields"]:
        p=field["p"];g=field["primitive_root"]
        for i in (0,1):
            A=(-3*pow(g,2*i,p))%p; B=pow(g,3*i,p)
            lines.append(f'E=ellinit([Mod({A},{p}),Mod({B},{p})]);print([ellcard(E),ellgroup(E)]);')
    generic_checks=iter(pari(lines))
    # At p61, f_pi=62; exact division-polynomial identity excludes conductor31.
    p=FIELDS[0][2]
    test=pari([f'p={p};E=ellinit([Mod(-3,p),Mod(1,p)]);F=elldivpol(E,31);print([poldegree(F),Mod(x,F)^p==Mod(x,F)]);'])
    assert test==[[480,1]]
    for field in result["fields"]:
        p=field["p"];g=field["primitive_root"]
        N,group=next(generic_checks)
        DK,f,_=fundamental_discriminant((p+1-N)**2-4*p)
        proof=None
        if field["id"]=="p61":
            assert f==62 and len(group)==2 and group[1]%2==0
            assert (p+1-N)//2%31==30 and p%31==1
            proof=("f_pi=62. Full rational 2-torsion proves (pi-1)/2 is an endomorphism, "
                   "excluding factor2 from f_E. For ell=31, trace/2=-1 mod31 and "
                   "x^p=x modulo the degree480 division polynomial psi_31, so "
                   "Frobenius acts as [-1] on E[31]. Thus (pi+1)/31 is an endomorphism, "
                   "excluding factor31. Consequently f_E=1.")
        elif field["id"]=="p83s":
            assert f==2 and len(group)==2 and group[1]%2==0
            proof="f_pi=2. Full rational 2-torsion gives (pi-1)/2 in End(E), so f_E=1."
        else:
            assert f==1
        base=make_generic(p,g,1,N,group,proof)
        Nt,groupt=next(generic_checks)
        assert N+Nt==2*p+2
        # Twists have isomorphic geometric endomorphism rings; all ordinary endomorphisms are F_p-defined.
        twist=make_generic(p,g,1,Nt,groupt,
                           "Quadratic twist of ordinary base with proven maximal End=O_K; twisting preserves this order.")
        twist.update({"index":1,"A":(-3*g*g)%p,"B":pow(g,3,p),
                      "curve_discriminant_mod_p":(1296*pow(g,6,p))%p,
                      "quadratic_twist_index":0})
        base["quadratic_twist_index"]=1
        field["families"].append({"family":"G","base_equation":"y^2=x^3-3x+1",
                                 "j_integer":2304,"j_mod_p":2304%p,
                                 "status":"ordinary","twist_class_count":2,"twists":[base,twist]})
    result["pari_crosschecked_CM_rows"]=len(cm_rows)
    # Deterministic modest search for a prime/low-cofactor example with a large twist subgroup.
    p=FIELDS[-1][2];g=result["fields"][-1]["primitive_root"]
    search=pari([f'p={p};',
      'for(b=1,1000,if(b==2,next);E=ellinit([Mod(-3,p),Mod(b,p)]);N=ellcard(E);found=0;for(h=1,4,if(N%h==0 && N/h>=2^79 && isprime(N/h),Nt=2*p+2-N;ft=factor(Nt);rt=ft[matsize(ft)[1],1];if(rt>=2^74,print([b,h,N,Nt,rt]);found=1;break)));if(found,break));'])
    assert len(search)==1
    B,h,N,Nt,rt=search[0]
    checks=pari([f'E=ellinit([Mod(-3,{p}),Mod({B},{p})]);print([ellcard(E),ellgroup(E)]);',
                 f'E=ellinit([Mod({(-3*g*g)%p},{p}),Mod({B*pow(g,3,p)%p},{p})]);print([ellcard(E),ellgroup(E)]);'])
    assert checks[0][0]==N and checks[1][0]==Nt
    base=make_generic(p,g,B,N,checks[0][1])
    assert base["Frobenius_order_conductor"]==1
    twist=make_generic(p,g,B,Nt,checks[1][1],
                       "Quadratic twist preserves maximal ordinary endomorphism order O_K.")
    twist.update({"index":1,"A":(-3*g*g)%p,"B":B*pow(g,3,p)%p,
                  "curve_discriminant_mod_p":base["curve_discriminant_mod_p"]*pow(g,6,p)%p,
                  "quadratic_twist_index":0})
    base["quadratic_twist_index"]=1
    assert twist["largest_prime_factor"]==rt
    result["searched_example"]={"field_id":"p83","p":p,"nonsquare_twist_parameter":g,
        "base_equation":f"y^2=x^3-3x+{B}","b":B,
        "search_rule":"Try b=1,2,...,1000 in y^2=x^3-3x+b; skip singular b=2. Accept first curve with N/h prime >=2^79 for h in {1,2,3,4}, and largest twist prime >=2^74.",
        "last_b_tested":B,"nonsingular_curves_tested":B-1,"twists":[base,twist]}
    print("searched_example",B,"N",N,"Nt",Nt,flush=True)

def main():
    result={"scope":"CM examples and all F_p-isomorphism twist classes at j=1728 and j=0",
            "certificate_method":"Lucas primality; CM candidate traces eliminated by exact point multiplications",
            "fields":[]}
    for key,form,p in FIELDS:
        certify_prime(p)
        g=int(primitive_root(p))
        assert all(pow(g,(p-1)//r,p)!=1 for r in factors(p-1))
        field={"id":key,"form":form,"p":p,"bits":p.bit_length(),
               "p_mod_12":p%12,"primitive_root":g,
               "p_minus_1_factors":{str(q):e for q,e in factors(p-1).items()},
               "families":[]}
        for family in [4,3]:
            ordinary=(p%4==1) if family==4 else (p%3==1)
            m=4 if family==4 else 6
            number=math.gcd(m,p-1)
            rep=None
            if ordinary:
                a,b=cornacchia(p,1 if family==4 else 3)
                rep={"a":a,"b":b,"d":1 if family==4 else 3,
                     "identity":f"p = {a}^2 + {1 if family==4 else 3}*{b}^2"}
                candidates=({2*a,-2*a,2*b,-2*b} if family==4 else
                            {2*a,-2*a,a+3*b,-a-3*b,a-3*b,-a+3*b})
            else:
                candidates={0}
            fam={"family":"E4" if family==4 else "E3",
                 "base_equation":"y^2 = x^3 - x" if family==4 else "y^2 = x^3 + 1",
                 "j_integer":1728 if family==4 else 0,"j_mod_p":1728%p if family==4 else 0,
                 "status":"ordinary" if ordinary else "supersingular",
                 "CM_representation":rep,"CM_candidate_traces":sorted(candidates),
                 "twist_class_count":number,"twists":[]}
            for i in range(number):
                A=(-pow(g,i,p))%p if family==4 else 0
                B=0 if family==4 else pow(g,i,p)
                t,witnesses=count_cm(p,A,B,candidates)
                N=p+1-t
                delta=(-16*(4*A**3+27*B**2))%p
                Dpi=t*t-4*p
                if ordinary:
                    DK=-4 if family==4 else -3
                    f=math.isqrt(Dpi//DK)
                    assert DK*f*f==Dpi
                    if family==4:
                        aa=t//2
                        bb=math.isqrt(p-aa*aa)
                        assert aa*aa+bb*bb==p
                        d1=math.gcd(aa-1,bb)
                    else:
                        assert (t+f)%2==0
                        d1=math.gcd((t+f)//2-1,f)
                    end_disc=DK
                    ring="Z[i]" if family==4 else "Z[(1+sqrt(-3))/2]"
                    geometric_ring=ring
                else:
                    DK=-p if p%4==3 else -4*p
                    f=2 if p%4==3 else 1
                    full_two=(family==4 and pow((-A)%p,(p-1)//2,p)==1)
                    d1=2 if full_two else 1
                    end_disc=-p if full_two else -4*p
                    ring="O_{Q(sqrt(-p))}" if full_two else "Z[sqrt(-p)]"
                    geometric_ring="maximal order in quaternion algebra B_{p,infinity}"
                assert N%(d1*d1)==0
                twist={"index":i,"A":A,"B":B,"coefficient_form":
                       f"A=-g^{i}" if family==4 else f"B=g^{i}",
                       "quadratic_twist_index":((i+2) % number if family==4 else (i+3) % number),
                       "curve_discriminant_mod_p":delta,"trace":t,"order":N,
                       "Frobenius_discriminant":Dpi,"quadratic_field_discriminant":DK,
                       "Frobenius_order_conductor":f,
                       "endomorphism_ring_over_Fp":ring,"endomorphism_discriminant_over_Fp":end_disc,
                       "geometric_endomorphism_ring":geometric_ring,
                       "group_invariant_factors":[d1,N//d1],
                       "count_elimination_witnesses":witnesses}
                twist.update(order_data(N,p))
                fam["twists"].append(twist)
                print(key,fam["family"],i,"t",t,"N",N,"r",twist["largest_prime_factor"],
                      "k",twist["embedding_degree"],flush=True)
            field["families"].append(fam)
        result["fields"].append(field)
        (OUT/"prime_data.json").write_text(json.dumps(result,indent=2)+"\n")
        (OUT/"prime_certificates.json").write_text(json.dumps(CERTS,indent=2)+"\n")
    enrich_pari(result)
    (OUT/"prime_data.json").write_text(json.dumps(result,indent=2)+"\n")
    (OUT/"prime_certificates.json").write_text(json.dumps(CERTS,indent=2)+"\n")
    return result

if __name__=="__main__":
    main()
