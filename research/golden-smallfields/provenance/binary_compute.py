#!/usr/bin/env python3
"""Exact arithmetic for fourteen binary elliptic-curve examples.

Requires SymPy. All integer-prime assertions used in the report are supported by
recursive Lucas certificates with completely factored p-1. No point enumeration
in fields of size 2^53 to 2^83 is attempted: traces use the exact Frobenius
recurrence from E/F_2, and general-coefficient counts use PARI/GP.
"""
import json
import math
import os
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parent
POLYS = {53: [53,6,2,1,0], 79: [79,9,0], 83: [83,7,4,2,0]}


def factors(n):
    f = {int(p):int(e) for p,e in sp.factorint(n).items()}
    assert math.prod(p**e for p,e in f.items()) == n
    return f


@lru_cache(None)
def prime_certificate(n):
    """Lucas criterion: full n-1 factorization, one witness for each prime factor."""
    if n == 2:
        return {"n":2,"method":"base case"}
    assert n > 2 and n % 2
    fs = factors(n-1)
    certs = {str(p):prime_certificate(p) for p in fs}
    witnesses = {}
    for p in fs:
        for a in range(2, 1000):
            if pow(a,n-1,n)==1 and math.gcd(pow(a,(n-1)//p,n)-1,n)==1:
                witnesses[str(p)] = a
                break
        else:
            raise AssertionError(f"No Lucas witness found for {n}, factor {p}")
    return {"n":n,"method":"Lucas full-factorization criterion", "n_minus_1_factors":fs,
            "witness_for_each_prime_factor":witnesses,"prime_factor_certificates":certs}


def gf2_mod(a,b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a


def gf2_square(a,f):
    sq=0
    for i in range(a.bit_length()):
        if (a>>i)&1:
            sq ^= 1 << (2*i)
    return gf2_mod(sq,f)


def gf2_gcd(a,b):
    while b:
        a,b=b,gf2_mod(a,b)
    return a


def irreducible(m,exps):
    f=sum(1<<i for i in exps)
    powers={0:2}
    x=2
    for j in range(1,m+1):
        x=gf2_square(x,f)
        powers[j]=x
    # Rabin irreducibility test.
    gcds={str(p):gf2_gcd(powers[m//p]^2,f) for p in factors(m)}
    assert powers[m]==2 and all(v==1 for v in gcds.values())
    return {"degree":m,"exponents":exps,"irreducible":True,
            "rabin_x_2_power_m_equals_x":True,"rabin_gcds":gcds}


def cm_coords_power(a,m):
    # w=(1+sqrt(-7))/2, w^2=w-2. Frobenius w-1 (a=0) or w (a=1).
    A,B=(-1,1) if a==0 else (0,1)
    u,v=1,0
    for _ in range(m):
        u,v=u*A-2*v*B,u*B+v*A+v*B
    return u,v


def trace(t0,m):
    x,y=2,t0
    if m==0: return x
    for _ in range(2,m+1):
        x,y=y,t0*y-2*x
    return y


def class_number_conductor(f):
    h=sp.Integer(f)
    for p in factors(f):
        h *= 1-sp.Rational(sp.kronecker_symbol(-7,p),p)
    assert h.q==1
    return int(h)


def embedding_degree(q,r):
    assert q%r
    k=int(sp.n_order(q,r))
    kf=factors(k)
    assert pow(q,k,r)==1
    checks={str(p):pow(q,k//p,r)!=1 for p in kf}
    assert all(checks.values())
    return {"k":k,"k_factors":kf,"q_power_k_mod_r":1,
            "minimality_checks":checks,"ord_r_2":int(sp.n_order(2,r))}


def row(a,m):
    q=2**m
    t0=3-(4 if a==0 else 2)
    t=trace(t0,m)
    n=q+1-t
    u,v=cm_coords_power(a,m)
    assert t==2*u+v and u*u+u*v+2*v*v==q
    assert t*t-4*q==-7*v*v
    d1=math.gcd(u-1,v)
    assert n==((u-1)*(u-1)+(u-1)*v+2*v*v)
    assert n%d1==0 and (n//d1)%d1==0
    fs=factors(n)
    for p in fs:
        prime_certificate(p)
    r=max(p for p in fs if p!=2)
    assert fs[r]==1
    k=embedding_degree(q,r)
    odd_embeddings={str(p):embedding_degree(q,p) for p in fs if p!=2}
    f=abs(v)
    conductor_options=[]
    for c in sorted(int(x) for x in sp.divisors(f)):
        conductor_options.append({"conductor":c,"order_discriminant":-7*c*c,
                                  "order_class_number":class_number_conductor(c)})
    # Frobenius of F_2 acts on the unique prime-order subgroup with order m.
    sqroots=sp.sqrt_mod(-7,r,all_roots=True)
    lambdas=[((t0+int(s))*pow(2,-1,r))%r for s in sqroots]
    valid=[x for x in lambdas if pow(x,m,r)==1]
    assert len(valid)==1 and valid[0]!=1 and sp.isprime(m)
    return {"label":f"K{m}-a{a}","field_degree":m,"q":q,"a":a,"b":1,
            "equation":f"y^2 + x*y = x^3 + {a}*x^2 + 1",
            "base_field":"F_2","base_trace":t0,"trace":t,"group_order":n,
            "group_order_factors":fs,"group_invariant_factors":[d1,n//d1],
            "group_cyclic":d1==1,"largest_prime_subgroup":r,
            "r_bit_length":r.bit_length(),"r_log2":math.log2(r),"cofactor":n//r,
            "generic_security_log2_sqrt_r":math.log2(r)/2,
            "embedding_degree":k["k"],"embedding_degree_details":k,
            "all_odd_subgroup_embedding_degrees":odd_embeddings,
            "curve_discriminant_in_field":1,"j_invariant":1,"ordinary":True,
            "CM_field_discriminant":-7,"CM_field":"Q(sqrt(-7))",
            "endomorphism_ring":"Z[(1 + sqrt(-7))/2]",
            "endomorphism_ring_conductor":1,"endomorphism_ring_discriminant":-7,
            "endomorphism_ring_class_number":1,
            "frobenius_discriminant":t*t-4*q,
            "frobenius_order_conductor":f,"frobenius_order_conductor_factors":factors(f),
            "frobenius_order_class_number":class_number_conductor(f),
            "possible_endomorphism_orders_in_isogeny_class":conductor_options,
            "frobenius_CM_coordinates_basis_1_w":[u,v],
            "twist_label":f"K{m}-a{1-a}",
            "F2_frobenius_scalar_on_largest_subgroup":valid[0],
            "F2_frobenius_order_on_largest_subgroup":m,
            "log2_sqrt_r_divided_by_sqrt_m":(math.log2(r)-math.log2(m))/2}



def gp_binary_controls():
    """PARI cardinalities and rational groups for B-style examples and search."""
    gp=os.environ.get('CURVE_REPORT_GP') or shutil.which('gp') or '/tmp/pari-local/usr/bin/gp'
    code=r"""
default(parisizemax,1000000000)
for(i=1,3,m=[53,79,83][i];f=[x^53+x^6+x^2+x+1,x^79+x^9+1,x^83+x^7+x^4+x^2+1][i];u=ffgen(Mod(f,2),'u);for(a=0,1,E=ellinit([1,a,0,0,u],u);N=ellcard(E);t=2^m+1-N;print([m,a,1,N,t,t^2-4*2^m,ellgroup(E)])))
f=x^83+x^7+x^4+x^2+1;u=ffgen(Mod(f,2),'u);q=2^83;
found=0; hits=0; tested=0;
for(i=1,256,tested=i;b=u^i;E0=ellinit([1,0,0,0,b],u);N0=ellcard(E0);for(a=0,1,N=if(a==0,N0,2*q+2-N0);h=if(N%4==0,4,2);r=N/h;if(r>=2^79 && isprime(r),hits++;Nt=2*q+2-N;Ft=factor(Nt);rt=Ft[matsize(Ft)[1],1];if(!found,found=[i,a,N,h,r,Nt,rt]);if(rt>=2^74,found=[i,a,N,h,r,Nt,rt];break(2)))));
print([found[1],found[2],hits,tested]);
if(found,for(a=0,1,b=u^found[1];E=ellinit([1,a,0,0,b],u);N=ellcard(E);t=q+1-N;print([83,a,found[1],N,t,t^2-4*q,ellgroup(E)])))
"""
    proc=subprocess.run([gp,'-fq','-s','256M'],input=code,text=True,capture_output=True,timeout=60,check=True)
    lines=[json.loads(l) for l in proc.stdout.splitlines() if l.startswith('[')]
    assert len(lines)==9
    generic=[binary_control_row(*x) for x in lines[:6]]
    searched=[binary_control_row(*x, searched=True) for x in lines[7:]]
    assert lines[6]==[251,1,1,256]
    metadata={"rule":"b=u^i, i=1,...,256; a=0,1; require h<=4 and r>=2^79; prefer twist largest prime >=2^74",
              "tested_exponents":256,"curve_candidates":512,"primary_hits":1,
              "selected_exponent":251,"selected_a":1,"preferred_twist_threshold_met":False,
              "randomness":"Deterministic search; these are not independent random coefficients. Frobenius-related exponents duplicate arithmetic.",
              "cardinality_method":"PARI ellcard, with second twist order 2q+2-N during screening; final ellcard both sides",
              "group_method":"PARI ellgroup"}
    return generic,searched,metadata


def binary_control_row(m,a,exponent,n,t,D,group, searched=False):
    q=2**m
    assert math.gcd(exponent,q-1)==1  # u^exponent and its j remain outside F_2
    assert n==q+1-t and D==t*t-4*q and t%2
    fs=factors(n)
    df=factors(abs(D))
    for p in set(fs)|set(df):
        prime_certificate(p)
    d=-math.prod(p for p,e in df.items() if e%2)
    dk=d if d%4==1 else 4*d
    ratio=D//dk
    conductor=math.isqrt(ratio)
    assert conductor**2==ratio
    r=max(p for p in fs if p!=2)
    k=embedding_degree(q,r)
    assert math.prod(group)==n
    invariants=[1,n] if len(group)==1 else sorted(group)
    label=f"S{m}-a{a}" if searched else f"B{m}-a{a}"
    b='u' if exponent==1 else f'u^{exponent}'
    return {"label":label,"field_degree":m,"q":q,"a":a,"b":b,"b_exponent":exponent,
            "equation":f"y^2 + x*y = x^3 + {a}*x^2 + {b}",
            "trace":t,"group_order":n,"group_order_factors":fs,
            "group_invariant_factors":invariants,"group_cyclic":len(group)==1,
            "largest_prime_subgroup":r,"r_bit_length":r.bit_length(),"r_log2":math.log2(r),
            "cofactor":n//r,"generic_security_log2_sqrt_r":math.log2(r)/2,
            "embedding_degree":k['k'],"embedding_degree_details":k,
            "all_odd_subgroup_embedding_degrees":{str(p):embedding_degree(q,p) for p in fs if p!=2},
            "curve_discriminant_in_field":b,"j_invariant":f'({b})^-1',"ordinary":True,
            "frobenius_discriminant":D,"frobenius_discriminant_abs_factors":df,
            "CM_field_discriminant":dk,"CM_field":f'Q(sqrt({dk}))',
            "frobenius_order_conductor":conductor,
            "endomorphism_ring":"maximal order O_K" if conductor==1 else "not fully determined; conductor divides Frobenius-order conductor",
            "endomorphism_ring_conductor":1 if conductor==1 else None,
            "endomorphism_ring_discriminant":dk if conductor==1 else None,
            "endomorphism_ring_class_number":None,
            "endomorphism_ring_class_number_status":"Not certified/computed for this report",
            "twist_label":f'S{m}-a{1-a}' if searched else f'B{m}-a{1-a}',
            "descends_to_F2":False,"field_of_definition_degree_from_j":m,
            "cardinality_method":"PARI ellcard","group_method":"PARI ellgroup"}


def main():
    fields=[irreducible(m,e) for m,e in POLYS.items()]
    curves=[row(a,m) for m in POLYS for a in [0,1]]
    generic,searched,screening=gp_binary_controls()
    data={"schema":"binary elliptic-curve computations v2", "fields":fields,"curves":curves,
          "generic_curves":generic,"searched_examples":searched,"search_metadata":screening,
          "verification":{"trace_method":"Koblitz: S0=2, S1=t0, Sn=t0*S_(n-1)-2*S_(n-2); B-style: PARI ellcard",
             "prime_certificates_file":"binary_certificates.json",
             "embedding_method":"Exact multiplicative order plus pow/minimality checks",
             "group_method":"Koblitz: Smith invariants of multiplication by pi^m-1 on maximal CM ring; B-style: PARI ellgroup"}}
    allcurves=curves+generic+searched
    cert_primes={p for c in allcurves for p in c['group_order_factors']}
    cert_primes|={p for c in generic+searched for p in c['frobenius_discriminant_abs_factors']}
    certs={str(p):prime_certificate(p) for p in sorted(cert_primes)}
    (OUT/'binary_certificates.json').write_text(json.dumps(certs,indent=2)+'\n')
    (OUT/'binary_data.json').write_text(json.dumps(data,indent=2)+'\n')
    for c in allcurves:
        print(c['label'], 't=',c['trace'],'#E=',c['group_order'],'factors=',c['group_order_factors'],
              'r=',c['largest_prime_subgroup'],'h=',c['cofactor'],'k=',c['embedding_degree'])
    print('All fields, recurrence/PARI counts, factorizations, Lucas certificates, and embedding orders verified.')

if __name__=='__main__':
    main()
