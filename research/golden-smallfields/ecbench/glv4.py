"""Verified supplementary GLV+GLS cases and honest variable-time timings."""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import statistics
import time

from .extension import Fp2,ExtensionCurve,babai_decompose
from .validate import certify,integer


def _load(cases=None):
    if cases is None:
        packaged=Path(__file__).resolve().parent/'data'/'glv4.json'
        source=packaged if packaged.is_file() else Path(__file__).resolve().parents[1]/'build'/'glv4.json'
        cases=json.loads(source.read_text())
    elif isinstance(cases,(str,Path)):
        cases=json.loads(Path(cases).read_text())
    if isinstance(cases,dict):
        return list(cases['cases']) if 'cases' in cases else [cases]
    if isinstance(cases,list):return cases
    raise ValueError('Expected a GLV+GLS case, case list, or JSON path')


def _require(condition,message):
    if not condition:raise ValueError('GLV+GLS verification failed: '+message)


def _point(record):
    return None if record is None else (integer(record['x']),integer(record['y']))


def _dot(a,b):return sum(x*y for x,y in zip(a,b))


def _gram(basis):
    stars=[];norms=[];mu=[]
    for row in basis:
        v=list(map(Fraction,row));coefficients=[]
        for star,norm in zip(stars,norms):
            coefficient=_dot(row,star)/norm;coefficients.append(coefficient)
            v=[x-coefficient*y for x,y in zip(v,star)]
        stars.append(v);norms.append(_dot(v,v));mu.append(coefficients)
    _require(all(norm>0 for norm in norms),'singular lattice')
    return stars,norms,mu


def _determinant(matrix):
    a=[list(map(Fraction,row)) for row in matrix];n=len(a);det=Fraction(1)
    _require(all(len(row)==n for row in a),'non-square lattice')
    for col in range(n):
        pivot=next((i for i in range(col,n) if a[i][col]),None)
        if pivot is None:return 0
        if pivot!=col:a[col],a[pivot]=a[pivot],a[col];det=-det
        value=a[col][col];det*=value
        for row in range(col+1,n):
            multiplier=a[row][col]/value
            for j in range(col,n):a[row][j]-=multiplier*a[col][j]
    _require(det.denominator==1,'nonintegral lattice determinant')
    return det.numerator


def _certificates(case):
    records=[]
    for key,proof in case['certificates'].items():
        n=integer(proof['n']);_require(str(n)==key,'certificate key mismatch')
        if proof['method']=='base':
            _require(n==2,'invalid base certificate')
            records.append({'n':str(n),'method':'base','n_minus_1_factorization':[]})
        else:
            _require(proof['method']=='Lucas','unsupported primality certificate')
            factors=[{'prime':row['prime'],'exponent':integer(row['exponent'])}
                     for row in proof['n_minus_one_factorization']]
            witness=proof['witness']
            records.append({'n':str(n),'method':'lucas_complete_n_minus_1',
                            'n_minus_1_factorization':factors,
                            'witnesses':{row['prime']:witness for row in factors}})
    return certify(records)


def _setup(case):
    p=integer(case['field']['base_prime']);d=integer(case['field']['nonresidue_d'])
    f=Fp2(p,d);E=ExtensionCurve(f,integer(case['curve']['a4']),integer(case['curve']['a6']))
    r=integer(case['curve']['subgroup_order']);G=_point(case['curve']['generator'])
    phi_record=case['endomorphisms']['phi'];psi_record=case['endomorphisms']['psi']
    zeta=integer(phi_record['zeta3']);cx=integer(psi_record['x_multiplier']);cy=integer(psi_record['y_multiplier'])
    def phi(P):return None if P is None else (f.mul(zeta,P[0]),P[1])
    def psi(P):return None if P is None else (f.mul(cx,f.frobenius(P[0])),f.mul(cy,f.frobenius(P[1])))
    return f,E,r,G,phi,psi


def _cm_base_count(case,field):
    """Exact j0 trace theorem plus witnesses excluding five wrong traces."""
    p=field.p;t=integer(case['base_curve']['trace'])
    a=integer(case['lattices']['glv4']['cm_a']);c=integer(case['lattices']['glv4']['cm_c'])
    _require(a*a-a*c+c*c==p,'CM norm representation')
    candidates={2*a-c,-2*a+c,a+c,-a-c,a-2*c,-a+2*c}
    _require(len(candidates)==6 and t in candidates,'CM trace candidates')
    _require(p%4==3 and p%3==1,'unsupported base-field CM counting case')
    E0=ExtensionCurve(field,0,integer(case['base_curve']['b']))
    witnesses=[];points=0
    for x in range(1000):
        rhs=(x**3+integer(case['base_curve']['b']))%p
        if rhs and pow(rhs,(p-1)//2,p)!=1:continue
        y=pow(rhs,(p+1)//4,p) if rhs else 0;y=min(y,(-y)%p)
        P=(x,y);points+=1
        for candidate in sorted(list(candidates)):
            if E0.mul(p+1-candidate,P) is not None:
                candidates.remove(candidate);witnesses.append(candidate)
        if len(candidates)==1:break
    _require(candidates=={t},'base curve point count')
    _require(integer(case['base_curve']['order'])==p+1-t,'base point count metadata')
    return {'method':'j0 CM trace theorem and exact scalar-multiplication elimination',
            'base_points_used':points,'incorrect_traces_eliminated':len(witnesses)}


def _verify_case(case):
    primes,proofs=_certificates(case)
    f,E,r,G,phi,psi=_setup(case);p=f.p;q=f.q
    _require(p in primes and r in primes,'uncertified characteristic or subgroup prime')
    _require(integer(case['field']['cardinality'])==q,'extension cardinality')
    _require(integer(case['field']['element_bits'])==q.bit_length(),'extension bit length')
    eta=integer(case['field']['eta']);_require(f.unpack(eta)==tuple(map(integer,case['field']['eta_coefficients'])),'eta coordinates')
    e0,e1=f.unpack(eta)
    _require(integer(case['field']['eta_norm'])==(e0*e0-f.d*e1*e1)%p,'eta norm')
    _require(f.power(eta,(q-1)//2)==p-1,'eta is not a nonsquare')
    b=integer(case['base_curve']['b']);t=integer(case['base_curve']['trace']);N=(p-1)**2+t*t
    _require(t%p!=0 and p%3==1,'ordinarity or j0 split condition')
    _require(E.A==0 and E.B==f.mul(b,f.power(eta,3)),'quadratic twist model')
    _require(f.unpack(E.B)==tuple(map(integer,case['curve']['a6_coefficients'])),'twist coefficient components')
    _require(integer(case['curve']['order'])==N==4*r,'twist point count formula')
    _require(integer(case['curve']['cofactor'])==4,'subgroup cofactor')
    _require(integer(case['curve']['subgroup_bits'])==r.bit_length(),'subgroup bit length')
    factors=case['curve']['order_factorization']
    _require(math.prod(integer(row['prime'])**integer(row['exponent']) for row in factors)==N,'order factorization')
    _require(all(integer(row['prime']) in primes for row in factors),'uncertified order factor')
    count_proof=_cm_base_count(case,f)
    _require(integer(case['curve']['trace'])==q+1-N,'extension trace')
    _require(integer(case['curve']['Frobenius_discriminant'])==(q+1-N)**2-4*q,'Frobenius discriminant')
    _require(integer(case['base_curve']['Frobenius_discriminant'])==t*t-4*p,'base Frobenius discriminant')
    conductor=integer(case['curve']['Frobenius_order_conductor'])
    _require(conductor>0,'Frobenius-order conductor positivity')
    _require((q+1-N)**2-4*q==-3*conductor*conductor,'Frobenius-order conductor')
    zeta=integer(case['endomorphisms']['phi']['zeta3'])
    g=integer(case['endomorphisms']['phi']['base_primitive_root'])
    _require(1<g<p and all(pow(g,(p-1)//integer(row['prime']),p)!=1 for row in proofs[p]['n_minus_1_factorization']),'base primitive root')
    _require(1<zeta<p and pow(zeta,3,p)==1 and zeta==pow(g,(p-1)//3,p),'primitive cube root')
    _require(integer(case['endomorphisms']['phi']['degree'])==1,'phi degree')
    _require(integer(case['endomorphisms']['psi']['degree'])==p,'psi degree')
    _require(integer(case['endomorphisms']['psi']['x_multiplier'])==f.power(eta,1-p),'psi x multiplier')
    _require(integer(case['endomorphisms']['psi']['y_multiplier'])==f.power(eta,3*(1-p)//2),'psi y multiplier')
    torsion=[(f.neg(eta),0),(f.neg(f.mul(eta,zeta)),0),(f.neg(f.mul(eta,zeta*zeta%p)),0)]
    _require(len(set(torsion))==3 and all(E.on_curve(P) and E.add(P,P) is None for P in torsion),'full rational two-torsion')
    _require(list(map(integer,case['curve']['group_invariant_factors']))==[2,2*r],'point-group structure')
    _require(G is not None and E.on_curve(G) and E.mul(r,G) is None,'subgroup generator exact order')
    source=_point(case['curve']['generator_source']);projector=integer(case['curve']['generator_projector'])
    _require(projector==2 and E.mul(projector,source)==G,'subgroup projection')
    order_checks=0
    cm_relation=case['endomorphisms']['psi']['actual_CM_linear_relation']
    cm_a,cm_c=integer(cm_relation['a']),integer(cm_relation['c'])
    _require(cm_a*cm_a-cm_a*cm_c+cm_c*cm_c==p and 2*cm_a-cm_c==t,'actual CM orientation')
    for name,record in case['named_points'].items():
        P=_point(record['point']);order=integer(record['order'])
        _require(E.on_curve(P) and E.mul(order,P) is None,'point membership/order '+name)
        _require(order>0 and N%order==0,'point order divisor '+name)
        for prime in (2,r):
            if order%prime==0:_require(E.mul(order//prime,P) is not None,'non-exact point order '+name)
        if P is None:_require(order==1,'infinity order')
        _require(E.on_curve(phi(P)) and E.on_curve(psi(P)),'endomorphism closure '+name)
        _require(phi(phi(phi(P)))==P and psi(psi(P))==E.neg(P),'endomorphism restrictions '+name)
        _require(phi(psi(P))==psi(phi(P)),'endomorphism commutation '+name)
        _require(E.add(E.mul(cm_a,P),E.mul(cm_c,phi(P)))==psi(P),'actual CM linear relation '+name)
        relation=E.add(E.add(psi(psi(P)),E.mul(-t,psi(P))),E.mul(p,P))
        _require(relation is None,'transported Frobenius relation '+name)
        order_checks+=1
    lp=integer(case['endomorphisms']['phi']['subgroup_eigenvalue']);ls=integer(case['endomorphisms']['psi']['subgroup_eigenvalue'])
    _require((lp*lp+lp+1)%r==0 and lp!=1,'phi eigenvalue polynomial')
    _require(ls==(p-1)*pow(t,-1,r)%r and ls*ls%r==r-1,'psi eigenvalue polynomial')
    _require(E.mul(lp,G)==phi(G) and E.mul(ls,G)==psi(G),'eigenvalue point identities')
    values=[1,lp,ls,lp*ls%r]
    _require(list(map(integer,case['endomorphisms']['four_term_eigenvalues']))==values and len(set(values))==4,'four distinct map eigenvalues')
    k=integer(case['curve']['embedding_degree'])
    _require(k>0 and (r-1)%k==0 and pow(q,k,r)==1,'embedding degree')
    for row in proofs[r]['n_minus_1_factorization']:
        ell=integer(row['prime'])
        if k%ell==0:_require(pow(q,k//ell,r)!=1,'nonminimal embedding degree')
    bases={}
    seed=[list(map(integer,row)) for row in case['lattices']['glv4']['balanced_seed']]
    _require(abs(_determinant(seed))==r and all(_dot(row,values)%r==0 for row in seed),'balanced lattice seed')
    for name,dimension,eigenvalues in [('glv2',2,values[:2]),('glv4',4,values)]:
        basis=[list(map(integer,row)) for row in case['lattices'][name]['basis']]
        _require(len(basis)==dimension and all(len(row)==dimension for row in basis),'lattice dimension')
        determinant=_determinant(basis)
        _require(abs(determinant)==r and determinant==integer(case['lattices'][name]['determinant']),'lattice determinant')
        _require(all(_dot(row,eigenvalues)%r==0 for row in basis),'lattice kernel')
        stars,norms,mu=_gram(basis);delta=Fraction(case['lattices']['LLL_delta'])
        _require(Fraction(1,4)<delta<1,'LLL delta')
        _require(all(abs(value)<=Fraction(1,2) for row in mu for value in row),'LLL size reduction')
        _require(all(norms[i]>=(delta-mu[i][i-1]**2)*norms[i-1] for i in range(1,dimension)),'LLL Lovasz condition')
        bounds=[sum(abs(row[j]) for row in stars)/2 for j in range(dimension)]
        recorded=case['lattices']['Babai_universal_bounds'][name]
        _require([Fraction(integer(v['numerator']),integer(v['denominator'])) for v in recorded['coordinate_absolute_bound_rationals']]==bounds,'Babai rational bounds')
        ceilings=[math.ceil(value) for value in bounds]
        _require(list(map(integer,recorded['coordinate_absolute_bound_ceilings']))==ceilings,'Babai ceiling bounds')
        _require(integer(recorded['guaranteed_component_bits'])==max(value.bit_length() for value in ceilings),'Babai bit bounds')
        bases[name]=basis
    basis_points=[G,phi(G),psi(G),phi(psi(G))]
    for vector in case['decomposition_vectors']:
        scalar=integer(vector['scalar']);expected=_point(vector['expected_point'])
        _require(E.mul(scalar,G)==expected,'ordinary scalar known answer')
        for name,points in [('glv2',basis_points[:2]),('glv4',basis_points)]:
            components=list(map(integer,vector[name+'_components']))
            _require(components==babai_decompose(scalar,r,bases[name]),'Babai decomposition '+name)
            _require(_dot(components,values[:len(components)])%r==scalar%r,'scalar reconstruction '+name)
            _require(E.joint_mul(components,points)==expected,'point reconstruction '+name)
            _require(max(abs(x).bit_length() for x in components)==integer(vector[name+'_max_component_bits']),'component bit count '+name)
    for vector in case['field_vectors']:
        a=integer(vector['input']);operation=vector['operation']
        answer={'square':f.square,'inverse':f.inverse,'frobenius':f.frobenius}[operation](a)
        _require(answer==integer(vector['expected']),'extension field vector '+operation)
    for vector in case['endomorphism_vectors']:
        P=_point(vector['point']);operation=vector['operation']
        functions={'phi':phi,'psi':psi,'phi_psi':lambda P:phi(psi(P)),
                   'psi_phi':lambda P:psi(phi(P)),'phi_cubed':lambda P:phi(phi(phi(P))),
                   'psi_squared':lambda P:psi(psi(P))}
        _require(functions[operation](P)==_point(vector['expected_point']),'endomorphism known answer '+operation)
    _require(integer(case['curve']['endomorphism_discriminant'])==-3,'known j0 CM discriminant metadata')
    return {'case_id':case['supplement_id'],'status':'pass','base_prime_bits':p.bit_length(),
            'field_bits':q.bit_length(),'subgroup_bits':r.bit_length(),'primality_certificates':len(primes),
            'named_point_orders':order_checks,'field_vectors':len(case['field_vectors']),
            'endomorphism_vectors':len(case['endomorphism_vectors']),
            'decomposition_vectors':len(case['decomposition_vectors']),
            'point_count_proof':count_proof,'point_group_proof':'twist point count, certified factorization, and full rational two-torsion',
            'endomorphism_proof_boundary':'special ordinary j=0 model and explicit cubic automorphism; no general endomorphism-ring algorithm',
            'guaranteed_glv4_component_bits':integer(case['lattices']['Babai_universal_bounds']['glv4']['guaranteed_component_bits']),
            'guaranteed_glv2_component_bits':integer(case['lattices']['Babai_universal_bounds']['glv2']['guaranteed_component_bits']),
            'source_pari_checks':integer(case['verification'].get('PARI_independent_checks','0')),
            'source_pari_recomputed':False}


def verify_cases(cases=None):
    """Check supplementary cases with the independent standard-library oracle."""
    start=time.perf_counter_ns();results=[_verify_case(case) for case in _load(cases)]
    return {'status':'pass','cases':results,'case_count':len(results),
            'elapsed_ms':(time.perf_counter_ns()-start)/1_000_000,
            'primary_dataset_modified':False}


def benchmark_cases(cases=None,samples=3,iterations=1):
    """Time full ordinary, GLV2, and GLV+GLS4 paths, including decomposition/maps.

    Timings compare this variable-time affine Python oracle, not optimized
    native cryptography. Component widths are not reported as speedups.
    """
    if type(samples) is not int or samples<=0 or type(iterations) is not int or iterations<=0:
        raise ValueError('samples and iterations must be positive integers')
    loaded=_load(cases);preflight=verify_cases(loaded);results=[]
    for case in loaded:
        f,E,r,G,phi,psi=_setup(case)
        bases={name:[list(map(integer,row)) for row in case['lattices'][name]['basis']] for name in ['glv2','glv4']}
        records=[]
        for index in range(samples):
            digest=int.from_bytes(hashlib.sha256(f"ecbench-glv-gls-bench-{index}".encode()).digest(),'big')
            scalar=r//2+digest%(r-r//2)
            timings={'ordinary_total_ns':0}
            for name in ['glv2','glv4']:
                for phase in ['decomposition_ns','maps_ns','joint_multiplication_ns','total_ns']:timings[name+'_'+phase]=0
            components_by_name={}
            for _ in range(iterations):
                start=time.perf_counter_ns();expected=E.mul(scalar,G)
                timings['ordinary_total_ns']+=time.perf_counter_ns()-start
                for name in ['glv2','glv4']:
                    start=time.perf_counter_ns();components=babai_decompose(scalar,r,bases[name]);decomposed=time.perf_counter_ns()
                    if name=='glv2':points=[G,phi(G)]
                    else:
                        ph=phi(G);ps=psi(G);points=[G,ph,ps,phi(ps)]
                    mapped=time.perf_counter_ns();answer=E.joint_mul(components,points);finished=time.perf_counter_ns()
                    _require(answer==expected,'benchmark output mismatch')
                    timings[name+'_decomposition_ns']+=decomposed-start
                    timings[name+'_maps_ns']+=mapped-decomposed
                    timings[name+'_joint_multiplication_ns']+=finished-mapped
                    timings[name+'_total_ns']+=finished-start
                    components_by_name[name]=components
            averaged={key:value/iterations for key,value in timings.items()}
            records.append({'scalar':str(scalar),'scalar_bits':scalar.bit_length(),'timings_ns_per_iteration':averaged,
                            'glv2_component_bits':max(abs(x).bit_length() for x in components_by_name['glv2']),
                            'glv4_component_bits':max(abs(x).bit_length() for x in components_by_name['glv4']),
                            'outputs_equal':True})
        medians={name:statistics.median(record['timings_ns_per_iteration'][name] for record in records)
                 for name in records[0]['timings_ns_per_iteration']}
        results.append({'case_id':case['supplement_id'],'samples':records,
                        'median_timings_ns_per_iteration':medians,
                        'measured_total_ratio_to_ordinary':{name:medians[name+'_total_ns']/medians['ordinary_total_ns'] for name in ['glv2','glv4']}})
    return {'status':'pass','backend':'independent variable-time affine Python Fp2 oracle',
            'samples_per_case':samples,'iterations_per_sample':iterations,'cases':results,
            'includes':['exact rational Babai decomposition','endomorphism map evaluation','signed interleaved binary multiplication'],
            'verification_preflight':preflight,'timing_limitations':'Python overhead and affine inversions are included; these measurements do not predict native implementation speed.'}
