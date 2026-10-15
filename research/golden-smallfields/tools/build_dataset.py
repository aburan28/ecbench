#!/usr/bin/env python3
"""Deterministically assemble the golden dataset from independently generated answers."""
import collections
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ecbench.reference import Field,Curve

DATA=ROOT/'ecbench'/'data'
SOURCE=ROOT/'provenance'
REPORT=ROOT.parent/'curve_report'

def read(path):return json.loads(Path(path).read_text())
def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,sort_keys=True,indent=2,ensure_ascii=True)+'\n')
def factors(f):return [{'prime':str(p),'exponent':int(e)} for p,e in sorted(f.items(),key=lambda x:int(x[0]))]
def decimalize(value):
    if type(value) is int:return str(value)
    if isinstance(value,list):return [decimalize(x) for x in value]
    if isinstance(value,dict):return {k:decimalize(v) for k,v in value.items()}
    return value

def main():
    SOURCE.mkdir(parents=True,exist_ok=True)
    for name in ['normalized_data.json','prime_data.json','binary_data.json','prime_certificates.json','binary_certificates.json',
                 'prime_compute.py','binary_compute.py','normalize.py','verify_report.py']:
        dest=SOURCE/name
        if not dest.exists():shutil.copyfile(REPORT/name,dest)
    normalized=read(SOURCE/'normalized_data.json')
    original_prime=read(SOURCE/'prime_data.json');original_binary=read(SOURCE/'binary_data.json')
    prime_raw=read(ROOT/'build'/'prime_vectors.json');binary_raw=read(ROOT/'build'/'binary_vectors.json')
    rawprime={c['id']:c for c in prime_raw['curves']}
    proof_raw=read(SOURCE/'prime_certificates.json')
    def collect(cert):
        n=str(cert['n'])
        if n not in proof_raw:
            witnesses=cert.get('witness_for_each_prime_factor',{})
            proof_raw[n]={'n':cert['n'],'method':'Lucas separate witnesses','n_minus_1_factors':cert.get('n_minus_1_factors',{}),'witnesses':witnesses}
        for child in cert.get('prime_factor_certificates',{}).values():collect(child)
    for cert in read(SOURCE/'binary_certificates.json').values():collect(cert)
    proofs=[]
    for n,c in sorted(proof_raw.items(),key=lambda x:int(x[0])):
        proofs.append({'n':n,'method':'base' if n=='2' else 'lucas_complete_n_minus_1',
          'n_minus_1_factorization':factors(c.get('n_minus_1_factors',{})),
          'witnesses':({p:str(c['a']) for p in c['n_minus_1_factors']} if 'a' in c else {p:str(a) for p,a in c.get('witnesses',{}).items()})})
    proof_by_n={p['n']:p for p in proofs}

    fields=[]
    for f in normalized['fields']:
        binary=f['kind']=='binary';bits=f['bits'];q=f['q']
        modulus=sum(1<<e for e in f['exponents']) if binary else q
        fields.append({'id':f['id'],'kind':f['kind'],'characteristic':str(2 if binary else q),'degree':bits if binary else 1,
          'cardinality':str(q),'element_bits':bits,'coordinate_bytes':(bits+7)//8,
          'representation':{'basis':'polynomial' if binary else 'prime_residue','modulus':hex(modulus),'modulus_integer':str(modulus),
           **({'polynomial_exponents':f['exponents'],'basis_generator':'0x'+format(2,f'0{2*((bits+7)//8)}x')} if binary else {'form':f['form'],'multiplicative_primitive_root':str(f['g'])})},
          'proof':{'method':'rabin_irreducibility' if binary else 'lucas_primality','certificate_n':None if binary else str(q)}})
    fobjects={f['id']:Field(f) for f in fields}
    psources={}
    for f in original_prime['fields']:
        for family in f['families']:
            for row in family['twists']:psources[f"{f['id']}-{family['family']}-{row['index']}"]=row
    for row in original_prime['searched_example']['twists']:psources[f"p83-S-{row['index']}"]=row
    bsources={row['label']:row for key in ['curves','generic_curves','searched_examples'] for row in original_binary[key]}

    curves=[]
    for c in normalized['curves']:
        f=fobjects[c['field']];binary=f.kind=='binary';raw=binary_raw['curves'][c['id']] if binary else rawprime[c['id']]
        def point(p):return None if p is None else {'x':f.hex(p[0]),'y':f.hex(p[1])}
        if binary:
            coefficients={k:f.hex(v) for k,v in raw['coefficients'].items()}
            named={k:raw[k] for k in ['P','Q','R']};named.update(G=raw['subgroup']['generator'],T=raw['two_torsion_point'],O=None)
            point_orders=raw['point_orders'];source=raw['subgroup']['source_point'];projector=raw['subgroup']['projection_scalar'];scan=raw['subgroup'];row=bsources[c['id']]
        else:
            coefficients={k:f.hex(v) for k,v in {'a1':0,'a2':0,'a3':0,'a4':c['A'],'a6':c['B']}.items()}
            named=raw['points'];point_orders=raw['point_orders'];source=raw['deterministic_scan']['G_source'];projector=raw['deterministic_scan']['projector'];scan=raw['deterministic_scan'];row=psources[c['id']]
        dk=c['DK'];de=c['Dend'];fe=math.isqrt(de//dk)
        if abs(dk) in (3,4,7):dkfacts={str(2 if abs(dk)==4 else abs(dk)):2 if abs(dk)==4 else 1}
        elif not c['ordinary']:dkfacts={str(c['q']):1}
        else:
            df=row.get('Frobenius_discriminant_absolute_factors',row.get('frobenius_discriminant_abs_factors'))
            dkfacts={str(p):int(e) for p,e in df.items()}
            conductor=c['fpi']
            for p in list(dkfacts):
                count=0
                while conductor%int(p)==0:count+=1;conductor//=int(p)
                dkfacts[p]-=2*count
                if dkfacts[p]==0:del dkfacts[p]
            assert conductor==1
        ring={'type':'imaginary_quadratic_order','field_discriminant':str(dk),'field_discriminant_factorization':factors(dkfacts),
          'order_conductor':str(fe),'order_discriminant':str(de),'maximal':fe==1,
          'class_number':str(1) if c['ordinary'] and abs(de) in (3,4,7) else None,
          'class_number_status':'known_exact' if c['ordinary'] and abs(de) in (3,4,7) else 'not_computed'}
        geometric=ring if c['ordinary'] else {'type':'maximal_quaternion_order','ramified_prime':str(c['q']),'ramified_at_infinity':True,
           'explicit_order_basis':None,'identification_status':'algebra_and_maximality_only'}
        fpi={'conductor':str(c['fpi']),'discriminant':str(c['Dpi']),
             'class_number':str(c['frob_classno']) if c.get('frob_classno') is not None else None}
        embeddings=[]
        for r in sorted(map(int,c['factors'])):
            if r==int(f.record['characteristic']):embeddings.append({'prime':str(r),'degree':None,'status':'characteristic_factor'});continue
            k=r-1
            for factor in proof_by_n[str(r)]['n_minus_1_factorization']:
                ell=int(factor['prime'])
                while k%ell==0 and pow(c['q'],k//ell,r)==1:k//=ell
            minimal=[v['prime'] for v in proof_by_n[str(r)]['n_minus_1_factorization'] if k%int(v['prime'])==0]
            embeddings.append({'prime':str(r),'degree':str(k),'status':'exact_multiplicative_order','minimality_prime_divisors':minimal})
        def factor_divisor(n):
            result={}
            for p in c['factors']:
                count=0
                while n%int(p)==0:count+=1;n//=int(p)
                if count:result[p]=count
            assert n==1
            return factors(result)
        rank=2 if c['group'][0]%c['r']==0 else 1
        curve={'id':c['id'],'field_id':c['field'],'family':c['family'],'model':'ordinary_binary' if binary else 'short_weierstrass',
          'coefficients':coefficients,'invariants':{'delta':f.hex(raw['delta'] if binary else c['delta']),
          'j':f.hex(raw['j'] if binary else c['j']),'trace':str(c['t']),'order':str(c['N']),
          'order_factorization':factors(c['factors']),'group_invariant_factors':[str(d) for d in c['group']],
          'group_exponent':str(c['group'][1]),'group_cyclic':c['group'][0]==1,'ordinary':c['ordinary'],
          'frobenius_polynomial':['1',str(-c['t']),str(c['q'])],'frobenius_discriminant':str(c['Dpi']),
          'endomorphism':{'geometric':geometric,'base_field':ring,'frobenius_order':fpi,
             'proof_method':c.get('endo_proof','Degree-two Koblitz Frobenius has fundamental -7 discriminant.' if c['family']=='K' else 'Fundamental Frobenius discriminant forces maximality.')},
          'subgroup_embedding_degrees':embeddings},
          'subgroup':{'order':str(c['r']),'cofactor':str(c['h']),'generator':point(raw['subgroup']['generator'] if binary else raw['subgroup_generator']),
            'embedding_degree':str(c['k']),'projection_scalar':str(projector),'projection_formula':'group_exponent / r',
            'source_point':point(source),'r_torsion_rank':rank,'unique_order_r_subgroup':rank==1,
            'cofactor_coprime_to_r':math.gcd(c['h'],c['r'])==1,
            'membership_rule':'on_curve_and_r_annihilation' if rank==1 else 'exact_chosen_cyclic_subgroup_bsgs',
            'generator_selection':decimalize(scan)},
          'named_points':{name:{'point':point(p),'order':str(point_orders[name]),'order_factorization':factor_divisor(point_orders[name])} for name,p in named.items()},
          'twist':{'partner_id':c['twist_id'],'kind':'quadratic','parameter':f.hex(1 if binary else (7 if c['family']=='S' else c['g']))},
          'provenance':{'source_label':c['id'],'computation_engine':'PARI/GP 2.17.2','scope':'experimental_small_field_benchmark',
              'legacy_source_file':'binary_data.json' if binary else 'prime_data.json'}}
        curves.append(curve)
    curve_by_id={c['id']:c for c in curves};cobjects={c['id']:Curve(c,fobjects[c['field_id']]) for c in curves}
    for c in curves:
        f=fobjects[c['field_id']]
        twist=c['twist']
        if f.kind=='binary':
            twist['partner_relation']='literal_standard_quadratic_twist'
            continue
        partner=curve_by_id[twist['partner_id']]
        d=f.parse(twist['parameter'])
        if c['family'] in ('E4','E3'):
            power=4 if c['family']=='E4' else 6
            shift=2 if power==4 else 3
            i=int(c['id'].rsplit('-',1)[1]);j=int(partner['id'].rsplit('-',1)[1])
            common=math.gcd(power,f.q-1);modulus=(f.q-1)//common
            exponent=((i+shift-j)//common)*pow(power//common,-1,modulus)%modulus
            scale=pow(d,exponent,f.q)
        else:scale=1 if c['id'].endswith('-0') else d
        twist['partner_relation']='base_field_isomorphic_representative_of_standard_twist'
        twist['representative_isomorphism_scale']=f.hex(scale)
        twist['representative_point_map']='direct twist (x,y) -> partner (x/scale^2,y/scale^3); infinity maps to infinity'

    vectors=[];counts=collections.Counter();seen_field=set()
    def add(domain,operation,scope,inputs,expected,purpose='known_answer'):
        index=counts[(domain,scope)];counts[(domain,scope)]+=1
        vectors.append({'id':f'{domain}:{scope}:{index:04d}','domain':domain,'operation':operation,
          ('field_id' if domain=='field' else 'curve_id'):scope,'inputs':inputs,'expected':expected,'purpose':purpose})
    def field_vector(fid,op,args,output=None,error=None):
        f=fobjects[fid];op={'subtract':'sub','multiply':'mul','divide':'div','reduce_modulus':'reduce','frobenius_period':'frobenius'}.get(op,op)
        if op=='reduce':inputs={'value':str(args[0])}
        elif op in ('pow','frobenius'):inputs={'a':f.hex(args[0]),'exponent':str(args[1])}
        else:inputs={'a':f.hex(args[0]),**({'b':f.hex(args[1])} if len(args)>1 else {})}
        expected={'error':error} if error else ({'integer':str(output)} if op=='trace' else {'element':f.hex(output)})
        key=json.dumps([fid,op,inputs,expected],sort_keys=True)
        if key not in seen_field:seen_field.add(key);add('field',op,fid,inputs,expected)
    for raw in prime_raw['curves']:
        for v in raw['field_vectors']:field_vector(raw['field_id'],v['op'],v['args'],v.get('result'),v.get('error'))
    for fid,raw in binary_raw['fields'].items():
        for v in raw['field_vectors']:field_vector(fid,v['operation'],v['inputs'],v.get('output'),v.get('error'))
    for fid,f in fobjects.items():
        for op,args,out in [('add',[0,1],1),('sub',[1,1],0),('mul',[0,f.q-1],0),('square',[0],0),('inverse',[1],1)]:field_vector(fid,op,args,out)
        field_vector(fid,'inverse',[0],error='inverse_zero')
        for value in (-1,f.q):add('field','inverse',fid,{'raw_a':str(value)},{'error':'noncanonical_element'},'input_rejection')
        add('field','mul',fid,{'a':f.hex(1),'raw_b':'-1'},{'error':'noncanonical_element'},'input_rejection')
        add('field','canonical',fid,{'value':str(f.q)},{'boolean':False},'input_rejection')
        if f.kind=='binary':add('field','reduce',fid,{'value':'-1'},{'error':'negative_polynomial'},'input_rejection')

    def point_for(cid,p):return cobjects[cid].point(None if p is None else tuple(p))
    def curve_vector(cid,op,args,out):
        inputs={'scalar':str(args[0]),'P':point_for(cid,args[1])} if op=='mul' else {'P':point_for(cid,args[0]),**({'Q':point_for(cid,args[1])} if op=='add' else {})}
        add('curve',op,cid,inputs,{'point':point_for(cid,out)})
    def encoding_pair(cid,p,compressed,uncompressed):
        pp=point_for(cid,p)
        for mode,encoded in [('compressed',compressed),('uncompressed',uncompressed)]:
            add('encoding','encode',cid,{'P':pp,'format':mode},{'octets':encoded})
            add('encoding','decode',cid,{'octets':encoded},{'point':pp})
    def decode_error(cid,encoded,error):
        if error=='invalid_encoding':error='invalid_length' if not encoded or int(encoded[:2],16) in (0,2,3,4) else 'invalid_prefix'
        error={'noncanonical_x0_sign':'invalid_y_selector','no_quadratic_solution':'point_not_on_curve'}.get(error,error)
        add('encoding','decode',cid,{'octets':encoded},{'error':error},'input_rejection')
    for raw in prime_raw['curves']:
        cid=raw['id']
        for v in raw['curve_vectors']:curve_vector(cid,v['op'],v['args'],v['result'])
        for v in raw['encoding_vectors']:encoding_pair(cid,v['point'],v['sec1_compressed_hex'],v['sec1_uncompressed_hex'])
        for v in raw['invalid_vectors']:
            if v['op']=='decode_point':decode_error(cid,v['encoding_hex'],v['error'])
            else:add('membership','validate_point',cid,{'raw_point':None if v['point'] is None else {'x':str(v['point'][0]),'y':str(v['point'][1])}},{'boolean':v['result']},'input_rejection')
        for v in raw['subgroup_validation_vectors']:
            add('membership',v['op'],cid,{'raw_point':None if v['point'] is None else {'x':str(v['point'][0]),'y':str(v['point'][1])}},{'boolean':v['result']},'subgroup_conformance')
    for cid,raw in binary_raw['curves'].items():
        for v in raw['curve_vectors']+raw['subgroup_vectors']:curve_vector(cid,v['operation'],v['inputs'],v['output'])
        for v in raw['sec1_vectors']['valid']:encoding_pair(cid,v['point'],v['compressed'],v['uncompressed'])
        for v in raw['sec1_vectors']['invalid']:decode_error(cid,v['encoding'],v['error'])
        c=cobjects[cid]
        for p,op in [(None,'subgroup_contains'),(None,'public_key_valid'),(c.G,'subgroup_contains'),(c.G,'public_key_valid'),(tuple(raw['P']),'is_r_torsion'),(tuple(raw['T']) if 'T' in raw else tuple(raw['two_torsion_point']),'public_key_valid')]:
            result={'subgroup_contains':c.subgroup_contains,'public_key_valid':c.public_key_valid,'is_r_torsion':c.r_torsion}[op](p)
            add('membership',op,cid,{'P':c.point(p)},{'boolean':result},'subgroup_conformance')

    for filename in ['prime_maps.json','binary_maps.json']:
        rawmaps=read(ROOT/'build'/filename)
        for v in (rawmaps['vectors'] if isinstance(rawmaps,dict) else rawmaps):
            cid=v['curve_id'];f=cobjects[cid].field
            inputs={'P':point_for(cid,v['point']),'target_coefficients':{k:f.hex(a) for k,a in v['target_coefficients'].items()}}
            if v['parameter'] is not None:inputs['parameter']=f.hex(v['parameter'])
            add('map',v['operation'],cid,inputs,{'point':point_for(cid,v['result'])},'explicit_map_known_answer')

    relations={'isogeny_classes':[],'geometric_classes_over_base_field':[],'ordinary_endomorphism_order_classes':[],'twist_pairs':[]}
    def grouping(key):
        groups=collections.defaultdict(list)
        for c in curves:groups[key(c)].append(c['id'])
        return groups
    for (fid,t),ids in grouping(lambda c:(c['field_id'],c['invariants']['trace'])).items():
        relations['isogeny_classes'].append({'field_id':fid,'trace':t,'curve_ids':ids,'evidence':'Tate trace criterion','explicit_maps_computed':False})
    for (fid,j),ids in grouping(lambda c:(c['field_id'],c['invariants']['j'])).items():relations['geometric_classes_over_base_field'].append({'field_id':fid,'j':j,'curve_ids':ids})
    groups=collections.defaultdict(list)
    for c in curves:
        if c['invariants']['ordinary']:
            ring=c['invariants']['endomorphism']['geometric'];groups[(ring['field_discriminant'],ring['order_conductor'])].append(c['id'])
    for (dk,fe),ids in groups.items():relations['ordinary_endomorphism_order_classes'].append({'field_discriminant':dk,'conductor':fe,'curve_ids':ids,'relation':'abstract_order_isomorphism'})
    seen=set()
    for c in curves:
        pair=tuple(sorted([c['id'],c['twist']['partner_id']]))
        if pair in seen:continue
        seen.add(pair);relations['twist_pairs'].append({'curve_ids':list(pair),'field_id':c['field_id'],'may_be_self_isomorphic':len(set(pair))==1,'evidence':'quadratic twist model and trace reversal'})
    relations['reference_maps']=[{'type':'separable_degree_two_isogeny','characteristic_condition':'p > 3','source_equation':'y^2=x^3-x','target_equation':'Y^2=X^3+4X',
      'map':{'X':'x-1/x','Y':'y*(1+1/x^2)'},'kernel':['infinity','(0,0)'],'exceptional_points_map_to':'infinity','target_is_auxiliary_model':True},
      {'type':'relative_frobenius_isogeny','characteristic_condition':'2','degree':'2','source_equation':'y^2+xy=x^3+a*x^2+b (a in F2)',
       'target_equation':'y^2+xy=x^3+a*x^2+b^2','map':{'X':'x^2','Y':'y^2'},'target_is_auxiliary_model':True},
      {'type':'gaussian_automorphism','characteristic_condition':'p > 3, p = 1 mod 4','source_equation':'y^2=x^3+A*x',
       'map':{'X':'-x','Y':'i*y; i^2=-1'},'target_is_auxiliary_model':False},
      {'type':'eisenstein_automorphism','characteristic_condition':'p > 3, p = 1 mod 3','source_equation':'y^2=x^3+B',
       'map':{'X':'zeta*x; zeta^3=1, zeta!=1','Y':'y'},'target_is_auxiliary_model':False}]

    payloads={'fields.json':fields,'curves.json':curves,'vectors.json':vectors,'relations.json':relations,'certificates.json':proofs}
    for name,obj in payloads.items():dump(DATA/name,obj)
    manifest={'dataset_id':'ecbench-golden-smallfields-v1','schema_version':'1.0.0','release_version':'1.0.0','created_date':'2026-10-12',
      'scope':'56 verified example curves over 7 fields; complete benchmark records, not all elliptic curves',
      'integer_encoding':'canonical decimal strings','field_element_encoding':'0x + fixed-width lowercase big-endian hex','point_at_infinity':None,
      'sec1_octet_encoding':'lowercase hex without 0x; compressed/uncompressed/infinity only',
      'counts':{'fields':len(fields),'curves':len(curves),'vectors':len(vectors),'primality_certificates':len(proofs),'named_points':sum(len(c['named_points']) for c in curves),
        'vectors_by_domain':dict(collections.Counter(v['domain'] for v in vectors)),'represented_isogeny_classes':len(relations['isogeny_classes'])},
      'files':{name:hashlib.sha256((DATA/name).read_bytes()).hexdigest() for name in payloads},
      'record_hashes':{c['id']:hashlib.sha256(json.dumps(c,sort_keys=True,separators=(',',':')).encode()).hexdigest() for c in curves},
      'proof_boundaries':['Actual ordinary rings and supersingular rational quadratic rings determined; explicit quaternion order bases not included.',
         'Default verifier checks certificates and all vectors; --with-pari independently recomputes point counts and group structure.',
         'Endomorphism-order metadata is source-proven; the core verifier checks discriminant consistency, not every model-specific ring proof.',
         'Large generic class numbers are not computed; retained class-number evidence is distinguished from core consistency checks.',
         'R generates the whole rational group only for cyclic rows; no full basis claimed for noncyclic groups.']}
    dump(DATA/'manifest.json',manifest)
    bundle={'manifest':manifest,**{name.removesuffix('.json'):obj for name,obj in payloads.items()}}
    dump(ROOT/'golden_ecbench.json',bundle)
    with (ROOT/'golden_ecbench.jsonl').open('w') as handle:
        for kind,rows in [('manifest',[manifest]),('field',fields),('curve',curves),('vector',vectors),('certificate',proofs),('relations',[relations])]:
            for row in rows:handle.write(json.dumps({'record_type':kind,'record':row},sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps(manifest['counts'],sort_keys=True))

if __name__=='__main__':main()
