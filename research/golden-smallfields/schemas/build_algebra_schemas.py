"""Strict contracts for algebra fixtures/results; does not modify primary data."""
from pathlib import Path
import copy
import json
import re
import jsonschema

ROOT = Path(__file__).resolve().parent.parent
BASE = json.loads((ROOT/'schemas/dataset.schema.json').read_text())['$defs']
CATALOGS = [json.loads((ROOT/'build'/name).read_text()) for name in ('algebra_cases.json','algebra_coverage_cases.json')]
RESULTS = [json.loads((ROOT/'build'/name).read_text()) for name in ('algebra_buchberger.json','algebra_f5b.json','algebra_coverage_buchberger.json','algebra_coverage_f5b.json')]

def ref(name): return {'$ref':'#/$defs/'+name}
def obj(props,required=None): return {'type':'object','properties':props,'required':list(props) if required is None else required,'additionalProperties':False}
def arr(items,minimum=0,maximum=None):
    schema={'type':'array','items':items,'minItems':minimum}
    if maximum is not None:schema['maxItems']=maximum
    return schema
def nullable(schema): return {'anyOf':[{'type':'null'},schema]}
def bounded(minimum=0,maximum=9007199254740991): return {'type':'integer','minimum':minimum,'maximum':maximum}
def number():return {'type':'number','minimum':0}
def dynamic(pattern,value):return {'type':'object','patternProperties':{pattern:value},'additionalProperties':False}
def text():return {'type':'string','minLength':1}

def infer(values,path,defs):
    key=path.rsplit('.',1)[-1]
    if key in ('curve','field'):return ref(key)
    if key in ('point','generator','P','Q') and all(v is None or isinstance(v,dict) and set(v)=={'x','y'} for v in values):return ref('point')
    nonnull=[v for v in values if v is not None]
    if len(nonnull)!=len(values):return {'type':'null'} if not nonnull else nullable(infer(nonnull,path,defs))
    if all(isinstance(v,dict) for v in values):
        if key=='points_by_x':return dynamic(r'^0x(?:[0-9a-f]{2})+$',arr(ref('finitePoint'),1,2))
        if key in ('x_encoding_masks','positive_scalar_by_x'):return dynamic(r'^0x(?:[0-9a-f]{2})+$',ref('unsigned'))
        if key=='domains':return dynamic(r'^[A-Za-z][A-Za-z0-9_]*$',arr(ref('unsigned'),1))
        if key=='source_hashes':return dynamic(r'^[A-Za-z0-9_.-]+$',ref('hash'))
        if key=='stage_times':return dynamic(r'^[a-z][a-z0-9_]*_seconds$',nullable(number()))
        if key=='cases_by_field':return dynamic(r'^(?:p61|p64|p83s?|b53|b79|b83)$',bounded())
        if key=='cohorts':return dynamic(r'^[A-Za-z][A-Za-z0-9_-]*$',ref('aggregate'))
        keys=sorted(set().union(*(set(v) for v in values)))
        common=set(values[0]).intersection(*(set(v) for v in values[1:]))
        return obj({k:infer([v[k] for v in values if k in v],path+'.'+k,defs) for k in keys},sorted(common))
    if all(isinstance(v,list) for v in values):
        if key in ('algebraic_assignments','root_assignments','successful_root_assignments'):return arr(arr(ref('unsigned'),1))
        if key=='ordered_relation_points':return arr(arr(ref('finitePoint'),2,2))
        children=[child for v in values for child in v]
        item=infer(children,path+'[]',defs) if children else {}
        result=arr(item)
        if key=='exponents':result=arr(bounded(0,4096),1,4096)
        if key in ('variables','variable_order'):result=arr({'type':'string','pattern':r'^[A-Za-z][A-Za-z0-9_]*$'},1,4096);result['uniqueItems']=True
        return result
    if all(type(v) is bool for v in values):return {'type':'boolean'}
    if all(type(v) is int for v in values):
        return bounded(-4096,4096) if key=='worker_exit_code' else bounded(0,4096) if key in ('coordinate','total_degree','max_total_degree','value','dimension','arity','summation_polynomial_arity') else bounded()
    if all(type(v) in (int,float) for v in values):
        result=number()
        if 'fraction' in key:result['maximum']=1
        return result
    if all(isinstance(v,str) for v in values):
        if key=='coefficient':return ref('unsigned')
        if key in ('curve_id',):return ref('curveId')
        if key in ('field_id',):return ref('fieldId')
        if key in ('case_sha256','hash_digest_hex','original_55_catalog_sha256'):return ref('hash')
        if all(v.startswith('0x') for v in values):return ref('element')
        if all(re.fullmatch(r'(?:0|-?[1-9][0-9]*)',v) for v in values):return ref('integer')
        if key in ('method','monomial_order','order','role','kind','basis_kind','encoding','model_variant'):
            return {'type':'string','enum':sorted(set(values))}
        return text()
    raise TypeError((path,{type(v).__name__ for v in values}))

def prune(defs,root):
    visited=set()
    def scan(value):
        if isinstance(value,dict):
            if '$ref' in value:
                name=value['$ref'].split('/')[-1]
                if name not in visited:visited.add(name);scan(defs[name])
            for k,v in value.items():
                if k!='$ref':scan(v)
        elif isinstance(value,list):
            for v in value:scan(v)
    scan(root)
    return {name:defs[name] for name in sorted(visited)}

def write(name,defs,root,examples,title):
    schema={'$schema':'https://json-schema.org/draft/2020-12/schema','$id':'urn:ecbench:'+name+':v1','title':title,
            'description':'Strict structural contract; mathematical equivalence, variable lengths, finite-field ranges, counts, root completeness, and timing claims require semantic verification.',
            '$defs':prune(defs,root),**root}
    jsonschema.Draft202012Validator.check_schema(schema)
    checker=jsonschema.Draft202012Validator(schema,format_checker=jsonschema.FormatChecker())
    for example in examples:
        errors=list(checker.iter_errors(example))
        if errors:
            for error in errors[:3]:print(name,list(error.absolute_path),error.message[:250])
            raise SystemExit(f'{name}: {len(errors)} errors')
    (ROOT/'schemas'/name).write_text(json.dumps(schema,indent=2,sort_keys=True)+'\n')
    print('Validated and wrote',name)

# Raw fixture contract, independent of the number of cases in a catalog.
cd=copy.deepcopy(BASE)
cd['case']=infer([case for catalog in CATALOGS for case in catalog['cases']],'cases[]',cd)
cd['case']['properties']['arity']={'const':2}
cd['case']['properties']['summation_polynomial_arity']={'const':3}
cd['case']['properties']['coefficient_modulus']=ref('positive')
cd['case']['properties']['expected']['properties']['algebraic_assignments']=arr(arr(ref('unsigned'),1))
cd['case']['properties']['target']['properties']['point']=ref('finitePoint')
cd['case']['properties']['target']['properties']['scalar']=nullable(ref('positive'))
cd['case']['properties']['polynomials']['items']['properties']['terms']['items']['properties']['coefficient']=ref('unsigned')
cd['caseCollection']=infer(CATALOGS,'case_collection',cd)
cd['caseCollection']['properties']['schema']={'const':'ecbench-semaev-two-point-cases-v1'}
cd['caseCollection']['properties']['cases']=arr(ref('case'),1)
cd['caseCollection']['properties']['primary_dataset_modified']={'const':False}
cd['caseCollection']['properties']['original_55_catalog_modified']={'const':False}
cd['caseCollection']['properties']['relation_arity']={'const':2}
cd['caseCollection']['properties']['summation_polynomial_arity']={'const':3}
write('algebra_cases.schema.json',cd,ref('caseCollection'),CATALOGS,'ECBench two-point Semaev fixture catalog')

# Result contract includes completed, censored, partial-diagnostic, and error rows.
rd=copy.deepcopy(BASE)
rows=[row for result in RESULTS for row in result['results']]
aggregates=[result['summary'] for result in RESULTS]+[v for result in RESULTS for v in result['cohorts'].values()]+[v for result in RESULTS for v in result['stratified_cohorts']]
rd['aggregate']=infer(aggregates,'aggregate',rd)
rd['aggregate']['properties'].pop('cohort',None)
rd['aggregate']['required']=[name for name in rd['aggregate']['required'] if name!='cohort']
rd['caseResult']=infer(rows,'result',rd)
rowprops=rd['caseResult']['properties']
rowprops['status']={'enum':['ok','timeout','error']}
rowprops['field_id']=nullable(ref('fieldId'))
rowprops['curve_id']=nullable(ref('curveId'))
rowprops['input']=nullable(rowprops['input'])
rowprops['groebner_basis']=nullable(rowprops['groebner_basis'])
rowprops['relations']=nullable(rowprops['relations'])
rowprops['error']=text()
rowprops['diagnostic_status']={'enum':['timeout','error']}
rowprops['backend']['properties']['characteristic']=ref('positive')
rowprops['backend']['required']=[k for k in rowprops['backend']['required'] if k!='characteristic']
# Backend peak-degree tracing is explicitly unavailable for these SymPy runs.
rowprops['backend']['properties']['solving_degree']={'type':'null'}
rowprops['backend']['properties']['method']={'enum':['buchberger','f5b']}
rowprops['backend']['properties']['monomial_order']={'enum':['lex','grlex','grevlex']}
rd['caseResult']['required']=[k for k in rd['caseResult']['required'] if k not in ('input','relations','groebner_basis','error','diagnostic_status')]
# A completed solve must retain its actual basis, roots, and verification stage.
rd['caseResult']['allOf']=[{'if':{'properties':{'status':{'const':'ok'}},'required':['status']},'then':{'required':['input','relations','groebner_basis'],'properties':{'input':{'type':'object'},'relations':{'type':'object'},'groebner_basis':{'type':'object'}}}}]
rd['resultCollection']=infer(RESULTS,'result_collection',rd)
rd['resultCollection']['properties']['schema_version']={'const':'ecbench-algebra-results-v1'}
rd['resultCollection']['properties']['results']=arr(ref('caseResult'),1)
rd['resultCollection']['properties']['summary']=ref('aggregate')
rd['resultCollection']['properties']['cohorts']=dynamic(r'^[A-Za-z][A-Za-z0-9_-]*$',ref('aggregate'))
rd['stratifiedAggregate']=copy.deepcopy(rd['aggregate'])
rd['stratifiedAggregate']['properties']['cohort']=infer([row['cohort'] for result in RESULTS for row in result['stratified_cohorts']],'cohort',rd)
rd['stratifiedAggregate']['required'].append('cohort')
rd['resultCollection']['properties']['stratified_cohorts']=arr(ref('stratifiedAggregate'))
rd['resultCollection']['properties']['timestamp_utc']={'type':'string','format':'date-time'}
rd['resultCollection']['properties']['configuration']['properties']['method']={'enum':['buchberger','f5b']}
rd['resultCollection']['properties']['configuration']['properties']['order']={'enum':['lex','grlex','grevlex']}
write('algebra_results.schema.json',rd,ref('resultCollection'),RESULTS,'ECBench measured algebraic relation-search results')
