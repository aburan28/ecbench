"""Golden fixtures, polynomial solvers, and repeatable EC benchmarks."""
import argparse
import collections
import datetime
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import sys
import time
from . import __version__
from .reference import Field,Curve
from .validate import load_dataset,verify


def catalog(name, data_dir=None, override=None):
    directory=Path(data_dir) if data_dir else Path(__file__).parent/'data'
    path=Path(override) if override else directory/(name+'.json')
    result=json.loads(path.read_text())
    if name=='algebra_cases' and not override:
        extra=directory/'algebra_coverage_cases.json'
        if extra.exists():
            coverage=json.loads(extra.read_text())
            cases=result['cases']+coverage['cases']
            if len({row['id'] for row in cases})!=len(cases):raise ValueError('Duplicate polynomial case IDs across catalogs')
            result={**result,'suite_id':'ecbench-all-field-polynomial-battery-v1','cases':cases,
                'component_suites':[{'name':'main','cases':len(result['cases'])},{'name':'all-field coverage','cases':len(coverage['cases'])}],
                'counts':{'cases':len(cases),'fields':len({row['field']['id'] for row in cases})}}
    return result


def select_algebra_cases(suite, case_ids=None, curve_ids=None, max_cases=None):
    selected=list(suite['cases'])
    known={row['id'] for row in selected}
    if set(case_ids or [])-known:raise ValueError('Unknown polynomial case ID')
    curves={row['curve_id'] for row in selected}
    if set(curve_ids or [])-curves:raise ValueError('No polynomial cases for selected curve ID')
    if case_ids:selected=[row for row in selected if row['id'] in case_ids]
    if curve_ids:selected=[row for row in selected if row['curve_id'] in curve_ids]
    if max_cases is not None:
        if max_cases<1:raise ValueError('max-cases must be positive')
        selected=selected[:max_cases]
    if not selected:raise ValueError('Selection contains no polynomial cases')
    return {**suite,'suite_id':suite.get('schema','ecbench-algebra-cases-v1'),'cases':selected}


def verify_catalog_hashes(data_dir=None):
    directory=Path(data_dir) if data_dir else Path(__file__).parent/'data'
    manifest=json.loads((directory/'research_manifest.json').read_text())
    for name,expected in manifest['sha256'].items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest()!=expected:
            raise ValueError(f'Research payload hash mismatch: {name}')
    return {'status':'pass','files_checked':len(manifest['sha256']),'counts':manifest['counts']}


def emit(value,output=None):
    text=json.dumps(value,sort_keys=True,indent=2)+'\n'
    if output:
        path=Path(output);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text)
    print(text,end='')


def benchmark(data,ids,samples,iterations):
    preflight=verify(data)
    fields={f['id']:Field(f) for f in data['fields']}
    rows=[]
    for record in data['curves']:
        if ids and record['id'] not in ids:continue
        f=fields[record['field_id']];c=Curve(record,f)
        p=c.parse_point(record['named_points']['R']['point']);q=c.double(p)
        scalar=(1<<(f.bits-1))+0x12345
        a=p[0] or 1;b=p[1] or 1
        operations=[('field_mul',lambda:f.mul(a,b),{'a':f.hex(a),'b':f.hex(b)}),
          ('field_inverse',lambda:f.inverse(a),{'a':f.hex(a)}),
          ('point_add_distinct',lambda:c.add(p,q),{'P':c.point(p),'Q':c.point(q)}),
          ('point_double',lambda:c.double(p),{'P':c.point(p)}),
          ('point_mul',lambda:c.mul(scalar,p),{'scalar':str(scalar),'P':c.point(p)})]
        seen_maps=set()
        for vector in data['vectors']:
            if vector['domain']!='map' or vector['curve_id']!=record['id'] or vector['operation'] in seen_maps:
                continue
            seen_maps.add(vector['operation']);operation_name=vector['operation']
            parameter=f.parse(vector['inputs']['parameter']) if 'parameter' in vector['inputs'] else None
            inputs={'P':c.point(p),'target_coefficients':vector['inputs']['target_coefficients']}
            if parameter is not None:inputs['parameter']=f.hex(parameter)
            operations.append(('map_'+operation_name,lambda name=operation_name,unit=parameter:c.map_point(name,p,unit),inputs))
        for name,operation,inputs in operations:
            operation() # Warm up outside the measured interval.
            elapsed=[]
            for _ in range(samples):
                started=time.perf_counter_ns()
                for _ in range(iterations):operation()
                elapsed.append((time.perf_counter_ns()-started)/iterations)
            ordered=sorted(elapsed)
            rows.append({'curve_id':record['id'],'field_id':record['field_id'],'operation':name,
              'inputs':inputs,'samples':samples,'iterations_per_sample':iterations,'unit':'nanoseconds_per_operation',
              'min':min(elapsed),'median':statistics.median(elapsed),'mean':statistics.mean(elapsed),
              'p95':ordered[max(0,math.ceil(.95*len(ordered))-1)],'sample_values':elapsed})
    return {'dataset_id':data['manifest']['dataset_id'],'release_version':data['manifest']['release_version'],
      'manifest_sha256':hashlib.sha256(json.dumps(data['manifest'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
      'implementation':'ecbench variable-time Python reference','implementation_version':__version__,
      'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'environment':{'python':platform.python_version(),'python_implementation':platform.python_implementation(),
        'platform':platform.platform(),'machine':platform.machine(),'processor':platform.processor()},
      'preflight':preflight,'timing_scope':'arithmetic call only; dataset load and preflight excluded',
      'results':rows}


def main(argv=None):
    parser=argparse.ArgumentParser(prog='ecbench',description='Golden EC fixtures and polynomial index-calculus benchmarks.')
    parser.add_argument('--data-dir',help='Directory containing the canonical dataset payloads.')
    parser.add_argument('--bundle',help='Use the combined golden_ecbench.json export.')
    sub=parser.add_subparsers(dest='command',required=True)
    listing=sub.add_parser('list',help='List the 56 benchmark curves.');listing.add_argument('--field');listing.add_argument('--json',action='store_true')
    show=sub.add_parser('show',help='Show one complete curve record.');show.add_argument('--curve',required=True)
    check=sub.add_parser('verify',help='Check hashes, certificates, parameters, point orders, and every vector.')
    check.add_argument('--with-pari',action='store_true');check.add_argument('--gp',help='Path to PARI/GP executable.');check.add_argument('--output')
    bench=sub.add_parser('bench',help='Measure the reference arithmetic with deterministic inputs.')
    bench.add_argument('--curve',action='append',help='Curve ID; repeat to select several. Default: all curves.')
    bench.add_argument('--samples',type=int,default=5);bench.add_argument('--iterations',type=int,default=1);bench.add_argument('--output')
    case_list=sub.add_parser('algebra-cases',help='List the polynomial relation-search fixtures.')
    case_list.add_argument('--cases',help='Alternative polynomial case JSON catalog.');case_list.add_argument('--curve',action='append');case_list.add_argument('--json',action='store_true')
    case_show=sub.add_parser('algebra-show',help='Show one exact sparse polynomial system and its golden relations.')
    case_show.add_argument('--cases');case_show.add_argument('--case',required=True)
    case_verify=sub.add_parser('algebra-verify',help='Exhaustively verify original polynomial fixtures and signed golden relations without a solver.')
    case_verify.add_argument('--cases');case_verify.add_argument('--case',action='append');case_verify.add_argument('--curve',action='append');case_verify.add_argument('--output')
    algebra=sub.add_parser('algebra-bench',aliases=['algebra'],help='Measure SymPy relation solving, basis degree, regularity, and solve time.')
    algebra.add_argument('--cases');algebra.add_argument('--case',action='append');algebra.add_argument('--curve',action='append')
    algebra.add_argument('--method',choices=['buchberger','f5b'],default='buchberger')
    algebra.add_argument('--order',choices=['lex','grlex','grevlex'],default='grevlex')
    algebra.add_argument('--timeout',type=float,default=5.0,help='Wall-time budget in seconds for each isolated case.')
    algebra.add_argument('--max-cases',type=int);algebra.add_argument('--no-diagnostics',action='store_true');algebra.add_argument('--output')
    export=sub.add_parser('backend-export',help='Export an executable Sage/Singular solver adapter for a case.')
    export.add_argument('--cases');export.add_argument('--case',required=True);export.add_argument('--destination',required=True)
    export.add_argument('--order',choices=['lex','grlex','grevlex'],default='grevlex')
    replay=sub.add_parser('backend-replay',help='Verify an external sparse-basis result and its exact point relations.')
    replay.add_argument('--cases');replay.add_argument('--case',required=True);replay.add_argument('--result',required=True);replay.add_argument('--output')
    endo=sub.add_parser('endomorphism-cases',help='Show supported GLV2/Frobenius maps and unsupported subgroup cases.')
    endo.add_argument('--curve',action='append');endo.add_argument('--supported',action='store_true');endo.add_argument('--catalog')
    endov=sub.add_parser('endomorphism-verify',help='Verify every GLV2/Frobenius eigenvalue and decomposition.')
    endov.add_argument('--catalog');endov.add_argument('--output')
    endob=sub.add_parser('endomorphism-bench',help='Measure GLV2/Frobenius decomposition and signed joint multiplication.')
    endob.add_argument('--catalog');endob.add_argument('--curve',action='append');endob.add_argument('--samples',type=int,default=3);endob.add_argument('--iterations',type=int,default=1);endob.add_argument('--output')
    for name,help_text in [('glv4-show','Show the supplementary GLV/GLS4 quadratic-field curve.'),('glv4-verify','Verify GLV/GLS4 maps, lattices, and golden answers.'),('glv4-bench','Measure two- and four-term joint scalar multiplication.')]:
        command=sub.add_parser(name,help=help_text);command.add_argument('--catalog');command.add_argument('--output')
        if name=='glv4-bench':command.add_argument('--samples',type=int,default=3);command.add_argument('--iterations',type=int,default=1)
    battery=sub.add_parser('battery-verify',help='Verify the core corpus, original polynomial fixtures, release hashes, and endomorphism catalogs.')
    battery.add_argument('--output')
    args=parser.parse_args(argv)
    try:
        data=load_dataset(args.data_dir,args.bundle)
        ids={c['id'] for c in data['curves']}
        if args.command=='list':
            rows=[{'id':c['id'],'field_id':c['field_id'],'model':c['model'],'ordinary':c['invariants']['ordinary'],
               'subgroup_bits':int(c['subgroup']['order']).bit_length(),'cofactor':c['subgroup']['cofactor']} for c in data['curves'] if not args.field or c['field_id']==args.field]
            if args.json:emit(rows)
            else:
                print(f"{'CURVE':16} {'FIELD':7} {'r BITS':>7} {'COFACTOR':>22} STATUS")
                for row in rows:print(f"{row['id']:16} {row['field_id']:7} {row['subgroup_bits']:7} {row['cofactor']:>22} {'ordinary' if row['ordinary'] else 'supersingular'}")
        elif args.command=='show':
            if args.curve not in ids:raise ValueError(f'Unknown curve: {args.curve}')
            emit(next(c for c in data['curves'] if c['id']==args.curve))
        elif args.command=='verify':emit(verify(data,args.with_pari,args.gp),args.output)
        elif args.command=='bench':
            if args.samples<1 or args.iterations<1:raise ValueError('samples and iterations must be positive')
            if set(args.curve or [])-ids:raise ValueError('Unknown curve ID in benchmark selection')
            emit(benchmark(data,set(args.curve or []),args.samples,args.iterations),args.output)
        elif args.command in {'algebra-cases','algebra-show','algebra-verify','algebra-bench','algebra','backend-export','backend-replay'}:
            suite=catalog('algebra_cases',args.data_dir,args.cases)
            if args.command=='algebra-cases':
                selected=select_algebra_cases(suite,curve_ids=args.curve)['cases']
                rows=[{'id':row['id'],'curve_id':row['curve_id'],'variables':len(row['variables']),
                    'polynomials':len(row['polynomials']),'target_kind':row['target']['kind'],
                    'factor_base_x':len(row['factor_base']['points_by_x']),
                    'golden_ordered_relations':row['expected']['ordered_count']} for row in selected]
                if args.json:emit(rows)
                else:
                    for row in rows:print(f"{row['id']:64} vars={row['variables']:2} base={row['factor_base_x']:3} relations={row['golden_ordered_relations']}")
            elif args.command=='algebra-verify':
                from .fixtures import verify_cases
                emit(verify_cases(select_algebra_cases(suite,args.case,args.curve)),args.output)
            elif args.command in {'algebra-show','backend-export','backend-replay'}:
                case=select_algebra_cases(suite,[args.case])['cases'][0]
                if args.command=='algebra-show':emit(case)
                elif args.command=='backend-export':
                    from .backends import export_sage_case
                    emit(export_sage_case(case,args.destination,args.order))
                else:
                    from .backends import replay_external_basis
                    result=replay_external_basis(case,json.loads(Path(args.result).read_text()))
                    emit(result,args.output)
                    if result['status']!='ok':return 1
            else:
                from .algebra import run_suite
                selected=select_algebra_cases(suite,args.case,args.curve,args.max_cases)
                def progress(row):print(f"{row['case_id']}: {row['status']}",file=sys.stderr,flush=True)
                result=run_suite(selected,args.timeout,args.method,args.order,diagnostics=not args.no_diagnostics,progress=progress)
                emit(result,args.output)
                if result['summary']['errors']:return 1
        elif args.command.startswith('endomorphism-'):
            from .endomorphisms import verify_cases,benchmark_cases
            cases=catalog('endomorphism_cases',args.data_dir,args.catalog)
            if args.command=='endomorphism-cases':
                if set(args.curve or [])-ids:raise ValueError('Unknown curve ID')
                emit([row for row in cases['cases'] if (not args.curve or row['curve_id'] in args.curve) and (not args.supported or row['status']=='supported')])
            elif args.command=='endomorphism-verify':emit(verify_cases(data,cases),args.output)
            else:emit(benchmark_cases(data,cases,args.curve,args.samples,args.iterations),args.output)
        elif args.command.startswith('glv4-'):
            cases=catalog('glv4',args.data_dir,args.catalog)
            if args.command=='glv4-show':emit(cases,args.output)
            else:
                from .glv4 import verify_cases,benchmark_cases
                emit(verify_cases(cases) if args.command=='glv4-verify' else benchmark_cases(cases,args.samples,args.iterations),args.output)
        elif args.command=='battery-verify':
            from .endomorphisms import verify_cases as verify_endo
            from .glv4 import verify_cases as verify_glv4
            from .fixtures import verify_cases as verify_fixtures
            emit({'core':verify(data),'research_payloads':verify_catalog_hashes(args.data_dir),
                'glv2_frobenius':verify_endo(data,catalog('endomorphism_cases',args.data_dir)),
                'glv4':verify_glv4(catalog('glv4',args.data_dir)),
                'polynomial_fixtures':verify_fixtures(catalog('algebra_cases',args.data_dir))},args.output)
    except (ValueError,AssertionError,KeyError,RuntimeError,OSError,ImportError) as error:
        print(json.dumps({'status':'fail','error':str(error)}),file=sys.stderr)
        return 1
    return 0
