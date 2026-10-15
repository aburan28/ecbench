"""Exact finite-domain fixture validation, independent of Groebner solving."""
import itertools
import time

from .algebra import _normalize_case, _record_input, _eval_sparse, _verify_relations


def verify_cases(suite, progress=None):
    cases=suite['cases'] if isinstance(suite,dict) else suite
    rows=[]
    started=time.perf_counter()
    for original in cases:
        case=_normalize_case(original)
        _,modulus,polynomials,domains,metrics=_record_input(case)
        case_started=time.perf_counter()
        roots=[list(values) for values in itertools.product(*domains)
               if all(_eval_sparse(polynomial.items(),values,modulus)==0 for polynomial in polynomials)]
        relations=_verify_relations(case,roots)
        if not relations['truth_matches']:
            raise ValueError(f"Fixture truth mismatch: {case['id']}")
        row={'case_id':case['id'],'status':'pass','domain_assignments_checked':metrics['domain_grid_size'],
             'polynomial_roots':len(roots),'ordered_relations':relations['ordered_relations'],
             'unordered_relations':relations['unordered_relations'],
             'elapsed_seconds':time.perf_counter()-case_started}
        rows.append(row)
        if progress:progress(row)
    return {'status':'pass','cases_checked':len(rows),
            'fields_checked':len({case['field']['id'] for case in cases}),
            'domain_assignments_checked':sum(row['domain_assignments_checked'] for row in rows),
            'scope':'Exhaustive evaluation of original sparse systems over declared finite domains, independent subgroup checks and signed point reconstruction; not a Groebner benchmark or ideal-equality certificate',
            'elapsed_seconds':time.perf_counter()-started,'results':rows}
