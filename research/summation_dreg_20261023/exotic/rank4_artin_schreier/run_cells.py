import sys, os, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import as_lib as L
from collections import Counter

def run_cell(n, l, m, model):
    t0 = time.time()
    cell = L.make_cell(n, l, m, model)
    br = L.brute(cell)
    res = dict(n=n, l=l, m=m, model=model, seed=cell['seed'], xR=cell['xR'], xs=cell['xs'], aux=cell['aux'],
               B_adapted=cell['Bad'], nsol_numeric=len(br['sols']), n_tuples=len(br['SFx']), n_BF=len(br['BF']), n_BFc=len(br['BFc']))
    chk = dict(BF_subset_SF=br['BF'] <= br['SFx'], BFnondegenerate_subset_SF=br['BFc'] <= br['SFx'], BF_missing_from_SF=len(br['BF'] - br['SFx']),
               planted_in_BF=tuple(cell['xs']) in br['BF'],
               planted_in_SF=tuple(cell['xs']) in br['SFx'],
               spurious_tuples=len(br['SFx'] - br['BF']))
    res['full'] = L.measure(cell['eqs'], cell['vars'], mac=True)
    res['full0_origcoords'] = L.measure(cell['eqs0'], cell['vars0'])
    chk['full_vdim_eq_numeric'] = (res['full']['singular']['vdim'] == len(br['sols'])) if res['full']['singular']['status'] == 'ok' else None
    res['checks'] = chk
    for arm in ('as', 'control'):
        cnt = Counter(L.branch_key(cell, s, arm) for s in br['sols'])
        planted_key = tuple(cell['planted'][k] for k in next(L.branch_assigns(cell, arm))[1].keys())
        branches = []; vsum = 0; allok = True
        for bits, asg in L.branch_assigns(cell, arm):
            eq = L.subst(cell['eqs'], asg); vs = [v for v in cell['vars'] if v not in asg]
            r = L.measure(eq, vs, mac=(bits == planted_key))
            r['bits'] = list(bits); r['planted'] = (bits == planted_key); r['numeric_count'] = cnt.get(bits, 0)
            if r['singular']['status'] == 'ok':
                r['vdim_match'] = (r['singular']['vdim'] == r['numeric_count'])
                vsum += r['singular']['vdim']; allok &= r['vdim_match']
            else: allok = False
            branches.append(r)
        # planted must satisfy planted branch system
        pb = [b for b in branches if b['planted']][0]
        agg = {}
        for eng in ('singular', 'msolve'):
            ds = [b[eng]['deg'] for b in branches]
            to = any(b[eng]['status'] != 'ok' for b in branches)
            valid = [d for d in ds if d is not None]
            agg[eng] = dict(max=max(valid) if valid else None, min=min(valid) if valid else None,
                            planted=pb[eng]['deg'], any_timeout=to, n_none=sum(d is None for d in ds))
        res[arm] = dict(nbranches=len(branches), branches=branches, vdim_sum=vsum, all_branch_vdim_match=allok,
                        sum_vdim_eq_numeric=(vsum == len(br['sols'])), agg=agg,
                        nvars_branch=branches[0]['nvars'], neqs_max=max(b['neqs'] for b in branches))
    res['wall_s'] = round(time.time() - t0, 1)
    return res

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--run'); ap.add_argument('--ns', default='8,9,10,11,12')
    ap.add_argument('--ls', default='2,3'); a = ap.parse_args()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw', a.run); os.makedirs(out, exist_ok=True)
    cells = []
    for n in map(int, a.ns.split(',')):
        for l in map(int, a.ls.split(',')):
            cells += [(n, l, 2, 'single'), (n, l, 3, 'unsplit'), (n, l, 3, 'split')]
    for (n, l, m, model) in cells:
        fn = os.path.join(out, f"cell_n{n}_l{l}_m{m}_{model}.json")
        if os.path.exists(fn): continue
        r = run_cell(n, l, m, model)
        json.dump(r, open(fn, 'w'), indent=1, default=list)
        f = r['full']; print(n, l, m, model, 'full', f['singular']['deg'], f['msolve']['deg'], 'AS', r['as']['agg']['singular']['max'],
              r['as']['agg']['msolve']['max'], 'CTL', r['control']['agg']['singular']['max'], r['control']['agg']['msolve']['max'],
              r['checks'], r['wall_s'], flush=True)
