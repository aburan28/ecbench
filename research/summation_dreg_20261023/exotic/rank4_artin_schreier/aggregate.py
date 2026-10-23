import json, glob, os, hashlib, sys
HERE = os.path.dirname(os.path.abspath(__file__))
def load(run):
    d = {}
    for f in sorted(glob.glob(f"{HERE}/raw/{run}/cell_*.json")):
        r = json.load(open(f)); d[(r['n'], r['l'], r['m'], r['model'])] = r
    return d
def strip(r):  # deterministic content only (no wall time)
    r = json.loads(json.dumps(r)); r.pop('wall_s', None); return r
def det_digest(r):
    return hashlib.sha256(json.dumps(strip(r), sort_keys=True).encode()).hexdigest()[:12]
def fmt(x): return '-' if x is None else str(x)
def row(r):
    f = r['full']; a = r['as']['agg']; c = r['control']['agg']
    def pair(eng, g): return fmt(g[eng]['max']) + ('*' if g[eng]['any_timeout'] else '')
    ff = f['singular']['deg'], f['msolve']['deg']
    fs = fmt(f['singular']['deg']) + ('*' if f['singular']['status'] != 'ok' else '')
    fm = fmt(f['msolve']['deg']) + ('*' if f['msolve']['status'] != 'ok' else '')
    mac = f.get('mac', {})
    pb = [b for b in r['as']['branches'] if b['planted']][0]
    return dict(n=r['n'], l=r['l'], m=r['m'], model=r['model'], full_vars=f['nvars'], full_eqs=f['neqs'], din=f['din'],
                fullS=fs, fullM=fm, fullS0=fmt(r['full0_origcoords']['singular']['deg']), fullM0=fmt(r['full0_origcoords']['msolve']['deg']),
                asS=pair('singular', a), asM=pair('msolve', a), ctlS=pair('singular', c), ctlM=pair('msolve', c),
                asS_pl=fmt(a['singular']['planted']), asM_pl=fmt(a['msolve']['planted']),
                br_vars=r['as']['nvars_branch'], nbr=r['as']['nbranches'], vdim=f['singular']['vdim'],
                ff_full=fmt(mac.get('ff')) + f"(<= {mac.get('Dmax')})", ff_as_planted=fmt(pb.get('mac', {}).get('ff')) + f"(<= {pb.get('mac', {}).get('Dmax')})",
                ok=(r['as']['all_branch_vdim_match'] and r['as']['sum_vdim_eq_numeric'] and r['control']['all_branch_vdim_match']
                    and r['checks']['BFnondegenerate_subset_SF'] and r['checks']['planted_in_SF']))
if __name__ == '__main__':
    r1, r2 = load('run1'), load('run2')
    diffs = [k for k in r1 if k not in r2 or det_digest(r1[k]) != det_digest(r2[k])]
    diffs += [k for k in r2 if k not in r1]
    print("cells run1/run2:", len(r1), len(r2), "differing:", diffs)
    rows = [row(r1[k]) for k in sorted(r1)]
    keys = ['n','l','m','model','full_vars','br_vars','nbr','din','vdim','fullS','asS','ctlS','fullM','asM','ctlM','fullS0','fullM0','asS_pl','asM_pl','ff_full','ff_as_planted','ok']
    print(" | ".join(keys))
    for x in rows: print(" | ".join(str(x[k]) for k in keys))
    json.dump(dict(rows=rows, diffs=[list(k) for k in diffs], ncells=(len(r1), len(r2))), open(f"{HERE}/raw/summary.json", 'w'), indent=1)
