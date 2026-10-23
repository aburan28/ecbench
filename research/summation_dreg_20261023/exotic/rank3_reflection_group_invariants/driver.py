import sys, json, time
sys.path.insert(0,'.')
from systems import *; from run_solvers import *
run=sys.argv[1]
INST=[(2,2,p,s) for p in (11,13,17,31,101) for s in (1,2)]+[(3,3,p,s) for p in (11,13,17,19,23) for s in (1,2,3)]
res=[]
for (m,n,p,seed) in INST:
    try: out,meta=build_systems(p,n,m,seed)
    except Exception as e:
        res.append(dict(m=m,n=n,p=p,seed=seed,error=repr(e))); print('ERR',m,p,seed,e,flush=True); continue
    row=dict(meta=meta,sys={})
    for k,s in out.items():
        t=time.time()
        sg=run_singular(s,p); ms=run_msolve(s,p)
        rec=dict(vars=len(s['vars']),neq=len(s['eqs']),eqdeg=[tdeg(e) for e in s['eqs']],
                 bezout=int(__import__('math').prod(max(tdeg(e),1) for e in s['eqs'])),hash=syshash(s,p),
                 sing=sg,msolve_dsolv=ms.get('dsolv'),msolve_status=ms['status'],
                 msolve_planted=(ms_planted(ms['gb'],s,p) if ms.get('gb') else None))
        row['sys'][k]=rec
        print(m,p,seed,k,rec['vars'],rec['eqdeg'],'sing',sg.get('dsolv'),'',sg.get('vdim'),sg.get('planted'),'ms',ms.get('dsolv'),rec['msolve_planted'],flush=True)
    res.append(row)
json.dump(res,open(f'raw/{run}.json','w'),indent=1)
