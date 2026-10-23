import sys, json
from multiprocessing import Pool
from one import measure
S={211:16,521:24,1031:32}
cfgs=[]
for p,s in S.items():
    for (k,d) in ((3,2),(2,3),(2,4)):
        fams=['kummer','generic','isogeny']+(['lattes'] if d==4 else [])
        for f in fams:
            for seed in (1,2,3): cfgs.append((p,f,k,d,s,seed))
def w(c):
    try: return measure(*c)
    except Exception as e: return dict(p=c[0],family=c[1],k=c[2],d=c[3],seed=c[5],error=repr(e))
if __name__=='__main__':
    run=sys.argv[1]
    with Pool(2) as P: res=P.map(w,cfgs,chunksize=1)
    json.dump(res,open(f'raw/{run}.json','w'),indent=1)
    print(run,'done',len(res))
