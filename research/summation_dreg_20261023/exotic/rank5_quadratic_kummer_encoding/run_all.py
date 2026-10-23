import sys,json,hashlib,os
from kummer import *
outdir=sys.argv[1]; os.makedirs(outdir,exist_ok=True)
CONFIGS=[(1009,4,5),(1009,6,6),(1009,8,1),(1009,16,2),(2017,8,3)]
res=[]
for (p,d,seed) in CONFIGS:
    inst=pick_instance(p,d,seed)
    planted={'x1':inst['planted'][0],'x2':inst['planted'][1],'x3':inst['planted'][2]}
    for m in (2,3):
        if p==2017 and m==3 and d==12: pass
        ins=inst
        if m==2:
            ins=dict(inst); P=[Q for Q in inst['pts'] if Q and Q[0] in (inst['planted'][0],inst['planted'][1])]
            # m=2 planted: x1,x2 -> target x(P1+P2) using the same y-choices as pick (first matching points)
            P1=next(Q for Q in inst['pts'] if Q and Q[0]==inst['planted'][0]); P2=next(Q for Q in inst['pts'] if Q and Q[0]==inst['planted'][1])
            R2=ec_add(P1,P2,inst['a'],inst['p']); ins['xR']=R2[0] if R2 else None
            if ins['xR'] is None or ins['xR'] in inst['V']: ins['xR']=None
            if ins['xR'] is None: continue
        sysd=build(ins,m)
        cor=correctness(ins,m) 
        for name in (['a','b','c1','c2'] if m==3 else ['a','c1','c2']):
            eqs=sysd[name]; vs=varsof(eqs)
            h=hashlib.sha256((str(vs)+str([str(e) for e in eqs])).encode()).hexdigest()[:16]
            bz,bzfull,ds=bezout(eqs,vs,p)
            ms,_=msolve_deg(eqs,vs,p)
            sg=singular(eqs,vs,p,planted)
            row=dict(p=p,d=d,seed=seed,m=m,sys=name,hash=h,nvars=len(vs),neqs=len(eqs),maxdeg=max(ds),
                     degs=ds,bezout_n=bz,bezout_full=bzfull,msolve_deg=ms,**sg,
                     a_b=inst['a'],b_b=inst['b'],c=inst['c'],xR=ins['xR'],planted=inst['planted'])
            if name=='a': row['correctness']=cor
            print(row['p'],d,m,name,'vars',len(vs),'maxdeg',max(ds),'ms',ms,'sing',sg['sing_deg'],'vdim',sg['vdim'],'pl',sg['planted_ok'],flush=True)
            res.append(row)
json.dump(res,open(f'{outdir}/results.json','w'),indent=1)
