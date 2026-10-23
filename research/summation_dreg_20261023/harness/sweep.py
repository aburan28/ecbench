import sys,json,time
H='/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness'
sys.path.insert(0,H)
import exp_chain as E
GRID=[
 # (n,l,m,auxV)
 (6,3,2,False),(8,4,2,False),(10,5,2,False),(12,6,2,False),
 (6,2,3,False),(8,2,3,False),(10,2,3,False),(8,3,3,False),(9,3,3,False),(12,2,3,False),
 (6,2,3,True),(8,2,3,True),(8,3,3,True),(9,3,3,True),
]
out=[]
for (n,l,m,auxV) in GRID:
    nvars_est = m*l + (0 if m==2 else (m-2)*(l if auxV else n))
    Dmax = 5 if nvars_est>12 else 6
    t=time.time()
    try:
        r=E.run(1,1,n,l,m,seed=1,Dmax=Dmax,aux_in_V=auxV,do_msolve=True)
    except Exception as ex:
        r=dict(n=n,l=l,m=m,aux_in_V=auxV,status='error',err=str(ex)[:200])
    r['wall']=round(time.time()-t,1)
    out.append(r)
    sg=r.get('singular'); 
    print(json.dumps({k:r.get(k) for k in ('n','l','m','aux_in_V','nvars','neqs','first_fall_degree','singular','wall','status')}),flush=True)
    json.dump(out,open(H+'/out/sweep_chain.json','w'),indent=1)
print("DONE")
