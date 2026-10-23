import sys,json,time
H='/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness'
sys.path.insert(0,H)
import exp_chain as E
n,l,m,auxV,Dmax=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]=='1',int(sys.argv[5])
t=time.time()
r=E.run(1,1,n,l,m,seed=1,Dmax=Dmax,aux_in_V=auxV,do_msolve=True)
print(json.dumps({k:r[k] for k in ('n','l','m','aux_in_V','nvars','neqs','first_fall_degree','singular','sys_sha','msolve','t_macaulay','t_singular','t_msolve') if k in r}))
print("TOTAL_WALL",round(time.time()-t,1))
