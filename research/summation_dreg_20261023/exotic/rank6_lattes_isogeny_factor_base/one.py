import time
from tower import *
def measure(p,family,k,d,s,seed):
    inst=make_instance(p,family,k,d,s,seed)
    if inst['plant'] is None: return dict(p=p,family=family,k=k,d=d,seed=seed,error='no plant',F=len(inst['F']))
    vars_,eqs,xR=build_system(inst)
    pv=planted_values(inst)
    rec=dict(p=p,family=family,k=k,d=d,s=inst['s'],seed=seed,nvars=len(vars_),neq=len(eqs),bezout=bezout(inst),
             Fsize=len(inst['F']),Fpts=len([x for x in inst['F'] if x in set(inst['ptsx'])]),
             hash=sys_hash(vars_,eqs,p),root_check=check_roots(inst,vars_,eqs))
    rng=random.Random(99); Rl=rng.sample(inst['ptsx'],min(60,len(inst['ptsx'])))
    rec['avg_decomp']=decomp_count(inst['F'],set(inst['ptsx']),inst['a'],inst['b'],p,Rl)
    # random-base control, same size, same seed-derived
    cr=random.Random(f"ctl-{p}-{seed}-{family}-{k}-{d}"); Fc=sorted(cr.sample(range(p),len(inst['F'])))
    rec['avg_decomp_randbase']=decomp_count(Fc,set(inst['ptsx']),inst['a'],inst['b'],p,Rl)
    rec['msolve']=run_msolve(vars_,eqs,p)
    rec['singular']=run_singular(vars_,eqs,p,pv)
    return rec
if __name__=='__main__':
    import sys
    p,fam,k,d,s,seed=sys.argv[1],sys.argv[2],*map(int,sys.argv[3:7])
    t=time.time(); print(json.dumps(measure(int(p),fam,k,d,s,seed),indent=1),time.time()-t)
