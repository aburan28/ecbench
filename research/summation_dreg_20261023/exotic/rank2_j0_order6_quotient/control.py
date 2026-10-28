import sys,os,json,re,subprocess,tempfile; sys.argv=['run.py',sys.argv[1]]
OUT=sys.argv[1]
exec(open('run.py').read().split("if __name__")[0])
def msolve_full(vars_,p,eqs):
    txt=",".join(vars_)+"\n"+str(p)+"\n"+",\n".join(eqs)+"\n"
    with tempfile.NamedTemporaryFile('w',suffix='.ms',delete=False) as f: f.write(txt); path=f.name
    try:
        r=subprocess.run(f'nice -n 10 timeout 300 prlimit --as=2500000000 msolve -t 1 -v 2 -g 2 -f {path} -o /dev/null',shell=True,capture_output=True,text=True,timeout=320); out=r.stdout+r.stderr
    finally: os.unlink(path)
    deg=[];mat=[]
    for l in out.splitlines():
        m=re.match(r'\s*(\d+)\s+\S+\s+\d+\s+(\d+)\s+x\s+(\d+)',l)
        if m: deg.append(int(m.group(1))); mat.append(int(m.group(2))*int(m.group(3)))
    return dict(deg=max(deg) if deg else None,maxmat=max(mat) if mat else None)
grid={61:[6,12,30,60],97:[6,12,24,48,96],193:[6,12,24,48,96],1009:[9,18,36,63,126],4093:[6,12,33,66],65521:[6,12,18,36,45,90]}
res=[]
for p,ds in grid.items():
  for d in ds:
    for seed in (0,1):
        inst=instance(p,d,seed); b,c=inst['b'],inst['c']; V=inst['V']; x0=inst['planted'][0]
        c2=pow(x0,d//3,p); W=[x for x in V if pow(x,d//3,p)==c2]; assert len(W)==d//3
        rnd=random.Random(f"ctl-{p}-{d}-{seed}")
        cand=[x for x in W if sqrt_mod(x**3+b,p) not in (None,0)]
        while True:
            x2=rnd.choice(cand); y0=sqrt_mod(x0**3+b,p); y2=sqrt_mod(x2**3+b,p)
            R=ec_add((x0,y0),(x2,y2),b,p)
            if R is not None and R[1]!=0: break
        xRc=R[0]
        Nctl=[(a,bb) for a in W for bb in W if S3n(a,bb,xRc,b,p)==0]
        assert (x0,x2) in Nctl
        ceqs=[f"x1^{d//3}-{c2}",f"x2^{d//3}-{c2}",poly_str_S3(p,b,xRc)]
        rad="LIB \"primdec.lib\"; ideal Rd=radical(I); ideal GR=std(Rd); print(\"RADVDIM=\"+string(vdim(GR)));\n"
        sg=run_singular(['x1','x2'],p,ceqs,rad); ms=msolve_full(['x1','x2'],p,ceqs)
        # original systems for proxy
        coef=Qnum_coeffs(inst['xR'],b,p)
        peq=[f"x1^{d}-{c}",f"x2^{d}-{c}",poly_str_S3(p,b,inst['xR'])]
        qeq=[f"u1^{d//3}-{c}",f"u2^{d//3}-{c}",poly_str_Q(coef,p)]
        Np=sum(1 for a in V for bb in V if S3n(a,bb,inst['xR'],b,p)==0)
        U=sorted({pow(x,3,p) for x in V}); Nq=sum(1 for a in U for bb in U if Qval(coef,a,bb,p)==0)
        r=dict(p=p,d=d,seed=seed,N_plain=Np,N_quot=Nq,N_quot_lifted=Np,N_ctrl=len(Nctl),
          ctrl=dict(hash=hashlib.sha256(("|".join(ceqs)+f"|{p}").encode()).hexdigest()[:16],singular=sg,msolve=ms,vdim=sg['vdim'],rad_ok=sg['rad']==len(Nctl)),
          plain=dict(singular=run_singular(['x1','x2'],p,peq),msolve=msolve_full(['x1','x2'],p,peq)),
          quot=dict(singular=run_singular(['u1','u2'],p,qeq),msolve=msolve_full(['u1','u2'],p,qeq)))
        res.append(r); print(p,d,seed,Np,Nq,len(Nctl),r['ctrl']['singular']['deg'],r['ctrl']['msolve']['deg'],flush=True)
os.makedirs(OUT,exist_ok=True); json.dump(res,open(os.path.join(OUT,'control.json'),'w'),indent=1,sort_keys=True)
