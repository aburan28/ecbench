import subprocess, re, os, tempfile, hashlib
from lib import poly_text, tdeg
PRE=['nice','-n','10','timeout','300','prlimit','--as=2500000000']
def sing_script(s,p):
    v=','.join(s['vars']); eqs=[poly_text(e,s['vars'],p) for e in s['eqs']]
    pl=','.join(f'{x}-{val}' for x,val in zip(s['vars'],s['planted']))
    # d_solv(Singular) := least degBound d for which std(I) (degBound=d) already has vdim == vdim of the full GB
    return (f'ring r={p},({v}),dp;\nideal I={",".join(eqs)};\ndegBound=0;\nideal G=std(I);\nint V=vdim(G);\n'
            f'print("VDIM "+string(V));print("DIM "+string(dim(G)));\n'
            f'ideal Mp={pl};ideal H=std(G+Mp);print("PLANTED "+string(H[1]!=1));\n'
            f'int d=1;int found=0;ideal Gd;\n'
            f'while(found==0 and d<=60){{degBound=d;Gd=std(I);if(vdim(Gd)==V){{found=d;}}d=d+1;}}\n'
            f'print("DSOLV "+string(found));\n$\n')
def ms_text(s,p):
    return ','.join(s['vars'])+f'\n{p}\n'+',\n'.join(poly_text(e,s['vars'],p) for e in s['eqs'])+'\n'
def syshash(s,p): return hashlib.sha256(ms_text(s,p).encode()).hexdigest()[:16]
def sing_degs(seg):
    """Singular prot: bare integers = degree of new element (printed on change); '[a:b]' = init data, '(n)' = #pairs,
    's' new element, '-' zero reduction. Adjacent degree changes concatenate ('78' = 7 then 8): split each digit run into
    the maximum number of non-decreasing tokens."""
    seg=re.sub(r'\[[^\]]*\]','',seg); seg=re.sub(r'\(\d+\)','',seg)
    prev=0; toks=[]
    for run in re.findall(r'\d+',seg):
        best=None
        def rec(i,last,acc):
            nonlocal best
            if i==len(run):
                if best is None or len(acc)>len(best): best=list(acc)
                return
            for j in range(i+1,len(run)+1):
                v=int(run[i:j])
                if run[i]=='0': break
                if v>=last: rec(j,v,acc+[v])
        rec(0,prev,[])
        if best is None: best=[int(run)]
        toks+=best; prev=best[-1]
    return toks

def run_singular(s,p):
    try:
        r=subprocess.run(PRE+['Singular','-q'],input=sing_script(s,p),capture_output=True,text=True,timeout=320)
    except subprocess.TimeoutExpired: return dict(status='timeout')
    o=r.stdout+r.stderr
    vd=re.search(r'VDIM (-?\d+)',o); pl=re.search(r'PLANTED (\d)',o); ds=re.search(r'DSOLV (\d+)',o)
    if not vd: return dict(status='fail',tail=o[-300:])
    return dict(status='ok',dsolv=(int(ds.group(1)) or None) if ds else None,vdim=int(vd.group(1)),planted=bool(int(pl.group(1))) if pl else None)

def run_msolve(s,p):
    with tempfile.NamedTemporaryFile('w',suffix='.ms',delete=False) as fh: fh.write(ms_text(s,p)); path=fh.name
    outp=path+'.out'
    try:
        r=subprocess.run(PRE+['msolve','-t','1','-v','2','-f',path,'-g','2','-o',outp],capture_output=True,text=True,timeout=320)
    except subprocess.TimeoutExpired: return dict(status='timeout')
    finally: pass
    o=r.stdout+r.stderr
    degs=[int(m.group(1)) for line in o.splitlines() for m in [re.match(r'\s*(\d+)\s+\d+\s+\d+\s+\d+ x \d+',line)] if m]
    gb=open(outp).read() if os.path.exists(outp) else ''
    os.unlink(path); 
    if os.path.exists(outp): os.unlink(outp)
    return dict(status='ok' if degs or gb else 'fail',dsolv=max(degs) if degs else None,gb=gb)

def ms_planted(gb,s,p):
    body=gb.split('#---\n')[-1].strip().rstrip(':').strip()
    if body.startswith('[') : body=body[1:]
    if body.endswith(']'): body=body[:-1]
    env=dict(zip(s['vars'],s['planted']))
    for poly in body.split(',\n'):
        poly=poly.strip().replace('^','**')
        if poly and eval(poly,{},env)%p!=0: return False
    return True
