import json,collections
R=json.load(open('run1/control.json'))
by=collections.defaultdict(list)
for r in R: by[r['d']].append(r)
def m(L,f): return sum(f(x) for x in L)/len(L)
print("rad_ok all",all(r['ctrl']['rad_ok'] for r in R))
print("d n | msolve P/Q/C | sing P/Q/C | N P(=Qlift)/C | vdim C | maxmat P/Q/C | dec/mat(x1e3) P/Q/C")
for d in sorted(by):
    L=by[d]
    f=lambda k,e:m(L,lambda x:x[k][e]['deg'])
    mat=lambda k:m(L,lambda x:x[k]['msolve']['maxmat'])
    NP=m(L,lambda x:x['N_plain']); NC=m(L,lambda x:x['N_ctrl'])
    print(d,len(L),'|',f('plain','msolve'),f('quot','msolve'),f('ctrl','msolve'),'|',round(f('plain','singular'),1),round(f('quot','singular'),1),round(f('ctrl','singular'),1),'|',round(NP,1),round(NC,1),'|',round(m(L,lambda x:x['ctrl']['vdim']),1),'|',round(mat('plain')),round(mat('quot')),round(mat('ctrl')),'|',round(1e3*NP/mat('plain'),2),round(1e3*NP/mat('quot'),2),round(1e3*NC/mat('ctrl'),2))
n=len(R)
print("N totals plain/ctrl",sum(r['N_plain'] for r in R),sum(r['N_ctrl'] for r in R))
for e in('msolve','singular'):
    print(e,"ctrl<quot",sum(r['ctrl'][e]['deg']<r['quot'][e]['deg'] for r in R),"ctrl>quot",sum(r['ctrl'][e]['deg']>r['quot'][e]['deg'] for r in R),"of",n,
      "mean ratio quot/ctrl",round(m(R,lambda x:x['quot'][e]['deg']/x['ctrl'][e]['deg']),2))
