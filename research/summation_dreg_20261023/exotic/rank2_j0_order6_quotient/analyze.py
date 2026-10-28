import json, collections
R=json.load(open('run1/results.json'))
rows=[]
for r in R:
    P,Q=r['plain'],r['quot']
    rows.append(dict(p=r['p'],d=r['d'],seed=r['seed'],Np=r['N_plain'],Nq=r['N_quot'],vdP=P['singular']['vdim'],vdQ=Q['singular']['vdim'],radP=P['singular']['rad'],radQ=Q['singular']['rad'],
      sP=P['singular']['deg'],sQ=Q['singular']['deg'],mP=P['msolve']['deg'],mQ=Q['msolve']['deg'],BciP=P['B_ci'],BciQ=Q['B_ci'],BaP=P['B_all'],BaQ=Q['B_all'],
      ok=r['expansion_ok'] and r['planted_ok'] and r['nine_factor_ok'] and P['singular']['rad']==r['N_plain'] and Q['singular']['rad']==r['N_quot']))
n=len(rows)
print("n",n,"all correctness ok:",all(x['ok'] for x in rows), "bad:",[ (x['p'],x['d'],x['seed']) for x in rows if not x['ok']])
print("p d s | Np Nq | vdimP vdimQ | sing P/Q | msolve P/Q | Bci P/Q (ratio)")
for x in rows: print(x['p'],x['d'],x['seed'],'|',x['Np'],x['Nq'],'|',x['vdP'],x['vdQ'],'|',x['sP'],x['sQ'],'|',x['mP'],x['mQ'],'|',x['BciP'],x['BciQ'],round(x['BciP']/x['BciQ'],2))
lower_m=sum(x['mQ']<x['mP'] for x in rows); lower_s=sum(x['sQ']<x['sP'] for x in rows)
ge_both=sum(x['mQ']>=x['mP'] and x['sQ']>=x['sP'] for x in rows)
b2=sum(x['BciP']/x['BciQ']>=2 for x in rows)
print("quot<plain: msolve",lower_m,"/",n," singular",lower_s,"/",n," ; quot>=plain on both:",ge_both,"; Bci ratio>=2:",b2)
by=collections.defaultdict(list)
for x in rows: by[x['d']].append(x)
print("by d: [d: n, msolve P->Q mean, singular P->Q mean]")
for d in sorted(by): 
    L=by[d]; print(d,len(L),sum(x['mP'] for x in L)/len(L),sum(x['mQ'] for x in L)/len(L),sum(x['sP'] for x in L)/len(L),sum(x['sQ'] for x in L)/len(L))
print("Np/Nq total",sum(x['Np'] for x in rows),sum(x['Nq'] for x in rows))
json.dump(rows,open('raw/summary_table.json','w'),indent=1)
# figure
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11,'text.color':'#0b0b0b','axes.labelcolor':'#0b0b0b','xtick.color':'#0b0b0b','ytick.color':'#0b0b0b'})
fig,ax=plt.subplots(1,2,figsize=(10,4.2))
for a,(kp,kq,t) in zip(ax,[('mP','mQ','msolve F4'),('sP','sQ','Singular std')]):
    a.scatter([x['d'] for x in rows],[x[kp] for x in rows],c='#2a78d6',label='plain (x1,x2)',s=28)
    a.scatter([x['d'] for x in rows],[x[kq] for x in rows],c='#eb6834',marker='s',label='quotient (u1,u2)',s=28)
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel('d (factor base size |V|)'); a.set_ylabel('max degree reached'); a.set_title(t); a.legend(frameon=False)
    a.grid(alpha=.25)
fig.suptitle('j=0 order-6 quotient, m=2: plain vs quotient solving degree',fontsize=12)
fig.tight_layout(); fig.savefig('fig_rank2_j0_order6_quotient.png',dpi=150); fig.savefig('fig_rank2_j0_order6_quotient.svg')
