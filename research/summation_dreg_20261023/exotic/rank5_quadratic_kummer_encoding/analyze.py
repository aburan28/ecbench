import json,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
r1=json.load(open('raw/run1/results.json')); r2=json.load(open('raw/run2/results.json'))
key=lambda x:(x['p'],x['d'],x['m'],x['sys'])
d2={key(x):x for x in r2}; diffs=[]
for x in r1:
    y=d2.get(key(x))
    for f in ('hash','msolve_deg','sing_deg','vdim','bezout_n','nvars'):
        if y is None or x[f]!=y[f]: diffs.append((key(x),f,x[f],None if y is None else y[f]))
print('diff count',len(diffs),diffs, 'rows',len(r1),len(r2))
L=['| p | d | m | sys | vars | eqs | max deg | Bezout(n smallest) | msolve | ms-d | Singular | sing-d | vdim | planted |','|'+'---|'*14]
for x in r1:
    f=lambda v:'n/r' if v is None else v
    L.append(f"| {x['p']} | {x['d']} | {x['m']} | {x['sys']} | {x['nvars']} | {x['neqs']} | {x['maxdeg']} | {x['bezout_n']} | {f(x['msolve_deg'])} | {f(None if x['msolve_deg'] is None else x['msolve_deg']-x['d'])} | {f(x['sing_deg'])} | {f(None if x['sing_deg'] is None else x['sing_deg']-x['d'])} | {f(x['vdim'])} | {x.get('planted_ok')} |")
open('table.md','w').write('\n'.join(L)); print('\n'.join(L))
# figure: offset (d_solv-d) m=3
cfg=sorted({(x['p'],x['d']) for x in r1 if x['m']==3}); sy=['a','b','c1','c2']; col=['#2a78d6','#eb6834','#1baf7a','#eda100']
fig,ax=plt.subplots(1,2,figsize=(11,4))
for k,(eng,fld) in enumerate([('msolve F4','msolve_deg'),('Singular degBound','sing_deg')]):
    for i,s in enumerate(sy):
        ys=[]
        for c in cfg:
            row=[x for x in r1 if (x['p'],x['d'])==c and x['m']==3 and x['sys']==s]
            v=row[0][fld] if row else None; ys.append(float('nan') if v is None else v-c[1])
        ax[k].bar([j+0.2*i-0.3 for j in range(len(cfg))],ys,0.2,color=col[i],label=s)
    ax[k].set_xticks(range(len(cfg))); ax[k].set_xticklabels([f'p={p},d={d}' for p,d in cfg],fontsize=8,color='#0b0b0b')
    ax[k].set_title(eng+' (m=3)',color='#0b0b0b'); ax[k].set_ylabel('d_solv - d (offset)',color='#0b0b0b')
ax[0].legend(title='system')
plt.tight_layout(); plt.savefig('fig_rank5_quadratic_kummer_encoding.png',dpi=140); plt.savefig('fig_rank5_quadratic_kummer_encoding.svg')
