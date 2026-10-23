import json, statistics as st, collections
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
r1=json.load(open('raw/run1.json')); r2=json.load(open('raw/run2.json'))
def key(row):
    out={}
    for k,s in row['sys'].items():
        out[k]=(s['hash'],s['sing'].get('dsolv'),s['sing'].get('vdim'),s['sing'].get('planted'),s['msolve_dsolv'],s['msolve_planted'],s['bezout'],tuple(s['eqdeg']))
    return (row['meta']['p'],row['meta']['seed'],row['meta']['m'],tuple(sorted(out.items())))
d1=[key(r) for r in r1 if 'sys' in r]; d2=[key(r) for r in r2 if 'sys' in r]
same=(d1==d2)
print('twice-run diff identical:',same,len(d1),len(d2),'errors',[r for r in r1 if 'error' in r])
SYS=['a','b','c1','c1p','c2']
agg=collections.defaultdict(lambda: collections.defaultdict(list))
for r in r1:
    if 'sys' not in r: continue
    m=r['meta']['m']
    for k in SYS:
        s=r['sys'][k]; eff_s=max(s['sing']['dsolv'] or 0,max(s['eqdeg']))
        agg[m][k].append(dict(sing=s['sing']['dsolv'],sing_eff=eff_s,ms=s['msolve_dsolv'],ms_eff=max(s['msolve_dsolv'] or 0,0),vdim=s['sing']['vdim'],bez=s['bezout'],vars=s['vars'],maxdeg=max(s['eqdeg']),
                                  rec=s['sing']['planted'] and s['msolve_planted']))
lines=[]
for m in (2,3):
    lines.append(f'm=n={m}: instances={len(agg[m]["a"])}')
    lines.append('| system | vars | max eq deg | Bezout | vdim | Singular d (min..max) | Singular d_eff=max(d,eqdeg) | msolve d (min..max) | planted recovered |')
    lines.append('|---|---|---|---|---|---|---|---|---|')
    for k in SYS:
        L=agg[m][k]
        rng=lambda f:(lambda v:f'{min(v)}..{max(v)}' if v else 'n/a')([x[f] for x in L if x[f] is not None])
        lines.append(f'| {k} | {L[0]["vars"]} | {rng("maxdeg")} | {rng("bez")} | {rng("vdim")} | {rng("sing")} | {rng("sing_eff")} | {rng("ms")} | {all(x["rec"] for x in L)} |')
    # pairwise per-instance comparison c vs b
    for c in ('c1','c1p'):
        for eng in ('sing','ms'):
            cmpv=[]
            for xb,xc in zip(agg[m]['b'],agg[m][c]):
                vb=xb[eng] if xb[eng] is not None else xb['maxdeg']; vc=xc[eng] if xc[eng] is not None else xc['maxdeg']
                cmpv.append('<' if vc<vb else ('=' if vc==vb else '>'))
            lines.append(f'  {c} vs b, {eng}: '+str(dict(collections.Counter(cmpv))))
open('raw/summary_table.txt','w').write('\n'.join(lines)); print('\n'.join(lines))
# expansion factors
vd={k:agg[3][k][0]['vdim'] for k in SYS}; print('vdim m=3',vd)
vd2={k:agg[2][k][0]['vdim'] for k in SYS}; print('vdim m=2',vd2)
# figure
cols=['#2a78d6','#eb6834','#1baf7a','#eda100']; txt='#0b0b0b'
fig,axs=plt.subplots(1,2,figsize=(10,4.2),sharey=False)
names=[('a','plain Y'),('b','S_m inv e_k(Y)'),('c1','D_m inv (u,pi), u_m=pi^2'),('c1p','D_m inv + relation')]
for ax,eng,lab in ((axs[0],'sing','Singular std (min degBound)'),(axs[1],'ms','msolve F4 (max round degree)')):
    w=0.2
    for i,(k,nm) in enumerate(names):
        mx=[];mn=[]
        vals=[]
        for x in agg[3][k]:
            vals.append(x[eng] if x[eng] is not None else x['maxdeg'])
        ax.bar(i,st.mean(vals),color=cols[i],width=0.7)
        ax.plot([i,i],[min(vals),max(vals)],color=txt,lw=1.5)
        ax.text(i,max(vals)+0.4,f'{st.mean(vals):.1f}',ha='center',color=txt,fontsize=9)
    ax.set_xticks(range(4)); ax.set_xticklabels([n for _,n in names],rotation=18,ha='right',color=txt,fontsize=8)
    ax.set_title(lab,color=txt,fontsize=10); ax.set_ylabel('solving degree, m=n=3 (mean, range)',color=txt,fontsize=9)
    ax.set_yscale('log'); ax.tick_params(colors=txt); 
    for s in ('top','right'): ax.spines[s].set_visible(False)
fig.suptitle('Rank 3: signed-permutation invariants, 3-summand decomposition over F_{p^3}',color=txt,fontsize=11)
fig.tight_layout(); fig.savefig('fig_rank3_reflection_group_invariants.png',dpi=140); fig.savefig('fig_rank3_reflection_group_invariants.svg')
