import json, statistics as st, collections
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
r1=json.load(open('raw/run1.json')); r2=json.load(open('raw/run2.json'))
def key(r): return (r['p'],r['family'],r['k'],r['d'],r['seed'])
d2={key(r):r for r in r2}
diff=[]
for r in r1:
    o=d2[key(r)]
    if 'error' in r or 'error' in o: 
        if r!=o: diff.append((key(r),'error'))
        continue
    a=(r['hash'],r['msolve'].get('deg'),r['singular'].get('deg'),r['singular'].get('vdim'),r['Fsize'],r['avg_decomp'])
    b=(o['hash'],o['msolve'].get('deg'),o['singular'].get('deg'),o['singular'].get('vdim'),o['Fsize'],o['avg_decomp'])
    if a!=b: diff.append((key(r),a,b))
out=[]
out.append(f"twice-run diff: {len(r1)} systems, mismatches = {len(diff)}  {diff[:3]}")
errs=[r for r in r1 if 'error' in r]; out.append(f"errors: {len(errs)} {[ (key(e),e['error']) for e in errs[:3]]}")
ok=[r for r in r1 if 'error' not in r]
bad=[key(r) for r in ok if not r['root_check'] or r['singular'].get('plant_not_eliminated') is not True]
out.append(f"planted-check failures: {bad}")
nt=[key(r) for r in ok if r['msolve']['status']!='ok' or r['singular']['status']!='ok']
out.append(f"timeouts/fails: {nt}")
G=collections.defaultdict(list)
for r in ok: G[(r['p'],r['k'],r['d'],r['family'])].append(r)
def med(x): x=[v for v in x if v is not None]; return st.median(x) if x else None
rows=[]
out.append("\n| p | k,d | family | vars | Bezout | s | |F| | decomp/randbase | msolve deg (seeds) | Singular deg (seeds) | vdim (seeds) |\n|---|---|---|---|---|---|---|---|---|---|---|")
for (p,k,d,f),L in sorted(G.items(),key=lambda t:(t[0][0],t[0][1],t[0][2],['generic','kummer','isogeny','lattes'].index(t[0][3]))):
    ms=[x['msolve']['deg'] for x in L]; sg=[x['singular']['deg'] for x in L]; vd=[x['singular'].get('vdim') for x in L]
    ratio=sum(x['avg_decomp'] for x in L)/max(1e-9,sum(x['avg_decomp_randbase'] for x in L))
    out.append(f"| {p} | {k},{d} | {f} | {L[0]['nvars']} | {L[0]['bezout']} | {L[0]['s']} | {med([x['Fsize'] for x in L])} | {ratio:.2f} | {ms} | {sg} | {vd} |")
# falsifier
out.append("\nFalsifier comparison (structured median vs generic median, per config):")
res=collections.Counter(); detail=[]
for (p,k,d,f),L in G.items():
    if f not in('lattes','isogeny'): continue
    g=G[(p,k,d,'generic')]
    for eng,get in (('msolve',lambda x:x['msolve']['deg']),('singular',lambda x:x['singular']['deg'])):
        s_=med([get(x) for x in L]); g_=med([get(x) for x in g]); gr=(min(get(x) for x in g),max(get(x) for x in g))
        v='equal' if s_==g_ else ('lower' if s_<g_ else 'higher')
        outside = s_<gr[0] or s_>gr[1]
        res[(f,eng,v,'outside_generic_range' if outside else 'within')]+=1
        detail.append((p,k,d,f,eng,s_,g_,gr,v,outside))
for t in detail: out.append(str(t))
out.append(str(dict(res)))
open('analysis.txt','w').write("\n".join(out)); print("\n".join(out))
# figure: Singular degree vs vdim for each family
pal={'generic':'#0b0b0b','kummer':'#eda100','isogeny':'#1baf7a','lattes':'#eb6834'}
fig,ax=plt.subplots(1,2,figsize=(11,4.4))
for fam,c in pal.items():
    xs=[x['singular'].get('vdim') for x in ok if x['family']==fam]; ys=[x['singular']['deg'] for x in ok if x['family']==fam]
    ax[0].scatter(xs,ys,c=c,label=fam,s=26,alpha=.8)
    ax[1].scatter([x['msolve']['deg'] for x in ok if x['family']==fam],[x['singular']['deg'] for x in ok if x['family']==fam],c=c,label=fam,s=26,alpha=.8)
ax[0].set_xlabel('vdim (Singular)');ax[0].set_ylabel('Singular d_solv');ax[0].set_title('Singular degree vs solution count',color='#0b0b0b')
ax[1].set_xlabel('msolve F4 max degree');ax[1].set_ylabel('Singular d_solv');ax[1].set_title('Engines compared',color='#0b0b0b')
ax[0].legend(); 
for a in ax: a.grid(alpha=.25)
plt.tight_layout(); plt.savefig('fig_rank6_lattes_isogeny.png',dpi=150); plt.savefig('fig_rank6_lattes_isogeny.svg')
