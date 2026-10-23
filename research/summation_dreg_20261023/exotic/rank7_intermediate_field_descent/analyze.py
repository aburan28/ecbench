import json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
r1=json.load(open('raw/run1.json')); r2=json.load(open('raw/run2.json'))
print("twice-run identical (hashes, degrees, vdim, everything):", r1==r2)
rows=[]
for r in r1:
    if 'status' in r: print("SKIP",r['n'],r['k'],r['m'],r['a6'],r['status']); continue
    rows.append(r)
print("| n | k | m | curve | brute #sol | A0 eqdeg / Bezout / vdim | A0 Sing [k:d]/label | A1 vdim | B vars | B msolve | B Sing [k:d]/label | B vdim | planted |")
for r in rows:
    A0,A1,Bs,Bm=r['A0'],r['A1'],r['B_singular'],r['B_msolve']
    print(f"| {r['n']} | {r['k']} | {r['m']} | {r['a6']} | {r['brute_count']} | {A0['eq_degrees']} / {A0['bezout']} / {A0['vdim']} | {A0['solving_degree']}/{A0['max_label_degree']} | {A1['vdim']} | {r['B_nvars']} | {Bm['solving_degree']} | {Bs['solving_degree']}/{Bs['max_label_degree']} | {Bs['vdim']} | {A0['planted_ok'] and A1['planted_ok'] and r['B_planted_ok']} |")
    assert A1['vdim']==Bs['vdim']==r['brute_count'] and r['numeric_eq_brute'], r
print("count gate OK")
# figure: Singular label degree vs k for m=2,3; A0 vs B(F2)
pal={'A0':'#2a78d6','B':'#eb6834','msolve':'#1baf7a'}
fig,axs=plt.subplots(1,2,figsize=(10,4.2))
for ax,m in zip(axs,(2,3)):
    R=[r for r in rows if r['m']==m and r['a6']=='generic']
    ks=[r['k'] for r in R]
    ax.plot(ks,[r['A0']['max_label_degree'] for r in R],'o-',color=pal['A0'],label='GF(q)-descent, Singular (no field eqs)')
    ax.plot(ks,[r['B_singular']['max_label_degree'] for r in R],'s-',color=pal['B'],label='F2-descent, Singular (k*m bits)')
    ax.plot(ks,[r['B_msolve']['solving_degree'] for r in R],'^--',color=pal['msolve'],label='F2-descent, msolve F4')
    ax.set_title(f'm = {m}  (S_{m+1}), generic a6',color='#0b0b0b'); ax.set_xlabel('k  (q = 2^k)',color='#0b0b0b'); ax.set_ylabel('max degree reached',color='#0b0b0b')
    ax.set_xticks(ks); ax.grid(alpha=.25); ax.set_ylim(0,None)
axs[0].legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig_rank7_intermediate_field_descent.png',dpi=150); plt.savefig('fig_rank7_intermediate_field_descent.svg')
