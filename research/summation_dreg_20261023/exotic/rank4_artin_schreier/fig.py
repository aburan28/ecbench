import json, os, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(f"{HERE}/raw/summary.json"))['rows']
def num(x):
    try: return int(str(x).rstrip('*'))
    except: return None
C = dict(full='#2a78d6', as_='#eb6834', ctl='#1baf7a'); T = '#0b0b0b'
plt.rcParams.update({'text.color': T, 'axes.labelcolor': T, 'xtick.color': T, 'ytick.color': T, 'font.size': 10})
cats = ['m2 single', 'm3 unsplit', 'm3 split']
fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey='row')
for ei, (eng, kf, ka, kc) in enumerate([('Singular std', 'fullS', 'asS', 'ctlS'), ('msolve F4', 'fullM', 'asM', 'ctlM')]):
    for ci, (m, model) in enumerate([(2, 'single'), (3, 'unsplit'), (3, 'split')]):
        ax = axes[ei][ci]
        rows = [r for r in S if r['m'] == m and r['model'] == model]
        labels = [f"n{r['n']}\nl{r['l']}" for r in rows]; xs = range(len(rows)); w = 0.27
        for k, (key, name, col) in enumerate([(kf, 'FULL', C['full']), (ka, 'AS residual (max branch)', C['as_']), (kc, 'CONTROL (non-trace fix)', C['ctl'])]):
            ys = [num(r[key]) or 0 for r in rows]
            ax.bar([x + (k - 1) * w for x in xs], ys, w, color=col, label=name if (ei == 0 and ci == 0) else None)
        ax.set_xticks(list(xs)); ax.set_xticklabels(labels, fontsize=7)
        ax.set_title(f"{eng}: {cats[ci]}", fontsize=10); ax.spines[['top', 'right']].set_visible(False)
        if ci == 0: ax.set_ylabel('solving degree d_solv')
fig.legend(loc='lower center', ncol=3, frameon=False)
fig.suptitle('AS trace-block elimination: residual lowers d_solv exactly as much as a non-trace control', fontsize=11)
fig.tight_layout(rect=(0, 0.04, 1, 0.96))
fig.savefig(f"{HERE}/fig_rank4_artin_schreier.png", dpi=150); fig.savefig(f"{HERE}/fig_rank4_artin_schreier.svg")
