"""Figures for the report. Validated dataviz palette (light mode, print).
Reads harness JSON outputs; writes SVG+PNG to out/fig/. Deterministic."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
H = os.path.dirname(os.path.abspath(__file__))
OUT = H + '/out/fig'; os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, AQUA, YELLOW, MAG, GREEN, VIOLET, RED = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'
INK, INK2, SURF = '#0b0b0b', '#52514e', '#fcfcfb'
plt.rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans', 'axes.edgecolor': INK2,
 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
 'axes.linewidth': 0.8, 'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF})

def style(ax):
    ax.grid(True, color='#e6e5e1', linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)

def load(name):
    p = f"{H}/out/{name}"
    return json.load(open(p)) if os.path.exists(p) else None

def fig_solving_vs_size():
    svu = load('split_vs_unsplit.json')
    d = load('digest_run1.json')
    if not svu: return
    us = {}; sp = {}
    for r in svu:
        if r.get('status'): continue
        if r['l'] != 2: continue
        us[r['n']] = r['unsplit']['singular']['solving_degree']
        sp[r['n']] = r['split']['singular']['solving_degree']
    ns = sorted(us)
    m2 = [(r['n'], r['solving_degree']) for r in (d or []) if r['kind'] == 'chain' and r.get('m') == 2 and r['vdim'] not in (0, None)]
    m2.sort()
    fig, ax = plt.subplots(figsize=(6.4, 4.0)); style(ax)
    ax.plot(ns, [sp[n] for n in ns], '-o', color=ORANGE, lw=2.2, ms=8,
            label=r'split $S_3$-chain ($m{=}3$, $3\ell{+}n$ vars)', zorder=4, mec=SURF, mew=1.3)
    ax.plot(ns, [us[n] for n in ns], '-o', color=AQUA, lw=2.2, ms=8,
            label=r'unsplit $S_4$ ($m{=}3$, $3\ell$ vars)', zorder=4, mec=SURF, mew=1.3)
    if m2:
        ax.plot([a for a, _ in m2], [b for _, b in m2], '--s', color=BLUE, lw=1.8, ms=6,
                label=r'$S_3$ PDP ($m{=}2$ baseline)', zorder=3, mec=SURF, mew=1.1)
    ax.set_xlabel(r'binary field degree  $n$  ($\ell{=}2$ for the $m{=}3$ lines)')
    ax.set_ylabel(r'solving degree  $d_{solv}$  (Singular F4)')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xticks(sorted(set(ns) | set(a for a, _ in m2)))
    ax.set_title('Splitting $S_4$ into an $S_3$-chain RAISES the solving degree\n'
                 r'binary Koblitz $y^2{+}xy{=}x^3{+}x^2{+}1$; matched $(n,\ell)$',
                 fontsize=10.5, color=INK, loc='left')
    ax.legend(frameon=False, fontsize=9, loc='center left', bbox_to_anchor=(0.02, 0.62))
    fig.tight_layout()
    fig.savefig(OUT + '/fig_solving_vs_size.svg'); fig.savefig(OUT + '/fig_solving_vs_size.png', dpi=160)
    plt.close(fig); print('wrote fig_solving_vs_size')

def fig_ff_gap():
    d = load('digest_run1.json')
    if not d: return
    rows = [r for r in d if r['kind'] in ('chain', 'split') and r.get('m') == 3 and r.get('auxV') in (False, None)
            and r.get('vdim') not in (0, None) and r.get('first_fall') is not None]
    seen = {}
    for r in rows: seen[r['sys_sha']] = r
    rows = sorted(seen.values(), key=lambda r: r['nvars'])
    labels = [f"n={r['n']}\n$\\ell$={r['l']}" for r in rows]
    ff = [r['first_fall'] for r in rows]; sd = [r['solving_degree'] for r in rows]
    x = np.arange(len(rows)); w = 0.38
    fig, ax = plt.subplots(figsize=(6.4, 4.0)); style(ax)
    b1 = ax.bar(x - w/2, ff, w, color=BLUE, label=r'first-fall degree $d_{ff}$', zorder=3)
    b2 = ax.bar(x + w/2, sd, w, color=RED, label=r'solving degree $d_{solv}$', zorder=3)
    ax.bar_label(b1, padding=2, fontsize=9, color=INK2); ax.bar_label(b2, padding=2, fontsize=9, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('degree'); ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim(top=max(sd) + 1.2)
    ax.set_title('First-fall degree stays flat while solving degree grows\n'
                 r'(split $S_3$-chain, $m{=}3$) — the Huang–Kosters–Yeo gap, reproduced',
                 fontsize=10.5, color=INK, loc='left')
    ax.legend(frameon=False, fontsize=9.5, loc='upper left')
    fig.tight_layout()
    fig.savefig(OUT + '/fig_ff_gap.svg'); fig.savefig(OUT + '/fig_ff_gap.png', dpi=160)
    plt.close(fig); print('wrote fig_ff_gap')

def fig_split_engines():
    """Split vs unsplit under BOTH engines: the direction of the split's degree effect
    flips with the engine. Data from engine_compare.json (n=8, ell=2, m=3)."""
    ec = load('engine_compare.json')
    if not ec: return
    def pick(model):
        for r in ec:
            if r['model'] == model and r['n'] == 8 and r['ell'] == 2:
                return r
        return None
    us, sp = pick('unsplit'), pick('split')
    if not us or not sp: return
    groups = ['unsplit $S_4$\n(3$\\ell$ vars)', 'split $S_3$-chain\n(3$\\ell$+n vars)']
    msolve = [us['msolve_f4'], sp['msolve_f4']]
    singular = [us['singular_std'], sp['singular_std']]
    x = np.arange(2); w = 0.36
    fig, ax = plt.subplots(figsize=(6.4, 4.0)); style(ax)
    b1 = ax.bar(x - w/2, msolve, w, color=BLUE, label='msolve F4 (fall reuse)', zorder=3)
    b2 = ax.bar(x + w/2, singular, w, color=ORANGE, label='Singular std', zorder=3)
    ax.bar_label(b1, padding=2, fontsize=10, color=INK2); ax.bar_label(b2, padding=2, fontsize=10, color=INK2)
    # arrows showing direction of change from unsplit->split
    ax.annotate('', xy=(1-w/2, msolve[1]), xytext=(0-w/2, msolve[0]),
                arrowprops=dict(arrowstyle='->', color=BLUE, lw=1.6, ls=':'))
    ax.annotate('', xy=(1+w/2, singular[1]), xytext=(0+w/2, singular[0]),
                arrowprops=dict(arrowstyle='->', color=ORANGE, lw=1.6, ls=':'))
    ax.text(0.5, msolve[0]+0.15, 'splitting LOWERS\nF4 degree (7→4 in survey C)', ha='center', fontsize=8, color=BLUE)
    ax.text(0.5, singular[1]+0.2, 'splitting RAISES\nSingular degree', ha='center', fontsize=8, color=ORANGE)
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=9.5)
    ax.set_ylabel('solving degree reached (same ideal)')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True)); ax.set_ylim(0, max(msolve+singular)+2)
    ax.set_title('The "solving degree" depends on the engine\n'
                 'binary Koblitz, n=8, $\\ell$=2, m=3 — the split helps F4 but not Singular std',
                 fontsize=10.5, color=INK, loc='left')
    ax.legend(frameon=False, fontsize=9.5, loc='upper right')
    fig.tight_layout()
    fig.savefig(OUT + '/fig_split_engines.svg'); fig.savefig(OUT + '/fig_split_engines.png', dpi=160)
    plt.close(fig); print('wrote fig_split_engines')

def fig_fb_distribution():
    """Exotic Rank-1: solving degree vs config, showing ZERO spread over hundreds of V."""
    cons = load('exotic_fb_consolidated.json')
    if not cons: return
    cons = sorted(cons, key=lambda s: (s['model'], s['ell']))
    labels = [f"{s['model']}\nn={s['n']}, $\\ell$={s['ell']}" for s in cons]
    x = np.arange(len(cons))
    mins = [s['deg_min'] for s in cons]; meds = [s['deg_med'] for s in cons]; maxs = [s['deg_max'] for s in cons]
    Ks = [s['K'] + s['n_frob'] for s in cons]
    cols = [ORANGE if s['model'] == 'split' else AQUA if s['model'] == 'unsplit' else BLUE for s in cons]
    fig, ax = plt.subplots(figsize=(6.6, 4.0)); style(ax)
    # draw min-max range (zero height => a dot); overlay median; annotate K and that min=max
    for i in range(len(cons)):
        ax.plot([x[i], x[i]], [mins[i], maxs[i]], color=cols[i], lw=8, solid_capstyle='round', zorder=3, alpha=0.85)
        ax.plot(x[i], meds[i], 'o', color=SURF, ms=5, zorder=4)
        ax.annotate(f"{Ks[i]} subspaces\nall $d_{{solv}}{{=}}{mins[i]}$", (x[i], maxs[i]),
                    xytext=(0, 8), textcoords='offset points', ha='center', fontsize=8, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel(r'solving degree $d_{solv}$  (min / median / max over $V$)')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim(min(mins) - 1, max(maxs) + 1.4)
    ax.set_title('Exotic lever #1 (factor-base subspace $V$): FALSIFIED\n'
                 r'min$_V$=median$_V$=max$_V$ at every config — $d_{solv}$ is intrinsic, not a design knob',
                 fontsize=10.3, color=INK, loc='left')
    fig.tight_layout()
    fig.savefig(OUT + '/fig_fb_distribution.svg'); fig.savefig(OUT + '/fig_fb_distribution.png', dpi=160)
    plt.close(fig); print('wrote fig_fb_distribution')

def fig_tiers(tier_counts=None):
    """Confidence-tier distribution. tier_counts: dict like {'T1':n,...}."""
    if tier_counts is None:
        p = f"{H}/out/tier_counts.json"
        if not os.path.exists(p): return
        tier_counts = json.load(open(p))
    order = ['T1', 'T2', 'T3', 'T4', 'T5']
    labels = {'T1': 'T1 certificate', 'T2': 'T2 replicated', 'T3': 'T3 single run',
              'T4': 'T4 prose/extrap.', 'T5': 'T5 proposal'}
    cols = {'T1': GREEN, 'T2': AQUA, 'T3': BLUE, 'T4': YELLOW, 'T5': INK2}
    vals = [tier_counts.get(t, 0) for t in order]
    fig, ax = plt.subplots(figsize=(6.4, 3.6)); style(ax)
    b = ax.barh([labels[t] for t in order][::-1], vals[::-1], color=[cols[t] for t in order][::-1], zorder=3)
    ax.bar_label(b, padding=3, fontsize=10, color=INK2)
    ax.set_xlabel('number of technique records')
    ax.set_title('Confidence-tier distribution of catalogued techniques', fontsize=11, color=INK, loc='left')
    fig.tight_layout()
    fig.savefig(OUT + '/fig_tiers.svg'); fig.savefig(OUT + '/fig_tiers.png', dpi=160)
    plt.close(fig); print('wrote fig_tiers')

if __name__ == '__main__':
    fig_solving_vs_size(); fig_ff_gap(); fig_fb_distribution(); fig_tiers()
    print("FIGS DONE")
