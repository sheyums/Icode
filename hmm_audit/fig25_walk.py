# Descriptive companion to PRESPEC 25 (not a test): each fly's cumulative net forward cycles,
# (#sleep->active wake) - (#active wake->sleep), against time, and the per-window net counts.
import json, numpy as np, sqw, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sqw15 import seqs
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
G = ['CS', 'CSDD', 'per0LD', 'per0DD', 'sss']; LD = {'CS', 'per0LD', 'sss'}
NAME = {'CS': 'CS, light–dark', 'CSDD': 'CS, constant dark', 'per0LD': '$per^0$, light–dark', 'per0DD': '$per^0$, constant dark', 'sss': '$sleepless$, light–dark'}
COL = {'CS': '#2a78d6', 'CSDD': '#6250d6', 'per0LD': '#1baf7a', 'per0DD': '#eda100', 'sss': '#6a3d9a'}
D = json.load(open(sqw.OUTD + 's25_windows.json'))
fig = plt.figure(figsize=(7.2, 5.4)); gs = fig.add_gridspec(2, 5, height_ratios=[1.5, 1], hspace=0.45, wspace=0.35)
stats = {}
for k, g in enumerate(G):
    ax = fig.add_subplot(gs[0, k]); tr = []
    for i in range(len(sqw.files(g))):
        y, _, _ = seqs(g, i); y = np.asarray(y); ok = (y[:-1] >= 0) & (y[1:] >= 0)
        step = np.where(ok & (y[:-1] == 2) & (y[1:] == 1), 1, 0) - np.where(ok & (y[:-1] == 1) & (y[1:] == 2), 1, 0)
        if np.mean(y >= 0) < 0.6: continue
        c = np.r_[0, np.cumsum(step)]; t = np.arange(len(c)) / 60; tr.append((t, c))
        ax.plot(t, c, color=COL[g], lw=0.5, alpha=0.35)
    T = min(len(c) for _, c in tr); M = np.array([c[:T] for _, c in tr]); tt = np.arange(T) / 60
    ax.plot(tt, M.mean(0), color='k', lw=1.6)
    if g in LD:
        for d in range(int(tt[-1] // 24) + 1): ax.axvspan(24 * d + 12, 24 * d + 24, color='0.9', lw=0, zorder=0)
    ax.axhline(0, color='0.5', lw=0.5); ax.set_xlim(0, tt[-1]); ax.set_xticks(np.arange(0, tt[-1] + 1, 24))
    ax.set_title(NAME[g], fontsize=7.5); ax.set_xlabel('hours')
    if k == 0: ax.set_ylabel('cumulative net forward cycles\n(one line per fly; black: mean)')
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ends = M[:, -1]; stats[g] = dict(n=len(M), hours=float(tt[-1]), end_mean=float(ends.mean()), frac_end_pos=float(np.mean(ends > 0)), frac_end_neg=float(np.mean(ends < 0)))
    ax2 = fig.add_subplot(gs[1, k]); net = []
    for f in D[f'{g}|6h']['flies']:
        net += [a - b for c, a, b in zip(f['cls'], f['nSW'], f['nWS']) if c != 9]
    net = np.clip(np.array(net), -6, 6); bins = np.arange(-6.5, 7.5, 1); h, _ = np.histogram(net, bins); h = h / h.sum(); x = np.arange(-6, 7)
    ax2.bar(x, h, width=0.85, color=[('#eb6834' if v < 0 else ('#2a78d6' if v > 0 else '0.75')) for v in x])
    ax2.set_xlabel('net forward cycles\nin one 6-h window'); ax2.set_xticks([-6, -3, 0, 3, 6]); ax2.set_xticklabels(['≤−6', '−3', '0', '3', '≥6'])
    if k == 0: ax2.set_ylabel('fraction of windows')
    ax2.text(0.97, 0.95, f"fwd {np.mean(net > 0):.0%}\nback {np.mean(net < 0):.0%}", transform=ax2.transAxes, ha='right', va='top', fontsize=6.5)
    for s in ('top', 'right'): ax2.spines[s].set_visible(False)
fig.text(0.01, 0.97, 'a', fontsize=10, fontweight='bold'); fig.text(0.01, 0.40, 'b', fontsize=10, fontweight='bold')
fig.text(0.5, -0.04, 'Grey shading: dark phase. Each fly takes small forward and backward steps; the forward bias appears only as the steps accumulate over days.', ha='center', fontsize=6.5)
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig25_cumulative_walk.{ext}', dpi=300, bbox_inches='tight')
print(json.dumps(stats, indent=0))
