# PRESPEC 23 figure: per-fly observed-chain J (cycle flies), day and night: fumin, sss, Sh^mns (CS for reference).
import json, numpy as np, sqw, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
O = json.load(open(sqw.OUTD + 's_observable.json')); S = json.load(open(sqw.OUTD + 's23_shmns.json')); rng = np.random.default_rng(23)
def J(g, p):
    if g == 'shmns': R = [d['main'][p] for d in S['flies'] if not d['excluded']]
    else: R = [d['main'][p] for d in O if d['g'] == g]
    return np.array([r['J_chain'] for r in R if r['ok'] and r['J_chain'] is not None and abs(r['J_chain']) >= 1e-12])
G = [('CS', 'CS', '#2a78d6'), ('fumin', '$fumin$', '#d43d51'), ('sss', '$sleepless$', '#6a3d9a'), ('shmns', '$Sh^{mns}$', '#1baf7a')]
fig, axs = plt.subplots(1, 2, figsize=(6.0, 2.8), sharey=True, constrained_layout=True)
for ax, p, ttl in ((axs[0], 0, 'light phase'), (axs[1], 1, 'dark phase')):
    for k, (g, lab, c) in enumerate(G):
        v = J(g, p); ax.plot(k + rng.uniform(-0.22, 0.22, len(v)), v, 'o', ms=2.4, mfc=c if p else 'white', mec=c, mew=0.6, alpha=0.8)
        bs = [np.median(rng.choice(v, len(v))) for _ in range(2000)]; lo, hi = np.percentile(bs, [2.5, 97.5]); m = np.median(v)
        ax.plot([k - 0.3, k + 0.3], [m, m], 'k', lw=1.3); ax.plot([k, k], [lo, hi], 'k', lw=0.8)
        ax.text(k, 1.0, f"{np.mean(v < 0):.0%}\nbackward", transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=6)
    ax.axhline(0, color='0.5', lw=0.6); ax.set_xticks(range(len(G))); ax.set_xticklabels([g[1] for g in G]); ax.set_title(ttl, fontsize=8.5, loc='left', pad=20)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
axs[0].set_ylabel('cycle flux $J$ (cycles h$^{-1}$)\n+ forward W→Q→S→W'); axs[0].set_ylim(-0.6, 1.0)
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig23_shmns.{ext}', dpi=300, bbox_inches='tight')
print('saved')
