# PRESPEC 29 figure: observed test-retest, J on days {0,2} vs days {1,3}, per fly.
import json, numpy as np, sqw, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sqw29 import halves
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
S = json.load(open(sqw.OUTD + 's29_hetero.json')); OB = S['obs_flies']; O = S['out']['obs']
P = [('sss', 0, '$sleepless$, light', '#6a3d9a'), ('shmns', 0, '$Sh^{mns}$, light', '#1baf7a'), ('CSDD', 0, 'CS, constant dark', '#6250d6'), ('per0DD', 0, '$per^0$, constant dark', '#eda100')]
fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.2), constrained_layout=True)
for ax, (g, p, ttl, c) in zip(axs, P):
    pr = [h for r in OB if r['g'] == g for h in [halves(r['Jd'], r['scored'], p)] if h is not None]; a, b = np.array(pr).T
    ax.plot(a, b, 'o', ms=3, mfc=c, mec='white', mew=0.3, alpha=0.85); lim = max(np.abs(np.r_[a, b]).max(), 0.1) * 1.05
    ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', lw=0.5); ax.plot([-lim, lim], [-lim, lim], color='0.4', lw=0.5, ls=':')
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal'); k = f"{g}|{'light' if g in ('sss', 'shmns') else 'whole'}"; t = O[k]['T1']
    ax.set_title(f"{ttl}\nretest ρ {t['rho']:+.2f}, p {t['p']:.3f}", fontsize=7.5, loc='left'); ax.set_xlabel('J, days 1 & 3 (cycles h$^{-1}$)' if False else 'J, days 0 & 2')
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
axs[0].set_ylabel('J, days 1 & 3 (cycles h$^{-1}$)')
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig29_retest.{ext}', dpi=300, bbox_inches='tight')
print('saved')
