# PRESPEC 24 figure: per-fly cycle flux J (count J, cycles/h), 23E10>Kir2.1 vs 23E10>Kir2.1-NC, light and dark.
import json, numpy as np, sqw, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
S = json.load(open(sqw.OUTD + 's24_23E10.json')); F = [d for d in S['flies'] if not d['excluded']]; rng = np.random.default_rng(24)
COL = {'nc': '#8a8a8a', 'kir': '#2a78d6'}; LAB = {'nc': '23E10 > Kir2.1-NC\n(control)', 'kir': '23E10 > Kir2.1\n(silenced)'}
fig, ax = plt.subplots(figsize=(3.6, 3.0), constrained_layout=True); x0 = 0; ticks = []
for p, ph in enumerate(['light', 'dark']):
    for k, g in enumerate(['nc', 'kir']):
        v = np.array([d['main'][p]['J'] for d in F if d['g'] == g and d['main'][p]['ok']]); x = x0 + 0.9 * k
        ax.plot(x + rng.uniform(-0.22, 0.22, len(v)), v, 'o', ms=2.4, mfc=COL[g] if p else 'white', mec=COL[g], mew=0.6, alpha=0.8)
        bs = [np.median(rng.choice(v, len(v))) for _ in range(2000)]; lo, hi = np.percentile(bs, [2.5, 97.5]); m = np.median(v)
        ax.plot([x - 0.3, x + 0.3], [m, m], 'k', lw=1.3); ax.plot([x, x], [lo, hi], 'k', lw=0.8)
    ticks.append(x0 + 0.45); x0 += 2.4
pr = S['results']['main']['primary']
ax.text(x0 - 2.4 + 0.45, 0.97, f"night: one-sided p {pr['p']:.3f}\n(not significant)", transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=6.5)
ax.axhline(0, color='0.6', lw=0.5); ax.set_xticks(ticks); ax.set_xticklabels(['light', 'dark']); ax.set_xlim(-0.6, x0 - 1.0)
ax.set_ylabel('cycle flux $J$ (cycles h$^{-1}$)\n+ forward W→Q→S→W')
for g in ['nc', 'kir']: ax.plot([], [], 's', color=COL[g], label=LAB[g].replace('\n', ' '))
ax.legend(fontsize=6, frameon=False, loc='lower left')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.set_title('Silencing 23E10 neurons: cycle flux per fly', fontsize=8, loc='left', pad=6)
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig24_23E10_J.{ext}', dpi=300, bbox_inches='tight')
print('saved')
