# PRESPEC 30 figure: (a) long-bout (>= 25 min) sleep per day, per fly and median, deprivation night marked;
# (b) per-fly change after deprivation (day 4 - days 0-2) against the baseline drift (day 2 - days 0-1).
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
S = json.load(open('sqw_out/s30_longbouts.json')); F = S['flies']; days = range(9)
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8), gridspec_kw=dict(width_ratios=[1.6, 1]), constrained_layout=True)
ax = axs[0]
for f in F: ax.plot(days, [f[str(d)]['L25'] / 60 for d in days], color='#2a78d6', lw=0.4, alpha=0.3)
for key, c, lab in (('L25', '#1f4e79', 'sleep in bouts ≥ 25 min'), ('all', '0.45', 'all sleep (bouts ≥ 5 min)')):
    ax.plot(days, [np.median([f[str(d)][key] for f in F]) / 60 for d in days], 'o-', color=c, lw=1.6, ms=3.5, label=lab)
ax.axvspan(2.6, 3.4, color='#eb6834', alpha=0.15, lw=0); ax.text(3, 0.02, 'deprivation\nZT15–24', transform=ax.get_xaxis_transform(), ha='center', fontsize=6.5, color='#b5491f')
ax.set_xticks(list(days)); ax.set_xlabel('day of recording'); ax.set_ylabel('sleep (h per day)'); ax.legend(fontsize=6.5, frameon=False, loc='lower right')
ax.set_title('Long-bout sleep climbs every day, deprived or not', fontsize=8, loc='left')
ax = axs[1]
D = [f['4']['L25'] - np.mean([f[str(d)]['L25'] for d in (0, 1, 2)]) for f in F]; Dp = [f['2']['L25'] - np.mean([f[str(d)]['L25'] for d in (0, 1)]) for f in F]
ax.plot(np.array(Dp) / 60, np.array(D) / 60, 'o', ms=3.5, mfc='#2a78d6', mec='white', mew=0.4); lim = max(np.abs(np.r_[D, Dp]).max() / 60, 1) * 1.05
ax.plot([-lim, lim], [-lim, lim], color='0.4', lw=0.6, ls=':'); ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', lw=0.5)
p = S['primary']; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
ax.set_xlabel('baseline drift, day 2 − days 0–1 (h)'); ax.set_ylabel('after deprivation, day 4 − days 0–2 (h)')
ax.set_title(f"above the diagonal = rebound beyond drift\nexcess median {p['excess_med'] / 60:+.2f} h, p {p['p_excess']:.2f}", fontsize=7, loc='left')
for a in axs:
    for s in ('top', 'right'): a.spines[s].set_visible(False)
for a, l in zip(axs, 'ab'): a.text(-0.12, 1.04, l, transform=a.transAxes, fontsize=10, fontweight='bold')
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig30_longbouts.{ext}', dpi=300, bbox_inches='tight')
print('saved')
