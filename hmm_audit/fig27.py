# PRESPEC 27 figure: (a,b) pooled 2-h bin cycle flux, full and with only one state's exits following the clock
# (CS DD, CS LD); (c) daily rhythm of each exit probability, relative amplitude, CS DD vs per0 DD.
import json, numpy as np, sqw, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
D = json.load(open(sqw.OUTD + 's27_decomp.json')); t = 2 * np.arange(12) + 1
C = {'Q': '#1baf7a', 'W': '#eb6834', 'S': '#6250d6'}; L = {'Q': 'only quiet-wake exits', 'W': 'only active-wake exits', 'S': 'only sleep exits'}
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), constrained_layout=True)
for ax, g, ttl in ((axs[0], 'CSDD', 'CS, constant darkness'), (axs[1], 'CS', 'CS, light–dark')):
    p = D[g]['point']; ax.plot(t, p['full']['J'], 'k-o', ms=3, lw=1.6, label='all exits follow the clock (data)')
    for s in 'QWS': ax.plot(t, p[s]['J_only'], '-', color=C[s], lw=1, label=L[s])
    if g == 'CS': ax.axvspan(12, 24, color='0.92', lw=0, zorder=0)
    ax.set_xlim(0, 24); ax.set_xticks([0, 6, 12, 18, 24]); ax.set_xlabel('clock time (h)'); ax.set_title(ttl, fontsize=8.5, loc='left')
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
axs[0].set_ylabel('cycle flux $J$ (cycles h$^{-1}$)'); axs[0].legend(fontsize=5.8, frameon=False, loc='upper left')
ax = axs[2]; keys = list(D['CSDD']['tp']); x = np.arange(len(keys))
for k, (g, col, off) in enumerate((('CSDD', '#6250d6', -0.2), ('per0DD', '#eda100', 0.2))):
    v = [D[g]['tp'][kk]['rel_amp'] for kk in keys]; lo = [D[g]['tp'][kk]['ci'][0] for kk in keys]; hi = [D[g]['tp'][kk]['ci'][1] for kk in keys]
    ax.bar(x + off, v, 0.38, color=col, label={'CSDD': 'CS DD', 'per0DD': '$per^0$ DD'}[g])
    ax.errorbar(x + off, v, yerr=[np.clip(np.array(v) - lo, 0, None), np.clip(np.array(hi) - v, 0, None)], fmt='none', ecolor='k', lw=0.6)
lab = {'Q': 'Q', 'W': 'W', 'S': 'S'}
ax.set_xticks(x); ax.set_xticklabels([k.replace('->', '→') for k in keys], fontsize=7); ax.set_ylabel('daily rhythm of the exit probability\n(amplitude / daily mean)')
ax.legend(fontsize=6.5, frameon=False); ax.set_title('which transitions follow the clock', fontsize=8.5, loc='left')
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
for a, l in zip(axs, 'abc'): a.text(-0.16, 1.04, l, transform=a.transAxes, fontsize=10, fontweight='bold')
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig27_clock_decomposition.{ext}', dpi=300, bbox_inches='tight')
print('saved')
