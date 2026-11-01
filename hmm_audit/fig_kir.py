# Fig C12 "Silencing sleep circuits makes sleep shallow without changing the cycle's drive" (PRESPEC 20y; y_kir.json).
# Per fly; light (open) and dark (filled); median with fly-bootstrap 95% CI; '*' = differs from control (Holm < 0.05,
# flags copied from the 20y ledger entry, not recomputed here).
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})
F = [d for d in json.load(open('y_kir.json')) if not d['dead']]; rng = np.random.default_rng(1)
G = ['empty', 'dfb', 'vnc']; COL = {'empty': '#8a8a8a', 'dfb': '#2a78d6', 'vnc': '#eb6834'}
LAB = {'empty': 'control\n(empty>Kir)', 'dfb': 'dFB\nsilenced', 'vnc': 'VNC\nsilenced'}
PANELS = [  # key, transform, y label, title, Holm-significant (group, phase) from the ledger
    ('sleep', lambda r: r['sleep'], 'sleep (min per hour)', 'Sleep amount', {('vnc', 0), ('dfb', 1), ('vnc', 1)}),
    ('QS', lambda r: r['P']['QS'], 'P(quiet wake → sleep) per min', 'Falling asleep', {('vnc', 0), ('dfb', 1), ('vnc', 1)}),
    ('SW', lambda r: r['P']['SW'], 'P(sleep → active wake) per min', 'Waking up', {('dfb', 0), ('vnc', 0), ('dfb', 1), ('vnc', 1)}),
    ('E31', lambda r: -r['E31'], 'depth of the sleep well,\n$E_Q - E_S$ (nats)', 'Sleep vs quiet wake (occupancy)', {('dfb', 0), ('vnc', 0), ('dfb', 1), ('vnc', 1)}),
    ('J', lambda r: r['J'], 'cycle flux $J$ (cycles h$^{-1}$)', 'Drive around W→Q→S→W', set())]
def mci(v):
    v = np.asarray(v, float); b = [np.median(rng.choice(v, len(v))) for _ in range(2000)]; return np.median(v), *np.percentile(b, [2.5, 97.5])
fig, axs = plt.subplots(1, 5, figsize=(7.2, 2.7), constrained_layout=True)
for ax, (key, f, yl, title, sig) in zip(axs, PANELS):
    x0 = 0; ticks = []
    for p in (0, 1):
        for k, g in enumerate(G):
            x = x0 + k * 0.9
            v = [f(d['main'][p]) for d in F if d['g'] == g and (key == 'sleep' or d['main'][p]['ok']) and (key not in ('J', 'E31') or d['main'][p]['cycle'])]
            ax.plot(x + rng.uniform(-0.2, 0.2, len(v)), v, 'o', ms=2.2, mfc=COL[g] if p else 'white', mec=COL[g], mew=0.5, alpha=0.85, zorder=2)
            m, lo, hi = mci(v); ax.plot([x - 0.3, x + 0.3], [m, m], color='k', lw=1.2, zorder=3); ax.plot([x, x], [lo, hi], color='k', lw=0.8, zorder=3)
            if (g, p) in sig: ax.text(x, 1.0, '*', transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=11)
        ticks.append(x0 + 0.9); x0 += 3.3
    ax.set_xticks(ticks); ax.set_xticklabels(['light', 'dark']); ax.set_xlim(-0.6, x0 - 0.9)
    ax.set_ylabel(yl); ax.set_title(title, fontsize=8.5)
    if key == 'J': ax.axhline(0, color='0.6', lw=0.5)
    if key == 'E31':
        ax.set_ylim(-0.6, 3.3); ax.text(0.02, 0.01, 'one control fly at 7.7\n(dark) not shown', transform=ax.transAxes, fontsize=5.5, va='bottom')
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
for g in G: axs[0].plot([], [], 's', ms=5, color=COL[g], label=LAB[g].replace('\n', ' '))
fig.legend(loc='outside upper center', ncol=3, frameon=False, fontsize=7.5)
fig.text(0.5, -0.06, 'Each point is one fly (open: light phase, filled: dark phase); bar: median with 95% CI; * differs from control.',
         ha='center', fontsize=6.5)
for ext in ('png', 'pdf'): fig.savefig(f'figures/figC12_kir_silencing.{ext}', dpi=300, bbox_inches='tight')
print('saved')
