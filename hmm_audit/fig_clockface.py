# Fig C10 "The clock drives the cycle" (PRESPEC 20w-B; data planB_20wB.json from planB-feasibility 803a3c4).
# (a) pooled cycle flux J per 2-h clock bin, fitted 24-h cosine, and the band a random clock phase would give
#     (mean +/- the null's 95th-percentile amplitude); (b) observed amplitude vs that null amplitude per genotype.
import json, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})
src = sys.argv[1] if len(sys.argv) > 1 else 'planB_20wB.json'
R = json.load(open(src))['real']
COL = {'CS': '#2a78d6', 'CSDD': '#6250d6', 'per0LD': '#1baf7a', 'per0DD': '#eda100'}
NAME = {'CS': 'CS, light–dark', 'CSDD': 'CS, constant dark', 'per0LD': '$per^0$, light–dark', 'per0DD': '$per^0$, constant dark'}
VERD = {'CS': 'clock-locked rhythm', 'CSDD': 'clock-locked rhythm', 'per0LD': 'not detected', 'per0DD': 'none (control)'}
G = ['CS', 'CSDD', 'per0LD', 'per0DD']; NB = 12; b = np.arange(NB); t = 2 * b + 1
X = np.column_stack([np.ones(NB), np.cos(2 * np.pi * b / NB), np.sin(2 * np.pi * b / NB)])
fig = plt.figure(figsize=(7.0, 4.6)); gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.75], hspace=0.55, wspace=0.55)
for k, g in enumerate(G):
    ax = fig.add_subplot(gs[k // 2, k % 2]); J = np.array(R[g]['J'], float); (m, a, c), *_ = np.linalg.lstsq(X, J, rcond=None)
    tt = np.linspace(0, 24, 200); fit = m + a * np.cos(2 * np.pi * (tt - 1) / 24) + c * np.sin(2 * np.pi * (tt - 1) / 24)
    if g in ('CS', 'per0LD'): ax.axvspan(12, 24, color='0.88', lw=0, zorder=0)
    else: ax.axvspan(0, 24, color='0.94', lw=0, zorder=0)
    A0 = R[g]['nullJ95']; ax.fill_between([0, 24], m - A0, m + A0, color='0.75', alpha=0.45, lw=0, zorder=1)
    ax.plot(tt, fit, color=COL[g], lw=2, zorder=3); ax.plot(t, J, 'o', ms=4, mfc=COL[g], mec='white', mew=0.6, zorder=4)
    ax.set_xlim(0, 24); ax.set_xticks([0, 6, 12, 18, 24]); ax.set_ylim(-0.02, 0.30)
    ax.set_xlabel('clock time (h)' + (' (ZT)' if g in ('CS', 'per0LD') else ' from start of each day'))
    if k % 2 == 0: ax.set_ylabel('cycle flux $J$ (cycles h$^{-1}$)\nW→Q→S→W, all flies pooled')
    ax.set_title(NAME[g], fontsize=8.5, loc='left'); top = g in ('CS', 'per0LD'); ax.text(0.03 if top else 0.98, 0.97 if top else 0.04, VERD[g], transform=ax.transAxes, ha='left' if top else 'right', va='top' if top else 'bottom', fontsize=7, style='italic')
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax = fig.add_subplot(gs[:, 2])
for k, g in enumerate(G):
    y = len(G) - 1 - k; ax.plot([0, R[g]['nullJ95']], [y, y], color='0.75', lw=6, solid_capstyle='butt')
    ax.plot(R[g]['AJ'], y, 'o', ms=7, mfc=COL[g], mec='white', mew=0.8)
    ax.text(max(R[g]['AJ'], R[g]['nullJ95']) + 0.008, y, f"p = {R[g]['pJ']:.3f}", va='center', fontsize=7)
SH = {'CS': 'CS LD', 'CSDD': 'CS DD', 'per0LD': '$per^0$ LD', 'per0DD': '$per^0$ DD'}
ax.set_yticks(range(len(G))); ax.set_yticklabels([SH[g] for g in G][::-1], fontsize=7.5); ax.set_xlim(0, 0.12); ax.set_ylim(-0.6, len(G) - 0.4)
ax.tick_params(axis='y', length=0, pad=6)
ax.set_xlabel('size of the daily rhythm in $J$\n(cosine amplitude, cycles h$^{-1}$)')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.set_title('grey bar: largest rhythm (95%) a random\nclock phase per fly gives from the same data', fontsize=6.5, loc='left')
fig.text(0.01, 0.97, 'a', fontsize=10, fontweight='bold'); fig.text(0.70, 0.97, 'b', fontsize=10, fontweight='bold')
fig.text(0.01, 0.005, 'Shading: dark phase (light–dark) or constant darkness. Grey band in (a): mean ± the null amplitude.',
         fontsize=6.5)
for ext in ('png', 'pdf'): fig.savefig(f'figures/figC10_clock_face.{ext}', dpi=300, bbox_inches='tight')
print('saved', {g: (round(R[g]['AJ'], 4), round(R[g]['phJ'], 1), R[g]['pJ']) for g in G})
