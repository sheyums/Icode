# Figures for PRESPEC 25 (s25_windows.json).
# fig25_raster_{6h,4h}: one row per fly, one column per window in time order; blue = measurably forward, orange =
#   measurably backward, light grey = no measurable cycle. Bar above: dark windows (LD) / constant dark (DD).
# fig25_counts: observed measurable backward and forward windows and the backward fraction against the clock-matched
#   null (grey: 95% of 500 simulated records per fly; tick: null mean), per genotype and window size.
import json, numpy as np, sqw, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})
D = json.load(open(sqw.OUTD + 's25_windows.json'))
G = ['CS', 'CSDD', 'per0LD', 'per0DD', 'sss']
NAME = {'CS': 'CS, light–dark', 'CSDD': 'CS, constant dark', 'per0LD': '$per^0$, light–dark', 'per0DD': '$per^0$, constant dark', 'sss': '$sleepless$, light–dark'}
SH = {'CS': 'CS LD', 'CSDD': 'CS DD', 'per0LD': '$per^0$ LD', 'per0DD': '$per^0$ DD', 'sss': '$sss$ LD'}
FWD, BWD, NONE = '#2a78d6', '#eb6834', '#ececec'; CM = ListedColormap([BWD, NONE, FWD, 'white'])
def raster(size):
    fig, axs = plt.subplots(1, len(G), figsize=(7.2, 3.6), gridspec_kw=dict(width_ratios=[len(D[f'{g}|{size}']['flies'][0]['cls']) for g in G]), constrained_layout=True)
    for ax, g in zip(axs, G):
        v = D[f'{g}|{size}']; fl = v['flies']; M = np.array([f['cls'] for f in fl]); order = np.lexsort((-(M == -1).sum(1), (M == 1).sum(1) - (M == -1).sum(1)))[::-1]
        M = M[order]; M = np.where(M == 9, 2, M); nw = M.shape[1]; per_day = 24 // int(size[:-1])
        ax.imshow(M + 1, cmap=CM, vmin=0, vmax=3, aspect='auto', interpolation='nearest', extent=(0, nw, M.shape[0], 0))
        for d in range(1, nw // per_day): ax.axvline(d * per_day, color='white', lw=1.2)
        dark = [(k % per_day) >= per_day // 2 for k in range(nw)] if g in ('CS', 'per0LD', 'sss') else [True] * nw
        for k, dk in enumerate(dark):
            if dk: ax.add_patch(plt.Rectangle((k, -1.6), 1, 1.0, color='0.25', clip_on=False, lw=0))
            else: ax.add_patch(plt.Rectangle((k, -1.6), 1, 1.0, facecolor='white', edgecolor='0.25', lw=0.4, clip_on=False))
        ax.set_ylim(M.shape[0], -1.8); ax.set_xticks([d * per_day + per_day / 2 for d in range(nw // per_day)]); ax.set_xticklabels([f'd{d + 1}' for d in range(nw // per_day)], fontsize=6.5)
        ax.set_yticks([]); ax.set_title(NAME[g], fontsize=7.5, pad=10)
        r = v['res']; c = r['clock']
        ax.set_xlabel(f"fwd {r['F']} | back {r['B']}\n(constant process: back ≈ {c['B_null_mean']:.0f})", fontsize=6.5)
        for s in ax.spines.values(): s.set_visible(False)
    axs[0].set_ylabel('flies (sorted by net direction)')
    h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (FWD, BWD, NONE)]
    fig.legend(h, ['measurably forward (W→Q→S→W)', 'measurably backward', 'no measurable cycle'], loc='outside upper center', ncol=3, frameon=False, fontsize=7)
    fig.text(0.5, -0.02, f'Each cell is one {size} window on one day. Bar above each panel: dark windows (filled) and light windows (open). Columns are grouped by day (d1, d2, ...).', ha='center', fontsize=6.5)
    for ext in ('png', 'pdf'): fig.savefig(f'figures/fig25_raster_{size}.{ext}', dpi=300, bbox_inches='tight')
def counts():
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9), constrained_layout=True, sharey=True)
    rows = [(g, s) for g in G for s in ('6h', '4h')]; ys = np.arange(len(rows))[::-1]
    for ax, what, xl in ((axs[0], 'B', 'measurably backward windows'), (axs[1], 'F', 'measurably forward windows'), (axs[2], 'frac', 'backward share of measurable windows')):
        for y, (g, s) in zip(ys, rows):
            r = D[f'{g}|{s}']['res']; c = r['clock']
            if what == 'frac':
                lo, hi = c['frac_B_null_95']; mid = c['frac_B_null_median']; obs = r['frac_B']
            else:
                lo, hi = c[f'{what}_null_95']; mid = c[f'{what}_null_mean']; obs = r[what]
            ax.plot([lo, hi], [y, y], color='0.78', lw=6, solid_capstyle='butt'); ax.plot([mid, mid], [y - 0.3, y + 0.3], color='0.35', lw=1)
            ax.plot(obs, y, 'o', ms=5, mfc='k' if s == '6h' else 'white', mec='k', mew=0.8)
        ax.set_xlabel(xl)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    axs[0].set_yticks(ys); axs[0].set_yticklabels([f'{SH[g]}, {s}' for g, s in rows], fontsize=7); axs[0].tick_params(axis='y', length=0)
    axs[2].axvline(0.5, color='0.6', lw=0.5, ls=':')
    for a, l in zip(axs, 'abc'): a.text(-0.08 if l != 'a' else -0.42, 1.02, l, transform=a.transAxes, fontsize=10, fontweight='bold')
    fig.text(0.5, -0.04, 'Dots: observed (filled 6-h windows, open 4-h). Grey bar: 95% range from 500 simulated records per fly of a constant process with each fly\'s own\n'
             'clock-time-matched transition rates; tick: its mean. Dotted line in (c): no preferred direction.', ha='center', fontsize=6.5)
    for ext in ('png', 'pdf'): fig.savefig(f'figures/fig25_counts.{ext}', dpi=300, bbox_inches='tight')
raster('6h'); raster('4h'); counts(); print('saved')
