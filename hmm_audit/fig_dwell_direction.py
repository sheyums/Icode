# DESCRIPTIVE figure (no new test; 2026-11-01): direction-dependent dwell times in CS at night (light-dark, ZT 12-24),
# the block with the largest waiting-time excess in PRESPEC 20v (descriptive there). Same definitions as 20v:
# 1-min symbols (sqw_ref.sqw_symbols, theta 0.5, light/dark references); runs of Q (or W) whose entry and exit
# neighbours are classified and lie, with the run, in the dark phase. Q runs: W->Q->S vs S->Q->W; W runs: Q->W->S vs
# S->W->Q. Shown as survival curves P(dwell >= t), flies pooled with fly-bootstrap 95% bands, and per-fly mean dwell.
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sqw, sqw_ref as R
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
Q, W, S = 0, 1, 2
def runs(y, ph):
    out = []; ch = np.r_[0, np.flatnonzero(np.diff(y)) + 1, len(y)]
    for k in range(1, len(ch) - 2):
        a, b = ch[k], ch[k+1]; s = y[a]
        if s not in (Q, W): continue
        e, x = y[a-1], y[b]
        if e < 0 or x < 0 or not (ph[a-1:b+1] == 1).all(): continue
        out.append((s, int(e), int(x), b - a))
    return out
F = []
for f in sqw.files('CS'):
    x = np.loadtxt(f).astype(int); y, _ = R.sqw_symbols(x, 0.5, light_dark=True); F.append(runs(np.asarray(y), R.phases(len(y), '12h')))
PAIRS = {Q: [((W, S), 'into sleep  W→Q→S'), ((S, W), 'out of sleep  S→Q→W')], W: [((Q, S), 'into sleep  Q→W→S'), ((S, Q), 'out of sleep  S→W→Q')]}
COLS = ['#2a78d6', '#eb6834']; T = np.arange(1, 31); rng = np.random.default_rng(3)
def surv(d): d = np.asarray(d); return np.array([(d >= t).mean() for t in T]) if len(d) else np.full(len(T), np.nan)
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6), gridspec_kw=dict(width_ratios=[1, 1, 0.8]), constrained_layout=True)
summ = {}
for ax, s, nm in ((axs[0], Q, 'quiet wake (Q)'), (axs[1], W, 'active wake (W)')):
    for (pair, lab), c in zip(PAIRS[s], COLS):
        per = [[r[3] for r in fl if r[0] == s and (r[1], r[2]) == pair] for fl in F]
        pool = [v for p in per for v in p]; sv = surv(pool)
        bs = np.array([surv([v for i in rng.integers(0, len(per), len(per)) for v in per[i]]) for _ in range(1000)])
        lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
        ax.fill_between(T, lo, hi, color=c, alpha=0.18, lw=0); ax.plot(T, sv, color=c, lw=2, label=f'{lab} (n {len(pool)})')
        summ[(s, pair)] = [np.mean(p) if len(p) >= 2 else np.nan for p in per]
    ax.set_yscale('log'); ax.set_xlim(1, 30); ax.set_ylim(5e-3, 1.05); ax.set_xlabel('time in state (min)')
    ax.set_ylabel('fraction of stays lasting ≥ t'); ax.set_title(nm + ', CS at night', fontsize=8.5, loc='left')
    ax.legend(frameon=False, fontsize=6.3, loc='upper right')
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax = axs[2]
for k, s in enumerate((Q, W)):
    a = np.array(summ[(s, PAIRS[s][0][0])]); b = np.array(summ[(s, PAIRS[s][1][0])]); m = np.isfinite(a) & np.isfinite(b)
    for ai, bi in zip(a[m], b[m]): ax.plot([2*k, 2*k + 0.8], [ai, bi], color='0.75', lw=0.6)
    ax.plot(np.full(m.sum(), 2*k), a[m], 'o', ms=3, color=COLS[0]); ax.plot(np.full(m.sum(), 2*k + 0.8), b[m], 'o', ms=3, color=COLS[1])
    print(f"{'QW'[s]}: flies {m.sum()} | mean dwell into sleep {np.median(a[m]):.2f} min vs out of sleep {np.median(b[m]):.2f} (medians over flies); "
          f"longer into sleep in {np.mean(a[m] > b[m]):.0%} of flies")
ax.set_xticks([0.4, 2.4]); ax.set_xticklabels(['Q', 'W']); ax.set_ylabel('mean stay per fly (min)')
ax.set_title('per fly: into vs out of sleep', fontsize=8.5, loc='left')
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.text(0.0, -0.06, 'DESCRIPTIVE. Blue: stays on the way into sleep; orange: stays on the way out of sleep. Bands: 95% bootstrap over flies. '
         'Grey lines join one fly\'s two means.', fontsize=6)
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig_dwell_direction_CS_night.{ext}', dpi=300, bbox_inches='tight')
