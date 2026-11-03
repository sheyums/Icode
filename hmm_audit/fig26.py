# Figures for PRESPEC 26 (s26_loopseq.json).
# (a) loop barcodes: the sequence of completed loops (blue forward, orange backward) for example flies, at their
#     completion times; (b) per-fly excess persistence E (same-direction fraction minus the clock-matched chain's);
# (c) observed vs simulated same-direction fraction per fly; (d) run lengths of same-direction loops, observed vs chain.
import json, numpy as np, sqw, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sqw15 import seqs
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
D = json.load(open(sqw.OUTD + 's26_loopseq.json')); FL = D['flies']; GR = D['genotypes']
G = ['CS', 'CSDD', 'per0LD', 'per0DD', 'fumin', 'sss']
SH = {'CS': 'CS LD', 'CSDD': 'CS DD', 'per0LD': '$per^0$ LD', 'per0DD': '$per^0$ DD', 'fumin': '$fumin$', 'sss': '$sss$'}
COL = {'CS': '#2a78d6', 'CSDD': '#6250d6', 'per0LD': '#1baf7a', 'per0DD': '#eda100', 'fumin': '#d43d51', 'sss': '#6a3d9a'}
FWD_C, BWD_C = '#2a78d6', '#eb6834'; rng = np.random.default_rng(26)
FWD = {(1, 0, 2), (0, 2, 1), (2, 1, 0)}
def loops_timed(y):
    out = []; stack = []; prev = None
    for t, x in enumerate(y.tolist()):
        if x < 0: stack = []; prev = None; continue
        if x == prev: continue
        prev = x
        if x in stack:
            k = stack.index(x); loop = stack[k:]
            if len(loop) == 3: out.append((t, 1 if tuple(loop) in FWD else -1))
            del stack[k + 1:]
        else: stack.append(x)
    return out
fig = plt.figure(figsize=(7.2, 6.4)); gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.25, 1.25], hspace=0.6, wspace=0.42)
# (a) barcodes: the fly closest to its genotype's median loop count, CS LD and per0 DD, plus CS DD
axa = fig.add_subplot(gs[0, :]); ex = []
for g in ('CS', 'CSDD', 'per0DD'):
    R = [r for r in FL if r['g'] == g and r['n_pairs'] >= 10]; m = np.median([r['n_loops'] for r in R]); ex.append(min(R, key=lambda r: abs(r['n_loops'] - m)))
for row, r in enumerate(ex):
    y, _, _ = seqs(r['g'], r['fly']); L = loops_timed(np.asarray(y)); yb = len(ex) - 1 - row
    if r['g'] == 'CS':
        for d in range(4): axa.axvspan(24 * d + 12, 24 * d + 24, ymin=(yb + 0.1) / len(ex), ymax=(yb + 0.9) / len(ex), color='0.9', lw=0)
    for t, d in L: axa.plot([t / 60, t / 60], [yb - 0.35, yb + 0.35], color=FWD_C if d == 1 else BWD_C, lw=1.1)
    axa.text(-1, yb, f"{SH[r['g']]} (fly {r['fly'] + 1})", ha='right', va='center', fontsize=7)
axa.set_yticks([]); axa.set_xlim(0, 96); axa.set_ylim(-0.6, len(ex) - 0.4); axa.set_xticks(range(0, 97, 24)); axa.set_xlabel('time (h)')
for s in ('top', 'right', 'left'): axa.spines[s].set_visible(False)
axa.set_title('each tick: one completed loop (blue forward, orange backward); grey: dark phase', fontsize=7, loc='left')
# (b) per-fly excess
axb = fig.add_subplot(gs[1, :2])
for k, g in enumerate(G):
    e = np.array([r['E_clock'] for r in FL if r['g'] == g and r['n_pairs'] >= 10])
    axb.plot(k + rng.uniform(-0.2, 0.2, len(e)), e, 'o', ms=2.6, mfc=COL[g], mec='white', mew=0.3, alpha=0.85)
    bs = [np.median(rng.choice(e, len(e))) for _ in range(2000)]; lo, hi = np.percentile(bs, [2.5, 97.5]); m = np.median(e)
    axb.plot([k - 0.3, k + 0.3], [m, m], color='k', lw=1.3); axb.plot([k, k], [lo, hi], color='k', lw=0.8)
    c = GR[g]['clock']; lab = f"p {c['p']:.2f}" + (f"\nHolm {c['holm']:.2f}" if 'holm' in c else '\n(descr.)')
    axb.text(k, 1.0, lab, transform=axb.get_xaxis_transform(), ha='center', va='bottom', fontsize=6)
axb.axhline(0, color='0.5', lw=0.6); axb.set_xticks(range(len(G))); axb.set_xticklabels([SH[g] for g in G], fontsize=7)
axb.set_ylabel('excess persistence $E$\n(same-direction fraction, data − chain)')
for s in ('top', 'right'): axb.spines[s].set_visible(False)
axb.set_title('above 0: runs of same-direction loops longer than the chain predicts (episodes); below 0: more alternation', fontsize=6.5, loc='left', pad=22)
# (c) observed vs null S
axc = fig.add_subplot(gs[1, 2])
for g in G:
    R = [r for r in FL if r['g'] == g and r['n_pairs'] >= 10]
    axc.plot([r['S_null_clock'] for r in R], [r['S_obs'] for r in R], 'o', ms=2.4, mfc=COL[g], mec='white', mew=0.3, alpha=0.85, label=SH[g])
axc.plot([0.3, 1], [0.3, 1], color='k', lw=0.6); axc.set_xlim(0.35, 1); axc.set_ylim(0.2, 1); axc.set_xlabel('same-direction fraction,\nclock-matched chain'); axc.set_ylabel('same-direction fraction, data')
axc.legend(fontsize=5.5, frameon=False, loc='lower right', handletextpad=0.1)
for s in ('top', 'right'): axc.spines[s].set_visible(False)
# (d) run-length survival, observed vs null, pooled per genotype
for j, g in enumerate(('CS', 'CSDD', 'per0DD')):
    ax = fig.add_subplot(gs[2, j]); R = [r for r in FL if r['g'] == g and r['n_pairs'] >= 10]
    for key, ls, lab in (('runs_obs', '-', 'data'), ('runs_null_clock', '--', 'chain')):
        v = np.concatenate([r[key] for r in R]); k = np.arange(1, 11); sv = [np.mean(v >= x) for x in k]
        ax.semilogy(k, sv, ls, color=COL[g] if key == 'runs_obs' else '0.3', lw=1.3, marker='o' if key == 'runs_obs' else None, ms=2.5, label=lab)
    ax.set_xlabel('run of same-direction loops (length ≥ k)'); ax.set_title(SH[g], fontsize=8, loc='left')
    if j == 0: ax.set_ylabel('fraction of runs')
    ax.legend(fontsize=6, frameon=False)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
for l, (x, y_) in zip('abcd', ((0.02, 0.965), (0.02, 0.62), (0.68, 0.62), (0.02, 0.31))): fig.text(x, y_, l, fontsize=10, fontweight='bold')
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig26_loop_sequence.{ext}', dpi=300, bbox_inches='tight')
print('saved')
