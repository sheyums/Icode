# Act 1 opening figure (descriptive): each CS fly's cumulative net forward cycles, (#sleep->active wake) -
# (#active wake->sleep), equal to completed net cycles +/- the loop in progress (checked: |difference| <= 1 for every
# minute of every fly). Black: mean over flies; band: mean +/- 1 SD across flies. Inset: completed loops per hour from
# cycle_times.json. Light-dark (left) and constant darkness (right).
import json, numpy as np, sqw, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sqw15 import seqs
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
CT = json.load(open(sqw.OUTD + 'cycle_times.json'))
P = [('CS', 'CS, light–dark', '#2a78d6', True), ('CSDD', 'CS, constant darkness', '#6250d6', False)]
fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.3), sharey=True, constrained_layout=True)
for ax, (g, title, col, ld) in zip(axs, P):
    M = []
    for i in range(len(sqw.files(g))):
        y, _, _ = seqs(g, i); y = np.asarray(y); ok = (y[:-1] >= 0) & (y[1:] >= 0)
        st = np.where(ok & (y[:-1] == 2) & (y[1:] == 1), 1, 0) - np.where(ok & (y[:-1] == 1) & (y[1:] == 2), 1, 0)
        M.append(np.r_[0, np.cumsum(st)])
    T = min(len(m) for m in M); M = np.array([m[:T] for m in M]); t = np.arange(T) / 60
    import matplotlib.transforms as mtr; tr = mtr.blended_transform_factory(ax.transData, ax.transAxes)
    for d in range(int(round(t[-1])) // 24 + 1):
        if d * 24 >= t[-1]: break
        if ld:
            ax.add_patch(plt.Rectangle((24 * d, 1.005), 12, 0.03, transform=tr, facecolor='white', edgecolor='k', lw=0.5, clip_on=False))
            ax.add_patch(plt.Rectangle((24 * d + 12, 1.005), 12, 0.03, transform=tr, facecolor='k', lw=0, clip_on=False))
        else:
            ax.add_patch(plt.Rectangle((24 * d, 1.005), 24, 0.03, transform=tr, facecolor='0.45', lw=0, clip_on=False))
    for m in M: ax.plot(t, m, color=col, lw=0.6, alpha=0.35, zorder=1)
    mu, sd = M.mean(0), M.std(0)
    ax.fill_between(t, mu - sd, mu + sd, color=col, alpha=0.18, lw=0, zorder=2)
    ax.plot(t, mu, color='k', lw=1.8, zorder=3); ax.axhline(0, color='0.5', lw=0.5)
    ax.set_xlim(0, t[-1]); ax.set_xticks(np.arange(0, t[-1] + 1, 24)); ax.set_xlabel('time (h)'); ax.set_title(title, fontsize=9, loc='left', pad=12)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    end = M[:, -1]
    if ld:
        d, n = CT['CS|0'], CT['CS|1']
        txt = (f"completed loops per hour\n  day:   {d['fwd_per_h']:.2f} fwd, {d['back_per_h']:.2f} back\n  night: {n['fwd_per_h']:.2f} fwd, {n['back_per_h']:.2f} back\n"
               f"one net cycle per ~{d['h_per_net_cycle']:.0f} h by day, ~{n['h_per_net_cycle']:.0f} h at night\none loop lasts ~{d['dur_fwd_median_min']:.0f}–{n['dur_fwd_median_min']:.0f} min (median)")
    else:
        c = CT['CSDD|0']
        txt = (f"completed loops per hour\n  {c['fwd_per_h']:.2f} fwd, {c['back_per_h']:.2f} back\n"
               f"one net cycle per ~{c['h_per_net_cycle']:.0f} h\none loop lasts ~{c['dur_fwd_median_min']:.0f} min (median)")
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va='top', ha='left', fontsize=6.5, bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='0.7', lw=0.5))
    ax.text(0.98, 0.02, f"{np.mean(end > 0):.0%} of {len(end)} flies end net forward\nmean {end.mean():+.1f} cycles in {t[-1]:.0f} h", transform=ax.transAxes, ha='right', va='bottom', fontsize=6.5)
axs[0].set_ylabel('net forward cycles completed\n(W→Q→S→W minus reverse)')
for a, l in zip(axs, 'ab'): a.text(-0.10, 1.03, l, transform=a.transAxes, fontsize=10, fontweight='bold')
fig.text(0.5, -0.06, 'Thin lines: one fly each. Black: mean; band: ±1 SD across flies. Bar on top: light (open) and dark (filled) phases; grey bar: constant darkness. The slope of the black line is the cycle flux J;\n'
         'single flies wander forward and back, and the forward drift stands clear of their spread only over days.', ha='center', fontsize=6.5)
for ext in ('png', 'pdf'): fig.savefig(f'figures/figAct1_cumulative_walk.{ext}', dpi=300, bbox_inches='tight')
print('saved')
