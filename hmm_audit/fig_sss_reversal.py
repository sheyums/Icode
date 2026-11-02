# Main-figure candidate (2026-11-02): sleepless runs the cycle BACKWARDS in the light phase. Per fly, cycle flux of the
# observed S/Q/W chain (exact chain J, defined for every fly) beside the HMM's flux (the check), CS for reference.
# Data: s_observable.json (20s-A). Descriptive figure; the pre-registered sign test is 20d (model-free current).
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
O = json.load(open('sqw_out/s_observable.json')); COL = {'CS':'#2a78d6','sss':'#e87ba4'}
BL = [('CS',0,'CS\nlight'),('CS',1,'CS\ndark'),('sss',0,'$sss$\nlight'),('sss',1,'$sss$\ndark')]
fig, ax = plt.subplots(figsize=(5.2, 3.0), constrained_layout=True); rng = np.random.default_rng(2); out = {}
for k, (g, p, lab) in enumerate(BL):
    R = [d['main'][p] for d in O if d['g'] == g and d['main'][p]['ok'] and not d['main'][p].get('hmm_deg', True) and d['main'][p].get('J_hmm') is not None]
    jo = np.array([r['J_chain'] for r in R]); jh = np.array([r['J_hmm'] for r in R]); x0 = 3*k
    jo = np.where(np.abs(jo) < 1e-12, 0.0, jo)   # no cycle in the observed chain (20s-A rule): exactly 0, not round-off
    cyc = jo != 0
    for a, b in zip(jo, jh): ax.plot([x0, x0+1], [a, b], color='0.82', lw=0.5, zorder=1)
    ax.scatter(x0 + rng.uniform(-0.12, 0.12, len(jo)), jo, s=9, facecolor=COL[g] if p else 'white', edgecolor=COL[g], lw=0.7, zorder=2)
    ax.scatter(x0 + 1 + rng.uniform(-0.12, 0.12, len(jh)), jh, s=9, marker='^', facecolor=COL[g] if p else 'white', edgecolor=COL[g], lw=0.7, zorder=2)
    for xx, v in ((x0, jo), (x0+1, jh)): ax.plot([xx-0.3, xx+0.3], [np.mean(v)]*2, color='k', lw=1.4, zorder=3)
    po = wilcoxon(jo[cyc]).pvalue; ph = wilcoxon(jh).pvalue
    ax.text(x0 + 0.5, 0.60, f'J < 0 in {np.mean(jo[cyc] < 0):.0%} | {np.mean(jh < 0):.0%}\np {po:.1g} | {ph:.1g}' + (f'\n({(~cyc).sum()} no cycle)' if (~cyc).any() else ''), ha='center', fontsize=6.0)
    out[f'{g}{p}'] = dict(n=len(jo), n_noCycle=int((~cyc).sum()), mean_obs=float(jo.mean()), mean_hmm=float(jh.mean()), neg_obs_cycleOnly=float(np.mean(jo[cyc]<0)), neg_hmm=float(np.mean(jh<0)), p_obs_cycleOnly=float(po), p_hmm=float(ph))
ax.axhline(0, color='0.5', lw=0.6); ax.set_ylim(-0.5, 0.75)
ax.set_xticks([3*k + 0.5 for k in range(4)]); ax.set_xticklabels([b[2] for b in BL])
ax.set_ylabel('cycle flux $J$ (cycles h$^{-1}$)\n> 0: W→Q→S→W,  < 0: reversed')
ax.scatter([], [], s=9, facecolor='0.4', edgecolor='0.4', label='observed S/Q/W chain'); ax.scatter([], [], s=9, marker='^', facecolor='0.4', edgecolor='0.4', label='hidden Markov model (check)')
ax.legend(frameon=False, fontsize=6.5, loc='lower left')
fig.text(0.0, -0.04, 'Per fly: circle = observed chain, triangle = hidden Markov model, joined by a line; bar = mean. Percentages: flies with J < 0 among flies whose observed chain has a cycle (observed | HMM); flies without a cycle (|J| < 1e-12) are plotted at 0 and left out of the % and the test.', fontsize=5.8)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.set_title('Sleepless flies run the cycle backwards on average by day', fontsize=9, loc='left')
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig_sss_reversal.{ext}', dpi=300, bbox_inches='tight')
json.dump(out, open('sqw_out/sss_reversal_fig.json', 'w'), indent=1); print(json.dumps(out, indent=1))
