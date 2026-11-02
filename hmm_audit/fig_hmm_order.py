# SI figure (2026-11-02): HMM order selection. Per fly, K = 2, 3, 4 fitted (20b, theta 0.5, 20 starts each); BIC with
# k = (K-1) + nP K(K-1) + 2K and n = minutes (primary) or number of symbol changes (sensitivity). Bars: fraction of
# flies whose BIC is lowest at each K. Data: sqw20b_fit2_theta0.5.jsonl.
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6})
L = [json.loads(l) for l in open('sqw20b_fit2_theta0.5.jsonl')]
G = ['CS','per0LD','per0DD','fumin','sss']; NAME = {'CS':'CS','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','fumin':'$fmn$','sss':'$sss$'}
KC = {2:'#c9d6e8', 3:'#1f4e79', 4:'#8fb6e0'}
fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.4), sharey=True, constrained_layout=True); rows = []
for ax, (nm, lab) in zip(axs, (('T', 'n = minutes (primary)'), ('nchg', 'n = symbol changes (sensitivity)'))):
    for i, g in enumerate(G):
        D = [d for d in L if d['g'] == g]; nP = 1 if g == 'per0DD' else 2; cnt = {2:0, 3:0, 4:0}
        for d in D:
            n = d['T'] if nm == 'T' else d['nchg']
            bic = {K: -2*d['F'][str(K)]['ll'] + ((K-1) + nP*K*(K-1) + 2*K)*np.log(n) for K in (2, 3, 4)}
            cnt[min(bic, key=bic.get)] += 1
        left = 0
        for K in (2, 3, 4):
            f = cnt[K]/len(D); ax.barh(i, f, left=left, color=KC[K], edgecolor='white', lw=1.5)
            if f >= 0.08: ax.text(left + f/2, i, f'{f:.0%}', ha='center', va='center', fontsize=6.5, color='white' if K == 3 else 'black')
            left += f
        rows.append((nm, g, len(D), cnt))
    ax.set_yticks(range(len(G))); ax.set_yticklabels([f'{NAME[g]} ({sum(1 for d in L if d["g"]==g)})' for g in G]); ax.invert_yaxis()
    ax.set_xlim(0, 1); ax.set_xlabel('fraction of flies, lowest BIC'); ax.set_title(lab, fontsize=8, loc='left')
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
import matplotlib.patches as mp
fig.legend(handles=[mp.Patch(color=KC[K], label=f'K = {K} hidden states') for K in (2, 3, 4)], frameon=False, fontsize=7, loc='outside lower center', ncol=3)
for ext in ('png', 'pdf'): fig.savefig(f'figures/figS_hmm_order_bic.{ext}', dpi=300, bbox_inches='tight')
for r in rows: print(r)
