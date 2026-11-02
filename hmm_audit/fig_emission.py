# Fig 1 panel (2026-11-02): the HMM's hidden states coincide with the observed symbols. Per genotype, the median over
# non-degenerate flies of the labelled 3x3 emission matrix B3 (rows: hidden H1/H2/H3 = quiet/active/sleep-like;
# columns: observed Q, W, S), from the stored 20b fits (sqw20b_fit2_theta0.5.jsonl, K = 3, theta 0.5).
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8})
L = [json.loads(l) for l in open('sqw20b_fit2_theta0.5.jsonl')]
G = ['CS', 'per0LD', 'per0DD', 'fumin', 'sss']; NAME = {'CS':'CS','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','fumin':'$fmn$','sss':'$sss$'}
fig, axs = plt.subplots(1, len(G), figsize=(7.0, 1.9), constrained_layout=True); out = {}
for ax, g in zip(axs, G):
    Bs = np.array([d['B3'] for d in L if d['g'] == g and not d.get('deg', False) and d.get('B3') is not None], float)
    M = np.median(Bs, 0); out[g] = dict(n=len(Bs), diag=[float(M[i, i]) for i in range(3)])
    ax.imshow(M, cmap='Greys', vmin=0, vmax=1)
    for i in range(3):
        for j in range(3): ax.text(j, i, f'{M[i,j]:.2f}', ha='center', va='center', fontsize=6.5, color='white' if M[i, j] > 0.55 else 'black')
    ax.set_xticks(range(3)); ax.set_xticklabels(['Q', 'W', 'S']); ax.set_yticks(range(3)); ax.set_yticklabels(['$H_1$', '$H_2$', '$H_3$'] if g == 'CS' else [])
    ax.set_title(f'{NAME[g]} (n {len(Bs)})', fontsize=7.5)
    if g == 'CS': ax.set_ylabel('hidden state')
fig.supxlabel('observed symbol (median emission probability over flies)', fontsize=7)
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig1_hmm_emissions.{ext}', dpi=300, bbox_inches='tight')
print(json.dumps(out))
