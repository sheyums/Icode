# PRESPEC 20m report and figure: 30-s vs 1-min.
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, spearmanr
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
O=json.load(open('sqw_out/m_resolution.json')); L=[]; COL={'CS':'#2a78d6','fumin':'#eb6834'}
keys=[('h1_dwell_min','H1 mean dwell (min, log)'),('h1_one_step','H1 runs of one step (fraction)'),('J','count $J_{cycle}$ (cycles/h)'),('EPR','EPR bound (nats/h)'),('c_ffbs','routing $c$ (posterior)'),('R1_D','R1 hidden $D$')]
deg={b:sum(1 for r in O if r[b]['deg']) for b in ('60','30')}; L.append(f"degenerate H1: 1 min {deg['60']} flies, 30 s {deg['30']} flies (of {len(O)})")
fig,axs=plt.subplots(2,3,figsize=(7.0,4.6)); axs=axs.ravel()
for k,(key,lab) in enumerate(keys):
    ax=axs[k]
    for g in ('CS','fumin'):
        for p,bn in ((0,'light'),(1,'dark')):
            a=[];b=[]
            for r in O:
                if r['g']!=g or r['60']['m'] is None or r['30']['m'] is None: continue
                m1,m3=r['60']['m'][p],r['30']['m'][p]
                if m1['noncycle'] or m3['noncycle'] or r['60']['refs'].get(bn) is None: continue
                if m1[key] is None or m3[key] is None: continue
                a.append(m1[key]); b.append(m3[key])
            a=np.array(a); b=np.array(b)
            ax.scatter(a,b,s=9,facecolor=COL[g] if p==1 else 'white',edgecolor=COL[g],lw=0.7)
            w=wilcoxon(b,a).pvalue if len(a)>5 and np.any(a!=b) else np.nan
            extra=f" | same sign {np.mean(np.sign(a)==np.sign(b)):.2f}" if key in ('c_ffbs','R1_D') else ''
            L.append(f"{key:13s} {g:5s} {bn:5s} n {len(a):2d} | 1 min median {np.median(a):+.3f} | 30 s median {np.median(b):+.3f} | ratio of medians {np.median(b)/np.median(a) if np.median(a)!=0 else np.nan:.2f} | paired Wilcoxon p {w:.2g} | Spearman {spearmanr(a,b).statistic:+.2f}{extra}")
    if key=='h1_dwell_min':
        ax.set_xscale('log'); ax.set_yscale('log'); lo,hi=0.8,200
    else:
        lo=min(ax.get_xlim()[0],ax.get_ylim()[0]); hi=max(ax.get_xlim()[1],ax.get_ylim()[1])
    ax.plot([lo,hi],[lo,hi],color='#888',lw=0.6,ls=(0,(3,2))); ax.set_xlim(lo,hi); ax.set_ylim(lo,hi)
    ax.set_xlabel('1-min steps'); ax.set_ylabel('30-s steps'); ax.set_title(lab,fontsize=8,loc='left')
    if key in ('c_ffbs','R1_D'): ax.axhline(0,color='#ccc',lw=0.5); ax.axvline(0,color='#ccc',lw=0.5)
axs[0].scatter([],[],s=9,color=COL['CS'],label='CS'); axs[0].scatter([],[],s=9,color=COL['fumin'],label='$fmn$'); axs[0].legend(frameon=False,fontsize=6.5)
fig.text(0.01,0.005,'Per fly and phase (open light, filled dark); dashed: identity. Same flies, same theta; 30-s model refitted (K = 3).',fontsize=6,color='#333')
fig.tight_layout(rect=(0,0.03,1,1),pad=0.4); fig.savefig('figures/fig17_resolution_30s_vs_1min.pdf'); fig.savefig('figures/fig17_resolution_30s_vs_1min.png',dpi=300)
# group-level CS dark R1 D
for b in ('60','30'):
    D=[r[b]['m'][1]['R1_D'] for r in O if r['g']=='CS' and r[b]['m'] is not None and r[b]['m'][1]['R1_D'] is not None]
    L.append(f"CS dark R1 hidden D at {b} s: mean {np.mean(D):+.3f}, flies < 0 {np.mean(np.array(D)<0):.2f} (n {len(D)}); Wilcoxon vs 0 p {wilcoxon(D).pvalue:.2g}")
open('sqw12.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
