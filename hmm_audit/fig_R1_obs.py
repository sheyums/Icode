# Main-text version of the sleep-length memory figure (C7), OBSERVED LEVEL ONLY (2026-11-01; the hidden panels of
# fig16_R1_memory_main move to the SI). (a) observed per-fly D with the within-fly null band; (b) CS dark: observed
# P(quiet exit) vs sleep-bout length, with the memoryless reference = the same flies' overall P(quiet exit), which a
# Markov chain predicts at every bout length (definition-consistent; replaces the HMM landscape line).
import json,numpy as np,sys,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker
sys.path.insert(0,'.'); import sqw_ref as R; from route_bins import fly_equal_bins
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
COL={'CS':'#2a78d6','per0LD':'#1baf7a','per0DD':'#eda100','fumin':'#eb6834','sss':'#e87ba4'}
NAME={'CS':'CS','fumin':'$fmn$','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','sss':'$sss$'}
res=json.load(open('sqw_out/k_route_duration.json')); S=json.load(open('sqw_out/k_summary.json'))
BL=[('CS','light'),('CS','dark'),('per0LD','light'),('per0LD','dark'),('per0DD','whole'),('fumin','light'),('fumin','dark'),('sss','light'),('sss','dark')]
rel='R1_exit_vs_sleep'; rng=np.random.default_rng(5)
fig=plt.figure(figsize=(6.2,2.9)); gs=fig.add_gridspec(1,2,width_ratios=[1.25,1])
ax=fig.add_subplot(gs[0,0]); ax.axhline(0,color='#888',lw=0.5)
for k,(g,bn) in enumerate(BL):
    p=1 if bn=='dark' else 0
    D=np.array([R.route_duration_D([x for x in r['obs'] if x['phase_on']==p],rel,p) for r in res if r['g']==g],float)
    D=np.clip(D[np.isfinite(D)],-1.5,1.5); o=S[f'{rel}|observed|{g}|{bn}']; sig=o['p_holm']<0.05
    ax.fill_between([k-0.33,k+0.33],o['null_lo'],o['null_hi'],color='#dcdcdc',lw=0,zorder=0)
    ax.scatter(k+rng.uniform(-0.18,0.18,len(D)),D,s=7,facecolor=COL[g] if bn!='light' else 'white',edgecolor=COL[g],lw=0.6,zorder=2,alpha=1 if sig else 0.55)
    ax.plot([k-0.28,k+0.28],[D.mean()]*2,color='#1a1a1a',lw=1.6,zorder=3)
    if sig: ax.text(k,1.38,'*',ha='center',fontsize=11)
ax.set_xticks(range(len(BL))); ax.set_xticklabels([NAME[g]+('' if bn=='whole' else (' L' if bn=='light' else ' D')) for g,bn in BL],rotation=60,ha='right',fontsize=6.5)
ax.set_ylim(-1.55,1.55); ax.set_title('(a) does the waking route depend on time asleep?',fontsize=8,loc='left')
ax.set_ylabel('$D$ = $\\langle\\ln t\\rangle_{\\rm quiet\\ exit}-\\langle\\ln t\\rangle_{\\rm active\\ exit}$')
ax=fig.add_subplot(gs[0,1])
xs,ys,lo,hi,_=fly_equal_bins(res,'CS',1,'q_out','sleep_min',[5,10,20,40,80,1e9],'observed')
ax.errorbar(xs,ys,yerr=[ys-lo,hi-ys],fmt='o',ls='-',color=COL['CS'],ms=4,lw=1.0,elinewidth=0.7,capsize=1.5,label='observed')
allp=[np.mean([bool(x['q_out']) for x in r['obs'] if x['phase_on']==1 and x['q_out'] is not None]) for r in res if r['g']=='CS']
flat=np.nanmean([v for v in allp if np.isfinite(v)])
ax.axhline(flat,color='#1a1a1a',lw=0.9,ls=(0,(4,2)),label='no memory (Markov): same at every length')
ax.set_xscale('log'); ax.set_xticks([5,10,20,40,80]); ax.set_xticklabels(['5','10','20','40','80']); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_ylim(0,0.5); ax.set_xlabel('sleep-bout length (min)'); ax.set_ylabel('P(wake through quiet wake)'); ax.set_title('(b) CS, night',fontsize=8,loc='left')
ax.legend(frameon=False,fontsize=6,loc='upper right')
fig.text(0.01,0.005,"(a) Points: flies (open light, filled dark/DD; faded where not significant); bar: mean; grey: 95% range of the mean if the route ignored time asleep\n(within-fly shuffle); *: Holm p < 0.05; |D| clipped at 1.5. (b) Mean over flies of each fly's proportion per bin (>= 10 flies), bootstrap 95% intervals.",fontsize=5.6,color='#333')
fig.tight_layout(rect=(0,0.08,1,1),pad=0.4,w_pad=1.2); fig.savefig('figures/figC7_sleep_memory_observed.pdf'); fig.savefig('figures/figC7_sleep_memory_observed.png',dpi=300)
print('saved; flat',round(flat,3))
