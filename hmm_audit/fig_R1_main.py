# Main-figure candidate: R1 (exit route vs sleep-bout length). (a) observed per-fly D, (b) hidden per-fly D,
# (c) CS dark: hidden P(quiet exit) vs sleep-bout length (flies weighted equally) with the flat landscape prediction.
import json,numpy as np,sys,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker
sys.path.insert(0,'.'); import sqw_ref as R; from route_bins import fly_equal_bins
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
COL={'CS':'#2a78d6','per0LD':'#1baf7a','per0DD':'#eda100','fumin':'#eb6834','sss':'#e87ba4'}
NAME={'CS':'CS','fumin':'$fmn$','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','sss':'$sss$'}
res=json.load(open('sqw_out/k_route_duration.json')); S=json.load(open('sqw_out/k_summary.json')); LS=json.load(open('sqw_out/l_landscape.json'))
BL=[('CS','light'),('CS','dark'),('per0LD','light'),('per0LD','dark'),('per0DD','whole'),('fumin','light'),('fumin','dark'),('sss','light'),('sss','dark')]
rel='R1_exit_vs_sleep'; rng=np.random.default_rng(5)
fig=plt.figure(figsize=(7.0,2.9)); gs=fig.add_gridspec(1,3,width_ratios=[1,1,0.95])
for c,level in enumerate(('observed','hidden')):
    ax=fig.add_subplot(gs[0,c]); ax.axhline(0,color='#888',lw=0.5)
    for k,(g,bn) in enumerate(BL):
        p=1 if bn=='dark' else 0
        D=[R.route_duration_D([x for x in r['obs'] if x['phase_on']==p],rel,p) for r in res if r['g']==g] if level=='observed' else [r['hid'].get(f'{rel}|{p}') for r in res if r['g']==g]
        D=np.array([d for d in D if d is not None and np.isfinite(d)]); D=np.clip(D,-1.5,1.5)
        o=S[f'{rel}|{level}|{g}|{bn}']; sig=o['p_holm']<0.05
        ax.fill_between([k-0.33,k+0.33],o['null_lo'],o['null_hi'],color='#dcdcdc',lw=0,zorder=0)
        ax.scatter(k+rng.uniform(-0.18,0.18,len(D)),D,s=7,facecolor=COL[g] if bn!='light' else 'white',edgecolor=COL[g],lw=0.6,zorder=2,alpha=1 if sig else 0.55)
        ax.plot([k-0.28,k+0.28],[D.mean()]*2,color='#1a1a1a',lw=1.6,zorder=3)
        if sig: ax.text(k,1.38,'*',ha='center',fontsize=11)
    ax.set_xticks(range(len(BL))); ax.set_xticklabels([NAME[g]+('' if bn=='whole' else (' L' if bn=='light' else ' D')) for g,bn in BL],rotation=60,ha='right',fontsize=6.5)
    ax.set_ylim(-1.55,1.55); ax.set_title(f"({'ab'[c]}) {'observed behaviour' if level=='observed' else 'hidden states'}",fontsize=8,loc='left')
    if c==0: ax.set_ylabel('$D$ = $\\langle\\ln t\\rangle_{\\rm quiet\\ exit}-\\langle\\ln t\\rangle_{\\rm active\\ exit}$')
ax=fig.add_subplot(gs[0,2])
for lev,mk,ls in (('hidden','o','-'),('observed','s',(0,(2,1.5)))):
    xs,ys,lo,hi,_=fly_equal_bins(res,'CS',1,'q_out','sleep_min',[5,10,20,40,80,1e9],lev)
    ax.errorbar(xs,ys,yerr=[ys-lo,hi-ys],fmt=mk,ls=ls,color=COL['CS'],ms=4 if lev=='hidden' else 3.5,lw=1.0,elinewidth=0.7,capsize=1.5,
                mfc=COL['CS'] if lev=='hidden' else 'white',label='hidden (H1 exit)' if lev=='hidden' else 'observed (Q exit)')
pred=np.mean([(r['F_SW']+1)/2 for r in LS if r['g']=='CS' and r['phase']==1])
ax.axhline(pred,color='#1a1a1a',lw=0.9,ls=(0,(4,2)),label='landscape (Markov) prediction')
ax.set_xscale('log'); ax.set_xticks([5,10,20,40,80]); ax.set_xticklabels(['5','10','20','40','80']); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_ylim(0,0.8); ax.set_xlabel('sleep-bout length (min)'); ax.set_ylabel('P(quiet exit)'); ax.set_title('(c) CS, dark',fontsize=8,loc='left')
ax.legend(frameon=False,fontsize=6,loc='upper right')
fig.text(0.01,0.005,'(a, b) Points: flies (open light, filled dark/DD; faded where not significant); bar: mean; grey: 95% range of the mean without dependence; *: Holm p < 0.05; |D| clipped at 1.5.\n(c) Mean over flies of each fly\'s proportion per bin (>= 10 flies), bootstrap 95% intervals; dashed: prediction from the fitted transition matrix, constant in bout length.',fontsize=5.6,color='#333')
fig.tight_layout(rect=(0,0.07,1,1),pad=0.4,w_pad=1.0); fig.savefig('figures/fig16_R1_memory_main.pdf'); fig.savefig('figures/fig16_R1_memory_main.png',dpi=300)
