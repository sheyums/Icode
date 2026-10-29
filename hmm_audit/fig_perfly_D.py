# Fig 13d: per-fly route-duration effect D (PRESPEC 20k statistic), flies as points, group mean, null range of the mean.
import json,numpy as np,sys,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'.'); import sqw_ref as R
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
COL={'CS':'#2a78d6','per0LD':'#1baf7a','per0DD':'#eda100','fumin':'#eb6834','sss':'#e87ba4'}
NAME={'CS':'CS','fumin':'$fmn$','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','sss':'$sss$'}
res=json.load(open('sqw_out/k_route_duration.json')); S=json.load(open('sqw_out/k_summary.json'))
BL=[('CS','light'),('CS','dark'),('per0LD','light'),('per0LD','dark'),('per0DD','whole'),('fumin','light'),('fumin','dark'),('sss','light'),('sss','dark')]
TIT={'R1_exit_vs_sleep':'R1 quiet exit vs sleep-bout length','R2_entry_vs_wake_before':'R2 quiet entry vs prior wake length','R3_entry_vs_sleep':'R3 quiet entry vs next sleep length'}
rng=np.random.default_rng(5); fig,axs=plt.subplots(2,3,figsize=(7.0,5.2),sharex=True,sharey='row'); L=[]
for col,rel in enumerate(R.RELATIONS):
    for row,level in enumerate(('observed','hidden')):
        ax=axs[row,col]; ax.axhline(0,color='#888',lw=0.5)
        for k,(g,bn) in enumerate(BL):
            p=1 if bn=='dark' else 0
            if level=='observed': D=[R.route_duration_D([x for x in r['obs'] if x['phase_on']==p],rel,p) for r in res if r['g']==g]
            else: D=[r['hid'].get(f'{rel}|{p}') for r in res if r['g']==g]
            D=np.array([d for d in D if d is not None and np.isfinite(d)])
            o=S[f'{rel}|{level}|{g}|{bn}']; filled=bn!='light'
            ax.fill_between([k-0.32,k+0.32],o['null_lo'],o['null_hi'],color='#d9d9d9',lw=0,zorder=0)
            ax.scatter(k+rng.uniform(-0.18,0.18,len(D)),D,s=7,facecolor=COL[g] if filled else 'white',edgecolor=COL[g],lw=0.6,zorder=2)
            ax.plot([k-0.28,k+0.28],[D.mean()]*2,color='#1a1a1a',lw=1.5,zorder=3)
            if o['p_holm']<0.05: ax.text(k,ax.get_ylim()[1] if False else 1.55,'*',ha='center',fontsize=10)
            L.append(f"{rel} {level} {g} {bn}: n {len(D)} mean {D.mean():+.3f} median {np.median(D):+.3f} flies D<0 {np.mean(D<0):.2f}")
        ax.set_xticks(range(len(BL))); ax.set_xticklabels([NAME[g]+('' if bn=='whole' else (' L' if bn=='light' else ' D')) for g,bn in BL],rotation=60,fontsize=6.5,ha='right')
        ax.set_ylim(-1.8,1.8)
        if row==0: ax.set_title(TIT[rel],fontsize=7.5)
        if col==0: ax.set_ylabel(f'{level}\n$D$ (per fly)')
fig.text(0.01,0.005,'Points: flies (open light, filled dark/DD). Bar: mean over flies. Grey: 95% range of the mean without dependence (observed: within-fly label permutation; hidden: Markov model simulated and decoded).\n*: Holm p < 0.05. D = mean ln dur(quiet route) - mean ln dur(active route); D < 0: quiet routes go with shorter bouts. Axis clipped at |D| = 1.8.',fontsize=5.6,color='#333')
fig.tight_layout(rect=(0,0.045,1,1),pad=0.4); fig.savefig('figures/fig13d_perfly_D.pdf'); fig.savefig('figures/fig13d_perfly_D.png',dpi=300)
open('fig13d.log','w').write('\n'.join(L)+'\n')
print(sum(1 for r in res for v in r['hid'].values() if v is not None and abs(v)>1.8),'hidden fly values beyond |1.8|')
