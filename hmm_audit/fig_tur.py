# Explainer figure for PRESPEC 20q (TUR). Panels a-b schematic; c-d data (q_tur.json, tau = 3 h, amendment D set).
import os,sys,json,numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
INK='#0b0b0b'; INK2='#52514e'; GRID='#e4e3df'; COL={'CSDD':'#2a78d6','per0DD':'#eda100'}; ACC='#eb6834'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2})
fig=plt.figure(figsize=(12,10.2),facecolor='white')
def clean(ax):
    for s in ('top','right'): ax.spines[s].set_visible(False)
# ---------- (a) the counter and its fan ----------
ax=fig.add_axes([0.07,0.60,0.38,0.31]); rng=np.random.default_rng(3); t=np.arange(0,48*60)
for k,(J,D,col,lab) in enumerate(((0.5,0.6,INK,'precise cycle (small D)'),(0.5,4.0,ACC,'sloppy cycle (large D)'))):
    for _ in range(6):
        steps=rng.normal(J/60,np.sqrt(2*D/60),len(t)); ax.plot(t/60,np.cumsum(steps),color=col,lw=0.8,alpha=0.55)
    th=t/60; ax.fill_between(th,J*th-2*np.sqrt(2*D*th),J*th+2*np.sqrt(2*D*th),color=col,alpha=0.07,lw=0)
    ax.plot([],[],color=col,lw=2,label=lab)
ax.plot(t/60,0.5*t/60,'--',color=INK2,lw=1); ax.text(47,-14,'dashed: mean slope = current J',fontsize=8,color=INK2,ha='right')
ax.set_xlabel('time (h)'); ax.set_ylabel('net cycles completed, N(t)\n(forward W→Q→S→W counts +1)'); ax.set_xlim(0,48); clean(ax); ax.grid(color=GRID,lw=0.6)
ax.legend(frameon=False,fontsize=8,loc='upper left',labelcolor=INK)
ax.set_title('a   A cycle counter: drift J, spread grows as √(2Dt)  (illustration)',loc='left',fontsize=10.5,fontweight='bold',color=INK,x=-0.14)
# ---------- (b) the price list ----------
ax=fig.add_axes([0.57,0.60,0.38,0.31]); x=np.logspace(-3,0,50)
ax.fill_between(x,x,1.5,color=ACC,alpha=0.10,lw=0); ax.fill_between(x,1e-3,x,color='#2a78d6',alpha=0.06,lw=0)
ax.plot(x,x,color=INK,lw=1.6); ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e-3,1); ax.set_ylim(1e-3,1)
ax.text(0.0016,0.3,'FORBIDDEN for any Markov process\nwith this dissipation:\nmore precise than it pays for\n→ hidden dissipation would be required',fontsize=8,color=INK,va='top')
ax.text(0.25,0.004,'ALLOWED\n(precision ≤ cost)',fontsize=8.5,color=INK,ha='center')
ax.text(0.13,0.16,'TUR: σ = J²/D',fontsize=8.5,color=INK,rotation=35,ha='center',va='bottom')
ax.plot(0.05,0.049,'s',ms=8,mfc='white',mec=INK,mew=1.6); ax.annotate('the observed S/Q/W chain\nsits ON the line\n(near-equilibrium cycle)',xy=(0.05,0.049),xytext=(0.075,0.012),fontsize=7.8,color=INK,arrowprops=dict(arrowstyle='-',color=INK2,lw=0.8))
ax.plot(0.05,0.12,'o',ms=8,mfc=ACC,mec=ACC); ax.annotate('data more regular than the chain\n→ violation',xy=(0.05,0.12),xytext=(0.0016,0.035),fontsize=7.8,color=INK,arrowprops=dict(arrowstyle='-',color=INK2,lw=0.8))
ax.set_xlabel('cost: entropy production σ (k$_B$/h)'); ax.set_ylabel('precision bought: J²/D (k$_B$/h)'); clean(ax)
ax.set_title('b   The thermodynamic price list  (schematic)',loc='left',fontsize=10.5,fontweight='bold',color=INK,x=-0.14)
# ---------- data ----------
O=json.load(open(sqw.OUTD+'q_tur.json'))
def rows(g,p,key):
    out=[]
    for d in O:
        if d['g']!=g or d['ref'][key] is None: continue
        e=d['ex'][p]; Bo=d['obs'][p]['B180']; Bs=np.array([s[p]['B180'] for s in d['simMC']],float); Bs=Bs[np.isfinite(Bs)]
        if not (np.isfinite(e['sig_mc']) and e['sig_mc']>1e-9 and np.isfinite(e['B_mc_ct']) and e['B_mc_ct']>1e-12 and np.isfinite(Bo) and len(Bs)>=100): continue
        kap=e['B_mc']/e['B_mc_ct']; Bc=Bo*e['B_mc']/np.median(Bs)/kap
        out.append(dict(Bc=Bc,smc=e['sig_mc'],Bmc=e['B_mc'],sh=e['sig_hmm'] if not d['deg'] else np.nan))
    return out
# ---------- (c) per fly on the price list ----------
ax=fig.add_axes([0.07,0.15,0.38,0.34]); x=np.logspace(-3.5,0.3,50)
ax.fill_between(x,x,3,color=ACC,alpha=0.10,lw=0); ax.plot(x,x,color=INK,lw=1.4)
for g,name in (('per0DD','per0 DD'),('CSDD','CS DD')):
    R=rows(g,0,'whole'); c=COL[g]
    for r in R:
        if np.isfinite(r['sh']) and r['sh']>0: ax.plot([r['smc'],r['sh']],[r['Bc'],r['Bc']],color=c,lw=0.7,alpha=0.45)
    ax.plot([r['smc'] for r in R],[r['Bc'] for r in R],'o',ms=6.5,mfc=c,mec='white',mew=1.2,label=f'{name}: data precision vs observed-chain cost (n {len(R)})')
    hs=[r for r in R if np.isfinite(r['sh']) and r['sh']>0]
    ax.plot([r['sh'] for r in hs],[r['Bc'] for r in hs],'o',ms=5.5,mfc='white',mec=c,mew=1.3)
ax.plot([],[],'o',ms=5.5,mfc='white',mec=INK2,mew=1.3,label='same flies vs the HMM\'s (larger) cost')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(3e-3,3); ax.set_ylim(3e-4,2)
ax.text(0.25,6e-4,'Single flies scatter onto both sides:\nD from 4 days is noisy, and no fly is\nbeyond its own simulated null (S1).\nThe test is the median over flies (panel d).',fontsize=7.6,color=INK,ha='left',va='bottom')
ax.set_xlabel('cost: entropy production σ (k$_B$/h), per fly'); ax.set_ylabel('precision bought: J²/D from the data\n(bias-corrected, k$_B$/h), per fly'); clean(ax)
ax.legend(frameon=False,fontsize=7.4,loc='upper left',labelcolor=INK)
ax.set_title('c   Each DD fly on the price list (20q data)',loc='left',fontsize=10.5,fontweight='bold',color=INK,x=-0.14)
# ---------- (d) block summary ----------
BL=[('per0 DD',('per0DD',0,'whole'),True),('CS DD',('CSDD',0,'whole'),True),('CS LD, light',('CS',0,'light'),False),('CS LD, dark',('CS',1,'dark'),False),
    ('per0 LD, light',('per0LD',0,'light'),False),('per0 LD, dark',('per0LD',1,'dark'),False),('fumin, light',('fumin',0,'light'),False),('fumin, dark',('fumin',1,'dark'),False),
    ('sss, light',('sss',0,'light'),False),('sss, dark',('sss',1,'dark'),False)]
ax=fig.add_axes([0.66,0.15,0.31,0.34]); rg=np.random.default_rng(1)
for k,(lab,a,prim) in enumerate(BL):
    v=np.array([r['Bc']/r['smc'] for r in rows(*a)]); m=np.median(v); b=[np.median(rg.choice(v,len(v))) for _ in range(2000)]; lo,hi=np.percentile(b,[2.5,97.5])
    y=len(BL)-1-k; ax.plot([lo,hi],[y,y],color=INK2,lw=1.6,solid_capstyle='round')
    ax.plot(m,y,'o',ms=8,mfc=INK if prim else 'white',mec=INK if prim else INK2,mew=1.6)
ax.axvline(1,color=ACC,lw=1.4); ax.text(1.02,-0.6,'TUR line',color=INK,fontsize=8,va='center')
ax.set_yticks(range(len(BL))); ax.set_yticklabels([b[0] for b in BL][::-1],fontsize=8.5); ax.set_xlim(0.2,1.5); ax.set_ylim(-0.9,len(BL)-0.4)
ax.set_xlabel('precision / observed-chain cost, B$_c$/σ$_{MC}$\n(median over flies of per-fly ratios, 95% CI)'); ax.grid(axis='x',color=GRID,lw=0.6); clean(ax)
ax.plot([],[],'o',ms=7,mfc=INK,mec=INK,label='primary (DD)'); ax.plot([],[],'o',ms=7,mfc='white',mec=INK2,mew=1.6,label='descriptive (LD)')
ax.legend(frameon=False,fontsize=7.6,loc='lower right',bbox_to_anchor=(1.03,1.0),ncol=2,labelcolor=INK,handletextpad=0.3,columnspacing=1.0)
ax.set_title('d   No block median lies above the line',loc='left',fontsize=10.5,fontweight='bold',color=INK,x=-0.40,pad=22)
fig.text(0.07,0.015,'Reading: the cycle is never more regular than the observed chain\'s own dissipation allows (no hidden dissipation required). per0 DD sits close to the line;\nCS DD and LD fall below it because slow modulation (free-running rhythm, day/phase drift) adds noise. Against the HMM\'s larger cost (open circles) the flies use only ~half.\nRule of three for one-way edges (fumin, sss) moves those LD ratios to 0.73-0.84, still below 1.',fontsize=7.8,color=INK2,va='bottom')
for ext in ('png','pdf'): fig.savefig(f'figures/fig19_tur.{ext}',dpi=200,facecolor='white')
print('saved')
