# Plan A Fig A4 "What lies beneath the observed cycle": (a) 20n hidden layer vs plain chain (per-fly dBIC), (b) 20o what the
# HMM's gain is made of (per-fly R_mem, R_drift), (c) 20p excess irreversibility from 3-min memory, (d) 20q TUR ratio.
# Schematic explainers stay in the supplement (Figs 18a-c, 19a-b). Data only; per fly; genotype order rule.
import os,sys,json,numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
COL={'CS':'#2a78d6','CSDD':'#6250d6','per0LD':'#1baf7a','per0DD':'#eda100','fumin':'#eb6834','sss':'#e87ba4'}
NAME={'CS':'CS','CSDD':'CS DD','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','fumin':'$fmn$','sss':'$sss$'}
rng=np.random.default_rng(1)
def mci(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]; b=[np.median(rng.choice(v,len(v))) for _ in range(2000)]; return np.median(v),*np.percentile(b,[2.5,97.5])
def clean(ax):
    for s in ('top','right'): ax.spines[s].set_visible(False)
fig=plt.figure(figsize=(7.0,5.6)); gs=fig.add_gridspec(2,2,hspace=0.55,wspace=0.55)
# (a) 20n
ax=fig.add_subplot(gs[0,0]); N=json.load(open(sqw.OUTD+'n_main.json'))
for k,g in enumerate(('CS','fumin')):
    v=np.array([d['dBIC_T'] for d in N if d['g']==g])   # dBIC_T = BIC(chain) - BIC(HMM) (sqw13); > 0: BIC prefers the HMM
    ax.plot(k+rng.uniform(-0.15,0.15,len(v)),v,'o',ms=3,mfc=COL[g],mec='white',mew=0.4)
    m,lo,hi=mci(v); ax.plot([k-0.28,k+0.28],[m,m],'k',lw=1.2); ax.plot([k,k],[lo,hi],'k',lw=0.8)
    ax.text(k,ax.get_ylim()[1] if False else 0,'',fontsize=6)
ax.axhline(0,color='0.5',lw=0.6); ax.set_xticks([0,1]); ax.set_xticklabels([NAME['CS'],NAME['fumin']]); ax.set_xlim(-0.6,1.6)
ax.set_ylabel('BIC(plain chain) $-$ BIC(HMM)\n(> 0: hidden layer preferred)'); clean(ax)
ax.text(0.5,-0.2,'hidden layer preferred by BIC in 28% (CS), 23% ($fmn$) of flies;\nlikelihood gain ~3 nats per 1000 min',transform=ax.transAxes,ha='center',va='top',fontsize=6.5)
ax.set_title('a   Hidden layer vs plain S/Q/W chain (20n)',loc='left',fontsize=8.5,fontweight='bold',x=-0.32)
# (b) 20o
ax=fig.add_subplot(gs[0,1]); O=json.load(open(sqw.OUTD+'o_memory_drift.json'))
for k,g in enumerate(('CS','fumin')):
    D=[d for d in O if d['g']==g and d['LL']['HMM_12']-d['LL']['MC1_12']>3]
    gain=np.array([d['LL']['HMM_12']-d['LL']['MC1_12'] for d in D])
    Rm=np.array([d['LL']['MC2_12']-d['LL']['MC1_12'] for d in D])/gain; Rd=1-np.array([d['LL']['HMM_6']-d['LL']['MC1_6'] for d in D])/gain
    for off,v,lab in ((-0.18,Rm,'one more minute of memory'),(0.18,Rd,'6-h drift of the rates')):
        x=k+off; vv=np.clip(v,-1.5,5); ax.plot(x+rng.uniform(-0.06,0.06,len(vv)),vv,'o',ms=2.6,mfc=COL[g] if off<0 else 'white',mec=COL[g],mew=0.6)
        m,lo,hi=mci(v); ax.plot([x-0.12,x+0.12],[m,m],'k',lw=1.2); ax.plot([x,x],[max(lo,-1.5),min(hi,5)],'k',lw=0.8)
ax.axhline(1,color='0.5',lw=0.6,ls='--'); ax.axhline(0,color='0.5',lw=0.6)
ax.set_xticks([0,1]); ax.set_xticklabels([NAME['CS'],NAME['fumin']]); ax.set_xlim(-0.6,1.6); ax.set_ylim(-1.6,5.2)
ax.set_ylabel("share of the HMM's likelihood gain\nmatched by ..."); clean(ax)
ax.plot([],[],'o',ms=3,mfc='0.3',mec='0.3',label='one more minute of memory (2nd-order chain)'); ax.plot([],[],'o',ms=3,mfc='white',mec='0.3',label='6-h drift of the rates')
ax.legend(frameon=False,fontsize=6.3,loc='upper center',bbox_to_anchor=(0.5,-0.12),ncol=1); ax.text(0.98,0.98,'values > 5 shown at 5',transform=ax.transAxes,ha='right',va='top',fontsize=6)
ax.set_title("b   What the HMM's gain is made of (20o)",loc='left',fontsize=8.5,fontweight='bold',x=-0.32)
# (c) 20p
ax=fig.add_subplot(gs[1,0]); P_=json.load(open(sqw.OUTD+'p_blockkl.json')); H=json.load(open(sqw.OUTD+'p_blockkl_ddhalves.json')); REF={(d['g'],d['fly']):d['ref'] for d in P_}
def E_of(src,g,p,rk,nul):
    out=[]
    for d in src:
        if d['g']!=g or REF[(g,d['fly'])][rk] is None: continue
        s=np.array([x[p]['D3'] for x in d['sims'][nul]]); out.append(d['obs'][p]['D3']-np.median(s))
    return out
rowsC=[('CS light',P_,'CS',0,'light','A','CS',False),('CS dark',P_,'CS',1,'dark','A','CS',False),('CS DD, ZT 0-12',H,'CSDD',0,'whole','A2','CSDD',False),
       ('CS DD, ZT 12-24',H,'CSDD',1,'whole','A2','CSDD',False),('$per^0$ LD light',P_,'per0LD',0,'light','A','per0LD',False),('$per^0$ LD dark',P_,'per0LD',1,'dark','A','per0LD',False),
       ('$per^0$ DD',P_,'per0DD',0,'whole','A','per0DD',True),('$fmn$ light',P_,'fumin',0,'light','A','fumin',False),('$fmn$ dark',P_,'fumin',1,'dark','A','fumin',True),
       ('$sss$ light',P_,'sss',0,'light','A','sss',False),('$sss$ dark',P_,'sss',1,'dark','A','sss',False)]
for k,(lab,src,g,p,rk,nul,c,q) in enumerate(rowsC):
    m,lo,hi=mci(E_of(src,g,p,rk,nul)); y=len(rowsC)-1-k; ax.plot([lo,hi],[y,y],color=COL[c],lw=1.4); ax.plot(m,y,'o',ms=5,mfc=COL[c] if q else 'white',mec=COL[c],mew=1.2)
ax.axvline(0,color='0.5',lw=0.6); ax.set_yticks(range(len(rowsC))); ax.set_yticklabels([r[0] for r in rowsC][::-1],fontsize=7)
ax.set_xlabel('extra irreversibility from 3-min memory, $E$ (nats h$^{-1}$)'); clean(ax)
ax.text(0.98,0.02,'filled: passes the\npre-registered rule',transform=ax.transAxes,ha='right',va='bottom',fontsize=6.3)
ax.set_title('c   Does the memory carry dissipation? (20p)',loc='left',fontsize=8.5,fontweight='bold',x=-0.32)
# (d) 20q
ax=fig.add_subplot(gs[1,1]); Q=json.load(open(sqw.OUTD+'q_tur.json'))
def ratio(g,p,key):
    out=[]
    for d in Q:
        if d['g']!=g or d['ref'][key] is None: continue
        e=d['ex'][p]; Bo=d['obs'][p]['B180']; Bs=np.array([s[p]['B180'] for s in d['simMC']],float); Bs=Bs[np.isfinite(Bs)]
        if not (np.isfinite(e['sig_mc']) and e['sig_mc']>1e-9 and np.isfinite(e['B_mc_ct']) and e['B_mc_ct']>1e-12 and np.isfinite(Bo) and len(Bs)>=100): continue
        out.append(Bo*e['B_mc']/np.median(Bs)/(e['B_mc']/e['B_mc_ct'])/e['sig_mc'])
    return out
rowsD=[('CS light','CS',0,'light',False),('CS dark','CS',1,'dark',False),('CS DD','CSDD',0,'whole',True),('$per^0$ LD light','per0LD',0,'light',False),('$per^0$ LD dark','per0LD',1,'dark',False),
       ('$per^0$ DD','per0DD',0,'whole',True),('$fmn$ light','fumin',0,'light',False),('$fmn$ dark','fumin',1,'dark',False),('$sss$ light','sss',0,'light',False),('$sss$ dark','sss',1,'dark',False)]
for k,(lab,g,p,key,prim) in enumerate(rowsD):
    m,lo,hi=mci(ratio(g,p,key)); y=len(rowsD)-1-k; ax.plot([lo,hi],[y,y],color=COL[g],lw=1.4); ax.plot(m,y,'o',ms=5,mfc=COL[g] if prim else 'white',mec=COL[g],mew=1.2)
ax.axvline(1,color='k',lw=0.9); ax.text(1.02,len(rowsD)-0.6,'TUR bound',fontsize=6.5)
ax.set_yticks(range(len(rowsD))); ax.set_yticklabels([r[0] for r in rowsD][::-1],fontsize=7); ax.set_xlim(0.2,1.5)
ax.set_xlabel('cycle precision / its observed dissipation, $B_c/\\sigma_{chain}$'); clean(ax)
ax.text(0.98,0.02,'filled: primary (DD);\nopen: LD, descriptive',transform=ax.transAxes,ha='right',va='bottom',fontsize=6.3)
ax.set_title('d   Is the cycle more precise than it pays for? (20q)',loc='left',fontsize=8.5,fontweight='bold',x=-0.32)
for ext in ('png','pdf'): fig.savefig(f'figures/figA4_beneath_the_surface.{ext}',dpi=300,bbox_inches='tight')
print('saved')
