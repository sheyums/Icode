# Plan A figures at the observable level (PRESPEC 20s-A data, s_observable.json). Genotype order rule: wild type -> clock
# mutant -> sleep mutants; light before dark. Fig A1: direction and flux (replaces 6/6b); Fig A2: EPR vs J and Onsager L
# (replaces 11/12); Fig A3: energy landscape (replaces 14/15). Per fly points; medians with fly-bootstrap 95% CIs.
import os,sys,json,numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
O=json.load(open(sqw.OUTD+'s_observable.json')); rng=np.random.default_rng(20)
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
COL={'CS':'#2a78d6','CSDD':'#6250d6','per0LD':'#1baf7a','per0DD':'#eda100','fumin':'#eb6834','sss':'#e87ba4'}
NAME={'CS':'CS','CSDD':'CS DD','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','fumin':'$fmn$','sss':'$sss$'}
GO=['CS','CSDD','per0LD','per0DD','fumin','sss']; DDg={'CSDD','per0DD'}
BL=[(g,p) for g in GO for p in ([0] if g in DDg else [0,1])]
def rows(g,p,cyc=False):
    out=[]
    for d in O:
        if d['g']!=g: continue
        m=d['main'][p]
        if not m['ok']: continue
        if cyc and (m['J_chain'] is None or abs(m['J_chain'])<1e-12): continue
        out.append(dict(m,fly=d['fly']))
    return out
def med_ci(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]; b=[np.median(rng.choice(v,len(v))) for _ in range(2000)]; return np.median(v),*np.percentile(b,[2.5,97.5])
def strip(ax,key,ylabel,hmm=None,ref0=True):
    x=0; ticks=[]; labs=[]
    for g,p in BL:
        R_=rows(g,p); v=np.array([r[key] for r in R_],float); fill=(p==1 or g in DDg)
        ax.plot(x+rng.uniform(-0.18,0.18,len(v)),v,'o',ms=2.6,mfc=COL[g] if fill else 'white',mec=COL[g],mew=0.6,alpha=0.8,zorder=2)
        m,lo,hi=med_ci(v); ax.plot([x-0.3,x+0.3],[m,m],color='k',lw=1.2,zorder=3); ax.plot([x,x],[lo,hi],color='k',lw=0.8,zorder=3)
        if hmm:
            hv=[r[hmm] for r in R_ if hmm in r and not r.get('hmm_deg',True)]
            if hv: ax.plot(x+0.36,np.median(hv),marker='<',ms=4,color='0.35',zorder=3)
        ticks.append(x); labs.append(NAME[g]+('' if g in DDg else ('\nlight' if p==0 else '\ndark'))); x+=1+(0.5 if (p==1 or g in DDg) else 0)
    if ref0: ax.axhline(0,color='0.6',lw=0.5,zorder=1)
    ax.set_xticks(ticks); ax.set_xticklabels(labs,fontsize=6.5); ax.set_ylabel(ylabel)
    for s in ('top','right'): ax.spines[s].set_visible(False)
# ---------------- Fig A1 ----------------
fig,axs=plt.subplots(3,1,figsize=(7.0,6.6),constrained_layout=True)
strip(axs[0],'c_obs','observable current $c_{obs}$\n(P(Q before sleep) $-$ P(Q after waking))')
strip(axs[1],'J','cycle flux $J$ (cycles h$^{-1}$)\n(forward W$\\to$Q$\\to$S$\\to$W $>0$)',hmm='J_hmm')
strip(axs[2],'EPR','entropy production (nats h$^{-1}$)\ncount estimate, rule of three',hmm='EPR_hmm',ref0=False)
for a,l in zip(axs,'abc'): a.text(-0.07,1.02,l,transform=a.transAxes,fontsize=10,fontweight='bold')
axs[0].plot([],[],'o',ms=3,mfc='white',mec='k',mew=0.6,label='light (one point per fly)'); axs[0].plot([],[],'o',ms=3,mfc='k',mec='k',label='dark / DD')
axs[0].plot([],[],'-',color='k',lw=1.2,label='median, 95% CI'); axs[0].plot([],[],'<',ms=4,color='0.35',ls='none',label='HMM (model-based upper bracket), median')
axs[0].legend(frameon=False,fontsize=6.5,ncol=4,loc='lower left',bbox_to_anchor=(0.04,1.03))
for ext in ('png','pdf'): fig.savefig(f'figures/figA1_direction_flux_observable.{ext}',dpi=300)
# ---------------- Fig A2 ----------------
fig=plt.figure(figsize=(7.0,5.0)); gs=fig.add_gridspec(2,4,width_ratios=[1,1,1,1.25],hspace=0.45,wspace=0.45); Ls={}
for k,g in enumerate(GO):
    ax=fig.add_subplot(gs[k//3,k%3]); pts=[]
    for p in ([0] if g in DDg else [0,1]):
        for r in rows(g,p,cyc=True):
            if np.isfinite(r['EPR']): pts.append((r['fly'],p,r['J'],r['EPR']))
    J=np.array([q[2] for q in pts]); E=np.array([q[3] for q in pts]); fill=np.array([q[1]==1 or g in DDg for q in pts])
    ax.plot(J[~fill],E[~fill],'o',ms=3,mfc='white',mec=COL[g],mew=0.7); ax.plot(J[fill],E[fill],'o',ms=3,mfc=COL[g],mec=COL[g])
    s=(J**2@E)/(J**2@J**2); xx=np.linspace(min(-0.05,J.min()),max(0.4,J.max()),100); ax.plot(xx,s*xx**2,color='k',lw=0.9)
    flies=sorted({q[0] for q in pts}); B=[]
    for _ in range(2000):
        pick=rng.choice(flies,len(flies)); qq=[q for f in pick for q in pts if q[0]==f]; Jb=np.array([q[2] for q in qq]); Eb=np.array([q[3] for q in qq]); B.append((Jb**2@Jb**2)/(Jb**2@Eb))
    Ls[g]=(1/s,*np.percentile(B,[2.5,97.5]))
    ax.set_title(NAME[g],fontsize=8,color=COL[g]); ax.set_xlabel('$J$ (cycles h$^{-1}$)'); ax.set_ylabel('EPR (nats h$^{-1}$)') if k%3==0 else None
    ax.text(0.04,0.92,f'$L$ = {1/s:.2f}',transform=ax.transAxes,fontsize=7)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax=fig.add_subplot(gs[:,3])
for k,g in enumerate(GO):
    L,lo,hi=Ls[g]; y=len(GO)-1-k; ax.plot([lo,hi],[y,y],color=COL[g],lw=1.6); ax.plot(L,y,'o',ms=5,color=COL[g])
ax.set_yticks(range(len(GO))); ax.set_yticklabels([NAME[g] for g in GO][::-1]); ax.set_xlabel('Onsager conductance $L$\n(cycles$^2$ nats$^{-1}$ h$^{-1}$; EPR $=J^2/L$)')
ax.axvline(Ls['CS'][0],color='0.6',lw=0.5,ls='--')
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.text(0.01,0.97,'a',fontsize=10,fontweight='bold'); fig.text(0.70,0.97,'b',fontsize=10,fontweight='bold')
for ext in ('png','pdf'): fig.savefig(f'figures/figA2_epr_j_onsager_observable.{ext}',dpi=300,bbox_inches='tight')
# ---------------- Fig A3 ----------------
fig,axs=plt.subplots(1,3,figsize=(7.0,2.8),gridspec_kw=dict(width_ratios=[1.3,1.3,1]),constrained_layout=True)
def strip_cyc(ax,key,ylabel):
    x=0; ticks=[]; labs=[]
    for g,p in BL:
        v=np.array([r[key] for r in rows(g,p,cyc=True)],float); fill=(p==1 or g in DDg)
        ax.plot(x+rng.uniform(-0.18,0.18,len(v)),v,'o',ms=2.4,mfc=COL[g] if fill else 'white',mec=COL[g],mew=0.6,alpha=0.8)
        m,lo,hi=med_ci(v); ax.plot([x-0.3,x+0.3],[m,m],color='k',lw=1.2); ax.plot([x,x],[lo,hi],color='k',lw=0.8)
        ticks.append(x); labs.append(NAME[g]+('' if g in DDg else (' L' if p==0 else ' D'))); x+=1+(0.4 if fill else 0)
    ax.axhline(0,color='0.6',lw=0.5); ax.set_xticks(ticks); ax.set_xticklabels(labs,fontsize=6,rotation=90); ax.set_ylabel(ylabel)
    for s in ('top','right'): ax.spines[s].set_visible(False)
strip_cyc(axs[0],'E21','$E_W - E_Q$ = ln($\\pi_Q/\\pi_W$) (nats)'); strip_cyc(axs[1],'E31','$E_S - E_Q$ = ln($\\pi_Q/\\pi_S$) (nats)')
ax=axs[2]
for g,p in BL:
    R_=rows(g,p,cyc=True); fill=(p==1 or g in DDg)
    for kp,kc,mk in (('F_SW','Xn_cnt','o'),('F_WS','Yn_cnt','s')):
        a=np.array([r[kp] for r in R_]); b=np.array([r[kc] for r in R_]); m=np.isfinite(a)&np.isfinite(b)
        ax.plot(a[m],b[m],mk,ms=2.2,mfc=COL[g] if fill else 'white',mec=COL[g],mew=0.5,alpha=0.7)
ax.plot([-1,1],[-1,1],color='k',lw=0.7); ax.set_xlim(-1,1); ax.set_ylim(-1,1); ax.set_aspect('equal')
ax.set_xlabel('landscape prediction'); ax.set_ylabel('counted routing (observed)'); ax.text(0.05,-0.88,'o S$\\to$W ($X_n$)\n□ W$\\to$S ($Y_n$)',fontsize=6)
for s in ('top','right'): ax.spines[s].set_visible(False)
for a,l in zip(axs,'abc'): a.text(-0.12,1.03,l,transform=a.transAxes,fontsize=10,fontweight='bold')
for ext in ('png','pdf'): fig.savefig(f'figures/figA3_landscape_observable.{ext}',dpi=300)
print('saved; L:',{g:tuple(round(x,3) for x in v) for g,v in Ls.items()})
