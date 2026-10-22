# Journal-style figures from PRESPEC 15-17 quantities. Fig 3 uses one DESCRIPTIVE post hoc quantity
# (exit depth, the mirror of the pre-registered approach depth), labelled as such.
import glob,json,numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from multiprocessing import Pool
import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit"); import taper as T
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in',
    'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':3,'ytick.major.size':3,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'pdf.fonttype':42,'savefig.dpi':300})
EXIT='#1f5f99'; APP='#d4572a'; INK='#1a1a1a'; GREY='#8a8a8a'
LAB={'CS':'Wild type (CS)','fumin':r'$\it{fumin}$','per0LD':r'$\it{per}^{0}$ LD','per0DD':r'$\it{per}^{0}$ DD','sss':r'$\it{sleepless}$'}
J=json.load(open('taper.json'))
def per_fly(g):
    X=[np.loadtxt(p).astype(int) for p in sorted(glob.glob(T.CD+g+'/**/*.txt',recursive=True))]
    with Pool(4) as pl: B=pl.map(T.fly_bouts,X)
    return B
B={g:per_fly(g) for g in T.GROUPS}
def fly_profiles(g,ph):
    out=[]
    for bs in B[g]:
        bb=[b for b in bs if ph is None or b['ph']==ph]
        if len(bb)>=3:
            mv=lambda k: np.mean([b[k] for b in bb],0)
            ref=np.mean([b['ref'] for b in bb]); out.append(dict(ex=mv('ex10'),ap=mv('ap10'),ref=ref,
                dep_in=np.mean([b['depth'] for b in bb]), dep_out=ref-np.mean([np.mean(b['ex10'][:2]) for b in bb])))
    return out
# ---------- Figure 1: per-fly asymmetry
blocks=[('CS','light',0),('CS','dark',1),('fumin','light',0),('fumin','dark',1),('per0LD','light',0),('per0LD','dark',1),('per0DD','whole',None),('sss','light',0),('sss','dark',1)]
fig,ax=plt.subplots(figsize=(7.0,2.6)); rng=np.random.default_rng(0); xs=[]
pos=[0,1,2.6,3.6,5.2,6.2,7.8,9.4,10.4]
for x,(g,bn,ph) in zip(pos,blocks):
    r=J[f'{g}|{bn}']; D=np.array(r['D'])*100
    ax.scatter(x+rng.uniform(-.18,.18,len(D)),D,s=6,color=GREY,lw=0,alpha=.8)
    ax.plot([x-.3,x+.3],[np.median(D)]*2,color=INK,lw=1.4)
    star='*' if r['sign_p']<0.0056 else ''
    ax.text(x,13.2,f"{100*np.mean(D>0):.0f}%{star}",ha='center',fontsize=7)
    xs.append((x,bn if g!='per0DD' else 'all'))
ax.axhline(0,color=INK,lw=.5)
ax.set_xticks([p for p,_ in xs]); ax.set_xticklabels([b for _,b in xs],fontsize=7)
for x0,g in [(0.5,'CS'),(3.1,'fumin'),(5.7,'per0LD'),(7.8,'per0DD'),(9.9,'sss')]: ax.text(x0,-0.19,LAB[g],ha='center',fontsize=8,transform=ax.get_xaxis_transform())
ax.set_ylim(-10,14); ax.set_ylabel(r'$\Delta = \langle a_{\rm exit}-a_{\rm approach}\rangle$ (%)')
ax.text(-0.9,13.2,'flies with $\\Delta>0$:',ha='right',fontsize=7)
fig.subplots_adjust(left=.13,right=.99,bottom=.24,top=.93); fig.savefig('figures/fig1_perfly_asymmetry.pdf'); fig.savefig('figures/fig1_perfly_asymmetry.png'); plt.close(fig)
# ---------- Figure 2: approach vs exit profiles
sel=[('CS','dark',1,'Wild type, dark'),('per0DD','whole',None,r'$\it{per}^{0}$, constant dark'),('fumin','light',0,r'$\it{fumin}$, light'),('sss','dark',1,r'$\it{sleepless}$, dark')]
fig,axs=plt.subplots(1,4,figsize=(7.0,2.2),sharey=True); t=np.arange(2,11)
for ax,(g,bn,ph,title) in zip(axs,sel):
    P=fly_profiles(g,ph); EX=np.array([p['ex'] for p in P])*100; AP=np.array([p['ap'] for p in P])*100; n=len(P)
    for M,c,lab,ls in [(EX,EXIT,'after waking','-'),(AP,APP,'before sleep','--')]:
        m=M.mean(0); s=M.std(0,ddof=1)/np.sqrt(n); ax.fill_between(t,m-s,m+s,color=c,alpha=.18,lw=0); ax.plot(t,m,color=c,lw=1.2,ls=ls,label=lab)
    ax.set_title(title,fontsize=8,loc='left'); ax.set_xlabel(r'$|t - t_{\rm transition}|$ (min)'); ax.set_xticks([2,4,6,8,10])
    ax.text(.97,.04,f"$n={n}$; $p={J[f'{g}|{bn}']['sign_p']:.1g}$" if False else f"$n={n}$ flies",transform=ax.transAxes,ha='right',fontsize=7)
axs[0].set_ylabel('active seconds (%)'); h,l=axs[0].get_legend_handles_labels(); fig.legend(h,l,loc='upper center',ncol=2,fontsize=7.5,handlelength=2.2,bbox_to_anchor=(0.55,1.0))
fig.subplots_adjust(left=.08,right=.99,bottom=.22,top=.80,wspace=.12); fig.savefig('figures/fig2_profiles.pdf'); fig.savefig('figures/fig2_profiles.png'); plt.close(fig)
# ---------- Figure 3: time-reversal symmetry plane (approach depth vs exit depth)
fig,ax=plt.subplots(figsize=(3.4,3.9)); mk={'CS':'o','fumin':'s','per0LD':'^','per0DD':'D','sss':'v'}; lim=[-7,10]
ax.plot(lim,lim,color=GREY,lw=.6,ls=':'); 
for g,bn,ph in blocks:
    P=fly_profiles(g,ph); di=np.array([p['dep_in'] for p in P])*100; do=np.array([p['dep_out'] for p in P])*100; n=len(P)
    mi,mo=di.mean(),do.mean(); si,so=di.std(ddof=1)/np.sqrt(n),do.std(ddof=1)/np.sqrt(n)
    fc=INK if bn in('dark','whole') else 'white'
    ax.errorbar(mi,mo,xerr=si,yerr=so,fmt=mk[g],ms=5,mfc=fc,mec=INK,ecolor=INK,elinewidth=.6,capsize=0,label=f"{LAB[g]}" if bn in ('dark','whole') else None)
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
ax.set_xlabel(r'quieting before sleep, $a_{\rm mid}-a_{\rm approach}$ (%)'); ax.set_ylabel(r'quieting after waking, $a_{\rm mid}-a_{\rm exit}$ (%)')
ax.text(-5.6,-4.6,'time-reversal symmetric',rotation=45,rotation_mode='anchor',ha='left',va='bottom',fontsize=7,color=GREY)
ax.legend(fontsize=7,loc='upper center',bbox_to_anchor=(0.5,-0.17),ncol=3,handletextpad=.3,columnspacing=.8); ax.text(-6.6,9.6,'filled: dark ($per^0$ DD: whole record)\nopen: light',fontsize=7,va='top')
fig.subplots_adjust(left=.16,right=.97,bottom=.27,top=.98); fig.savefig('figures/fig3_symmetry_plane.pdf'); fig.savefig('figures/fig3_symmetry_plane.png'); plt.close(fig)
# ---------- Figure 4 (DESCRIPTIVE; replaces the invalid fits of PRESPEC 17 C): average relaxation 2-20 min
fig,axs=plt.subplots(1,5,figsize=(7.0,2.0),sharey=True); t20=np.arange(2,21); rows4=[]
for ax,(g,bn,ph) in zip(axs,[('CS','dark',1),('per0DD','whole',None),('per0LD','dark',1),('fumin','light',0),('sss','dark',1)]):
    EX=[];AP=[];RF=[]
    for bs in B[g]:
        bb=[b for b in bs if (ph is None or b['ph']==ph) and b['ex20'] is not None]
        if bb: EX.append(np.mean([b['ex20'] for b in bb],0)); AP.append(np.mean([b['ap20'] for b in bb],0)); RF.append(np.mean([b['ref'] for b in bb]))
    EX=np.array(EX)*100; AP=np.array(AP)*100; n=len(EX); rf=np.mean(RF)*100
    for M,c,ls,lab in [(EX,EXIT,'-','after waking'),(AP,APP,'--','before sleep')]:
        m=M.mean(0); se=M.std(0,ddof=1)/np.sqrt(n); ax.fill_between(t20,m-se,m+se,color=c,alpha=.18,lw=0); ax.plot(t20,m,color=c,ls=ls,lw=1.1,label=lab)
    ax.axhline(rf,color=GREY,lw=.6,ls=':'); ax.set_title({'CS':'Wild type, dark','per0DD':'$per^0$, DD','per0LD':'$per^0$ LD, dark','fumin':'$fumin$, light','sss':'$sleepless$, dark'}[g],fontsize=7.5,loc='left')
    ax.text(.97,.04,f"$n={n}$",transform=ax.transAxes,ha='right',fontsize=7); ax.set_xticks([2,10,20]); ax.set_xlabel('min from transition')
    rows4.append((g,bn,n))
axs[0].set_ylabel('active seconds (%)'); h,l=axs[0].get_legend_handles_labels(); h.append(plt.Line2D([],[],color=GREY,ls=':',lw=.6)); l.append('mid-wake level')
fig.legend(h,l,loc='upper center',ncol=3,fontsize=7,bbox_to_anchor=(0.55,1.0)); fig.subplots_adjust(left=.08,right=.99,bottom=.22,top=.76,wspace=.1)
fig.savefig('figures/fig4_relaxation.pdf'); fig.savefig('figures/fig4_relaxation.png'); plt.close(fig); print('fig4 flies (bouts >= 40 min):',rows4)
# numbers for the record
for g,bn,ph in blocks:
    P=fly_profiles(g,ph); print(f"{g:7s} {bn:5s} n={len(P)} depth_in {100*np.mean([p['dep_in'] for p in P]):+.2f}% depth_out {100*np.mean([p['dep_out'] for p in P]):+.2f}% (mean over flies)")
