# Figures 1-3 for PRESPEC 18 V2b (bridging (b), seconds 61-300). DESCRIPTIVE window (chosen post hoc).
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in',
    'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':3,'ytick.major.size':3,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'pdf.fonttype':42,'savefig.dpi':300})
EXIT='#1f5f99'; APP='#d4572a'; INK='#1a1a1a'; GREY='#8a8a8a'
LAB={'CS':'Wild type (CS)','fumin':r'$\it{fumin}$','per0LD':r'$\it{per}^{0}$ LD','per0DD':r'$\it{per}^{0}$ DD','sss':r'$\it{sleepless}$'}
R=json.load(open('v18b.json')); V='V2b'
blocks=[('CS','light'),('CS','dark'),('fumin','light'),('fumin','dark'),('per0LD','light'),('per0LD','dark'),('per0DD','whole'),('sss','light'),('sss','dark')]
# Fig 1
fig,ax=plt.subplots(figsize=(7.0,2.6)); rng=np.random.default_rng(0); pos=[0,1,2.6,3.6,5.2,6.2,7.8,9.4,10.4]
for x,(g,bn) in zip(pos,blocks):
    r=R[f'{V}|{g}|{bn}']; D=np.array(r['D'])*100
    ax.scatter(x+rng.uniform(-.18,.18,len(D)),D,s=6,color=GREY,lw=0,alpha=.8); ax.plot([x-.3,x+.3],[np.median(D)]*2,color=INK,lw=1.4)
    ax.text(x,1.03,f"{100*np.mean(D>0):.0f}%{'*' if r['sign_p']<0.0056 else ''}",ha='center',fontsize=7,transform=ax.get_xaxis_transform())
ax.axhline(0,color=INK,lw=.5); ax.set_xticks(pos); ax.set_xticklabels([b if g!='per0DD' else 'all' for g,b in blocks],fontsize=7)
for x0,g in [(0.5,'CS'),(3.1,'fumin'),(5.7,'per0LD'),(7.8,'per0DD'),(9.9,'sss')]: ax.text(x0,-0.19,LAB[g],ha='center',fontsize=8,transform=ax.get_xaxis_transform())
ax.set_ylim(-15,20); ax.set_ylabel(r'$\Delta_{5} = \langle a_{\rm exit}-a_{\rm approach}\rangle_{1-5\,{\rm min}}$ (%)'); ax.text(-0.45,1.03,'flies with $\\Delta_5>0$:',ha='right',fontsize=7,transform=ax.get_xaxis_transform())
fig.subplots_adjust(left=.13,right=.99,bottom=.24,top=.88); fig.savefig('figures/fig1_pm5_bridged.pdf'); fig.savefig('figures/fig1_pm5_bridged.png'); plt.close(fig)
# Fig 2
sel=[('CS','dark','Wild type, dark'),('per0DD','whole',r'$\it{per}^{0}$, constant dark'),('fumin','light',r'$\it{fumin}$, light'),('sss','dark',r'$\it{sleepless}$, dark')]
fig,axs=plt.subplots(1,4,figsize=(7.0,2.2),sharey=True); t=(60+30*np.arange(8)+15)/60
for ax,(g,bn,title) in zip(axs,sel):
    r=R[f'{V}|{g}|{bn}']; EX=np.array(r['EX'])*100; AP=np.array(r['AP'])*100; n=len(EX)
    for M,c,lab,ls in [(EX,EXIT,'after waking','-'),(AP,APP,'before sleep','--')]:
        m=M.mean(0); s=M.std(0,ddof=1)/np.sqrt(n); ax.fill_between(t,m-s,m+s,color=c,alpha=.18,lw=0); ax.plot(t,m,color=c,lw=1.2,ls=ls,marker='o',ms=2,label=lab)
    ax.set_title(title,fontsize=8,loc='left'); ax.set_xlabel(r'$|t - t_{\rm transition}|$ (min)'); ax.set_xticks([1,2,3,4,5]); ax.set_xlim(1,5)
    ax.text(.97,.04,f"$n={n}$ flies",transform=ax.transAxes,ha='right',fontsize=7)
axs[0].set_ylabel('active seconds (%)'); h,l=axs[0].get_legend_handles_labels(); fig.legend(h,l,loc='upper center',ncol=2,fontsize=7.5,handlelength=2.2,bbox_to_anchor=(0.55,1.0))
fig.subplots_adjust(left=.08,right=.99,bottom=.22,top=.80,wspace=.12); fig.savefig('figures/fig2_pm5_bridged.pdf'); fig.savefig('figures/fig2_pm5_bridged.png'); plt.close(fig)
# Fig 3
fig,ax=plt.subplots(figsize=(3.4,3.9)); mk={'CS':'o','fumin':'s','per0LD':'^','per0DD':'D','sss':'v'}; lim=[-7,10]
ax.plot(lim,lim,color=GREY,lw=.6,ls=':')
for g,bn in blocks:
    r=R[f'{V}|{g}|{bn}']; xi=np.array(r['x'])*100; yo=np.array(r['y'])*100; n=len(xi)
    ax.errorbar(xi.mean(),yo.mean(),xerr=xi.std(ddof=1)/np.sqrt(n),yerr=yo.std(ddof=1)/np.sqrt(n),fmt=mk[g],ms=5,mfc=INK if bn in('dark','whole') else 'white',mec=INK,ecolor=INK,elinewidth=.6,capsize=0,label=LAB[g] if bn in('dark','whole') else None)
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
ax.set_xlabel(r'quieting before sleep, $a_{\rm mid}-a_{\rm approach}$ (%)'); ax.set_ylabel(r'quieting after waking, $a_{\rm mid}-a_{\rm exit}$ (%)')
ax.text(-5.6,-4.6,'time-reversal symmetric',rotation=45,rotation_mode='anchor',ha='left',va='bottom',fontsize=7,color=GREY)
ax.legend(fontsize=7,loc='upper center',bbox_to_anchor=(0.5,-0.17),ncol=3,handletextpad=.3,columnspacing=.8); ax.text(-6.6,9.6,'filled: dark ($per^0$ DD: whole record)\nopen: light\nwindow: 1-5 min from transition',fontsize=7,va='top')
fig.subplots_adjust(left=.16,right=.97,bottom=.27,top=.98); fig.savefig('figures/fig3_pm5_bridged.pdf'); fig.savefig('figures/fig3_pm5_bridged.png'); plt.close(fig)
print('ok')
