# Figure 5: signed cycle current (PRESPEC 19; descriptive window).
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in',
    'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':3,'ytick.major.size':3,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'pdf.fonttype':42,'savefig.dpi':300})
INK='#1a1a1a'; GREY='#9a9a9a'; FWD='#1f5f99'; REV='#d4572a'
R=json.load(open('cycle.json'))
blocks=[('CS','light'),('CS','dark'),('fumin','light'),('fumin','dark'),('per0LD','light'),('per0LD','dark'),('per0DD','whole'),('sss','light'),('sss','dark')]
pos=[0,1,2.6,3.6,5.2,6.2,7.8,9.4,10.4]
LAB={'CS':'Wild type','fumin':r'$\it{fumin}$','per0LD':r'$\it{per}^{0}$ LD','per0DD':r'$\it{per}^{0}$ DD','sss':r'$\it{sleepless}$'}
fig=plt.figure(figsize=(7.0,4.4)); gs=fig.add_gridspec(2,2,width_ratios=[1,3.2],height_ratios=[1,1],wspace=.18,hspace=.35)
# (a) cycle diagram
ax=fig.add_subplot(gs[:,0]); ax.set_xlim(-1.3,1.3); ax.set_ylim(-1.8,1.5); ax.set_aspect('equal'); ax.axis('off')
P={'S':(0,1.0),'Q':(-0.87,-0.5),'A':(0.87,-0.5)}; NM={'S':'sleep','Q':'quiet\nwake','A':'active\nwake'}
for k,(x,y) in P.items(): ax.add_patch(Circle((x,y),0.32,fc='white',ec=INK,lw=.8)); ax.text(x,y,k,ha='center',va='center',fontsize=11); ax.text(x,y-0.62 if k!='S' else y+0.45,NM[k],ha='center',va='center',fontsize=7)
def arr(a,b,c,rad):
    (x1,y1),(x2,y2)=P[a],P[b]; ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),connectionstyle=f"arc3,rad={rad}",arrowstyle='-|>',mutation_scale=9,lw=1.1,color=c,shrinkA=13,shrinkB=13))
for a,b in [('A','Q'),('Q','S'),('S','A')]: arr(a,b,FWD,-0.25)
ax.text(0,-0.05,'forward\n$c>0$',ha='center',va='center',fontsize=7,color=FWD)
ax.text(0,-1.62,r'$c = P(Q\,|\,{\rm onset}) - P(Q\,|\,{\rm waking})$',ha='center',fontsize=7.5)
ax.text(-1.25,1.45,'(a)',fontsize=9,fontweight='bold')
# (b) c, (c) J
rng=np.random.default_rng(1)
for row,(key,ylab,tag) in enumerate([('c',r'cycle bias $c$','(b)'),('J',r'net current $J$ (cycles h$^{-1}$)','(c)')]):
    ax=fig.add_subplot(gs[row,1])
    for x,(g,bn) in zip(pos,blocks):
        r=R[f'0.5|{g}|{bn}']; v=np.array(r[key]); ci=r[f'{key}_ci']; m=v.mean(); col=FWD if ci[0]>0 else (REV if ci[1]<0 else INK)
        ax.scatter(x+rng.uniform(-.2,.2,len(v)),v,s=4,color=GREY,lw=0,alpha=.6)
        ax.errorbar(x,m,yerr=[[m-ci[0]],[ci[1]-m]],fmt='o',ms=4.5,color=col,elinewidth=1.1,capsize=2)
        for th,dx,mk in [(0.33,-0.33,'v'),(0.75,0.33,'^')]:
            ax.plot(x+dx,np.mean(R[f'{th}|{g}|{bn}'][key]),mk,ms=3,mfc='white',mec=GREY,mew=.6)
    ax.axhline(0,color=INK,lw=.5); ax.set_ylabel(ylab); ax.set_xticks(pos)
    ax.set_xticklabels([b if g!='per0DD' else 'all' for g,b in blocks],fontsize=7)
    lim=(-0.55,0.75) if key=='c' else (-0.25,0.35); ax.set_ylim(*lim)
    ax.text(-0.09,1.04,tag,fontsize=9,fontweight='bold',va='bottom',transform=ax.transAxes)
    if row==1:
        for x0,g in [(0.5,'CS'),(3.1,'fumin'),(5.7,'per0LD'),(7.8,'per0DD'),(9.9,'sss')]: ax.text(x0,-0.26,LAB[g],ha='center',fontsize=8,transform=ax.get_xaxis_transform())
    if row==0:
        ax.plot([],[],'o',color=FWD,ms=4,label='forward (95% CI > 0)'); ax.plot([],[],'o',color=REV,ms=4,label='reverse (95% CI < 0)'); ax.plot([],[],'o',color=INK,ms=4,label='CI includes 0')
        ax.plot([],[],'v',mfc='white',mec=GREY,ms=3,label=r'$\theta=0.33$'); ax.plot([],[],'^',mfc='white',mec=GREY,ms=3,label=r'$\theta=0.75$')
        ax.legend(fontsize=6.5,ncol=5,loc='upper center',bbox_to_anchor=(0.5,1.2),handletextpad=.2,columnspacing=.8)
fig.subplots_adjust(left=.02,right=.99,bottom=.13,top=.92); fig.savefig('figures/fig5_cycle_current.pdf'); fig.savefig('figures/fig5_cycle_current.png'); print('ok')
