# EPR vs J_cycle, genotype x phase medians (S/Q/W model, theta 0.5, 20b fits). Descriptive figure (post hoc request).
# J: count-based, forward-positive (H3->H2->H1->H3), cycles/h. EPR: rule-of-3 95% lower bound, per posterior path,
# averaged over paths (Jensen). Error bars: 95% bootstrap CI of the median over flies.
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sqw_ref as R
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,
                     'xtick.direction':'in','ytick.direction':'in','xtick.major.width':0.6,'ytick.major.width':0.6})
INK='#1a1a1a'; GREY='#9a9a9a'
COL={'CS':'#2a78d6','fumin':'#eb6834','per0LD':'#1baf7a','per0DD':'#eda100','sss':'#e87ba4'}
NAME={'CS':'CS','fumin':'$fmn$','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','sss':'$sss$'}
F={(d['g'],d['fly']):d for d in map(json.loads,open('sqw_out/fit2_0.5.jsonl'))}; RT=json.load(open('sqw_out/route2_0.5.json'))
rng=np.random.default_rng(20); rows=[]
# blocks where the count-based EPR overestimates the exact model EPR by > 50% on simulated flies (PRESPEC 20f)
WEAK={('CS','light'),('per0LD','light'),('sss','light'),('sss','dark')}
for g in COL:
    for p,bn in enumerate(['whole'] if g=='per0DD' else ['light','dark']):
        J=[];E=[]
        for k,d in F.items():
            if k[0]!=g or d['deg']: continue
            rt=RT[f'{g}|{k[1]}']; h=rt['hours'][p]; pp=np.array(rt['pairs_path'])[:,p]
            J.append(np.mean([R.j_cycle_count(n,h) for n in pp])); E.append(np.mean([R.epr_count_bound(n,h)[0] for n in pp]))
        J=np.array(J);E=np.array(E)
        bJ=[np.median(rng.choice(J,len(J))) for _ in range(2000)]; bE=[np.median(rng.choice(E,len(E))) for _ in range(2000)]
        Em=np.array([R.model_energetics(np.array(d['P3'])[p])['EPR'] for k,d in F.items() if k[0]==g and not d['deg']])
        rows.append(dict(Em=np.median(Em),weak=(g,bn) in WEAK,g=g,bn=bn,n=len(J),J=np.median(J),Jci=np.percentile(bJ,[2.5,97.5]),E=np.median(E),Eci=np.percentile(bE,[2.5,97.5]),Jf=J,Ef=E))
fig,ax=plt.subplots(figsize=(3.4,2.9))
ax.axhspan(0.255,0.327,color=GREY,alpha=0.18,lw=0)
ax.text(-0.075,0.29,'reversible-null 95th pct.\n(single fly)',fontsize=6,color='#555555',va='center')
jj=np.linspace(0,0.3,50)
for A in (0.5,1,2):
    ax.plot(jj,A*jj,color=GREY,lw=0.5,ls=(0,(3,2))); ax.plot(-jj,A*jj,color=GREY,lw=0.5,ls=(0,(3,2)))
ax.axvline(0,color=INK,lw=0.5)
for A,(xa,ya) in {2:(0.192,0.405),1:(0.272,0.262),0.5:(0.268,0.118)}.items(): ax.text(xa,ya,f'$A$ = {A} nats/cycle',fontsize=5.5,color='#555555',ha='right')
lab_xy={('CS','light'):(0.005,0.095),('CS','dark'):(0.212,0.312),('fumin','light'):(0.104,0.085),('fumin','dark'):(0.116,0.150),
        ('per0LD','light'):(0.075,0.012),('per0LD','dark'):(0.205,0.160),('per0DD','whole'):(0.240,0.215),('sss','light'):(-0.082,0.010),('sss','dark'):(-0.082,0.078)}
lead={('CS','light'),('per0LD','light'),('sss','light'),('sss','dark')}
for r in rows:
    c=COL[r['g']]; filled=r['bn']!='light'
    al=0.35 if r['weak'] else 1.0
    ax.errorbar(r['J'],r['E'],xerr=[[r['J']-r['Jci'][0]],[r['Jci'][1]-r['J']]],yerr=[[r['E']-r['Eci'][0]],[r['Eci'][1]-r['E']]],
                fmt='none',ecolor=c,elinewidth=0.9,capsize=0,alpha=al)
    ax.plot(r['J'],r['E'],'o',ms=5.5,mfc=c if filled else 'white',mec=c,mew=1.2,zorder=3,alpha=al)
    lab=NAME[r['g']]+('' if r['bn']=='whole' else (' L' if r['bn']=='light' else ' D'))
    key=(r['g'],r['bn'])
    ax.annotate(lab,(r['J'],r['E']),xytext=lab_xy[key],textcoords='data',fontsize=6.5,color=INK,
                arrowprops=dict(arrowstyle='-',lw=0.4,color='#555555',shrinkA=1,shrinkB=4) if key in lead else None)
ax.plot([],[],'o',ms=5,mfc='white',mec=INK,mew=1,label='light (ZT 0–12)'); ax.plot([],[],'o',ms=5,mfc=INK,mec=INK,label='dark (ZT 12–24) / DD')
ax.legend(frameon=False,fontsize=6.5,loc='upper left',handletextpad=0.3,borderaxespad=0.2)
ax.axhspan(0,0.025,color='#d9d9d9',alpha=0.6,lw=0,zorder=0)
ax.text(0.03,0.585,'bottom band: hidden-estimator offset\n(+0.01–0.025 nats h$^{-1}$; PRESPEC 20f)\nfaded: weak circulation, EPR ~ offset',transform=ax.transAxes,ha='left',va='top',fontsize=5.5,color='#555555')
ax.set_xlabel(r'$J_{\mathrm{cycle}}$ (cycles h$^{-1}$; $+$ = S$\to$W$\to$Reg$\to$S)'); ax.set_ylabel(r'EPR lower bound (nats h$^{-1}$)')
ax.set_xlim(-0.09,0.3); ax.set_ylim(0,0.42)
fig.tight_layout(pad=0.3); fig.savefig('figures/fig6_epr_vs_j.pdf'); fig.savefig('figures/fig6_epr_vs_j.png',dpi=300)
with open('fig6_epr_vs_j.log','w') as fh:
    for r in rows: fh.write(f"{r['g']:7s} {r['bn']:5s} n {r['n']:2d} | model-exact EPR median {r['Em']:.3f} | weak {r['weak']} | J median {r['J']:+.3f} [{r['Jci'][0]:+.3f},{r['Jci'][1]:+.3f}] | EPR median {r['E']:.3f} [{r['Eci'][0]:.3f},{r['Eci'][1]:.3f}] | EPR/|J| {r['E']/max(abs(r['J']),1e-9):.2f}\n")
print(open('fig6_epr_vs_j.log').read())
