# Onsager conductance per genotype (post hoc, user request): EPR = J^2 / L, OLS through the origin, as in the user's
# Methods (Eq. 25): slope s = 1/L, SE by OLS (df = n - 1), 95% CI of s by t, inverted for L (bounds swap); SE_L = L^2 SE_s.
# Also a fly bootstrap CI (resampling flies with both phases), and the variant with the calibrated estimator offset
# (0.02 nats/h, PRESPEC 20f) subtracted. Wald z-test of s against CS, Benjamini-Hochberg across the 4 contrasts.
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import t as tdist, norm
P=json.load(open('sqw_out/j_points.json')); G=['CS','per0LD','per0DD','fumin','sss']; rng=np.random.default_rng(20)
def ols0(J,E):
    x=np.asarray(J)**2; y=np.asarray(E); s=(x@y)/(x@x); n=len(x); res=y-s*x
    se=np.sqrt((res@res)/(n-1)/(x@x)); r2=1-(res@res)/(y@y); return s,se,n,r2
OUT={}; L=[]
for ph_sel,lab in ((None,'both phases'),(0,'light'),(1,'dark')):
    for g in G:
        if g=='per0DD' and ph_sel==0: continue
        rows=[r for r in P[g] if ph_sel is None or r[1]==ph_sel or (g=='per0DD' and ph_sel==1)]
        J=np.array([r[2] for r in rows]); E=np.array([r[3] for r in rows]); flies=sorted({r[0] for r in rows})
        s,se,n,r2=ols0(J,E); tq=tdist.ppf(0.975,n-1); lo,hi=s-tq*se,s+tq*se
        Lv=1/s; seL=Lv**2*se; Lci=(1/hi,1/lo if lo>0 else np.inf)
        bs=[]
        for _ in range(2000):
            pick=rng.choice(flies,len(flies)); rr=[r for f in pick for r in rows if r[0]==f]
            bs.append(1/ols0([r[2] for r in rr],[r[3] for r in rr])[0])
        s2=ols0(J,np.maximum(E-0.02,0))[0]
        OUT[(g,lab)]=dict(L=Lv,seL=seL,ci=Lci,boot=np.percentile(bs,[2.5,97.5]),s=s,se=se,n=n,r2=r2,L_offset=1/s2)
        L.append(f"{lab:11s} {g:7s} n {n:3d} | 1/L = {s:.2f} +- {se:.2f} | L = {Lv:.3f} +- {seL:.3f} (95% CI {Lci[0]:.3f}-{Lci[1]:.3f}; fly bootstrap {np.percentile(bs,2.5):.3f}-{np.percentile(bs,97.5):.3f}) | R2(origin) {r2:.2f} | L with offset subtracted {1/s2:.3f}")
    # Wald vs CS on 1/L, BH
    ps=[]; comp=[g for g in G if g!='CS' and (g,lab) in OUT]
    for g in comp:
        a,b=OUT[(g,lab)],OUT[('CS',lab)] if ('CS',lab) in OUT else None
        if b is None: continue
        z=(a['s']-b['s'])/np.hypot(a['se'],b['se']); ps.append((g,z,2*norm.sf(abs(z))))
    if ps:
        p=np.array([x[2] for x in ps]); o=np.argsort(p); m=len(p); adj=np.empty(m); run=1
        for k in range(m-1,-1,-1): run=min(run,p[o[k]]*m/(k+1)); adj[o[k]]=run
        for (g,z,pp),a in zip(ps,adj): L.append(f"   {lab}: {g} vs CS Wald z {z:+.2f} p {pp:.2g} BH {a:.2g}")
open('onsager_L.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
# ---- figure ----
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
COL={'CS':'#2a78d6','fumin':'#eb6834','per0LD':'#1baf7a','per0DD':'#eda100','sss':'#e87ba4'}
NAME={'CS':'CS','fumin':'$fmn$','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','sss':'$sss$'}
fig,axs=plt.subplots(1,2,figsize=(7.0,2.8),gridspec_kw=dict(width_ratios=[1,1.4]))
ax=axs[0]
for k,g in enumerate(G):
    o=OUT[(g,'both phases')]
    ax.errorbar(k,o['L'],yerr=[[o['L']-o['ci'][0]],[o['ci'][1]-o['L']]],fmt='o',color=COL[g],ms=5,elinewidth=1.1,capsize=2.5)
    ax.plot([k+0.18]*2,o['boot'],color=COL[g],lw=0.7,ls=(0,(1.5,1)))
ax.set_xticks(range(len(G))); ax.set_xticklabels([NAME[g] for g in G]); ax.set_xlim(-0.5,len(G)-0.5)
ax.set_ylabel('Onsager conductance $L$ (cycles$^2$ nats$^{-1}$ h$^{-1}$)'); ax.set_title('(a) light and dark pooled',fontsize=8,loc='left')
ax.set_ylim(0,None)
ax=axs[1]
for k,g in enumerate(G):
    for lab,dx,filled in (('light',-0.12,False),('dark',0.12,True)):
        if (g,lab) not in OUT: continue
        o=OUT[(g,lab)]; x=k+(0 if g=='per0DD' else dx)
        ax.errorbar(x,o['L'],yerr=[[o['L']-o['ci'][0]],[min(o['ci'][1],10)-o['L']]],fmt='o',color=COL[g],mfc=COL[g] if filled else 'white',mec=COL[g],ms=5,elinewidth=1.1,capsize=2.5)
ax.set_xticks(range(len(G))); ax.set_xticklabels([NAME[g] for g in G]); ax.set_xlim(-0.5,len(G)-0.5); ax.set_ylim(0,None)
ax.plot([],[],'o',mfc='white',mec='#1a1a1a',ms=5,label='light'); ax.plot([],[],'o',mfc='#1a1a1a',mec='#1a1a1a',ms=5,label='dark / DD'); ax.legend(frameon=False,fontsize=6.5)
ax.set_title('(b) by phase',fontsize=8,loc='left')
fig.text(0.01,0.005,'EPR = $J^2/L$ fitted through the origin per genotype (per-fly points, hidden states). Bars: 95% CI from the fit; dotted: fly bootstrap.',fontsize=6,color='#333')
fig.tight_layout(rect=(0,0.04,1,1),pad=0.4); fig.savefig('figures/fig12_onsager_L.pdf'); fig.savefig('figures/fig12_onsager_L.png',dpi=300)
# fly-bootstrap comparison of L against CS (pooled phases): difference in L, 95% CI and two-sided bootstrap p
def bootL(g,B=4000,seed=21):
    r=np.random.default_rng(seed); rows=P[g]; flies=sorted({x[0] for x in rows}); out=[]
    for _ in range(B):
        pick=r.choice(flies,len(flies)); rr=[x for f in pick for x in rows if x[0]==f]; out.append(1/ols0([x[2] for x in rr],[x[3] for x in rr])[0])
    return np.array(out)
bc=bootL('CS',seed=30); L2=['fly-bootstrap comparison with CS (pooled phases), L difference and 95% CI, two-sided bootstrap p:']; pr=[]
for k,g in enumerate(['per0LD','per0DD','fumin','sss']):
    d=bootL(g,seed=31+k)-bc; p=min(1,2*min(np.mean(d<=0),np.mean(d>=0))); pr.append(p)
    L2.append(f"   {g:7s} dL {OUT[(g,'both phases')]['L']-OUT[('CS','both phases')]['L']:+.3f} [{np.percentile(d,2.5):+.3f},{np.percentile(d,97.5):+.3f}] p {p:.3g}")
d=bootL('per0DD',seed=40)-bootL('per0LD',seed=41); L2.append(f"   per0DD - per0LD dL [{np.percentile(d,2.5):+.3f},{np.percentile(d,97.5):+.3f}] p {min(1,2*min(np.mean(d<=0),np.mean(d>=0))):.3g}")
open('onsager_L.log','a').write('\n'.join(L2)+'\n'); print('\n'.join(L2))
