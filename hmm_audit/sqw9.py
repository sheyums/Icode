# PRESPEC 20j: per-fly EPR vs J (hidden), linear vs quadratic, slopes by genotype; Fig. 11.
import json,numpy as np,sys,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'.'); import sqw_ref as R, sqw
F={(d['g'],d['fly']):d for d in map(json.loads,open('sqw_out/fit2_0.5.jsonl'))}; RT=json.load(open('sqw_out/route2_0.5.json'))
pts={}   # g -> list of (fly, phase, J, EPR)
for (g,i),d in F.items():
    if d['deg']: continue
    rt=RT[f'{g}|{i}']; P=np.array(d['P3'])
    for p in range(len(rt['hours'])):
        bn='whole' if g=='per0DD' else ['light','dark'][p]
        if abs(R.model_energetics(P[p])['J'])<1e-12: continue
        r=d['ref'].get(bn); 
        if r is None or not np.isfinite(r): continue
        pp=np.array(rt['pairs_path'])[:,p]; h=rt['hours'][p]
        pts.setdefault(g,[]).append((i,p,np.mean([R.j_cycle_count(n,h) for n in pp]),np.mean([R.epr_count_bound(n,h)[0] for n in pp])))
def fit(x,y):   # least squares y = b + a x; returns a, b, rss
    A=np.c_[np.ones_like(x),x]; coef,res,_,_=np.linalg.lstsq(A,y,rcond=None); return coef[1],coef[0],float(((A@coef-y)**2).sum())
def models(rows):
    J=np.array([r[2] for r in rows]); E=np.array([r[3] for r in rows])
    a,b1,r1=fit(np.abs(J),E); c,b2,r2=fit(J**2,E)
    m=(np.abs(J)>=0.05)&(E>=0.05); n=np.polyfit(np.log(np.abs(J[m])),np.log(E[m]-0.02),1)[0] if m.sum()>=5 else np.nan
    return dict(a=a,b1=b1,rss1=r1,c=c,b2=b2,rss2=r2,n=n,n_pts=len(J),n_exp=int(m.sum()))
rng=np.random.default_rng(20); L=[]; RES={}
for g in ['CS','fumin','per0LD','per0DD','sss']:
    for ph in ([None,0,1] if g!='per0DD' else [None]):
        rows=[r for r in pts[g] if ph is None or r[1]==ph]; flies=sorted({r[0] for r in rows}); est=models(rows); B=[]
        for _ in range(2000):
            pick=rng.choice(flies,len(flies)); rr=[r for f in pick for r in rows if r[0]==f]; B.append(models(rr))
        ci=lambda k: np.nanpercentile([b[k] for b in B],[2.5,97.5])
        pref=np.mean([b['rss1']<b['rss2'] for b in B])
        RES[(g,ph)]=(est,B)
        lab='both phases' if ph is None else ('whole' if g=='per0DD' else ['light','dark'][ph])
        L.append(f"{g:7s} {lab:11s} n {est['n_pts']:3d} | M1 EPR = {est['b1']:+.3f} + {est['a']:.2f}|J|  a CI [{ci('a')[0]:.2f},{ci('a')[1]:.2f}] RSS {est['rss1']:.3f} | M2 EPR = {est['b2']:+.3f} + {est['c']:.2f}J^2  c CI [{ci('c')[0]:.2f},{ci('c')[1]:.2f}] RSS {est['rss2']:.3f}"
                 f" | better: {'linear (M1)' if est['rss1']<est['rss2'] else 'quadratic (M2)'} (bootstrap share M1 better {pref:.2f}) | exponent n {est['n']:.2f} [{ci('n')[0]:.2f},{ci('n')[1]:.2f}] (n pts {est['n_exp']})")
L.append('slope a (M1, both phases) minus CS, bootstrap 95% CI:')
for g in ['fumin','per0LD','per0DD','sss']:
    d=np.array([b['a'] for b in RES[(g,None)][1]])-np.array([b['a'] for b in RES[('CS',None)][1]])
    L.append(f"  {g:7s} {RES[(g,None)][0]['a']-RES[('CS',None)][0]['a']:+.2f} [{np.percentile(d,2.5):+.2f},{np.percentile(d,97.5):+.2f}]")
open('sqw9.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
json.dump({g:v for g,v in pts.items()},open('sqw_out/j_points.json','w'))
# ---- Fig 11 ----
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
COL={'CS':'#2a78d6','fumin':'#eb6834','per0LD':'#1baf7a','per0DD':'#eda100','sss':'#e87ba4'}
NAME={'CS':'CS','fumin':'$fmn$','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','sss':'$sss$'}
fig,axs=plt.subplots(2,3,figsize=(7.0,4.8)); axs=axs.ravel()
ax=axs[0]
for g in COL:
    for (i,p,J,E) in pts[g]:
        filled=not (g!='per0DD' and p==0)
        ax.plot(J,E,'o',ms=3.2,mfc=COL[g] if filled else 'white',mec=COL[g],mew=0.7,alpha=0.85)
ax.axvline(0,color='#555',lw=0.4); ax.set_title('(a) all flies',fontsize=8,loc='left'); ax.set_xlabel('$J_{\\mathrm{cycle}}$ (cycles h$^{-1}$)'); ax.set_ylabel('EPR (nats h$^{-1}$)')
for k,g in enumerate(COL):
    ax=axs[k+1]; est=RES[(g,None)][0]
    for (i,p,J,E) in pts[g]:
        filled=not (g!='per0DD' and p==0); ax.plot(J,E,'o',ms=3.2,mfc=COL[g] if filled else 'white',mec=COL[g],mew=0.7,alpha=0.9)
    jj=np.linspace(min(r[2] for r in pts[g]),max(r[2] for r in pts[g]),200)
    ax.plot(jj,est['b1']+est['a']*np.abs(jj),color='#1a1a1a',lw=0.9,label=f"$b+a|J|$, $a$={est['a']:.2f}")
    ax.plot(jj,est['b2']+est['c']*jj**2,color='#1a1a1a',lw=0.9,ls=(0,(3,2)),label=f"$b+cJ^2$, $c$={est['c']:.2f}")
    ax.axvline(0,color='#555',lw=0.4); ax.set_title(f"({'bcdef'[k]}) {NAME[g]}  ($n$ = {est['n_pts']}; exponent {est['n']:.1f})",fontsize=8,loc='left')
    ax.legend(frameon=False,fontsize=6,loc='upper left'); ax.set_xlabel('$J_{\\mathrm{cycle}}$ (cycles h$^{-1}$)')
    if k in (2,): ax.set_ylabel('EPR (nats h$^{-1}$)')
fig.text(0.01,0.003,'Per fly and phase (open: light, filled: dark / DD). EPR = rule-of-three bound from posterior hidden paths. Fits pool both phases.',fontsize=6,color='#333')
fig.tight_layout(rect=(0,0.02,1,1),pad=0.4); fig.savefig('figures/fig11_epr_vs_j_perfly.pdf'); fig.savefig('figures/fig11_epr_vs_j_perfly.png',dpi=300)
