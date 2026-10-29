# PRESPEC 20l: (1) state energies, (2) landscape vs counted routing (Viterbi, FFBS), (3) Fig 13 hidden with landscape lines.
import json,numpy as np,sys,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker
from scipy.stats import spearmanr
sys.path.insert(0,'.'); import sqw_ref as R
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
COL={'CS':'#2a78d6','per0LD':'#1baf7a','per0DD':'#eda100','fumin':'#eb6834','sss':'#e87ba4'}
NAME={'CS':'CS','fumin':'$fmn$','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','sss':'$sss$'}; G=list(COL)
F={(d['g'],d['fly']):d for d in map(json.loads,open('sqw_out/fit2_0.5.jsonl'))}; RT=json.load(open('sqw_out/route2_0.5.json'))
rng=np.random.default_rng(20); L=[]
rows=[]   # per fly x phase
for (g,i),d in F.items():
    if d['deg']: continue
    P=np.array(d['P3']); rt=RT[f'{g}|{i}']
    for p in range(P.shape[0]):
        bn='whole' if g=='per0DD' else ['light','dark'][p]; ref=d['ref'].get(bn)
        if ref is None or not np.isfinite(ref) or abs(R.model_energetics(P[p])['J'])<1e-12: continue
        la=R.landscape(P[p]); v=np.array(rt['vit'][p]); Xv,Yv=R.xy_counts(v)
        rp=np.array(rt['routes_path'])[:,p]; X,Y=R.xy_counts(rp)
        Xf=np.nanmean(X) if np.isfinite(X).mean()>=0.5 else np.nan; Yf=np.nanmean(Y) if np.isfinite(Y).mean()>=0.5 else np.nan
        rows.append(dict(g=g,fly=i,phase=p,bn=bn,**{k:float(la[k]) for k in ('E2_minus_E1','E3_minus_E1','F_SW','F_WS')},Xv=float(Xv),Yv=float(Yv),Xf=float(Xf),Yf=float(Yf)))
json.dump(rows,open('sqw_out/l_landscape.json','w'))
def ci_med(v):
    v=np.asarray(v); return np.median(v),*np.percentile([np.median(rng.choice(v,len(v))) for _ in range(2000)],[2.5,97.5])
# ---- (1) energies ----
L.append('(1) state energies relative to H1 (nats), median [fly-bootstrap 95% CI]')
fig,axs=plt.subplots(1,2,figsize=(7.0,2.9))
for ax,key,lab in ((axs[0],'E2_minus_E1','$E_2-E_1$ (wake H2)'),(axs[1],'E3_minus_E1','$E_3-E_1$ (sleep H3)')):
    for k,g in enumerate(G):
        for p,dx,filled in ((0,-0.17,False),(1,0.17,True)):
            if g=='per0DD' and p==1: continue
            v=np.array([r[key] for r in rows if r['g']==g and r['phase']==p])
            x=k+(0 if g=='per0DD' else dx)
            ax.scatter(x+rng.uniform(-0.07,0.07,len(v)),v,s=10,facecolor=COL[g] if (filled or g=='per0DD') else 'white',edgecolor=COL[g],lw=0.7,zorder=2)
            m,lo,hi=ci_med(v); ax.plot([x-0.11,x+0.11],[m,m],color='#1a1a1a',lw=1.6,zorder=3)
            L.append(f"  {key} {g:7s} {'whole' if g=='per0DD' else ['light','dark'][p]:5s} n {len(v):2d} median {m:+.2f} [{lo:+.2f},{hi:+.2f}]")
    ax.axhline(0,color='#888',lw=0.5); ax.set_xticks(range(len(G))); ax.set_xticklabels([NAME[g] for g in G]); ax.set_ylabel(lab+' (nats)')
axs[0].scatter([],[],s=10,facecolor='white',edgecolor='#1a1a1a',label='light'); axs[0].scatter([],[],s=10,facecolor='#1a1a1a',edgecolor='#1a1a1a',label='dark / DD')
axs[0].legend(frameon=False,fontsize=6.5,loc='best')
fig.text(0.01,0.005,'Per fly and phase, $E_i=-\\ln\\pi_i$ of the phase matrix, relative to the regulator H1 ($E_1=0$); bar: median. Negative: state more occupied than H1.',fontsize=6,color='#333')
fig.tight_layout(rect=(0,0.04,1,1),pad=0.4); fig.savefig('figures/fig14_state_energies.pdf'); fig.savefig('figures/fig14_state_energies.png',dpi=300); plt.close(fig)
# ---- (2) landscape vs counted ----
for tag,xk,yk in (('viterbi','Xv','Yv'),('ffbs','Xf','Yf')):
    L.append(f'(2) landscape prediction vs counted routing ({tag}): median(counted - predicted), Spearman')
    fig,axs=plt.subplots(2,2,figsize=(6.0,6.0),sharex=True,sharey=True)
    for c,pp in enumerate((0,1)):
        for r_,(pk,ck,lab) in enumerate((('F_SW',xk,'$F_{S\\to W}$ (sleep$\\to$wake)'),('F_WS',yk,'$F_{W\\to S}$ (wake$\\to$sleep)'))):
            ax=axs[r_,c]; ax.plot([-1,1],[-1,1],color='#888',lw=0.6,ls=(0,(3,2)))
            for g in G:
                sel=[r for r in rows if r['g']==g and (r['phase']==pp or (g=='per0DD' and pp==1)) and not (g=='per0DD' and pp==0)]
                x=np.array([r[pk] for r in sel]); y=np.array([r[ck] for r in sel]); m=np.isfinite(x)&np.isfinite(y)
                ax.scatter(x[m],y[m],s=9,facecolor=COL[g] if pp==1 else 'white',edgecolor=COL[g],lw=0.7)
                if m.sum()>4:
                    L.append(f"  {tag} {lab[:14]:14s} {g:7s} {['light','dark/DD'][pp]:7s} n {m.sum():2d} median diff {np.median(y[m]-x[m]):+.3f} Spearman {spearmanr(x[m],y[m]).statistic:+.2f}")
            ax.set_xlim(-1.05,1.05); ax.set_ylim(-1.05,1.05); ax.set_aspect('equal')
            if r_==0: ax.set_title(['light (ZT 0-12)','dark (ZT 12-24) / DD'][c],fontsize=8)
            if c==0: ax.set_ylabel(f"counted ({'Viterbi' if tag=='viterbi' else 'posterior'})\n"+lab)
            if r_==1: ax.set_xlabel('landscape prediction tanh[($E_{direct}-E^{c}_{indirect}$)/2]')
    for g in G: axs[0,1].scatter([],[],s=9,color=COL[g],label=NAME[g])
    axs[0,1].legend(frameon=False,fontsize=6,loc='lower right')
    fig.tight_layout(pad=0.4); fig.savefig(f'figures/fig15_landscape_vs_counted_{tag}.pdf'); fig.savefig(f'figures/fig15_landscape_vs_counted_{tag}.png',dpi=300); plt.close(fig)
# ---- (3) Fig 13 hidden panels with landscape lines ----
res=json.load(open('sqw_out/k_route_duration.json'))
BS=[5,10,20,40,80,1e9]; BW=[10,20,40,80,160,1e9]
def wilson(k,n):
    p=k/n; z=1.96; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return p,c-h,c+h
fig,axs=plt.subplots(2,3,figsize=(7.0,4.6),sharey=True)
for col,rel in enumerate(R.RELATIONS):
    lab,dur=R.RELATIONS[rel]; edges=BW if dur=='wake_before_min' else BS; pk='F_SW' if lab=='q_out' else 'F_WS'
    for row,pp in enumerate((0,1)):
        ax=axs[row,col]
        for g in G:
            if (g=='per0DD' and pp==0): continue
            recs=[x for r in res if r['g']==g for path in r['hidrec'] for x in path if x['phase_on']==(0 if g=='per0DD' else pp) and x[lab] is not None]
            xs=[];ys=[];lo=[];hi=[]
            for a,b in zip(edges[:-1],edges[1:]):
                sel=[x for x in recs if a<=x[dur]<b]; n=len(sel)/20; k=sum(1 for x in sel if x[lab])/20
                if n<10: continue
                p_,l_,h_=wilson(k,n); xs.append(np.sqrt(a*min(b,a*2))); ys.append(p_); lo.append(l_); hi.append(h_)
            ax.errorbar(xs,ys,yerr=[np.array(ys)-np.array(lo),np.array(hi)-np.array(ys)],fmt='o-',ms=3,lw=0.9,elinewidth=0.6,color=COL[g],mfc=COL[g] if pp==1 else 'white',mec=COL[g])
            pred=np.mean([(r[pk]+1)/2 for r in rows if r['g']==g and r['phase']==(0 if g=='per0DD' else pp)])
            ax.axhline(pred,color=COL[g],lw=0.8,ls=(0,(4,2)),alpha=0.9)
        ax.set_xscale('log'); ax.set_ylim(0,1); tk=[10,20,40,80,160] if dur=='wake_before_min' else [5,10,20,40,80]
        ax.set_xticks(tk); ax.set_xticklabels([str(t) for t in tk]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        if row==0: ax.set_title({'R1_exit_vs_sleep':'R1 quiet exit vs sleep bout','R2_entry_vs_wake_before':'R2 quiet entry vs prior wake','R3_entry_vs_sleep':'R3 quiet entry vs next sleep bout'}[rel],fontsize=8)
        if col==0: ax.set_ylabel(['light','dark / DD'][row]+'\nP(quiet route), hidden')
        if row==1: ax.set_xlabel(('wake' if dur=='wake_before_min' else 'sleep')+' duration (min)')
for g in G: axs[0,2].plot([],[],'o-',color=COL[g],ms=3,label=NAME[g])
axs[0,2].plot([],[],color='#1a1a1a',ls=(0,(4,2)),lw=0.8,label='landscape prediction'); axs[0,2].legend(frameon=False,fontsize=6,loc='upper right')
fig.text(0.01,0.003,'Hidden level (posterior paths, matched definition). Dashed: landscape prediction (F+1)/2 from the phase matrix, mean over flies; it is constant in duration by construction (Markov).\nA slope of the counted curve against its flat line is memory. Bouts pooled over flies (descriptive; tests in Fig. 13b are within-fly).',fontsize=5.8,color='#333')
fig.tight_layout(rect=(0,0.045,1,1),pad=0.4); fig.savefig('figures/fig13c_route_duration_landscape.pdf'); fig.savefig('figures/fig13c_route_duration_landscape.png',dpi=300)
open('sqw11.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
