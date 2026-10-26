# PRESPEC 20k statistics and figure.
import json,numpy as np,sys,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker
sys.path.insert(0,'.'); import sqw_ref as R
res=json.load(open('sqw_out/k_route_duration.json')); REL=list(R.RELATIONS); rng=np.random.default_rng(20)
G=['CS','fumin','per0LD','per0DD','sss']; L=[]; OUT={}
def holm(p):
    p=np.asarray(p,float); o=np.argsort(p); m=len(p); adj=np.empty(m); run=0
    for k,j in enumerate(o): run=max(run,(m-k)*p[j]); adj[j]=min(1,run)
    return adj
for rel in REL:
    for level in ('observed','hidden'):
        rows=[]
        for g in G:
            for p in range(1 if g=='per0DD' else 2):
                bn='whole' if g=='per0DD' else ['light','dark'][p]; Rg=[r for r in res if r['g']==g]
                if level=='observed':
                    recs=[[x for x in r['obs'] if x['phase_on']==p] for r in Rg]
                    D=np.array([R.route_duration_D(rc,rel,p) for rc in recs]); ok=np.isfinite(D); real=D[ok].mean()
                    lab,dur=R.RELATIONS[rel]; null=[]
                    for _ in range(2000):
                        tot=[]
                        for rc,o in zip(recs,ok):
                            if not o: continue
                            rr=[x for x in rc if x[lab] is not None]; labs=rng.permutation([x[lab] for x in rr])
                            q=[np.log(x[dur]) for x,l in zip(rr,labs) if l]; a=[np.log(x[dur]) for x,l in zip(rr,labs) if not l]
                            tot.append(np.mean(q)-np.mean(a))
                        null.append(np.mean(tot))
                    null=np.array(null); pv=min(1,2*min(np.mean(null>=real),np.mean(null<=real))); nullm=null.mean(); n=int(ok.sum()); nlo,nhi=np.percentile(null,[2.5,97.5])
                else:
                    D=np.array([r['hid'].get(f'{rel}|{p}') if r['hid'].get(f'{rel}|{p}') is not None else np.nan for r in Rg],float)
                    ok=np.isfinite(D); real=D[ok].mean(); n=int(ok.sum())
                    S=np.array([[s.get(f'{rel}|{p}') if s.get(f'{rel}|{p}') is not None else np.nan for s in r['sims']] for r in Rg],float)  # fly x sim
                    simmeans=np.nanmean(S,0); nullm=simmeans.mean()
                    pv=min(1,2*min((1+np.sum(simmeans>=real))/(1+len(simmeans)),(1+np.sum(simmeans<=real))/(1+len(simmeans)))); nlo,nhi=np.percentile(simmeans,[2.5,97.5])
                rows.append([g,bn,n,real,nullm,pv,nlo,nhi])
        adj=holm([r[5] for r in rows])
        for r,a in zip(rows,adj):
            OUT[f'{rel}|{level}|{r[0]}|{r[1]}']=dict(n=r[2],D=r[3],null_mean=r[4],p=r[5],p_holm=float(a),null_lo=float(r[6]),null_hi=float(r[7]))
            L.append(f"{rel:24s} {level:8s} {r[0]:7s} {r[1]:5s} flies {r[2]:2d} | D {r[3]:+.3f} | null mean {r[4]:+.3f} | p {r[5]:.3g} Holm {a:.3g}{'  *' if a<0.05 else ''}")
json.dump(OUT,open('sqw_out/k_summary.json','w')); open('sqw10.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
# ---- figure: P(quiet route) vs duration bin, observed and hidden ----
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,'xtick.direction':'in','ytick.direction':'in'})
COL={'CS':'#2a78d6','fumin':'#eb6834','per0LD':'#1baf7a','per0DD':'#eda100','sss':'#e87ba4'}
NAME={'CS':'CS','fumin':'$fmn$','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','sss':'$sss$'}
BS=[5,10,20,40,80,1e9]; BW=[10,20,40,80,160,1e9]
def wilson(k,n):
    if n==0: return np.nan,np.nan,np.nan
    p=k/n; z=1.96; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return p,c-h,c+h
fig,axs=plt.subplots(2,3,figsize=(7.0,4.6),sharey=True)
for col,rel in enumerate(REL):
    lab,dur=R.RELATIONS[rel]; edges=BW if dur=='wake_before_min' else BS
    for row,level in enumerate(('observed','hidden')):
        ax=axs[row,col]
        for g in G:
            for p in range(1 if g=='per0DD' else 2):
                if level=='observed': recs=[x for r in res if r['g']==g for x in r['obs'] if x['phase_on']==p and x[lab] is not None]; w=None
                else:
                    recs=[x for r in res if r['g']==g for path in r['hidrec'] for x in path if x['phase_on']==p and x[lab] is not None]
                xs=[];ys=[];lo=[];hi=[]
                for a,b in zip(edges[:-1],edges[1:]):
                    sel=[x for x in recs if a<=x[dur]<b]; k=sum(1 for x in sel if x[lab]); n=len(sel)
                    if level=='hidden': n=n/20; k=k/20      # 20 paths kept per fly: average counts per path
                    if n<10: continue
                    pp,l,h=wilson(k,n); xs.append(np.sqrt(a*min(b,a*2))); ys.append(pp); lo.append(l); hi.append(h)
                filled=g=='per0DD' or p==1
                ax.errorbar(xs,ys,yerr=[np.array(ys)-np.array(lo),np.array(hi)-np.array(ys)],fmt='o-',ms=3,lw=0.8,elinewidth=0.6,color=COL[g],mfc=COL[g] if filled else 'white',mec=COL[g])
        ax.set_xscale('log'); ax.set_ylim(0,1)
        tk=[10,20,40,80,160] if dur=='wake_before_min' else [5,10,20,40,80]
        ax.set_xticks(tk); ax.set_xticklabels([str(t) for t in tk]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        if row==0: ax.set_title({'R1_exit_vs_sleep':'R1 quiet exit vs sleep bout','R2_entry_vs_wake_before':'R2 quiet entry vs prior wake','R3_entry_vs_sleep':'R3 quiet entry vs next sleep bout'}[rel],fontsize=8)
        if col==0: ax.set_ylabel(f"{level}\nP(quiet route)")
        if row==1: ax.set_xlabel(('wake' if dur=='wake_before_min' else 'sleep')+' duration (min)')
for g in G: axs[0,2].plot([],[],'o-',color=COL[g],ms=3,label=NAME[g])
axs[0,2].legend(frameon=False,fontsize=6,loc='upper right')
fig.text(0.01,0.003,'Open: light, filled: dark / DD. Bouts pooled over flies (descriptive; the tests in Fig. 13b are within-fly). Bins with < 10 bouts omitted; Wilson 95% intervals\n(hidden: per-path average counts, 20 paths per fly).',fontsize=6,color='#333')
fig.tight_layout(rect=(0,0.04,1,1),pad=0.4); fig.savefig('figures/fig13_route_vs_duration.pdf'); fig.savefig('figures/fig13_route_vs_duration.png',dpi=300)

# ---- Fig 13b: effect D vs null interval, per relation, level, genotype x phase ----
BL=[('CS','light'),('CS','dark'),('per0LD','light'),('per0LD','dark'),('per0DD','whole'),('fumin','light'),('fumin','dark'),('sss','light'),('sss','dark')]
TIT={'R1_exit_vs_sleep':'R1 quiet exit vs sleep-bout length','R2_entry_vs_wake_before':'R2 quiet entry vs prior wake length','R3_entry_vs_sleep':'R3 quiet entry vs next sleep length'}
fig,axs=plt.subplots(2,3,figsize=(7.0,4.8),sharey=True,sharex=True)
yy=np.arange(len(BL))[::-1]
for col,rel in enumerate(REL):
    for row,level in enumerate(('observed','hidden')):
        ax=axs[row,col]; ax.axvline(0,color='#888',lw=0.5)
        for y,(g,bn) in zip(yy,BL):
            o=OUT[f'{rel}|{level}|{g}|{bn}']
            ax.plot([o['null_lo'],o['null_hi']],[y,y],color='#cfcfcf',lw=4,solid_capstyle='butt',zorder=1)
            filled=bn!='light'
            ax.plot(o['D'],y,'o',ms=4.5,mfc=COL[g] if filled else 'white',mec=COL[g],mew=1.1,zorder=3)
            if o['p_holm']<0.05: ax.text(o['D'],y+0.28,'*',ha='center',va='center',fontsize=9,color='#1a1a1a')
        if row==0: ax.set_title(TIT[rel],fontsize=7.5)
        if col==0: ax.set_ylabel(level)
        if row==1: ax.set_xlabel('D = mean ln dur(quiet) - mean ln dur(active)',fontsize=7)
axs[0,0].set_yticks(yy); axs[0,0].set_yticklabels([NAME[g]+('' if bn=='whole' else (' L' if bn=='light' else ' D')) for g,bn in BL])
for ax in axs[1]: ax.tick_params(labelsize=7)
fig.text(0.01,0.003,'Points: mean D over flies. Grey: 95% range without dependence (observed: within-fly label permutation; hidden: fitted Markov model, simulated and decoded).\n*: Holm p < 0.05. D < 0: quiet routes go with shorter bouts.',fontsize=5.8,color='#333')
fig.tight_layout(rect=(0,0.045,1,1),pad=0.4); fig.savefig('figures/fig13b_route_duration_effects.pdf'); fig.savefig('figures/fig13b_route_duration_effects.png',dpi=300)
