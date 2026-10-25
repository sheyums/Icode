# Hidden vs observable EPR-J (post hoc request). Hidden: counts on posterior hidden paths (as fig_epr_j.py).
# Observable: the same counts on the run-compressed S/Q/W minute symbols, with Q<->H1, W<->H2, S<->H3
# (symbol codes 0,1,2 share the state indices), so J_obs = (N_SW - N_WS)/h and the same rule-of-3 EPR bound.
# Transitions into or out of a missing minute (-1) are dropped. Coarse-graining: EPR_obs is a lower bound.
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sqw_ref as R, sqw2
exec(open('fig_epr_j.py').read().split("fig,ax=plt.subplots")[0])     # styles, colours, hidden-level rows
def obs_pairs(y,ph,nph):
    y=np.asarray(y); keep=y>=0; t=np.flatnonzero(keep); ys=y[t]; ps=np.asarray(ph)[t]
    ch=np.flatnonzero(np.diff(ys))+1; ok=(t[ch]-t[ch-1])==1          # consecutive real minutes only
    N=np.zeros((nph,3,3)); np.add.at(N,(ps[ch[ok]],ys[ch-1][ok],ys[ch][ok]),1); return N
obs=[]
for g in COL:
    nph=1 if g=='per0DD' else 2; per={}
    for k,d in F.items():
        if k[0]!=g or d['deg']: continue
        y,ph,_=sqw2.seqs(g,k[1],0.5); N=obs_pairs(y,ph,nph)
        for p in range(nph):
            h=np.sum(ph==p)/60; per.setdefault(p,[]).append((R.j_cycle_count(N[p],h),R.epr_count_bound(N[p],h)[0]))
    for p,bn in enumerate(['whole'] if g=='per0DD' else ['light','dark']):
        a=np.array(per[p]); J,E=a[:,0],a[:,1]
        bJ=[np.median(rng.choice(J,len(J))) for _ in range(2000)]; bE=[np.median(rng.choice(E,len(E))) for _ in range(2000)]
        obs.append(dict(g=g,bn=bn,n=len(J),J=np.median(J),Jci=np.percentile(bJ,[2.5,97.5]),E=np.median(E),Eci=np.percentile(bE,[2.5,97.5])))
fig,axs=plt.subplots(1,2,figsize=(7.0,3.0),sharey=True)
for ax,rs,title in ((axs[0],rows,'(a) hidden states (posterior paths)'),(axs[1],obs,'(b) observable S/Q/W symbols')):
    jj=np.linspace(0,0.45,50)
    for A in (0.5,1,2): ax.plot(jj,A*jj,color=GREY,lw=0.5,ls=(0,(3,2)))
    ax.axvline(0,color=INK,lw=0.5)
    for r in rs:
        c=COL[r['g']]; filled=r['bn']!='light'
        al=0.35 if (r['g'],r['bn']) in WEAK else 1.0
        ax.errorbar(r['J'],r['E'],xerr=[[r['J']-r['Jci'][0]],[r['Jci'][1]-r['J']]],yerr=[[r['E']-r['Eci'][0]],[r['Eci'][1]-r['E']]],fmt='none',ecolor=c,elinewidth=0.9,alpha=al)
        ax.plot(r['J'],r['E'],'o',ms=5.5,mfc=c if filled else 'white',mec=c,mew=1.2,zorder=3,alpha=al)
    ax.set_title(title,fontsize=8,loc='left'); ax.set_xlim(-0.07,0.32)
    ax.set_xlabel(r'$J_{\mathrm{cycle}}$ (cycles h$^{-1}$; $+$ = S$\to$W$\to$Q/Reg$\to$S)')
axs[0].set_ylabel(r'EPR lower bound (nats h$^{-1}$)'); axs[0].set_ylim(0,0.42)
for A,(xa,ya) in {2:(0.205,0.415),1:(0.315,0.33),0.5:(0.315,0.17)}.items():
    for ax in axs: ax.text(xa,ya,f'$A$ = {A}',fontsize=5.5,color='#555555',ha='right',va='top')
for g in COL: axs[1].plot([],[],'o',ms=5,mfc=COL[g],mec=COL[g],label=NAME[g])
axs[1].plot([],[],'o',ms=5,mfc='white',mec=INK,mew=1,label='light'); axs[1].plot([],[],'o',ms=5,mfc=INK,mec=INK,label='dark / DD')
axs[0].axhspan(0,0.025,color='#d9d9d9',alpha=0.6,lw=0,zorder=0)
axs[0].text(0.31,0.004,'hidden-estimator offset +0.01–0.025 nats h$^{-1}$ (20f); faded: weak circulation',ha='right',va='bottom',fontsize=5.5,color='#555555')
axs[1].legend(frameon=False,fontsize=6.5,loc='upper left',ncol=2,handletextpad=0.3,columnspacing=0.8)
fig.tight_layout(pad=0.4); fig.savefig('figures/fig6b_epr_vs_j_hidden_vs_obs.pdf'); fig.savefig('figures/fig6b_epr_vs_j_hidden_vs_obs.png',dpi=300)
with open('fig6b_epr_vs_j.log','w') as fh:
    fh.write('block | hidden J, EPR | observable J, EPR | ratio hidden/observable (J, EPR)\n')
    for h,o in zip(rows,obs):
        fh.write(f"{h['g']:7s} {h['bn']:5s} | {h['J']:+.3f} {h['E']:.3f} | {o['J']:+.3f} [{o['Jci'][0]:+.3f},{o['Jci'][1]:+.3f}] {o['E']:.3f} [{o['Eci'][0]:.3f},{o['Eci'][1]:.3f}] | J {h['J']/o['J'] if abs(o['J'])>1e-3 else float('nan'):.2f}, EPR {h['E']/o['E']:.2f}\n")
print(open('fig6b_epr_vs_j.log').read())
