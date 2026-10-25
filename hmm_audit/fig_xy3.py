# Three-way Xn-Yn comparison per fly, whole record (post hoc figure, user request; all values from existing outputs).
# (a) observable, model-free (PRESPEC 19): Xn_obs = 2 P(Q|waking) - 1, Yn_obs = 2 P(Q|onset) - 1, theta 0.5.
# (b) S/Q/W HMM (20b, theta 0.5): posterior route counts over the whole record (light + dark), per-path Xn, Yn
#     averaged over paths (Jensen); non-degenerate flies.
# (c) legacy binary per-window HMM, Viterbi counts summed over accepted windows (analysis_C outputs C_*.npy).
# Above the diagonal (Yn > Xn): forward circulation (c = (Yn - Xn)/2 > 0).
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sqw, sqw_ref as R
from cycle import fly as cyc_fly, stats as cyc_stats
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':8,'axes.linewidth':0.6,
                     'xtick.direction':'in','ytick.direction':'in'})
INK='#1a1a1a'; GREY='#9a9a9a'; COL={'CS':'#2a78d6','fumin':'#eb6834','per0LD':'#1baf7a','per0DD':'#eda100','sss':'#e87ba4','per0':'#52514e'}
F={(d['g'],d['fly']):d for d in map(json.loads,open('sqw_out/fit2_0.5.jsonl'))}; RT=json.load(open('sqw_out/route2_0.5.json'))
def obs(g):
    out=[]
    for i,p in enumerate(sqw.files(g)):
        x=np.loadtxt(p).astype(int); ev,Ln=cyc_fly(x); s=cyc_stats(ev,None,0.5,Ln/3600)
        if s and s['nSQ']+s['nSA']>0 and s['nQS']+s['nAS']>0:
            out.append((2*s['nSQ']/(s['nSQ']+s['nSA'])-1,2*s['nQS']/(s['nQS']+s['nAS'])-1))
    return np.array(out)
def hmm(g):
    out=[]
    for k,d in F.items():
        if k[0]!=g or d['deg']: continue
        rp=np.array(RT[f'{g}|{k[1]}']['routes_path']).sum(1)          # paths x 4, both phases pooled
        X,Y=R.xy_counts(rp)
        if np.isfinite(X).mean()>=0.5 and np.isfinite(Y).mean()>=0.5: out.append((np.nanmean(X),np.nanmean(Y)))
    return np.array(out)
def legacy(g):
    C=np.load(f'C_{g}.npy'); out=[]
    for f in np.unique(C[:,0]):
        V=C[C[:,0]==f][:,3:7].sum(0)    # N32, N312, N23, N213
        if V[0]+V[1]>0 and V[2]+V[3]>0: out.append(((V[1]-V[0])/(V[1]+V[0]),(V[3]-V[2])/(V[3]+V[2])))
    return np.array(out)
ROWS=[('CS',['CS'],'CS'),('per0',['per0LD','per0DD'],'per0'),('fumin',['fumin'],'fumin'),('sss',['sss'],'sss')]
NAME={'CS':'CS','per0LD':'$per^0$ LD','per0DD':'$per^0$ DD','fumin':'$fmn$','sss':'$sss$','per0':'$per^0$ (LD+DD)'}
fig,axs=plt.subplots(4,3,figsize=(7.0,9.2),sharex=True,sharey=True); log=[]
titles=['(a) observable, model-free','(b) S/Q/W HMM, posterior','(c) binary HMM, Viterbi (legacy)']
for r,(row,newg,legg) in enumerate(ROWS):
    for c in range(3):
        ax=axs[r,c]; ax.plot([-1,1],[-1,1],color=GREY,lw=0.6,ls=(0,(3,2))); ax.axhline(0,color='#d9d9d9',lw=0.4); ax.axvline(0,color='#d9d9d9',lw=0.4)
        sets=[(g,(obs if c==0 else hmm)(g)) for g in newg] if c<2 else [(legg,legacy(legg))]
        txt=[]
        for g,P in sets:
            ax.scatter(P[:,0],P[:,1],s=9,color=COL[g],edgecolor='white',lw=0.3,label=NAME[g] if (c==0 and len(newg)>1) or (c==2 and g=='per0') else None,zorder=3)
            fwd=np.mean(P[:,1]>P[:,0]); txt.append(f"{NAME[g]}: n {len(P)}, above diag. {fwd:.0%}")
            log.append(f"{titles[c]} | {g} | n {len(P)} | mean Xn {P[:,0].mean():+.3f} Yn {P[:,1].mean():+.3f} | mean c {(P[:,1]-P[:,0]).mean()/2:+.3f} | above diagonal {fwd:.2f}")
        ax.text(0.03,0.03,'\n'.join(txt),transform=ax.transAxes,fontsize=5.5,color='#333333',va='bottom')
        if (c==0 and len(newg)>1): ax.legend(frameon=False,fontsize=6,loc='upper left',handletextpad=0.1)
        if r==0: ax.set_title(titles[c],fontsize=8,loc='left')
        ax.set_xlim(-1.05,1.05); ax.set_ylim(-1.05,1.05); ax.set_aspect('equal'); ax.set_xticks([-1,0,1]); ax.set_yticks([-1,0,1])
        if c==0: ax.set_ylabel(f"{NAME[newg[0]] if len(newg)==1 else '$per^0$'}\n$Y_n$ (wake$\\to$sleep)")
        if r==3: ax.set_xlabel('$X_n$ (sleep$\\to$wake)')
fig.text(0.5,0.005,'Points: flies, whole record. Dashed: $Y_n=X_n$ (no net circulation); above it: forward cycle S$\\to$W$\\to$Q/Reg$\\to$S.',ha='center',fontsize=6.5,color='#333333')
fig.tight_layout(rect=(0,0.015,1,1),h_pad=0.4,w_pad=0.2)
fig.savefig('figures/fig7_xn_yn_three_way.pdf'); fig.savefig('figures/fig7_xn_yn_three_way.png',dpi=300)
open('fig7_xn_yn.log','w').write('\n'.join(log)+'\n'); print('\n'.join(log))
