import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; GRID='#e4e3df'; FFBS='#2a78d6'; VIT='#eb6834'
fig,axs=plt.subplots(2,2,figsize=(10,10.4),facecolor=SURF)
rows=[]
for ax,g in zip(axs.ravel(),['CS','per0','fumin','sss']):
    B=np.load(f'B_{g}.npy'); C=np.load(f'C_{g}.npy')
    key={(int(r[0]),int(r[1])):r for r in C}
    pts=[]
    for fly in np.unique(B[:,0]).astype(int):
        b=B[B[:,0]==fly]; c=np.array([key[(fly,int(w))] for w in b[:,1]])
        V=c[:,3:7].sum(0); F=c[:,7:11].sum(0)
        ph=lambda d,i,j: (d[j]-d[i])/(d[j]+d[i]) if d[j]+d[i]>0 else np.nan
        fx,fy,vx,vy=ph(F,0,1),ph(F,2,3),ph(V,0,1),ph(V,2,3)
        xr=(np.nanmean(b[:,3]),np.nanmean(b[:,4])); yr=(np.nanmean(b[:,6]),np.nanmean(b[:,7]))
        if np.isfinite(fx) and np.isfinite(fy): pts.append((fly,len(b),fx,fy,vx,vy,*xr,*yr))
    P=np.array(pts); rows+= [(g,*p) for p in P]
    ax.set_facecolor(SURF)
    for p in P:
        ax.plot([p[6],p[7]],[p[3],p[3]],color=FFBS,alpha=.18,lw=1.5,solid_capstyle='round',zorder=1)
        ax.plot([p[2],p[2]],[p[8],p[9]],color=FFBS,alpha=.18,lw=1.5,solid_capstyle='round',zorder=1)
    ax.scatter(P[:,4],P[:,5],s=34,facecolor='none',edgecolor=VIT,lw=1.4,zorder=2,label='Viterbi counts (legacy)')
    ax.scatter(P[:,2],P[:,3],s=40,color=FFBS,edgecolor=SURF,lw=1.2,zorder=3,label='FFBS posterior counts')
    ax.plot([],[],color=FFBS,alpha=.35,lw=2,label='equally good fits (identified range)')
    ax.axhline(0,color=GRID,lw=1,zorder=0); ax.axvline(0,color=GRID,lw=1,zorder=0)
    ax.set_xlim(-1.05,1.05); ax.set_ylim(-1.05,1.05); ax.set_aspect('equal')
    ax.set_xticks([-1,-.5,0,.5,1]); ax.set_yticks([-1,-.5,0,.5,1])
    for s in ax.spines.values(): s.set_color(GRID)
    ax.tick_params(colors=INK2,labelsize=9)
    ax.set_title(f"{g}   n = {len(P)} flies, {int(P[:,1].sum())} windows",color=INK,fontsize=11,loc='left')
    ax.set_xlabel('Xn = (N312 − N32)/(N312 + N32)',color=INK2,fontsize=9.5)
    ax.set_ylabel('Yn = (N213 − N23)/(N213 + N23)',color=INK2,fontsize=9.5)
axs[0,0].legend(loc='lower left',fontsize=8.5,frameon=False,labelcolor=INK2)
fig.suptitle('Per-fly routing coordinates, windows where H1 and H2 emit identically',color=INK,fontsize=13,x=0.06,ha='left')
fig.text(0.06,0.012,'Bars: range of window-mean model Xn/Yn over transition matrices with identical likelihood (inner bound, grid search).\nPoints are where the fit happened to stop; any point inside a fly\'s bars fits the data equally well.',color=INK2,fontsize=8.5)
fig.tight_layout(rect=(0,0.04,1,0.97))
fig.savefig('fly_Xn_Yn_scatter.png',dpi=160,facecolor=SURF)
import csv
with open('fly_Xn_Yn_table.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['genotype','fly','n_windows','Xn_FFBS','Yn_FFBS','Xn_Viterbi','Yn_Viterbi','Xn_range_lo','Xn_range_hi','Yn_range_lo','Yn_range_hi'])
    for r in rows: w.writerow([r[0],int(r[1]),int(r[2])]+[f'{v:.4f}' for v in r[3:]])
for g in ['CS','per0','fumin','sss']:
    a=np.array([r[3:] for r in rows if r[0]==g],float)
    print(f"{g:6s} flies {len(a):3d} | FFBS Xn median {np.median(a[:,0]):+.2f} Yn {np.median(a[:,1]):+.2f} | Viterbi Xn {np.nanmedian(a[:,2]):+.2f} Yn {np.nanmedian(a[:,3]):+.2f} | median range Xn [{np.median(a[:,4]):+.2f},{np.median(a[:,5]):+.2f}] Yn [{np.median(a[:,6]):+.2f},{np.median(a[:,7]):+.2f}]")
