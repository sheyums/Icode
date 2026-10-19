# UNPLANNED (added after seeing A): identified set for windows where H1 is an emission twin of H3.
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from common import *
g_=np.linspace(-3,4,701); A,Bg=np.meshgrid(g_,g_,indexing='ij'); A=A.ravel(); Bg=Bg.ravel()
keep=np.abs(A*(1-Bg)-(1-A)*Bg)>1e-9; A,Bg=A[keep],Bg[keep]
S=np.zeros((A.size,3,3)); S[:,0,0]=A; S[:,0,2]=1-A; S[:,2,0]=Bg; S[:,2,2]=1-Bg; S[:,1,1]=1
Si=np.linalg.inv(S)
for g in ['CS','fumin']:
    rows=[]
    for w in windows(g):
        if np.abs(w['B'][0]-w['B'][2]).sum()>=1e-6: continue
        Pp=S@w['P']@Si; ok=Pp.min(axis=(1,2))>=0
        x=np.r_[xn(w['P']),xn(Pp[ok])]; y=np.r_[yn(w['P']),yn(Pp[ok])]; x=x[np.isfinite(x)]; y=y[np.isfinite(y)]
        mm=lambda v:(v.min(),v.max()) if v.size else (np.nan,np.nan)
        rows.append((*mm(x),*mm(y)))
    r=np.array(rows)
    if not len(r): print(g,'no H1==H3 windows'); continue
    print(f"{g}: H1==H3 windows {len(r)}; Xn undetermined {np.mean((r[:,0]<=0)&(r[:,1]>=0)):.3f}, Yn undetermined {np.mean((r[:,2]<=0)&(r[:,3]>=0)):.3f}; median widths X {np.nanmedian(r[:,1]-r[:,0]):.2f} Y {np.nanmedian(r[:,3]-r[:,2]):.2f}")
