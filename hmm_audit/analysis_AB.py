import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from common import *
g_=np.linspace(-3,4,701); A,Bg=np.meshgrid(g_,g_,indexing='ij'); A=A.ravel(); Bg=Bg.ravel()
det=A*(1-Bg)-(1-A)*Bg; keep=np.abs(det)>1e-9; A,Bg,det=A[keep],Bg[keep],det[keep]
S=np.zeros((A.size,3,3)); S[:,0,0]=A; S[:,0,1]=1-A; S[:,1,0]=Bg; S[:,1,1]=1-Bg; S[:,2,2]=1
Si=np.linalg.inv(S)
def loglik(P,B,y,pi0):
    f=pi0*B[:,y[0]]; ll=np.log(f.sum()); f/=f.sum()
    for t in range(1,len(y)):
        f=(f@P)*B[:,y[t]]; s=f.sum(); ll+=np.log(s); f/=s
    return ll
for g in (sys.argv[1:] or ['CS','fumin']):
    W=windows(g); L1=np.array([np.abs(w['B'][0]-w['B'][1]).sum() for w in W])
    print(f"\n=== {g}: {len(W)} accepted windows, {len({w['fly'] for w in W})} flies")
    for th in [1e-6,0.05,0.2]: print(f"  A: L1 < {th:g}: {np.mean(L1<th):.3f}  ({np.sum(L1<th)})")
    print("  A: median B rows:", np.round(np.median([w['B'] for w in W],0),4).tolist())
    rows=[]; spot=[]
    for w in W:
        if np.abs(w['B'][0]-w['B'][1]).sum()>=1e-6: continue
        Pp=S@w['P']@Si; ok=(Pp.min(axis=(1,2))>=0)
        with np.errstate(all='ignore'):
            x=np.r_[xn(w['P']),xn(Pp[ok])]; yv=np.r_[yn(w['P']),yn(Pp[ok])]
        x=x[np.isfinite(x)]; yv=yv[np.isfinite(yv)]
        mm=lambda v:(v.min(),v.max()) if v.size else (np.nan,np.nan)
        rows.append((w['fly'],w['win'],xn(w['P']),*mm(x),yn(w['P']),*mm(yv)))
        if w['win']==min(v['win'] for v in W if v['fly']==w['fly']) and len(spot)<3:
            xv=np.nan_to_num(xn(Pp[ok]),nan=-9); k=np.flatnonzero(ok)[np.argmax(xv)]
            pi0=np.ones(3)/3
            spot.append((w['fly'],w['win'],loglik(w['P'],w['B'],w['y'],pi0),
                         loglik(Pp[k],w['B'],w['y'],pi0@Si[k]),xn(w['P']),xn(Pp[k])))
    r=np.array(rows)
    print(f"  B: exactly-degenerate windows analysed: {len(r)}")
    print(f"  B: fitted Xn_pred  median {np.nanmedian(r[:,2]):+.3f};  identified interval median [{np.nanmedian(r[:,3]):+.3f}, {np.nanmedian(r[:,4]):+.3f}]")
    print(f"  B: fitted Yn_pred  median {np.nanmedian(r[:,5]):+.3f};  identified interval median [{np.nanmedian(r[:,6]):+.3f}, {np.nanmedian(r[:,7]):+.3f}]")
    print(f"  B: Xn undefined in all models: {np.mean(np.isnan(r[:,3])):.3f}; Yn: {np.mean(np.isnan(r[:,6])):.3f}")
    print(f"  B: Xn direction UNDETERMINED (interval contains 0): {np.mean((r[:,3]<=0)&(r[:,4]>=0)):.3f}")
    print(f"  B: Yn direction UNDETERMINED (interval contains 0): {np.mean((r[:,6]<=0)&(r[:,7]>=0)):.3f}")
    print(f"  B: Xn interval width median {np.nanmedian(r[:,4]-r[:,3]):.3f}; Yn {np.nanmedian(r[:,7]-r[:,6]):.3f}")
    for s in spot: print(f"  spot fly{s[0]} win{s[1]}: logL P {s[2]:.6f} | P' {s[3]:.6f} | Xn_pred {s[4]:+.3f} -> {s[5]:+.3f}")
    np.save(f'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit/B_{g}.npy',r)
