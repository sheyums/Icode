# Synthetic check: route counts from Viterbi vs FFBS posterior, KNOWN truth, TRUE parameters.
import numpy as np, sys
rng = np.random.default_rng(1)
T = 14400; REPS = 200; M = 400
P = np.array([[1-1/60, 0.7/60, 0.3/60],
              [0.5/600, 1-1/600, 0.5/600],
              [0.4/1200, 0.6/1200, 1-1/1200]])
def kX(P): return P[0,1]/(P[0,1]+P[0,2])
Xn_pred = (P[2,0]*kX(P)-P[2,1])/(P[2,0]*kX(P)+P[2,1])
def simulate(B):
    z = np.empty(T, int); z[0] = rng.integers(3)
    u = rng.random(T); cP = P.cumsum(1)
    for t in range(1,T): z[t] = np.searchsorted(cP[z[t-1]], u[t])
    y = (rng.random(T) < B[z,1]).astype(int)
    return z, y
def counts(z):
    c = z[np.r_[0, np.flatnonzero(np.diff(z))+1]]
    n32 = np.sum((c[:-1]==2)&(c[1:]==1))
    n312 = np.sum((c[:-2]==2)&(c[1:-1]==0)&(c[2:]==1)) if len(c)>2 else 0
    return n32, n312
def viterbi(y,B):
    lP, lB = np.log(P), np.log(B); d = np.log(1/3)+lB[:,y[0]]; psi = np.zeros((T,3),int)
    for t in range(1,T):
        s = d[:,None]+lP; psi[t] = s.argmax(0); d = s.max(0)+lB[:,y[t]]
    z = np.empty(T,int); z[-1] = d.argmax()
    for t in range(T-1,0,-1): z[t-1] = psi[t,z[t]]
    return z
def ffbs(y,B):
    F = np.empty((T,3)); f = B[:,y[0]]/3; F[0] = f/f.sum()
    for t in range(1,T):
        f = (F[t-1]@P)*B[:,y[t]]; F[t] = f/f.sum()
    Z = np.empty((M,T),int); cF = F[-1].cumsum(); Z[:,-1] = np.searchsorted(cF, rng.random(M))
    for t in range(T-2,-1,-1):
        Q = F[t][:,None]*P[:,Z[:,t+1]]; Q /= Q.sum(0)
        Z[:,t] = (rng.random(M)[None,:] > Q.cumsum(0)).sum(0)
    return np.array([counts(z) for z in Z])
for label,B in [("DEGENERATE  B1=B2=[.4 .6]", np.array([[.4,.6],[.4,.6],[.98,.02]])),
                ("SEPARATED   B1=[.85 .15], B2=[.2 .8]", np.array([[.85,.15],[.2,.8],[.98,.02]]))]:
    rows=[]
    for r in range(REPS):
        z,y = simulate(B); tr = counts(z); vi = counts(viterbi(y,B)); D = ffbs(y,B)
        lo,hi = np.quantile(D,[.025,.975],axis=0)
        cov = (lo<=tr)&(tr<=hi)
        xd = D[:,0]+D[:,1]; xn = np.where(xd>0,(D[:,1]-D[:,0])/np.maximum(xd,1),np.nan)
        rows.append([*tr,*vi,*D.mean(0),*cov,np.nanmean(xn)])
    a = np.array(rows,float)
    def xn(n32,n312): s=n32+n312; return (n312-n32)/s
    print(f"\n== {label}  ({REPS} reps x {T} s, true params; Xn_pred={Xn_pred:+.3f})")
    print(f"{'':10s}{'N32':>8s}{'N312':>8s}{'Xn(pooled)':>12s}")
    for nm,i in [("truth",0),("Viterbi",2),("FFBS mean",4)]:
        print(f"{nm:10s}{a[:,i].mean():8.2f}{a[:,i+1].mean():8.2f}{xn(a[:,i].sum(),a[:,i+1].sum()):12.3f}")
    print(f"FFBS 95% interval coverage of TRUE count: N32 {a[:,6].mean():.3f}, N312 {a[:,7].mean():.3f}")
    print(f"mean per-window posterior Xn: {np.nanmean(a[:,8]):+.3f}  vs Xn_pred {Xn_pred:+.3f}")
    sys.stdout.flush()
