# PRESPEC 20: shared HMM fitting for 3-symbol (Q=0, W=1, S=2) minute sequences.
import numpy as np, warnings, logging
from hmmlearn.hmm import CategoricalHMM
logging.getLogger('hmmlearn').setLevel(logging.ERROR); warnings.filterwarnings('ignore')
Q,W,S=0,1,2; NSTEP=240
def kpar(K): return (K-1)+K*(K-1)+2*K
def nchanges(y): return max(2,int(np.sum(np.diff(y)!=0)))
def fit(y,K,starts,seed):
    """Best-of-`starts` EM fit; returns (logL, startprob, P, B)."""
    X=np.asarray(y).reshape(-1,1); best=None
    for s in range(starts):
        m=CategoricalHMM(n_components=K,n_features=3,n_iter=500,tol=1e-6,random_state=seed*1000+s,init_params='ste',params='ste')
        try: m.fit(X); ll=m.score(X)
        except Exception: continue
        if not np.isfinite(ll): continue
        if best is None or ll>best[0]: best=(ll,m.startprob_.copy(),m.transmat_.copy(),m.emissionprob_.copy())
    return best
def bic(ll,K,n): return -2*ll+kpar(K)*np.log(n)
def label3(P,B,pi=None):
    """Reorder a K=3 fit to (H1, H2, H3): H3 = max P(S); of the rest H2 = max P(W)."""
    h3=int(np.argmax(B[:,S])); rest=[i for i in range(3) if i!=h3]
    h2=rest[int(np.argmax(B[rest,W]))]; h1=[i for i in rest if i!=h2][0]; o=[h1,h2,h3]
    P2=P[np.ix_(o,o)]; B2=B[o]; return (P2,B2,None if pi is None else pi[o])
def tv(a,b): return 0.5*np.abs(a-b).sum()
def degenerate(B,thr=0.05): return min(tv(B[0],B[1]),tv(B[0],B[2]))<thr
def xn(P):
    with np.errstate(all='ignore'):
        kX=P[0,1]/(P[0,1]+P[0,2]); a=P[2,0]*kX; return (a-P[2,1])/(a+P[2,1])
def yn(P):
    with np.errstate(all='ignore'):
        kY=P[0,2]/(P[0,1]+P[0,2]); a=P[1,0]*kY; return (a-P[1,2])/(a+P[1,2])
def simulate(pi,P,B,n,rng):
    K=len(pi); z=np.empty(n,int); y=np.empty(n,int); z[0]=rng.choice(K,p=pi)
    for t in range(1,n): z[t]=rng.choice(K,p=P[z[t-1]])
    for t in range(n): y[t]=rng.choice(3,p=B[z[t]])
    return z,y
def stationary(P):
    w,v=np.linalg.eig(P.T); p=np.real(v[:,np.argmin(np.abs(w-1))]); return p/p.sum()
