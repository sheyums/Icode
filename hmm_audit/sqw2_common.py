# PRESPEC 20b: HMM with shared emissions and phase-specific transition matrices.
# y: symbols (0=Q,1=W,2=S, -1 missing); ph: phase of each minute (transition INTO minute t uses P[ph[t]]).
import numpy as np
from sqw_common import label3, degenerate, xn, yn, stationary, Q, W, S
def kpar2(K,nP): return (K-1)+nP*K*(K-1)+2*K
def _E(B,y):                       # (T,S,K)
    E=B[:,:,np.clip(y,0,2)].transpose(2,0,1).copy(); E[y<0]=1.0; return E
def _forward(pi,P,E,ph):
    T,Sn,K=E.shape; A=np.empty((T,Sn,K)); c=np.empty((T,Sn))
    a=pi*E[0]; c[0]=a.sum(1); A[0]=a/c[0][:,None]
    for t in range(1,T):
        a=np.einsum('sk,skl->sl',A[t-1],P[:,ph[t]])*E[t]; c[t]=a.sum(1); A[t]=a/c[t][:,None]
    return A,c
def em(y,ph,nP,K,starts,seed,maxit=1000,tol=1e-7):
    """Batched Baum-Welch over `starts` random initialisations. Returns best (ll, pi, P[nP,K,K], B[K,3])."""
    rng=np.random.default_rng(seed); Sn=starts; y=np.asarray(y); ph=np.asarray(ph)
    pi=rng.dirichlet(np.ones(K),Sn); P=rng.dirichlet(np.ones(K),(Sn,nP,K))
    half=Sn//2; P[:half]=0.5*P[:half]+0.5*np.eye(K)           # half the starts diagonal-heavy
    B=rng.dirichlet(np.ones(3),(Sn,K)); ll_old=np.full(Sn,-np.inf); T=len(y)
    idx=[np.flatnonzero(ph[1:]==p)+1 for p in range(nP)]; onehot=np.stack([(y==v) for v in range(3)],1).astype(float)
    for it in range(maxit):
        E=_E(B,y); A,c=_forward(pi,P,E,ph); ll=np.log(c).sum(0)
        Bk=np.empty_like(A); Bk[-1]=1.0
        for t in range(T-2,-1,-1):
            Bk[t]=np.einsum('skl,sl->sk',P[:,ph[t+1]],E[t+1]*Bk[t+1])/c[t+1][:,None]
        G=A*Bk; G/=G.sum(2,keepdims=True); R=E*Bk/c[:,:,None]
        for p in range(nP):
            Xi=np.einsum('tsk,tsl->skl',A[idx[p]-1],R[idx[p]])*P[:,p]
            P[:,p]=Xi/np.maximum(Xi.sum(2,keepdims=True),1e-300)
        Bn=np.einsum('tsk,tv->skv',G,onehot); B=Bn/np.maximum(Bn.sum(2,keepdims=True),1e-300)
        pi=G[0]
        if np.all(ll-ll_old<tol*np.abs(ll)): break
        ll_old=ll
    E=_E(B,y); _,c=_forward(pi,P,E,ph); ll=np.log(c).sum(0); ll[~np.isfinite(ll)]=-np.inf
    s=int(np.argmax(ll)); return float(ll[s]),pi[s],P[s],B[s],it+1
def label3b(P,B,pi):
    """Label-align (H1,H2,H3) by emissions; apply the same permutation to every phase's P."""
    _,_,_=None,None,None
    h3=int(np.argmax(B[:,S])); rest=[i for i in range(3) if i!=h3]
    h2=rest[int(np.argmax(B[rest,W]))]; h1=[i for i in rest if i!=h2][0]; o=[h1,h2,h3]
    return P[:,o][:,:,o],B[o],pi[o]
def simulate2(pi,Pt,B,rng):
    """Pt: (T,K,K) transition into each minute (Pt[0] unused)."""
    T=len(Pt); K=len(pi); z=np.empty(T,int); z[0]=rng.choice(K,p=pi); u=rng.random(T)
    for t in range(1,T): z[t]=min(np.searchsorted(np.cumsum(Pt[t][z[t-1]]),u[t]),K-1)
    cb=np.cumsum(B,1); v=rng.random(T); y=np.array([min(np.searchsorted(cb[z[t]],v[t]),2) for t in range(T)])
    return z,y
def phases(T,offset=0):              # minute t at ZT ((t+offset) mod 1440)/60; light = ZT 0-12
    return (((np.arange(T)+offset)%1440)>=720).astype(int)
