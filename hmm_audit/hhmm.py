# Hierarchical (coarse x fine) categorical HMM with exact structured M-step. PRESPEC 12.
import os; os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
from hmmlearn.hmm import CategoricalHMM
from sklearn.utils import check_random_state
class HHMM(CategoricalHMM):
    def __init__(self, Nc=3, Nf=2, n_iter=500, tol=0.01, random_state=None):
        super().__init__(n_components=Nc*Nf, n_features=4, n_iter=n_iter, tol=tol, random_state=random_state,
                         init_params='se', params='ste')
        self.Nc, self.Nf = Nc, Nf
    def blk(self,c): return slice(c*self.Nf,(c+1)*self.Nf)
    def build(self, eps, A, G, Bin):
        Nc,Nf=self.Nc,self.Nf; K=Nc*Nf; P=np.zeros((K,K))
        for c in range(Nc):
            for f in range(Nf):
                s=c*Nf+f; P[s,self.blk(c)]=(1-eps[s])*A[c][f]
                for c2 in range(Nc):
                    if c2!=c: P[s,self.blk(c2)]=eps[s]*G[c,c2]*Bin[c2]
        return P
    def _init(self, X, lengths=None):
        super()._init(X, lengths)
        rs=check_random_state(self.random_state); Nc,Nf=self.Nc,self.Nf
        eps=10**rs.uniform(-3.3,-1.7,Nc*Nf)
        A=[]
        for c in range(Nc):
            a=np.full((Nf,Nf),0.0); stay=rs.uniform(0.85,0.97,Nf)
            for f in range(Nf): a[f]=(1-stay[f])/(Nf-1); a[f,f]=stay[f]
            A.append(a)
        G=rs.dirichlet(np.ones(Nc-1),Nc); Gm=np.zeros((Nc,Nc))
        for c in range(Nc): Gm[c,[k for k in range(Nc) if k!=c]]=G[c]
        Bin=[rs.dirichlet(np.ones(Nf)) for _ in range(Nc)]
        self.transmat_=self.build(eps,A,Gm,Bin)
    def _do_mstep(self, stats):
        super()._do_mstep(stats)       # start, emissions (and an unstructured trans we overwrite)
        N=stats['trans']; Nc,Nf=self.Nc,self.Nf; K=Nc*Nf; tiny=1e-300
        eps=np.zeros(K); A=[]; Gc=np.zeros((Nc,Nc)); Bc=np.zeros((Nc,Nf))
        for c in range(Nc):
            stay=N[self.blk(c)][:,self.blk(c)]
            A.append(stay/np.maximum(stay.sum(1,keepdims=True),tiny))
            for f in range(Nf):
                s=c*Nf+f; tot=N[s].sum(); lv=tot-stay[f].sum(); eps[s]=lv/tot if tot>0 else 0
            for c2 in range(Nc):
                if c2!=c:
                    L=N[self.blk(c)][:,self.blk(c2)]; Gc[c,c2]+=L.sum(); Bc[c2]+=L.sum(0)
        G=Gc/np.maximum(Gc.sum(1,keepdims=True),tiny); Bin=Bc/np.maximum(Bc.sum(1,keepdims=True),tiny)
        self.transmat_=self.build(eps,A,G,Bin); self._struct=(eps,A,G,Bin)
    def struct(self):
        return self._struct
def stationary(P):
    w,V=np.linalg.eig(P.T); v=np.abs(V[:,np.argmin(np.abs(w-1))].real); return v/v.sum()
def coarse_summary(m):
    eps,A,G,Bin=m.struct(); Nc,Nf=m.Nc,m.Nf; P=m.transmat_; pi=stationary(P); E=m.emissionprob_
    mf=np.array([(pi[m.blk(c)]*(1-E[m.blk(c),0])).sum()/pi[m.blk(c)].sum() for c in range(Nc)])
    dwell=[]
    for c in range(Nc):
        Q=P[m.blk(c)][:,m.blk(c)]; t=np.linalg.solve(np.eye(Nf)-Q,np.ones(Nf)); dwell.append(float(Bin[c]@t))
    o=np.argsort(mf); out={'mf':mf[o].tolist(),'dwell':[dwell[i] for i in o],'occ':[float(pi[m.blk(c)].sum()) for c in o],
       'fine_emis':[E[m.blk(c)].round(4).tolist() for c in o],'G':G[np.ix_(o,o)].round(6).tolist()}
    if Nc==3:
        S,T,W=o; out['Xn_c']=float((G[S,T]*G[T,W]-G[S,W])/(G[S,T]*G[T,W]+G[S,W])); out['Yn_c']=float((G[W,T]*G[T,S]-G[W,S])/(G[W,T]*G[T,S]+G[W,S]))
    return out
def fit(S,Nc,seed,Nf=2):
    m=HHMM(Nc=Nc,Nf=Nf,random_state=seed); X=np.concatenate(S).reshape(-1,1); m.fit(X,[len(s) for s in S]); return m
def score(m,S): return m.score(np.concatenate(S).reshape(-1,1),[len(s) for s in S])
