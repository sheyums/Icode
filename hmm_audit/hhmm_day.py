# Hierarchical count HMM with a DAY covariate (PRESPEC 13 amendment): coarse exit probabilities
# eps_{c,f,d} = eps_{c,f} * kappa_{c,d}, sum_d log kappa_{c,d} = 0. Fine dynamics A_c, emissions,
# coarse jump matrix G and entry distributions Bin are shared across days (routing assumed day-invariant).
import numpy as np
from scipy.optimize import minimize
from hmmlearn import _utils
from hmmlearn.hmm import CategoricalHMM
from hhmm import HHMM, stationary
class HHMMDay(HHMM):
    def __init__(self, D, Nc=3, Nf=2, n_iter=300, tol=0.01, random_state=None):
        super().__init__(Nc=Nc, Nf=Nf, n_iter=n_iter, tol=tol, random_state=random_state); self.D=D
    def set_days(self, days): self.days=np.asarray(days,int)
    def _init(self, X, lengths=None):
        super()._init(X, lengths)
        # recover struct of the random structured init
        self._unpack_init(); self.kappa=np.ones((self.Nc,self.D)); self._make_days()
    def _unpack_init(self):
        Nc,Nf=self.Nc,self.Nf; P=self.transmat_; K=Nc*Nf
        eps=np.array([1-P[s,self.blk(s//Nf)].sum() for s in range(K)])
        A=[P[self.blk(c)][:,self.blk(c)]/P[self.blk(c)][:,self.blk(c)].sum(1,keepdims=True) for c in range(Nc)]
        G=np.zeros((Nc,Nc)); Bin=[]
        for c2 in range(Nc):
            src=[c for c in range(Nc) if c!=c2][0]; L=P[self.blk(src)][:,self.blk(c2)]; Bin.append(L[0]/L[0].sum())
        for c in range(Nc):
            for c2 in range(Nc):
                if c2!=c: G[c,c2]=P[c*Nf,self.blk(c2)].sum()
            G[c]/=G[c].sum()
        self.base=(eps,A,G,Bin); self._struct=self.base
    def _make_days(self):
        eps,A,G,Bin=self.base; Nf=self.Nf
        self.Pday=[self.build(np.clip(eps*np.repeat(self.kappa[:,d],Nf),1e-12,0.5),A,G,Bin) for d in range(self.D)]
        self.transmat_=np.mean(self.Pday,0)
    def _do_estep(self, X, lengths):
        impl={"scaling":self._fit_scaling,"log":self._fit_log}[self.implementation]
        stats=self._initialize_sufficient_statistics(); self._estep_begin(); lp=0
        self._tday=np.zeros((self.D,)+self.transmat_.shape)
        for sub_X,d in zip(_utils.split_X_lengths(X,lengths),self.days):
            self.transmat_=self.Pday[d]; before=stats['trans'].copy()
            lattice,logprob,post,fwd,bwd=impl(sub_X)
            self._accumulate_sufficient_statistics(stats,sub_X,lattice,post,fwd,bwd)
            self._tday[d]+=stats['trans']-before; lp+=logprob
        self.transmat_=np.mean(self.Pday,0); return stats, lp
    def _do_mstep(self, stats):
        CategoricalHMM._do_mstep(self, stats)            # start + emissions
        Nc,Nf,D=self.Nc,self.Nf,self.D; T=self._tday; Tsum=T.sum(0); tiny=1e-300
        A=[]; Gc=np.zeros((Nc,Nc)); Bc=np.zeros((Nc,Nf)); eps=np.zeros(Nc*Nf); kappa=np.ones((Nc,D))
        for c in range(Nc):
            b=self.blk(c); stay=Tsum[b][:,b]; A.append(stay/np.maximum(stay.sum(1,keepdims=True),tiny))
            for c2 in range(Nc):
                if c2!=c: L=Tsum[b][:,self.blk(c2)]; Gc[c,c2]+=L.sum(); Bc[c2]+=L.sum(0)
            Lv=np.array([[T[d][c*Nf+f].sum()-T[d][c*Nf+f,b].sum() for d in range(D)] for f in range(Nf)])  # leave counts f x d
            Sv=np.array([[T[d][c*Nf+f,b].sum() for d in range(D)] for f in range(Nf)])
            def nll(z):
                le=z[:Nf]; lk=np.r_[z[Nf:],-z[Nf:].sum()]; e=np.clip(np.exp(le[:,None]+lk[None,:]),1e-15,0.5)
                return -(Lv*np.log(e)+Sv*np.log1p(-e)).sum()
            e0=np.clip(Lv.sum(1)/np.maximum(Lv.sum(1)+Sv.sum(1),tiny),1e-9,0.4)
            r=minimize(nll,np.r_[np.log(e0),np.zeros(D-1)],method='L-BFGS-B')
            eps[c*Nf:(c+1)*Nf]=np.exp(r.x[:Nf]); kappa[c]=np.exp(np.r_[r.x[Nf:],-r.x[Nf:].sum()])
        G=Gc/np.maximum(Gc.sum(1,keepdims=True),tiny); Bin=Bc/np.maximum(Bc.sum(1,keepdims=True),tiny)
        self.base=(eps,A,G,list(Bin)); self._struct=self.base; self.kappa=kappa; self._make_days()
def fit_day(S, days, D, Nc, seed, Nf=2, n_iter=300, init_from=None):
    m=HHMMDay(D,Nc=Nc,Nf=Nf,n_iter=n_iter,random_state=seed); m.set_days(days)
    if init_from is not None:
        m.init_params=''; m.startprob_=init_from.startprob_.copy(); m.emissionprob_=init_from.emissionprob_.copy()
        m.base=init_from.base if hasattr(init_from,'base') else init_from._struct; m._struct=m.base
        m.kappa=init_from.kappa.copy() if hasattr(init_from,'kappa') else np.ones((Nc,D)); m._make_days()
    m.fit(np.concatenate(S).reshape(-1,1),[len(s) for s in S]); return m
def score_day(m, S, days):
    tot=0.
    for s,d in zip(S,days):
        m.transmat_=m.Pday[d]; tot+=CategoricalHMM.score(m,s.reshape(-1,1))
    m.transmat_=np.mean(m.Pday,0); return tot
def fit_pooled_day(S, days, D, Nc, seed, n_sub=12, n_start=6, Nf=2):
    rs=np.random.RandomState(seed); idx=rs.choice(len(S),min(n_sub,len(S)),replace=False)
    sub=[S[i] for i in idx]; sd=[days[i] for i in idx]
    cands=[fit_day(sub,sd,D,Nc,seed*100+r,Nf) for r in range(n_start)]
    best=max(cands,key=lambda m:score_day(m,S,days))
    return fit_day(S,days,D,Nc,seed,Nf,init_from=best)
def coarse_summary_day(m):
    eps,A,G,Bin=m.base; Nc,Nf=m.Nc,m.Nf; P=np.mean(m.Pday,0); pi=stationary(P); E=m.emissionprob_
    mf=np.array([(pi[m.blk(c)]*(1-E[m.blk(c),0])).sum()/pi[m.blk(c)].sum() for c in range(Nc)]); o=np.argsort(mf)
    dwell=np.zeros((Nc,m.D))
    for d in range(m.D):
        for c in range(Nc):
            Q=m.Pday[d][m.blk(c)][:,m.blk(c)]; dwell[c,d]=Bin[c]@np.linalg.solve(np.eye(Nf)-Q,np.ones(Nf))
    out={'mf':mf[o].tolist(),'dwell_by_day':dwell[o].tolist(),'kappa':m.kappa[o].tolist(),'G':G[np.ix_(o,o)].tolist(),'fine_emis':[E[m.blk(c)].round(4).tolist() for c in o]}
    if Nc==3:
        S_,T_,W_=o; out['Xn_c']=float((G[S_,T_]*G[T_,W_]-G[S_,W_])/(G[S_,T_]*G[T_,W_]+G[S_,W_])); out['Yn_c']=float((G[W_,T_]*G[T_,S_]-G[W_,S_])/(G[W_,T_]*G[T_,S_]+G[W_,S_]))
    return out
