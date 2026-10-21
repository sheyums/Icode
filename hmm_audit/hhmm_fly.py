# Partial pooling hierarchical count HMM (PRESPEC 14): shared routing, fly-specific texture.
import numpy as np
from scipy.optimize import minimize
from hmmlearn import _utils
from hmmlearn.hmm import CategoricalHMM
from hhmm import stationary
from hhmm_day import HHMMDay
class HHMMFly(HHMMDay):
    def __init__(self, F, D, Nc=3, Nf=2, n_iter=300, tol=0.01, random_state=None):
        super().__init__(D, Nc=Nc, Nf=Nf, n_iter=n_iter, tol=tol, random_state=random_state); self.F=F
    def set_groups(self, flies, days): self.flies=np.asarray(flies,int); self.days=np.asarray(days,int)
    def _init(self, X, lengths=None):
        super()._init(X, lengths)
        eps=self.base[0]; self.eps_f=np.tile(eps,(self.F,1)); self.E_f=np.tile(self.emissionprob_,(self.F,1,1)); self._make_P()
    def _make_P(self):
        _,A,G,Bin=self.base; Nf=self.Nf
        self.Pfd={(i,d):self.build(np.clip(self.eps_f[i]*np.repeat(self.kappa[:,d],Nf),1e-12,0.5),A,G,Bin) for i in range(self.F) for d in range(self.D)}
    def _do_estep(self, X, lengths):
        impl={"scaling":self._fit_scaling,"log":self._fit_log}[self.implementation]
        stats=self._initialize_sufficient_statistics(); self._estep_begin(); lp=0; K=self.n_components
        self._t=np.zeros((self.F,self.D,K,K)); self._obs=np.zeros((self.F,K,self.n_features))
        E0=self.emissionprob_.copy(); T0=self.transmat_.copy()
        for sub_X,i,d in zip(_utils.split_X_lengths(X,lengths),self.flies,self.days):
            self.transmat_=self.Pfd[(i,d)]; self.emissionprob_=self.E_f[i]
            bt=stats['trans'].copy(); bo=stats['obs'].copy()
            lattice,logprob,post,fwd,bwd=impl(sub_X); self._accumulate_sufficient_statistics(stats,sub_X,lattice,post,fwd,bwd)
            self._t[i,d]+=stats['trans']-bt; self._obs[i]+=stats['obs']-bo; lp+=logprob
        self.emissionprob_=E0; self.transmat_=T0; return stats, lp
    def _do_mstep(self, stats):
        Nc,Nf,D,F=self.Nc,self.Nf,self.D,self.F; tiny=1e-300; T=self._t; Tsum=T.sum((0,1))
        if 's' in self.params:
            sp=np.maximum(stats['start'],0); self.startprob_=sp/sp.sum()
        self.E_f=self._obs/np.maximum(self._obs.sum(2,keepdims=True),tiny)
        self.emissionprob_=self.E_f.mean(0)
        A=[]; Gc=np.zeros((Nc,Nc)); Bc=np.zeros((Nc,Nf))
        for c in range(Nc):
            b=self.blk(c); stay=Tsum[b][:,b]; A.append(stay/np.maximum(stay.sum(1,keepdims=True),tiny))
            for c2 in range(Nc):
                if c2!=c: L=Tsum[b][:,self.blk(c2)]; Gc[c,c2]+=L.sum(); Bc[c2]+=L.sum(0)
            Lv=np.array([[[T[i,d][c*Nf+f].sum()-T[i,d][c*Nf+f,b].sum() for d in range(D)] for f in range(Nf)] for i in range(F)])
            Sv=np.array([[[T[i,d][c*Nf+f,b].sum() for d in range(D)] for f in range(Nf)] for i in range(F)])
            def nll(z):
                le=z[:F*Nf].reshape(F,Nf); lk=np.r_[z[F*Nf:],-z[F*Nf:].sum()]
                e=np.clip(np.exp(le[:,:,None]+lk[None,None,:]),1e-15,0.5); return -(Lv*np.log(e)+Sv*np.log1p(-e)).sum()
            z0=np.r_[np.log(np.clip(self.eps_f[:,c*Nf:(c+1)*Nf],1e-12,.5)).ravel(),np.log(self.kappa[c,:-1])]
            r=minimize(nll,z0,method='L-BFGS-B')
            self.eps_f[:,c*Nf:(c+1)*Nf]=np.exp(r.x[:F*Nf].reshape(F,Nf)); self.kappa[c]=np.exp(np.r_[r.x[F*Nf:],-r.x[F*Nf:].sum()])
        G=Gc/np.maximum(Gc.sum(1,keepdims=True),tiny); Bin=Bc/np.maximum(Bc.sum(1,keepdims=True),tiny)
        self.base=(np.exp(np.log(np.clip(self.eps_f,1e-15,1)).mean(0)),A,G,list(Bin)); self._struct=self.base
        self._make_P(); self.Pday=[np.mean([self.Pfd[(i,d)] for i in range(F)],0) for d in range(D)]; self.transmat_=np.mean(self.Pday,0)
def fit_fly(S, flies, days, F, D, Nc, seed, Nf=2, n_iter=300):
    m=HHMMFly(F,D,Nc=Nc,Nf=Nf,n_iter=n_iter,random_state=seed); m.set_groups(flies,days)
    m.fit(np.concatenate(S).reshape(-1,1),[len(s) for s in S]); return m
def summary_fly(m):
    eps,A,G,Bin=m.base; Nc=m.Nc; P=np.mean(m.Pday,0); pi=stationary(P); E=m.E_f.mean(0)
    mf=np.array([(pi[m.blk(c)]*(1-E[m.blk(c),0])).sum()/pi[m.blk(c)].sum() for c in range(Nc)]); o=np.argsort(mf)
    out={'mf':mf[o].tolist(),'kappa':m.kappa[o].tolist()}
    if Nc==3:
        S_,T_,W_=o; out['Xn_c']=float((G[S_,T_]*G[T_,W_]-G[S_,W_])/(G[S_,T_]*G[T_,W_]+G[S_,W_])); out['Yn_c']=float((G[W_,T_]*G[T_,S_]-G[W_,S_])/(G[W_,T_]*G[T_,S_]+G[W_,S_]))
    return out
