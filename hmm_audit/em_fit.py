# EM for a mixture of geometrics on excess e = d - (nmin-1) >= 1: the same likelihood as the
# discrete left-truncated hyperexponential (observed weights q), since truncation of a geometric
# is memoryless.
import numpy as np
def ll_mix(e,q,p):
    return np.log((q[None,:]*p[None,:]*(1-p[None,:])**(e[:,None]-1)).sum(1)).sum()
def em2(e,rng,starts=10,iters=3000,tol=1e-10):
    best=-np.inf; bq=bp=None; m=e.mean()
    for s in range(starts):
        p=np.sort(1/(m*np.exp(rng.uniform(-2,2,2))))[::-1].clip(1e-9,1-1e-12); q=rng.dirichlet([1,1]); old=-np.inf
        for _ in range(iters):
            lw=np.log(q)[None,:]+np.log(p)[None,:]+(e[:,None]-1)*np.log1p(-p)[None,:]
            mx=lw.max(1,keepdims=True); r=np.exp(lw-mx); tot=r.sum(1,keepdims=True); r/=tot
            L=(np.log(tot[:,0])+mx[:,0]).sum()
            q=r.mean(0).clip(1e-300); p=(r.sum(0)/(r*e[:,None]).sum(0)).clip(1e-12,1-1e-12)
            if L-old<tol: break
            old=L
        L=ll_mix(e,q,p)
        if L>best: best,bq,bp=L,q,p
    return best,bq,bp
def ll1(e):
    p=min(1/e.mean(),1-1e-12); return (np.log(p)+(e-1)*np.log1p(-p)).sum()
from scipy.optimize import minimize
def fit2(e,rng,starts=30):
    # EM multi-start, then one L-BFGS polish from the EM optimum in (log-rate, logit-q) space.
    L,q,p=em2(e,rng,starts=starts,iters=20000)
    f=lambda z: -ll_mix(e,np.array([1/(1+np.exp(-z[2])),1-1/(1+np.exp(-z[2]))]),1/(1+np.exp(-z[:2])))
    z0=np.r_[np.log(p/(1-p)), np.log(q[0]/q[1])]
    r=minimize(f,z0,method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-12,maxiter=4000))
    return max(L,-r.fun)
def fit2g(e,G=40,polish=3):
    # Profile grid over component means (m1<m2), weight q profiled by EM on q (concave in q),
    # then Nelder-Mead polish from the best `polish` grid points. Deterministic.
    ms=np.exp(np.linspace(np.log(1.0001),np.log(max(e.max(),2.0)*3),G)); P=1/ms
    lf=np.log(P)[None,:]+(e[:,None]-1)*np.log1p(-np.minimum(P,1-1e-12))[None,:]      # n x G
    i,j=np.triu_indices(G,1); A=lf[:,i]; Bm=lf[:,j]; mx=np.maximum(A,Bm); a=np.exp(A-mx); b=np.exp(Bm-mx)
    q=np.full(len(i),0.5)
    for _ in range(200):
        den=q*a+(1-q)*b; q=np.clip((q*a/den).mean(0),1e-12,1-1e-12)
    L=(np.log(q*a+(1-q)*b)+mx).sum(0); top=np.argsort(L)[::-1][:polish]; best=L[top[0]]; bq=np.array([q[top[0]],1-q[top[0]]]); bp=np.array([P[i[top[0]]],P[j[top[0]]]])
    f=lambda z: -ll_mix(e,np.array([1/(1+np.exp(-z[2])),1-1/(1+np.exp(-z[2]))]),1/(1+np.exp(-z[:2])))
    for t in top:
        p0=np.array([P[i[t]],P[j[t]]]).clip(1e-12,1-1e-9); z0=np.r_[np.log(p0/(1-p0)),np.log(q[t]/(1-q[t]))]
        r=minimize(f,z0,method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-12,maxiter=4000))
        if -r.fun>best: best=-r.fun; s_=1/(1+np.exp(-r.x[2])); bq=np.array([s_,1-s_]); bp=1/(1+np.exp(-r.x[:2]))
    return best,bq,bp
