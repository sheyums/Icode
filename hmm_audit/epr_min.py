# PRESPEC 15 step 1.
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from common import *
from scipy.optimize import minimize
from multiprocessing import Pool
g_=np.linspace(-3,4,201); A,Bg=np.meshgrid(g_,g_,indexing='ij'); A=np.r_[A.ravel(),1.]; Bg=np.r_[Bg.ravel(),0.]
keep=np.abs(A*(1-Bg)-(1-A)*Bg)>1e-9; A,Bg=A[keep],Bg[keep]
def Smat(a,b,blk):
    i,j=blk; o=3-i-j; S=np.zeros(np.shape(a)+(3,3)); S[...,i,i]=a; S[...,i,j]=1-a; S[...,j,i]=b; S[...,j,j]=1-b; S[...,o,o]=1; return S
def stat(P):
    M=np.swapaxes(P,-1,-2)-np.eye(3); M[...,2,:]=1; b=np.zeros(P.shape[:-1]); b[...,2]=1
    return np.linalg.solve(M,b[...,None])[...,0]
def epr(P):
    pi=stat(P); out=np.zeros(P.shape[:-2]); inf=np.zeros(P.shape[:-2],bool)
    for a in range(3):
        for b in range(a+1,3):
            f=pi[...,a]*P[...,a,b]; r=pi[...,b]*P[...,b,a]
            z=(f<=0)&(r<=0); one=(f<=0)^(r<=0); inf|=one
            with np.errstate(all='ignore'): t=np.where(z|one,0.,(f-r)*np.log(f/r))
            out+=t
    return np.where(inf,np.inf,out)*3600
def flux(P): pi=stat(P); return (pi[...,2]*P[...,2,1]-pi[...,1]*P[...,1,2])*3600
def job(w):
    blk=w['blk']; P=w['P']; S=Smat(A,Bg,blk); Si=np.linalg.inv(S); Pp=S@P@Si; ok=Pp.min(axis=(1,2))>=-1e-12
    Pp=np.clip(Pp[ok],0,None); Pp/=Pp.sum(-1,keepdims=True); e=epr(Pp); J=flux(Pp)
    fin=np.isfinite(e); a0=A[ok]; b0=Bg[ok]
    def obj(z):
        S1=Smat(z[0],z[1],blk); Q=S1@P@np.linalg.inv(S1)
        if Q.min()<-1e-12: return 1e6+1e3*(-Q.min())
        Q=np.clip(Q,0,None); Q/=Q.sum(-1,keepdims=True); v=epr(Q); return v if np.isfinite(v) else 1e6
    best=np.inf; bz=None
    for k in np.argsort(np.where(fin,e,np.inf))[:5]:
        r=minimize(obj,[a0[k],b0[k]],method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-12,maxiter=2000))
        if r.fun<best: best,bz=r.fun,r.x
    gmin=e[fin].min() if fin.any() else np.inf
    if gmin<best: best=gmin; k=int(np.argmin(np.where(fin,e,np.inf))); bz=[a0[k],b0[k]]
    S1=Smat(bz[0],bz[1],blk); Q=np.clip(S1@P@np.linalg.inv(S1),0,None); Q/=Q.sum(-1,keepdims=True)
    with np.errstate(all='ignore'):
        ef=float(epr(P)); xq=float(xn(Q)); yq=float(yn(Q))
    return dict(g=w['g'],blk=blk,epr_fit=ef,epr_min=float(best),epr_max=float(e[fin].max()) if fin.any() else np.nan,n_inf=int((~fin).sum()),
                J_lo=float(J.min()),J_hi=float(J.max()),J_fit=float(flux(P)),Xn_at_min=xq,Yn_at_min=yq)
if __name__=='__main__':
    W=[]
    for g in ['CS','per0','fumin','sss']:
        for w in windows(g):
            for blk in [(0,1),(0,2)]:
                if np.abs(w['B'][blk[0]]-w['B'][blk[1]]).sum()<1e-6: W.append(dict(g=g,blk=blk,P=w['P']))
    with Pool(4) as p: R=p.map(job,W,chunksize=8)
    import json; json.dump(R,open('epr_min.json','w'))
    for g in ['CS','per0','fumin','sss']:
        for blk,nm in [((0,1),'H1==H2'),((0,2),'H1==H3')]:
            r=[x for x in R if x['g']==g and tuple(x['blk'])==blk]
            if not r: continue
            em=np.array([x['epr_min'] for x in r]); ef=np.array([x['epr_fit'] for x in r]); fin=np.isfinite(ef)
            print(f"{g:6s} {nm} windows {len(r):4d} | EPR_min <= 1e-3 nats/h: {np.mean(em<=1e-3):.3f} | median EPR_min {np.median(em):.2e} | median fitted EPR {np.median(ef[fin]) if fin.any() else np.nan:.3g} nats/h (fitted infinite in {np.sum(~fin)}) | J range contains 0: {np.mean([(x['J_lo']<=0)&(x['J_hi']>=0) for x in r]):.3f} | median Xn,Yn at min-EPR {np.nanmedian([x['Xn_at_min'] for x in r]):+.2f},{np.nanmedian([x['Yn_at_min'] for x in r]):+.2f}",flush=True)
