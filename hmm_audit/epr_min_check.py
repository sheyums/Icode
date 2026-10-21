# POST HOC robustness of PRESPEC 15 step 1 (written after seeing results): the search is finite, so EPR_min
# is an UPPER bound on the true minimum. Re-search 40 fumin and 40 sss H1==H2 windows with EPR_min > 1e-3
# on a much wider grid ([-30,30]^2, 401 x 401) and 25 Nelder-Mead starts.
import sys,json; sys.path.insert(0,'.')
import epr_min as E
from common import *
from multiprocessing import Pool
g2=np.linspace(-30,30,401); A,B=np.meshgrid(g2,g2,indexing='ij'); A=A.ravel(); B=B.ravel(); k=np.abs(A*(1-B)-(1-A)*B)>1e-9
E.A,E.Bg=A[k],B[k]
def job(w):
    blk=w['blk']; P=w['P']; S=E.Smat(E.A,E.Bg,blk); Pp=S@P@np.linalg.inv(S); ok=Pp.min(axis=(1,2))>=-1e-12
    Pp=np.clip(Pp[ok],0,None); Pp/=Pp.sum(-1,keepdims=True); e=E.epr(Pp); a0=E.A[ok]; b0=E.Bg[ok]
    def obj(z):
        S1=E.Smat(z[0],z[1],blk); Q=S1@P@np.linalg.inv(S1)
        if Q.min()<-1e-12: return 1e6+1e3*(-Q.min())
        Q=np.clip(Q,0,None); Q/=Q.sum(-1,keepdims=True); v=E.epr(Q); return v if np.isfinite(v) else 1e6
    from scipy.optimize import minimize
    best=np.nanmin(np.where(np.isfinite(e),e,np.inf))
    for kk in np.argsort(np.where(np.isfinite(e),e,np.inf))[:25]:
        r=minimize(obj,[a0[kk],b0[kk]],method='Nelder-Mead',options=dict(xatol=1e-12,fatol=1e-14,maxiter=4000)); best=min(best,r.fun)
    return w['g'],w['old'],float(best)
if __name__=='__main__':
    R=json.load(open('epr_min.json')); W=[]
    for g in ['fumin','sss']:
        ws=[w for w in windows(g) if np.abs(w['B'][0]-w['B'][1]).sum()<1e-6]
        rr=[x for x in R if x['g']==g and tuple(x['blk'])==(0,1)]
        cand=[(w,x) for w,x in zip(ws,rr) if x['epr_min']>1e-3][:40]
        W+=[dict(g=g,blk=(0,1),P=w['P'],old=x['epr_min']) for w,x in cand]
    with Pool(4) as p: out=p.map(job,W)
    for g in ['fumin','sss']:
        o=np.array([(a,b) for gg,a,b in out if gg==g])
        print(f"{g}: {len(o)} windows | old median EPR_min {np.median(o[:,0]):.3f} -> wide search {np.median(o[:,1]):.3f} nats/h | now <= 1e-3: {np.sum(o[:,1]<=1e-3)} | wide/old ratio median {np.median(o[:,1]/o[:,0]):.3f}")
