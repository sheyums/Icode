# PRESPEC 10 (ab86161): count data — C1 model-free transitions, C2 inactivity runs, C3 count HMM.
import os; os.environ['OMP_NUM_THREADS']='1'
import sys, glob; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
import numpy as np
from multiprocessing import Pool
from scipy.stats import wilcoxon
from pooled import fit, loglik
from scipy.special import logsumexp
from analysis_C import counts as route_counts
from common import xn, yn
from hmmlearn.hmm import CategoricalHMM
CD='/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/counts/'
GEN=['CS','fumin']; DAY=86400; HALF=43200
def load(g): return [np.loadtxt(f).astype(int) for f in sorted(glob.glob(CD+g+'/**/*.txt',recursive=True))]
def zruns(x):
    z=(x==0).astype(int); ch=np.r_[0,np.flatnonzero(np.diff(z))+1,len(z)]
    return [(a,b-a) for a,b in zip(ch[:-1],ch[1:]) if z[a]==1]
phase=lambda t: int((t%DAY)>=HALF)
# ---------------- C1 ----------------
def C1(g,X):
    rows=[]
    for fi,x in enumerate(X):
        mv=(x>0).astype(float); sl=[(a,n) for a,n in zruns(x) if n>=300]
        for (a1,n1),(a2,_) in zip(sl[:-1],sl[1:]):
            w0=a1+n1; w1=a2                      # wake bout = [w0, w1)
            if w1-w0<1200: continue
            mid=(w0+w1)//2; ref=mv[mid-150:mid+150].mean()
            rows.append((fi,'pre',phase(w1),mv[w1-360:w1-60].mean()-ref))
            rows.append((fi,'post',phase(w0),mv[w0+60:w0+360].mean()-ref))
    out=[]
    for kind in ['pre','post']:
        for p in [0,1]:
            per=[np.mean([r[3] for r in rows if r[0]==fi and r[1]==kind and r[2]==p]) for fi in range(len(X)) if sum(1 for r in rows if r[0]==fi and r[1]==kind and r[2]==p)>=3]
            per=np.array(per); st=wilcoxon(per) if len(per)>=5 else None
            out.append(f"{g:6s} {kind:4s} {'light' if p==0 else 'dark '} flies {len(per):2d} median contrast {np.median(per):+.4f} (IQR {np.percentile(per,25):+.4f}..{np.percentile(per,75):+.4f}) Wilcoxon p {st.pvalue if st else float('nan'):.3g}")
    return out
# ---------------- C2 ----------------
def C2job(a):
    g,p,K,fold,d,dt=a; r=fit(d,K,1,np.random.default_rng([10,GEN.index(g),p,K,fold+1]),starts=30)
    return (g,p,K,fold,r.x,-r.fun, loglik(r.x,K,dt,1) if dt is not None else np.nan)
def C2data(X):
    D=[]
    for fi,x in enumerate(X):
        rr=zruns(x)[1:-1]
        D+= [(fi,phase(a),n) for a,n in rr]
    return np.array(D)
def comps(z,K):
    lam=np.exp(z[:K]); lw=np.r_[z[K:],0.] if K>1 else np.zeros(1); lw=lw-logsumexp(lw)
    mean=1/(-np.expm1(-lam)); o=np.argsort(mean); return np.exp(lw)[o],mean[o],lam[o]
# ---------------- C3 ----------------
def sym(x): return np.minimum(x,3).reshape(-1,1)
def seqs(x,p):
    nd=len(x)//DAY + (1 if len(x)%DAY>=HALF+(HALF if p else 0) else 0)
    out=[]
    for d in range(int(np.ceil(len(x)/DAY))):
        s=x[d*DAY+p*HALF: d*DAY+(p+1)*HALF]
        if len(s)>=HALF*0.9: out.append(s)
    return out
def hfit(S,K,seed):
    m=CategoricalHMM(n_components=K,n_features=4,n_iter=500,tol=0.01,random_state=seed,init_params='ste',params='ste')
    X=np.concatenate([sym(s) for s in S]); L=[len(s) for s in S]; m.fit(X,L); return m, m.score(X,L)
def relabel(m):
    o=np.argsort(m.emissionprob_[:,0]); H3,H2=o[-1],o[0]; H1=[k for k in range(3) if k not in (H2,H3)][0]
    order=[H1,H2,H3]; P=m.transmat_[np.ix_(order,order)]; B=m.emissionprob_[order]; pi=m.startprob_[order]; return P,B,pi
def ffbs(y,P,B,pi,rng,M=100):
    T=len(y); F=np.empty((T,3)); f=pi*B[:,y[0]]; F[0]=f/f.sum()
    for t in range(1,T):
        f=(F[t-1]@P)*B[:,y[t]]; F[t]=f/f.sum()
    Z=np.empty((M,T),np.int8); Z[:,-1]=np.searchsorted(F[-1].cumsum(),rng.random(M)*F[-1].sum())
    for t in range(T-2,-1,-1):
        Q=F[t][:,None]*P[:,Z[:,t+1]]; Q/=Q.sum(0)
        Z[:,t]=np.minimum((rng.random(M)[None,:]>Q.cumsum(0)).sum(0),2)
    return np.array([route_counts(z.astype(int)) for z in Z])   # N32,N312,N23,N213 per draw
def C3job(a):
    g,fi,p,x=a; S=seqs(x,p); res={'g':g,'fly':fi+1,'phase':p,'nseq':len(S)}
    fits=[hfit(S,3,s) for s in range(8)]; lls=np.array([f[1] for f in fits]); best=int(np.argmax(lls))
    near=[k for k in range(8) if lls[k]>=lls[best]-0.5]
    xs=[];ys=[]
    for k in near:
        P,B,pi=relabel(fits[k][0]); xs.append(xn(P)); ys.append(yn(P))
    P,B,pi=relabel(fits[best][0])
    res.update(logL=lls[best],n_near=len(near),Xn_pred=xn(P),Yn_pred=yn(P),Xn_range=np.nanmax(xs)-np.nanmin(xs),Yn_range=np.nanmax(ys)-np.nanmin(ys),
               L1_12=np.abs(B[0]-B[1]).sum(),L1_13=np.abs(B[0]-B[2]).sum(),B=B.round(4).tolist(),P=P.round(6).tolist())
    rng=np.random.default_rng([10,GEN.index(g),fi,p]); C=sum(ffbs(s.astype(int).clip(max=3),P,B,pi,rng) for s in S)
    mC=C.mean(0); xr=lambda c:(c[1]-c[0])/(c[1]+c[0]) if c[1]+c[0]>0 else np.nan; yr=lambda c:(c[3]-c[2])/(c[3]+c[2]) if c[3]+c[2]>0 else np.nan
    dx=np.array([xr(c) for c in C]); dy=np.array([yr(c) for c in C])
    res.update(N32=mC[0],N312=mC[1],N23=mC[2],N213=mC[3],Xn=xr(mC),Yn=yr(mC),Xn_ci=np.nanpercentile(dx,[2.5,97.5]).tolist(),Yn_ci=np.nanpercentile(dy,[2.5,97.5]).tolist())
    if len(S)>=2:
        ho={}
        for K in [2,3,4]:
            b=max((hfit(S[:-1],K,s) for s in range(2)),key=lambda f:f[1])[0]; ho[K]=b.score(sym(S[-1]))/len(S[-1])
        res['heldout_ll_per_s']=ho
    return res
if __name__=='__main__':
    data={g:load(g) for g in GEN}
    print('=== C1: movement fraction around transitions minus mid-wake reference ===',flush=True)
    for g in GEN:
        for line in C1(g,data[g]): print(line,flush=True)
    print('\n=== C2: inactivity runs, all lengths ===',flush=True)
    jobs=[]
    for g in GEN:
        D=C2data(data[g])
        for p in [0,1]:
            m=D[:,1]==p
            for K in range(1,6):
                jobs.append((g,p,K,-1,D[m,2],None))
                for fold in range(5):
                    tr=m&((D[:,0]%5)!=fold); te=m&((D[:,0]%5)==fold); jobs.append((g,p,K,fold,D[tr,2],D[te,2]))
    with Pool(4) as pl: res=pl.map(C2job,jobs,chunksize=1)
    for g in GEN:
        D=C2data(data[g])
        for p in [0,1]:
            d=D[D[:,1]==p,2]; cv={K:sum(r[6] for r in res if r[:3]==(g,p,K) and r[3]>=0) for K in range(1,6)}; Kb=max(cv,key=cv.get)
            z=[r for r in res if r[:4]==(g,p,Kb,-1)][0][4]; w,mu,lam=comps(z,Kb)
            big=d[d>=300]; lp=np.log(w)[None,:]+np.log(-np.expm1(-lam))[None,:]-lam[None,:]*(big[:,None]-1); r_=np.exp(lp-logsumexp(lp,axis=1,keepdims=True))
            frac=r_[:,mu<300].sum()/len(big) if len(big) else np.nan
            print(f"{g:6s} {'light' if p==0 else 'dark '} runs {len(d)} (>=300 s: {len(big)}) | CV K={Kb} | CV logL {[round(cv[K],1) for K in range(1,6)]}\n       weights {np.round(w,4).tolist()} means(s) {np.round(mu,1).tolist()} | runs>=300 s assigned to components with mean<300 s: {frac:.3f}",flush=True)
    print('\n=== C3: 3-state count HMM per fly x phase ===',flush=True)
    jobs=[(g,fi,p,x) for g in GEN for fi,x in enumerate(data[g]) for p in [0,1]]
    with Pool(4) as pl: R3=pl.map(C3job,jobs,chunksize=1)
    import json
    with open('C3_results.json','w') as f: json.dump(R3,f,default=float)
    for g in GEN:
        for p in [0,1]:
            r=[x for x in R3 if x['g']==g and x['phase']==p]
            deg=np.mean([min(x['L1_12'],x['L1_13'])<0.05 for x in r]); unid=np.mean([(x['Xn_range']>0.2) or (x['Yn_range']>0.2) for x in r])
            print(f"{g:6s} {'light' if p==0 else 'dark '} flies {len(r)} | H1 degenerate (min L1<0.05) {deg:.2f} | median L1(H1,H2) {np.median([x['L1_12'] for x in r]):.3f} L1(H1,H3) {np.median([x['L1_13'] for x in r]):.3f} | route-unidentified (restart range>0.2) {unid:.2f} | median restarts within 0.5 nats {np.median([x['n_near'] for x in r]):.0f}")
            print(f"        median Xn {np.nanmedian([x['Xn'] for x in r]):+.3f} Yn {np.nanmedian([x['Yn'] for x in r]):+.3f} | median Xn_range {np.nanmedian([x['Xn_range'] for x in r]):.3f} Yn_range {np.nanmedian([x['Yn_range'] for x in r]):.3f}")
            if 'heldout_ll_per_s' in r[0]:
                ho=np.array([[x['heldout_ll_per_s'][K] for K in [2,3,4]] for x in r if 'heldout_ll_per_s' in x])
                print(f"        held-out log-lik per s, median over flies: K2 {np.median(ho[:,0]):.5f} K3 {np.median(ho[:,1]):.5f} K4 {np.median(ho[:,2]):.5f}; flies where K3>K2 {np.mean(ho[:,1]>ho[:,0]):.2f}, K4>K3 {np.mean(ho[:,2]>ho[:,1]):.2f}",flush=True)
