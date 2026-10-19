# PRESPEC 5 (eb4076b): pooled bout analysis P1 (phases per class) and P2 (consecutive-bout coupling).
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from common import *
from scipy.optimize import minimize
from scipy.special import logsumexp
from multiprocessing import Pool
XMIN={0:300,1:2}; GEN=['CS','per0','fumin','sss']
def bouts(g):
    f,ck,rk=SETS[g]; R=sio.loadmat(f)[rk]; out=[]   # rows: fly, window, segment, cls, dur, order
    seg=0
    for i in range(R.shape[0]):
        runs=[]
        def flush():
            nonlocal seg
            for r in runs[1:-1]: out.append((i+1,r[0],seg,r[1],r[2],r[3]))
            seg+=1
        for j in range(R.shape[1]):
            y=np.asarray(R[i,j]).ravel().astype(int)-1
            if y.size==0:
                if runs: flush(); runs=[]
                continue
            ch=np.r_[0,np.flatnonzero(np.diff(y))+1,y.size]
            for a,b in zip(ch[:-1],ch[1:]): runs.append((j+1,y[a],b-a,len(runs)))
        if runs: flush()
    return np.array(out)
# ---------- P1: discrete left-truncated hyperexponential ----------
def unpack(z,K): lam=np.exp(z[:K]); lw=z[K:]-logsumexp(z[K:]) if K>1 else np.zeros(1); return lam,lw
def nll(z,K,vals,cnt,nmin):
    lam,lw=unpack(z,K)
    if K>1: lw=np.r_[z[K:],0.]; lw=lw-logsumexp(lw)
    l1m=np.log(-np.expm1(-lam))
    lp=logsumexp(lw[None,:]-lam[None,:]*(vals[:,None]-1)+l1m[None,:],axis=1)
    lS=logsumexp(lw-lam*(nmin-1))
    return -(cnt*(lp-lS)).sum()
def fit(d,K,nmin,rng,starts=30):
    vals,cnt=np.unique(d,return_counts=True); best=None
    lo,hi=np.log(1/(d.max()-nmin+2)),np.log(1/1.5)
    for s in range(starts):
        z=np.r_[np.sort(rng.uniform(lo,hi,K)), rng.normal(0,1,K-1)]
        r=minimize(nll,z,args=(K,vals,cnt,nmin),method='L-BFGS-B')
        if best is None or r.fun<best.fun: best=r
    return best
def loglik(z,K,d,nmin):
    vals,cnt=np.unique(d,return_counts=True); return -nll(z,K,vals,cnt,nmin)
def report(z,K,nmin):
    lam=np.exp(z[:K]); lw=np.r_[z[K:],0.] if K>1 else np.zeros(1); lw=lw-logsumexp(lw)
    lq=lw-lam*(nmin-1); q=np.exp(lq-logsumexp(lq)); mean=(nmin-1)+1/(-np.expm1(-lam))
    o=np.argsort(mean); return q[o],mean[o]
def p1_job(a):
    g,c,K,fold,d,dtest=a; rng=np.random.default_rng([5,GEN.index(g),c,K,fold+1])
    r=fit(d,K,XMIN[c],rng); return (g,c,K,fold,r.x,-r.fun, loglik(r.x,K,dtest,XMIN[c]) if dtest is not None else np.nan)
# ---------- P2 ----------
def centred_rank(B):
    u=np.full(len(B),np.nan)
    for key in set(zip(B[:,0],B[:,1],B[:,3])):
        m=(B[:,0]==key[0])&(B[:,1]==key[1])&(B[:,3]==key[2]); n=m.sum()
        if n>=3:
            r=np.argsort(np.argsort(B[m,4]+1e-6*np.random.default_rng(0).random(n)))+1
            u[m]=(r-(n+1)/2)/n
    return u
def pairs(B,u):
    sw=[];ws=[]
    for k in range(len(B)-1):
        if B[k,2]==B[k+1,2] and B[k+1,5]==B[k,5]+1 and np.isfinite(u[k]) and np.isfinite(u[k+1]):
            (sw if B[k,3]==0 else ws).append((u[k],u[k+1],B[k,0]))
    return np.array(sw),np.array(ws)
corr=lambda P: np.corrcoef(P[:,0],P[:,1])[0,1]
if __name__=='__main__':
    data={g:bouts(g) for g in GEN}
    for g in GEN:
        B=data[g]; print(f"{g}: flies {len(np.unique(B[:,0]))}, sleep bouts {np.sum(B[:,3]==0)}, wake bouts {np.sum(B[:,3]==1)}",flush=True)
    jobs=[]
    for g in GEN:
        B=data[g]
        for c in [0,1]:
            m=B[:,3]==c
            for K in [1,2,3,4]:
                jobs.append((g,c,K,-1,B[m,4],None))
                for fold in range(5):
                    tr=m&(((B[:,0]-1)%5)!=fold); te=m&(((B[:,0]-1)%5)==fold)
                    jobs.append((g,c,K,fold,B[tr,4],B[te,4]))
    with Pool(4) as p: res=p.map(p1_job,jobs,chunksize=1)
    np.save('P1_res.npy',np.array(res,dtype=object),allow_pickle=True)
    print("\n=== P1: phases per bout class (discrete, left-truncated hyperexponential) ===")
    for g in GEN:
        for c,nm in [(0,'sleep'),(1,'wake ')]:
            n=np.sum((data[g][:,3]==c)); rows=[]
            for K in [1,2,3,4]:
                full=[r for r in res if r[:4]==(g,c,K,-1)][0]; cv=sum(r[6] for r in res if r[0]==g and r[1]==c and r[2]==K and r[3]>=0)
                k=2*K-1; L=full[5]; aicc=-2*L+2*k+2*k*(k+1)/max(n-k-1,1); bic=-2*L+k*np.log(n)
                rows.append((K,L,cv,aicc,bic,full[4]))
            cvbest=max(rows,key=lambda r:r[2])[0]
            print(f"\n{g} {nm} n={n}  CV-selected K={cvbest}  (AICc-best K={min(rows,key=lambda r:r[3])[0]}, BIC-best K={min(rows,key=lambda r:r[4])[0]})")
            for K,L,cv,aicc,bic,z in rows:
                q,mu=report(z,K,XMIN[c])
                print(f"   K={K} logL {L:11.2f}  CV logL {cv:11.2f}  AICc {aicc:10.1f}  BIC {bic:10.1f} | q {np.round(q,3).tolist()} means(s) {np.round(mu,0).tolist()} min branch {q.min()*n:.1f}")
    print("\n=== P2: consecutive-bout coupling (within-window centred ranks) ===")
    rng=np.random.default_rng(5); boot={}
    for g in GEN:
        B=data[g]; u=centred_rank(B); sw,ws=pairs(B,u); rsw,rws=corr(sw),corr(ws)
        # within-window permutation null
        keys=[np.flatnonzero((B[:,0]==a)&(B[:,1]==b)&(B[:,3]==c)) for a,b,c in set(zip(B[:,0],B[:,1],B[:,3]))]
        nsw=[];nws=[]
        for _ in range(2000):
            up=u.copy()
            for ix in keys: up[ix]=u[rng.permutation(ix)]
            s,w=pairs(B,up); nsw.append(corr(s)); nws.append(corr(w))
        nsw,nws=np.array(nsw),np.array(nws)
        psw=np.mean(np.abs(nsw)>=abs(rsw)); pws=np.mean(np.abs(nws)>=abs(rws))
        flies=np.unique(B[:,0]); bs=[]
        for _ in range(2000):
            pick=rng.choice(flies,len(flies))
            S=np.vstack([sw[sw[:,2]==f] for f in pick]); W=np.vstack([ws[ws[:,2]==f] for f in pick])
            bs.append((corr(S),corr(W)))
        bs=np.array(bs); boot[g]=bs; A=bs[:,0]-bs[:,1]
        ci=lambda v: np.percentile(v,[2.5,97.5])
        print(f"{g:6s} pairs SW {len(sw)} WS {len(ws)} | rho_SW {rsw:+.3f} CI {np.round(ci(bs[:,0]),3)} perm p {psw:.4f} | rho_WS {rws:+.3f} CI {np.round(ci(bs[:,1]),3)} perm p {pws:.4f} | A {rsw-rws:+.3f} CI {np.round(ci(A),3)}",flush=True)
    print("\nGenotype differences in A (bootstrap CI of difference):")
    for i,g1 in enumerate(GEN):
        for g2 in GEN[i+1:]:
            d=(boot[g1][:,0]-boot[g1][:,1])-(boot[g2][:,0]-boot[g2][:,1]); print(f"   A[{g1}] - A[{g2}]: {np.median(d):+.3f} CI {np.round(np.percentile(d,[2.5,97.5]),3)}")
