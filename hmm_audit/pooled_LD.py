# PRESPEC 5b (fab217d): P1 and P2 within LIGHT (ZT0-12) and DARK (ZT12-24).
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from pooled import *
def bouts_t(g):
    f,ck,rk=SETS[g]; R=sio.loadmat(f)[rk]; out=[]; seg=0
    for i in range(R.shape[0]):
        runs=[]; t0=0
        def flush():
            nonlocal seg
            for r in runs[1:-1]: out.append((i+1,r[0],seg,r[1],r[2],r[3],r[4]))
            seg+=1
        for j in range(R.shape[1]):
            y=np.asarray(R[i,j]).ravel().astype(int)-1
            if y.size==0:
                if runs: flush(); runs=[]
                t0+=14400; continue
            ch=np.r_[0,np.flatnonzero(np.diff(y))+1,y.size]
            for a,b in zip(ch[:-1],ch[1:]): runs.append((j+1,y[a],b-a,len(runs),t0+a))
            t0+=y.size
        if runs: flush()
    return np.array(out)   # fly, win, seg, cls, dur, order, start_frame
def p1_job_ld(a):
    g,c,K,fold,d,dtest=a; rng=np.random.default_rng([5,GEN.index(g),c,K,fold+20])
    r=fit(d,K,XMIN[c],rng); return (g,c,K,fold,r.x,-r.fun, loglik(r.x,K,dtest,XMIN[c]) if dtest is not None else np.nan)
def phase(t): return ((t/3600)%24>=12).astype(int)   # 0 light, 1 dark
if __name__=='__main__':
    data={g:bouts_t(g) for g in GEN}
    print("=== crossing fraction (bouts whose span crosses ZT0/ZT12) ===")
    for g in GEN:
        B=data[g]; ph=phase(B[:,6]); ph_end=phase(B[:,6]+B[:,4]-1)
        cross=(ph!=ph_end)|(B[:,4]>=43200)
        for p,pn in [(0,'light'),(1,'dark ')]:
            for c,cn in [(0,'sleep'),(1,'wake')]:
                m=(ph==p)&(B[:,3]==c); print(f"  {g:6s} {pn} {cn:5s} n={m.sum():5d} crossing {cross[m].mean():.3f}")
    jobs=[]
    for g in GEN:
        B=data[g]; ph=phase(B[:,6])
        for p in [0,1]:
            for c in [0,1]:
                m=(B[:,3]==c)&(ph==p)
                for K in [1,2,3,4]:
                    jobs.append((g,c,K,-1-10*p,B[m,4],None))
                    for fold in range(5):
                        tr=m&(((B[:,0]-1)%5)!=fold); te=m&(((B[:,0]-1)%5)==fold)
                        jobs.append((g,c,K,fold+10*p,B[tr,4],B[te,4]))
    with Pool(4) as pl: res=pl.map(p1_job_ld,jobs,chunksize=1)
    np.save('P1_LD_res.npy',np.array(res,dtype=object),allow_pickle=True)
    print("\n=== P1 by phase ===")
    for g in GEN:
        for p,pn in [(0,'light'),(1,'dark')]:
            for c,cn in [(0,'sleep'),(1,'wake ')]:
                n=np.sum((data[g][:,3]==c)&(phase(data[g][:,6])==p)); rows=[]
                for K in [1,2,3,4]:
                    full=[r for r in res if r[:4]==(g,c,K,-1-10*p)][0]
                    cv=sum(r[6] for r in res if r[0]==g and r[1]==c and r[2]==K and 10*p<=r[3]<10*p+5)
                    k=2*K-1; L=full[5]; aicc=-2*L+2*k+2*k*(k+1)/max(n-k-1,1); bic=-2*L+k*np.log(n)
                    rows.append((K,L,cv,aicc,bic,full[4]))
                best=max(rows,key=lambda r:r[2]); q,mu=report(best[5],best[0],XMIN[c])
                print(f"{g:6s} {pn:5s} {cn} n={n:5d} CV K={best[0]} (AICc {min(rows,key=lambda r:r[3])[0]}, BIC {min(rows,key=lambda r:r[4])[0]}) | CV logL by K {[round(r[2],1) for r in rows]} | q {np.round(q,3).tolist()} means {np.round(mu,0).tolist()} min branch {q.min()*n:.1f}")
    print("\n=== P2 by phase ===")
    rng=np.random.default_rng(5)
    for g in GEN:
        B=data[g]; u=centred_rank(B[:,:6]); ph=phase(B[:,6])
        keys=[np.flatnonzero((B[:,0]==a)&(B[:,1]==b)&(B[:,3]==c)) for a,b,c in set(zip(B[:,0],B[:,1],B[:,3]))]
        def pp(uu):
            o={}
            for p in [0,1]:
                sw=[];ws=[]
                for k in range(len(B)-1):
                    if ph[k]==p and B[k,2]==B[k+1,2] and B[k+1,5]==B[k,5]+1 and np.isfinite(uu[k]) and np.isfinite(uu[k+1]):
                        (sw if B[k,3]==0 else ws).append((uu[k],uu[k+1],B[k,0]))
                o[p]=(np.array(sw),np.array(ws))
            return o
        obs=pp(u); nulls={0:[],1:[]}
        for _ in range(2000):
            up=u.copy()
            for ix in keys: up[ix]=u[rng.permutation(ix)]
            o=pp(up)
            for p in [0,1]: nulls[p].append((corr(o[p][0]),corr(o[p][1])))
        flies=np.unique(B[:,0]); bs={0:[],1:[]}
        for _ in range(2000):
            pick=rng.choice(flies,len(flies))
            for p in [0,1]:
                sw,ws=obs[p]; S=np.vstack([sw[sw[:,2]==f] for f in pick if np.any(sw[:,2]==f)]); W=np.vstack([ws[ws[:,2]==f] for f in pick if np.any(ws[:,2]==f)])
                bs[p].append((corr(S),corr(W)))
        for p,pn in [(0,'light'),(1,'dark ')]:
            sw,ws=obs[p]; rsw,rws=corr(sw),corr(ws); nl=np.array(nulls[p]); b=np.array(bs[p]); A=b[:,0]-b[:,1]
            ci=lambda v: np.round(np.percentile(v,[2.5,97.5]),3)
            print(f"{g:6s} {pn} pairs SW {len(sw)} WS {len(ws)} | rho_SW {rsw:+.3f} CI {ci(b[:,0])} p {np.mean(np.abs(nl[:,0])>=abs(rsw)):.4f} | rho_WS {rws:+.3f} CI {ci(b[:,1])} p {np.mean(np.abs(nl[:,1])>=abs(rws)):.4f} | A {rsw-rws:+.3f} CI {ci(A)}",flush=True)
        d=(np.array(bs[0])[:,0]-np.array(bs[0])[:,1])-(np.array(bs[1])[:,0]-np.array(bs[1])[:,1])
        print(f"{g:6s} A[light]-A[dark] {np.median(d):+.3f} CI {np.round(np.percentile(d,[2.5,97.5]),3)}",flush=True)
