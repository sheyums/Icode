# PRESPEC 9: edge-free within-fly two-phase test at 4 h (A) and 12 h (B) blocks.
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from pooled import *
from em_fit import fit2g, ll1
from math import comb, log
from scipy.stats import chi2
NMIN=300
def bouts_edgefree(g):
    f,ck,rk=SETS[g]; R=sio.loadmat(f)[rk]; out=[]
    for i in range(R.shape[0]):
        runs=[]; t0=0
        def flush():
            for r in runs[1:-1]:
                if not r[5]: out.append((i+1,r[0],r[1],r[2],r[4]))
        for j in range(R.shape[1]):
            y=np.asarray(R[i,j]).ravel().astype(int)-1
            if y.size==0:
                if runs: flush(); runs=[]
                t0+=14400; continue
            ch=np.r_[0,np.flatnonzero(np.diff(y))+1,y.size]; nr=len(ch)-1
            for k,(a,b) in enumerate(zip(ch[:-1],ch[1:])): runs.append((j+1,y[a],b-a,len(runs),t0+a,k==nr-1))
            t0+=y.size
        if runs: flush()
    return np.array(out)   # fly, win, cls, dur, start
def LR(d):
    e=d-(NMIN-1); L2,q,p=fit2g(e); return 2*(L2-ll1(e)),q,p
def geo(rng,e,n): return (NMIN-1)+rng.geometric(min(1/e.mean(),1-1e-12),n)
def sim(rng,d,groups):
    e=d-(NMIN-1); return np.concatenate([geo(rng,e[groups==u],np.sum(groups==u)) for u in np.unique(groups)]).astype(float)
def job(a):
    tag,g,gi,fly,b,d,day,sub=a; rng=np.random.default_rng(np.random.SeedSequence([9,tag,gi,int(fly),b]))
    obs,q,p=LR(d); one=np.zeros(len(d),int)
    nulls={'N1':one,'N2':day} if tag==0 else {'N1':one,'N2':day,'N3':day*10+sub}
    ps={k:(1+sum(LR(sim(rng,d,grp))[0]>=obs for _ in range(100)))/101 for k,grp in nulls.items()}
    o=np.argsort(1/p); return (tag,g,int(fly),b,len(d),obs,ps,q[o][0],(NMIN-1)+1/p[o][0],(NMIN-1)+1/p[o][1])
def binom_p(k,n,p=0.05): return sum(comb(n,j)*p**j*(1-p)**(n-j) for j in range(k,n+1))
if __name__=='__main__':
    jobs=[]; nblk={}
    for gi,g in enumerate(GEN):
        B=bouts_edgefree(g); B=B[B[:,2]==0]; zt=(B[:,4]/3600)%24; day=(B[:,4]//86400).astype(int); b4=(zt//4).astype(int); ph=(zt>=12).astype(int)
        for fly in np.unique(B[:,0]):
            for b in range(6):
                m=(B[:,0]==fly)&(b4==b)
                if m.sum()>=20: jobs.append((0,g,gi,fly,b,B[m,3].astype(float),day[m],b4[m]))
            for b in range(2):
                m=(B[:,0]==fly)&(ph==b)
                if m.sum()>=20: jobs.append((1,g,gi,fly,b,B[m,3].astype(float),day[m],b4[m]))
    for tag in [0,1]: print(('A 4h' if tag==0 else 'B 12h'),'blocks:',{g:sum(j[0]==tag and j[1]==g for j in jobs) for g in GEN},flush=True)
    with Pool(4) as pl: res=pl.map(job,jobs,chunksize=1)
    np.save('within_fly9_res.npy',np.array(res,dtype=object),allow_pickle=True)
    for tag,lab,prim in [(0,'A: 4 h blocks','N2'),(1,'B: 12 h blocks','N3')]:
        print(f"\n=== {lab} (primary null {prim}) ===")
        for g in GEN:
            r=[x for x in res if x[0]==tag and x[1]==g]; n=len(r)
            if not n: print(f"{g:6s} no blocks"); continue
            pp=np.array([x[6][prim] for x in r]); F=-2*np.log(pp).sum(); fp=chi2.sf(F,2*n)
            cnt=' | '.join(f"{k}: {sum(x[6][k]<0.05 for x in r)}/{n} (binom p {binom_p(sum(x[6][k]<0.05 for x in r),n):.2g})" for k in r[0][6])
            sig=[x for x in r if x[6][prim]<0.05]
            ex=f" | sig: q_short {np.median([x[7] for x in sig]):.2f}, means {np.median([x[8] for x in sig]):.0f}/{np.median([x[9] for x in sig]):.0f} s" if sig else ""
            print(f"{g:6s} blocks {n:3d} | FISHER ({prim}) p {fp:.3g} | {cnt}{ex}",flush=True)
