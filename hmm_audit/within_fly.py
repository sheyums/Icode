# PRESPEC 7 (7028df0): K=1 vs K=2 sleep bouts within fly x ZT-bin, parametric bootstrap under N1 and N2.
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from pooled import *
from pooled_LD import bouts_t
from math import comb
from em_fit import fit2g, ll1
NMIN=300
def LR(d,rng):
    e=d-(NMIN-1); L2,q,p=fit2g(e); return 2*(L2-ll1(e)), q, p
def geo(rng,mean_e,n): return (NMIN-1)+rng.geometric(1/mean_e,n)
def block_job(a):
    g,gi,fly,b,d,days=a; rng=np.random.default_rng(np.random.SeedSequence([7,gi,int(fly),b]))
    obs,q2,p2=LR(d,rng); e=d-(NMIN-1)
    n1=[LR(geo(rng,e.mean(),len(d)).astype(float),rng)[0] for _ in range(100)]
    ud=np.unique(days); n2=[]
    for _ in range(100):
        sim=np.concatenate([geo(rng,e[days==u].mean(),np.sum(days==u)) for u in ud]).astype(float)
        n2.append(LR(sim,rng)[0])
    o=np.argsort(1/p2); q=q2[o]; mu=(NMIN-1)+1/p2[o]
    return (g,int(fly),b,len(d),obs,(1+np.sum(np.array(n1)>=obs))/101,(1+np.sum(np.array(n2)>=obs))/101,q[0],mu[0],mu[1])
def binom_p(k,n,p=0.05): return sum(comb(n,j)*p**j*(1-p)**(n-j) for j in range(k,n+1))
if __name__=='__main__':
    jobs=[]; excl={}
    for gi,g in enumerate(GEN):
        B=bouts_t(g); B=B[B[:,3]==0]; zb=((B[:,6]/3600)%24//4).astype(int); day=(B[:,6]//86400).astype(int)+1; excl[g]=0
        for fly in np.unique(B[:,0]):
            for b in range(6):
                m=(B[:,0]==fly)&(zb==b)
                if m.sum()>=20: jobs.append((g,gi,fly,b,B[m,4].astype(float),day[m]))
                elif m.sum()>0: excl[g]+=1
    print("blocks analysed:",{g:sum(j[0]==g for j in jobs) for g in GEN},"excluded (<20 bouts, >0):",excl,flush=True)
    with Pool(4) as p: res=p.map(block_job,jobs,chunksize=1)
    np.save('within_fly_res.npy',np.array(res,dtype=object),allow_pickle=True)
    for g in GEN:
        r=[x for x in res if x[0]==g]; n=len(r)
        if not n: print(g,'no blocks'); continue
        k1=sum(x[5]<0.05 for x in r); k2=sum(x[6]<0.05 for x in r)
        sig=[x for x in r if x[6]<0.05]
        extra=f" | sig-N2 blocks: median q_short {np.median([x[7] for x in sig]):.2f}, means {np.median([x[8] for x in sig]):.0f} s / {np.median([x[9] for x in sig]):.0f} s" if sig else ""
        print(f"{g:6s} blocks {n:3d} | N1: {k1}/{n} = {k1/n:.2f} (binom p {binom_p(k1,n):.2g}) | N2 (PRIMARY): {k2}/{n} = {k2/n:.2f} (binom p {binom_p(k2,n):.2g}){extra}",flush=True)
    print("\nBy ZT bin (all genotypes' blocks, descriptive): fraction p<0.05 under N2")
    for b in range(6):
        r=[x for x in res if x[2]==b]
        if r: print(f"  ZT {4*b:2d}-{4*b+4:2d}: {sum(x[6]<0.05 for x in r)}/{len(r)}")
