# PRESPEC 6: heterogeneity-only null for the pooled multi-phase result.
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from pooled import *
def T1(d,nmin,rng):
    r1=fit(d,1,nmin,rng,starts=3); r2=fit(d,2,nmin,rng,starts=10); return r1.fun-r2.fun
def T2(groups,nmin):
    num=den=0.
    for d in groups:
        e=d-(nmin-1); n=len(e); m=e.mean(); s2=e.var(ddof=1)
        num+=(n-1)*(s2/m**2-(1-1/m)); den+=n-1
    return num/den
def job(a):
    g,c,groups,seed=a; rng=np.random.default_rng(seed); nmin=XMIN[c]
    if seed is None: return None
    sim=[ (nmin-1)+rng.geometric(1/(d-(nmin-1)).mean(),len(d)) for d in groups]
    return T1(np.concatenate(sim),nmin,rng), T2(sim,nmin)
if __name__=='__main__':
    for c,cn in [(0,'SLEEP'),(1,'WAKE')]:
        print(f"\n=== {cn} ===",flush=True)
        for gi,g in enumerate(GEN):
            B=bouts(g); m=B[:,3]==c; groups=[]; drop=0
            for key in sorted(set(zip(B[m,0],B[m,1]))):
                d=B[m&(B[:,0]==key[0])&(B[:,1]==key[1]),4]
                if len(d)>=3: groups.append(d.astype(float))
                else: drop+=1
            nmin=XMIN[c]; obs1=T1(np.concatenate(groups),nmin,np.random.default_rng([6,gi,c])); obs2=T2(groups,nmin)
            ss=np.random.SeedSequence([6,gi,c]).spawn(500)
            with Pool(4) as p: nul=np.array(p.map(job,[(g,c,groups,s) for s in ss]))
            for k,(lab,o) in enumerate([('T1 pooled dL(K2-K1)',obs1),('T2 within-window overdisp',obs2)]):
                v=nul[:,k]; pv=(1+np.sum(v>=o))/501
                extra=f"  share explained by heterogeneity {v.mean()/o:.2f}" if k==0 else ""
                print(f"{g:6s} windows {len(groups)} (dropped {drop}), bouts {sum(map(len,groups))} | {lab}: obs {o:9.3f}  null mean {v.mean():9.3f} [95%: {np.percentile(v,2.5):.3f}, {np.percentile(v,97.5):.3f}]  p {pv:.4f}{extra}",flush=True)
