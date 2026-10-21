import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from synth_test import P3, sample, fit, score, coarse_summary, np
from multiprocessing import Pool
T=coarse_summary(P3)
def job(seed):
    S=sample(P3,seed); fits=[fit(S,3,r) for r in range(6)]; best=max(fits,key=lambda m:score(m,S)); cs=coarse_summary(best)
    return cs['Xn_c']-T['Xn_c'], cs['Yn_c']-T['Yn_c'], *[a/b for a,b in zip(cs['dwell'],T['dwell'])]
if __name__=='__main__':
    with Pool(4) as p: R=np.array(p.map(job,range(101,113)))
    for k,n in enumerate(['Xn_c error','Yn_c error','dwell ratio S','dwell ratio T','dwell ratio W']):
        v=R[:,k]; print(f"{n:15s} mean {v.mean():+.3f}  SD {v.std(ddof=1):.3f}  95% CI of mean [{v.mean()-1.96*v.std(ddof=1)/np.sqrt(len(v)):+.3f}, {v.mean()+1.96*v.std(ddof=1)/np.sqrt(len(v)):+.3f}]  values {np.round(v,2).tolist()}")
    ok=all(abs(R[:,k].mean())<=0.10 and abs(R[:,k].mean())<=1.96*R[:,k].std(ddof=1)/np.sqrt(len(R)) for k in (0,1))
    print('REPLICATE VALIDATION', 'PASS' if ok else 'FAIL')
