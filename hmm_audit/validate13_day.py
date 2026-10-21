# PRESPEC 13 amendment 2: gates V5, V6 for the day-covariate model. Checkpointed.
import os,sys,json; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
import numpy as np
from hhmm_day import *
from hhmm import coarse_summary
from synth_test import P3, P2, truth_model
from validate13 import perturb
from multiprocessing import Pool
NF,ND,L=25,4,14400; KS=np.exp([-.4,-.13,.13,.4])
def make(tm,seed,het):
    rs=np.random.RandomState(seed); S=[]; days=[]
    for f in range(NF):
        base=perturb(tm,rs) if het else tm; eps,A,G,Bin=base._struct
        for d in range(ND):
            e=eps.copy(); e[:2]*=KS[d]; m=truth_model(tm.Nc,e,A,G,Bin,base.emissionprob_)
            S.append(m.sample(L,random_state=rs)[0].ravel()); days.append(d)
    return S,days
def job(a):
    kind,seed=a; out={'kind':kind,'seed':seed}; T=coarse_summary(P3)
    if kind=='V5':
        S,days=make(P3,seed,False); m=fit_pooled_day(S,days,ND,3,seed); cs=coarse_summary_day(m)
        out.update(dXn=cs['Xn_c']-T['Xn_c'],dYn=cs['Yn_c']-T['Yn_c'],kappaS_ratio=(np.array(cs['kappa'][0])/KS).tolist())
    else:
        S,days=make(P2,seed,True); rs=np.random.default_rng(seed); ho=set(rs.choice(NF,NF//5,replace=False).tolist())
        tr=[i for i in range(len(S)) if i//ND not in ho]; te=[i for i in range(len(S)) if i//ND in ho]
        m2=fit_pooled_day([S[i] for i in tr],[days[i] for i in tr],ND,2,seed); m3=fit_pooled_day([S[i] for i in tr],[days[i] for i in tr],ND,3,seed)
        g=[]
        for f in sorted(ho):
            idx=[i for i in te if i//ND==f]; s=[S[i] for i in idx]; dd=[days[i] for i in idx]
            g.append((score_day(m3,s,dd)-score_day(m2,s,dd))/sum(len(x) for x in s)*3600)
        out.update(gain_median=float(np.median(g)),gain_all=g)
    return out
if __name__=='__main__':
    path='validate13_day.jsonl'; done=set()
    if os.path.exists(path): done={(json.loads(l)['kind'],json.loads(l)['seed']) for l in open(path)}
    jobs=[j for j in [('V6',601),('V6',602),('V5',501),('V5',502),('V5',503)] if j not in done]; print('remaining',len(jobs),flush=True)
    with Pool(4) as p:
        for r in p.imap_unordered(job,jobs): open(path,'a').write(json.dumps(r,default=float)+'\n'); print('done',r['kind'],r['seed'],flush=True)
    R=[json.loads(l) for l in open(path)]; v5=[r for r in R if r['kind']=='V5']; ok5=True
    for k in ('dXn','dYn'):
        v=np.array([r[k] for r in v5]); print(f"V5 {k}: mean {v.mean():+.3f} SD {v.std(ddof=1):.3f} values {np.round(v,3).tolist()}"); ok5&=abs(v.mean())<=0.10 and v.std(ddof=1)<=0.15
    kr=np.array([r['kappaS_ratio'] for r in v5]).mean(0); print('V5 mean fitted/true kappa_S by day',np.round(kr,3).tolist()); ok5&=all(0.75<=x<=1.25 for x in kr)
    v6=[r for r in R if r['kind']=='V6']; ok6=all(r['gain_median']<=5 for r in v6); print('V6 median held-out gain per fly:',[round(r['gain_median'],2) for r in v6])
    print('V5',ok5,'V6',ok6,'=> DAY-MODEL VALIDATION','PASS' if ok5 and ok6 else 'FAIL')
