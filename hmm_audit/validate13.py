# PRESPEC 13 Part V. Checkpointed per job to validate13.jsonl.
import os,sys,json; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from hhmm import *
from synth_test import P3, P2, truth_model
from multiprocessing import Pool
NF,ND,L=25,4,14400
def perturb(tm,rs):
    eps,A,G,Bin=tm._struct; E=tm.emissionprob_
    E2=0.9*E+0.1*rs.dirichlet(np.ones(4),E.shape[0]); return truth_model(tm.Nc,np.clip(eps*np.exp(rs.normal(0,.3,len(eps))),1e-6,.5),A,G,Bin,E2)
def make(tm,seed,het):
    rs=np.random.RandomState(seed); S=[]; models=[]
    for f in range(NF):
        m=perturb(tm,rs) if het else tm; models.append(m)
        for d in range(ND): S.append(m.sample(L,random_state=rs)[0].ravel())
    return S, models   # fly-major
def job(a):
    kind,seed=a; T=coarse_summary(P3); out={'kind':kind,'seed':seed}
    if kind in ('V1','V3'):
        S,models=make(P3,seed,kind=='V3'); m=fit_pooled(S,3,seed); cs=coarse_summary(m)
        llfit=score(m,S); lltrue=sum(models[i//ND].score(s.reshape(-1,1)) for i,s in enumerate(S))
        out.update(dLL=float(llfit-lltrue))
        out.update(dXn=cs['Xn_c']-T['Xn_c'],dYn=cs['Yn_c']-T['Yn_c'],dwell_ratio=[a/b for a,b in zip(cs['dwell'],T['dwell'])],mf=cs['mf'])
    else:
        S,_=make(P2,seed,True); tr=[s for i,s in enumerate(S) if i%ND!=ND-1]; te=[s for i,s in enumerate(S) if i%ND==ND-1]
        m2=fit_pooled(tr,2,seed); m3=fit_pooled(tr,3,seed)
        g=[(m3.score(s.reshape(-1,1))-m2.score(s.reshape(-1,1)))/len(s)*3600 for s in te]; out.update(gain_median=float(np.median(g)),gain_all=g)
    return out
if __name__=='__main__':
    path='validate13.jsonl'; done=set()
    if os.path.exists(path): done={(json.loads(l)['kind'],json.loads(l)['seed']) for l in open(path)}
    jobs=[j for j in [('V2',301),('V2',302)]+[('V1',s) for s in range(201,207)]+[('V3',401),('V3',402)] if j not in done]
    print('remaining',len(jobs),flush=True)
    with Pool(4) as p:
        for r in p.imap_unordered(job,jobs):
            open(path,'a').write(json.dumps(r,default=float)+'\n'); print('done',r['kind'],r['seed'],flush=True)
    R=[json.loads(l) for l in open(path)]
    v1=[r for r in R if r['kind']=='V1']; ok1=True
    for k in ('dXn','dYn'):
        v=np.array([r[k] for r in v1]); print(f"V1 {k}: mean {v.mean():+.3f} SD {v.std(ddof=1):.3f} values {np.round(v,3).tolist()}"); ok1&=abs(v.mean())<=0.10 and v.std(ddof=1)<=0.15
    dr=np.array([r['dwell_ratio'] for r in v1]).mean(0); print('V1 mean dwell ratio S,T,W',np.round(dr,3).tolist()); ok1&=all(0.8<=x<=1.25 for x in dr)
    v2=[r for r in R if r['kind']=='V2']; ok2=all(r['gain_median']<=5 for r in v2); print('V2 median held-out gain per fly (nats/h):',[round(r['gain_median'],2) for r in v2])
    for r in [r for r in R if r['kind'] in ('V1','V3')]: print(f"{r['kind']} seed {r['seed']}: dLL(fit - truth) {r['dLL']:+.1f} nats -> {'OPTIMIZER FAILURE' if r['dLL']<-1 else 'fit >= truth (error is statistical)'}")
    for r in [r for r in R if r['kind']=='V3']: print(f"V3 seed {r['seed']}: dXn {r['dXn']:+.3f} dYn {r['dYn']:+.3f} dwell ratio {np.round(r['dwell_ratio'],2).tolist()}")
    print('V1',ok1,'V2',ok2,'=> VALIDATION','PASS' if ok1 and ok2 else 'FAIL')
