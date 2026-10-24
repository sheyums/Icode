# PRESPEC 20c step 3 summary and PRESPEC 20e amendment B group-level tests (6 h blocks, theta 0.5).
import json,numpy as np, sqw_ref as R
from scipy.stats import wilcoxon, friedmanchisquare
O='sqw_out/b6/'; F={(d['g'],d['fly']):d for d in map(json.loads,open(O+'fit2_0.5.jsonl'))}; RT=json.load(open(O+'route2_0.5.json'))
usable={l.split(':')[0]:('USABLE=True' in l) for l in open('sqw3_realistic.log') if 'USABLE=' in l}
rng=np.random.default_rng(20); L=[]; BL=['ZT0-6','ZT6-12','ZT12-18','ZT18-24']
def ci(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]; return v.mean(),*np.percentile([rng.choice(v,len(v)).mean() for _ in range(2000)],[2.5,97.5])
def holm(p):
    p=np.asarray(p); o=np.argsort(p); m=len(p); adj=np.empty(m); run=0
    for r,i in enumerate(o): run=max(run,(m-r)*p[i]); adj[i]=min(1,run)
    return adj
for g in ['CS','per0LD','fumin']:
    D=[F[k] for k in sorted(F) if k[0]==g]; best=[]
    for d in D:
        b={int(K):R.bic(v['ll'],int(K),4,d['T']) for K,v in d['F'].items()}; best.append(min(b,key=b.get))
    nd=[d for d in D if not d['deg']]
    L.append(f"{g} ({'USABLE' if usable[g] else 'NOT usable: descriptive'}): flies {len(D)} | best K "+' '.join(f"{K}:{np.mean([x==K for x in best]):.2f}" for K in (2,3,4))+f" | H1 degenerate {1-len(nd)/len(D):.2f}")
    C=[];Jc=[];nocyc=0
    for d in nd:
        P=np.array(d['P3']); e=[R.model_energetics(P[b]) for b in range(4)]
        if any(abs(x['J'])<1e-12 for x in e): nocyc+=1; continue
        C.append([R.c_pred(P[b]) for b in range(4)])
        m=np.array(RT[f"{g}|{d['fly']}"]['pairs_ffbs']); h=RT[f"{g}|{d['fly']}"]['hours']; Jc.append([R.j_cycle_count(m[b],h[b]) for b in range(4)])
    C=np.array(C); Jc=np.array(Jc)
    L.append(f"   flies used {len(C)} (no-cycle excluded {nocyc})")
    for b in range(4):
        L.append(f"   {BL[b]:8s} c_pred {ci(C[:,b])[0]:+.3f} [{ci(C[:,b])[1]:+.3f},{ci(C[:,b])[2]:+.3f}] (c>0 in {np.mean(C[:,b]>0):.2f}) | count J {ci(Jc[:,b])[0]:+.3f} [{ci(Jc[:,b])[1]:+.3f},{ci(Jc[:,b])[2]:+.3f}] /h")
    p1=[wilcoxon(C[:,b]).pvalue for b in range(4)]; a1=holm(p1)
    L.append("   one-sample Wilcoxon c vs 0 (Holm): "+' | '.join(f"{BL[b]} p {p1[b]:.2g} adj {a1[b]:.2g}" for b in range(4)))
    fr=friedmanchisquare(*[C[:,b] for b in range(4)])
    pairs=[(0,1),(1,2),(2,3),(3,0)]; p2=[wilcoxon(C[:,j]-C[:,i]).pvalue for i,j in pairs]; a2=holm(p2)
    L.append(f"   Friedman chi2 {fr.statistic:.2f} p {fr.pvalue:.2g} | adjacent within-fly Wilcoxon (Holm): "+' | '.join(f"{BL[i]}->{BL[j]} d {np.median(C[:,j]-C[:,i]):+.3f} p {p:.2g} adj {a:.2g}" for (i,j),p,a in zip(pairs,p2,a2)))
    sig=[b for b in range(4) if a1[b]<0.05]; claim=False
    for i in sig:
        for j in sig:
            if i<j and np.sign(np.median(C[:,i]))!=np.sign(np.median(C[:,j])):
                pij=wilcoxon(C[:,j]-C[:,i]).pvalue; claim|=pij*len(pairs)<0.05
    L.append(f"   'group direction changes across the day' (amendment B rule): {claim}")
open('sqw3_report.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
