# PRESPEC 20b section 4: synthetic validation V1b, V2b, V3b (20 flies x 5760 min each).
import numpy as np, json, time
from multiprocessing import Pool
from sqw2_common import *
from sqw_common import bic
from sqw_synth import P1 as PL1, B1
T=5760; NF=20; ST=20
PD1=np.array([[0.67,0.165,0.165],[0.04345,0.93,0.02655],[0.03,0.01,0.96]])   # dark: Xn +0.2, Yn -0.1
TRUTH={'light':(xn(PL1),yn(PL1)),'dark':(xn(PD1),yn(PD1))}
PL2=np.array([[0.93,0.07],[0.04,0.96]]); PD2=np.array([[0.95,0.05],[0.03,0.97]]); B2=np.array([[0.30,0.70,0.00],[0.03,0.00,0.97]])
def drift(P0,s):   # exit rates scaled by 2/3..4/3 (factor 2) across the phase, s in [0,1)
    f=2/3+(2/3)*s; P=P0*f; np.fill_diagonal(P,0); np.fill_diagonal(P,1-P.sum(1)); return P
def make(case,rng):
    ph=phases(T)
    if case=='V1b': Pt=np.where(ph[:,None,None]==0,PL1,PD1); B=B1; pi=stationary(PL1)
    else:
        B=B2; pi=stationary(PL2)
        if case=='V2b': Pt=np.where(ph[:,None,None]==0,PL2,PD2)
        else:
            s=(np.arange(T)%720)/720; Pt=np.stack([drift(PL2 if p==0 else PD2,si) for p,si in zip(ph,s)])
    return simulate2(pi,Pt,B,rng)[1],ph
def one(a):
    case,f=a; rng=np.random.default_rng(3000+100*['V1b','V2b','V3b'].index(case)+f); y,ph=make(case,rng)
    F={K:em(y,ph,2,K,ST,20+f) for K in (2,3,4)}; ne=max(2,int(np.sum(np.diff(y)!=0)))
    r=dict(case=case,f=f,nchg=ne,bic={K:-2*F[K][0]+kpar2(K,2)*np.log(T) for K in F},bic_eff={K:-2*F[K][0]+kpar2(K,2)*np.log(ne) for K in F},
           iters={K:F[K][4] for K in F})
    P,B,pi=label3b(F[3][2],F[3][3],F[3][1])
    r.update(deg=bool(degenerate(B)),xl=float(xn(P[0])),yl=float(yn(P[0])),xd=float(xn(P[1])),yd=float(yn(P[1])),B=B.tolist())
    return r
if __name__=='__main__':
    t0=time.time(); jobs=[(c,f) for c in ('V1b','V2b','V3b') for f in range(NF)]
    with Pool(4) as p: res=p.map(one,jobs,chunksize=1)
    json.dump(res,open('sqw2_synth.json','w'),default=lambda o:o.item() if hasattr(o,'item') else str(o))
    L=[f"truth light Xn {TRUTH['light'][0]:+.3f} Yn {TRUTH['light'][1]:+.3f}; dark Xn {TRUTH['dark'][0]:+.3f} Yn {TRUTH['dark'][1]:+.3f}"]; fr={}
    for c in ('V1b','V2b','V3b'):
        R=[r for r in res if r['case']==c]; best=[min(r['bic'],key=lambda k:r['bic'][k]) for r in R]; beste=[min(r['bic_eff'],key=lambda k:r['bic_eff'][k]) for r in R]
        fr[c]={K:np.mean([b==K for b in best]) for K in (2,3,4)}; fe={K:np.mean([b==K for b in beste]) for K in (2,3,4)}
        L.append(f"{c}: best K (n=T) {fr[c]}; (n_eff) {fe}; disagree {np.mean([a!=b for a,b in zip(best,beste)]):.2f}; H1 degenerate {np.mean([r['deg'] for r in R]):.2f}; median changes {np.median([r['nchg'] for r in R])}; max EM iters {max(max(r['iters'].values()) for r in R)}")
        q={k:np.array([r[k] for r in R]) for k in ('xl','yl','xd','yd')}
        L.append('   '+'; '.join(f"{k} mean {np.nanmean(v):+.3f} SD {np.nanstd(v):.3f}" for k,v in q.items()))
        if c=='V1b':
            tv={'xl':TRUTH['light'][0],'yl':TRUTH['light'][1],'xd':TRUTH['dark'][0],'yd':TRUTH['dark'][1]}
            ok=fr[c][3]>=0.7 and all(abs(np.nanmean(q[k])-tv[k])<0.1 and np.nanstd(q[k])<0.3 for k in q)
            L.append(f"   V1b PASS={ok}")
    L.append(f"   V2b PASS={fr['V2b'][2]>=0.8}")
    L.append(f"   V3b PASS={fr['V3b'][2]>=0.8 and fr['V3b'][3]-fr['V2b'][3]<=0.10}")
    L.append(f"elapsed {time.time()-t0:.0f} s"); open('sqw2_synth.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
