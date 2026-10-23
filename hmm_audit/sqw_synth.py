# PRESPEC 20 section 7: synthetic validation (V1 positive, V2 negative control).
import numpy as np, json, sys, time
from multiprocessing import Pool
from sqw_common import *
NW=50; NBOOT=50; STARTS=20
# V1 truth: Xn=-0.3, Yn=+0.4 (minute steps), distinct H1 emissions (order H1,H2,H3; symbols Q,W,S)
P1=np.array([[0.67,0.165,0.165],[0.05765,0.93,0.01235],[0.02074,0.01926,0.96]])
B1=np.array([[0.80,0.15,0.05],[0.15,0.85,0.00],[0.03,0.00,0.97]])
# V2 truth: 2 states (wake, sleep); wake emits Q at a fixed rate, no H1
P2=np.array([[0.93,0.07],[0.04,0.96]]); B2=np.array([[0.30,0.70,0.00],[0.03,0.00,0.97]])
def one(args):
    case,w=args; rng=np.random.default_rng(2000+100*(case=='V2')+w)
    P,B=(P1,B1) if case=='V1' else (P2,B2); pi=stationary(P)
    _,y=simulate(pi,P,B,NSTEP,rng)
    fits={K:fit(y,K,STARTS,20+w) for K in (2,3,4)}
    ne=nchanges(y)
    r=dict(case=case,w=w,nchg=int(np.sum(np.diff(y)!=0)),
           bic={K:bic(fits[K][0],K,NSTEP) for K in fits},bic_eff={K:bic(fits[K][0],K,ne) for K in fits})
    _,pi3,P3,Bf=fits[3]; P3,Bf,pi3=label3(P3,Bf,pi3)
    r.update(xn=float(xn(P3)),yn=float(yn(P3)),deg=bool(degenerate(Bf)),B3=Bf.tolist())
    if case=='V2':   # parametric bootstrap of the K=3 routing
        xs=[];ys=[]
        for b in range(NBOOT):
            _,yb=simulate(pi3,P3,Bf,NSTEP,rng); f=fit(yb,3,STARTS,9000+w*100+b)
            if f is None: continue
            Pb,Bb,_=label3(f[2],f[3]); xs.append(xn(Pb)); ys.append(yn(Pb))
        r.update(sdx=float(np.nanstd(xs)),sdy=float(np.nanstd(ys)))
    return r
if __name__=='__main__':
    t0=time.time(); jobs=[('V1',w) for w in range(NW)]+[('V2',w) for w in range(NW)]
    with Pool(4) as p: res=p.map(one,jobs,chunksize=1)
    json.dump(res,open('sqw_synth.json','w'),default=lambda o:o.item() if hasattr(o,'item') else str(o))
    out=[]
    for case in ('V1','V2'):
        R=[r for r in res if r['case']==case]
        best=[min(r['bic'],key=r['bic'].get) for r in R]; beste=[min(r['bic_eff'],key=r['bic_eff'].get) for r in R]
        fr={K:np.mean([b==K for b in best]) for K in (2,3,4)}; fre={K:np.mean([b==K for b in beste]) for K in (2,3,4)}
        x=np.array([r['xn'] for r in R]); y=np.array([r['yn'] for r in R])
        out.append(f"{case}: BIC(n=240) best K fractions {fr}; BIC(n_eff) {fre}; K differs in {np.mean([a!=b for a,b in zip(best,beste)]):.2f}")
        out.append(f"  median symbol changes {np.median([r['nchg'] for r in R])}; H1 degenerate in {np.mean([r['deg'] for r in R]):.2f}")
        out.append(f"  K=3 Xn mean {np.nanmean(x):+.3f} median {np.nanmedian(x):+.3f} SD {np.nanstd(x):.3f} (nan {np.isnan(x).sum()});"
                   f" Yn mean {np.nanmean(y):+.3f} median {np.nanmedian(y):+.3f} SD {np.nanstd(y):.3f}")
        if case=='V1':
            ok=fr[3]>=0.7 and abs(np.nanmean(x)+0.3)<0.1 and abs(np.nanmean(y)-0.4)<0.1
            out.append(f"  V1 PASS={ok} (need K=3 >= 0.70, |bias| < 0.1; truth Xn -0.3, Yn +0.4)")
        else:
            uninf=np.mean([r['deg'] or max(r['sdx'],r['sdy'])>=0.2 for r in R])
            ok=fr[2]>=0.8 and uninf>=0.8
            out.append(f"  bootstrap SD median x {np.median([r['sdx'] for r in R]):.3f} y {np.median([r['sdy'] for r in R]):.3f}; degenerate-or-SD>=0.2 in {uninf:.2f}")
            out.append(f"  V2 PASS={ok} (need K=2 >= 0.80 and uninformative K=3 routing >= 0.80)")
    out.append(f"elapsed {time.time()-t0:.0f} s")
    open('sqw_synth.log','w').write('\n'.join(out)+'\n'); print('\n'.join(out))
