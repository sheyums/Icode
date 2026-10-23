# Design check after the PRESPEC 20 V1 failure (synthetic only, no fly data): steps per fit needed
# for BIC to choose K=3 and for routing to be estimable, under the V1 truth.
import numpy as np, json
from multiprocessing import Pool
from sqw_common import *
from sqw_synth import P1,B1
NS=[240,960,2880,5760]; REP=30
def one(a):
    n,r=a; rng=np.random.default_rng(7000+n+r); _,y=simulate(stationary(P1),P1,B1,n,rng)
    F={K:fit(y,K,10,r) for K in (2,3,4)}; b={K:bic(F[K][0],K,n) for K in F}
    P,B,_=label3(F[3][2],F[3][3]); return dict(n=n,best=min(b,key=b.get),xn=float(xn(P)),yn=float(yn(P)),deg=bool(degenerate(B)))
if __name__=='__main__':
    with Pool(4) as p: res=p.map(one,[(n,r) for n in NS for r in range(REP)],chunksize=1)
    L=[]
    for n in NS:
        R=[x for x in res if x['n']==n]; x=np.array([r['xn'] for r in R]); y=np.array([r['yn'] for r in R])
        L.append(f"n {n:5d} ({n/240:.0f} windows): K=3 chosen {np.mean([r['best']==3 for r in R]):.2f}, K=2 {np.mean([r['best']==2 for r in R]):.2f}, K=4 {np.mean([r['best']==4 for r in R]):.2f}; "
                 f"Xn mean {np.nanmean(x):+.3f} SD {np.nanstd(x):.3f}; Yn mean {np.nanmean(y):+.3f} SD {np.nanstd(y):.3f}; degenerate {np.mean([r['deg'] for r in R]):.2f}")
    open('sqw_power.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
