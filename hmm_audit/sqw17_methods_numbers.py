# 20s-A, numbers the Methods text needs at the observable level (same procedures as 20j/onsager_L sensitivity lines):
# (1) per-fly Spearman of count J vs exact observed-chain J; (2) Onsager L through the origin vs with a free intercept.
import os,sys,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
from scipy.stats import spearmanr
O=json.load(open(sqw.OUTD+'s_observable.json'))
BL=[('CS',0),('CS',1),('per0LD',0),('per0LD',1),('per0DD',0),('fumin',0),('fumin',1),('sss',0),('sss',1),('CSDD',0)]
rh=[]
for g,p in BL:
    R=[d['main'][p] for d in O if d['g']==g and d['main'][p]['ok'] and d['main'][p]['J_chain'] is not None and abs(d['main'][p]['J_chain'])>=1e-12]
    r=spearmanr([x['J'] for x in R],[x['J_chain'] for x in R])[0]; rh.append(r); print(f'{g:6s} phase {p}: Spearman(count J, chain J) {r:.2f} (n {len(R)})')
print(f'range {min(rh):.2f}-{max(rh):.2f}')
for g in ('CS','CSDD','per0LD','per0DD','fumin','sss'):
    pts=[d['main'][p] for d in O if d['g']==g for p in range(len(d['main'])) if d['main'][p]['ok'] and d['main'][p]['J_chain'] is not None and abs(d['main'][p]['J_chain'])>=1e-12 and np.isfinite(d['main'][p]['EPR'])]
    x=np.array([q['J'] for q in pts])**2; y=np.array([q['EPR'] for q in pts]); s0=(x@y)/(x@x)
    A=np.c_[np.ones_like(x),x]; b,s1=np.linalg.lstsq(A,y,rcond=None)[0]
    print(f'{g:6s} L origin {1/s0:.3f} | free intercept: slope change {100*(s1-s0)/s0:+.1f}%, intercept {b:+.3f} nats/h')
