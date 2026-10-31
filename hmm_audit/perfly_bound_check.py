# POST HOC, descriptive (user question 2026-10-31): per fly-phase comparison of the HMM's exact EPR with the observed
# chain's exact EPR (sigma_MC) and with the triplet bound sigma_MC + E (E = that fly's 20p excess over its own null A).
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
from sqw13 import mc_fit
from sqw15 import seqs, DD
from scipy.stats import wilcoxon
O=json.load(open(sqw.OUTD+'p_blockkl.json')); F=sqw2.load(0.5)
def one(d):
    g,i=d['g'],d['fly']; y,ph,_=seqs(g,i); nP=1 if g in DD else 2; _,_,Pm=mc_fit(y,ph,nP)
    Ph=np.array(d['hmm']['P']) if d['hmm'] else np.array(F[(g,i)]['P3']); out=[]
    for p in range(nP):
        key='whole' if g in DD else ['light','dark'][p]
        if d['ref'][key] is None or d['deg']: out.append(None); continue
        s=np.array([x[p]['D3'] for x in d['sims']['A']]); E=d['obs'][p]['D3']-np.median(s)
        out.append(dict(mc=float(R.model_energetics(Pm[p])['EPR']),hmm=float(R.model_energetics(Ph[p])['EPR']),E=float(E)))
    return dict(g=g,fly=i,v=out)
if __name__=='__main__':
    with Pool(4) as pl: res=pl.map(one,O,chunksize=4)
    json.dump(res,open(sqw.OUTD+'perfly_bound.json','w'))
    print('per fly-phase; EPR in k_B/h; bound = sigma_MC + E (E from 20p, that fly vs its own null). Non-degenerate flies with a Q/W reference.')
    for g,p,lab in [('CS',0,'light'),('CS',1,'dark'),('CSDD',0,'whole'),('per0LD',0,'light'),('per0LD',1,'dark'),('per0DD',0,'whole'),('fumin',0,'light'),('fumin',1,'dark'),('sss',0,'light'),('sss',1,'dark')]:
        V=[r['v'][p] for r in res if r['g']==g and r['v'][p] is not None]
        mc=np.array([v['mc'] for v in V]); hm=np.array([v['hmm'] for v in V]); b=mc+np.array([v['E'] for v in V])
        fin=np.isfinite(hm)&np.isfinite(mc); inf=np.mean(~np.isfinite(hm))
        ok=hm[fin]>=b[fin]; rat=hm[fin]/np.maximum(b[fin],1e-9)
        pw=wilcoxon(hm[fin]-mc[fin]).pvalue if fin.sum()>5 else np.nan
        print(f"{g:6s} {lab:5s} n {len(V):2d} (HMM EPR infinite {inf:.2f}) | sigma_MC median {np.median(mc[fin]):.3f} | bound {np.median(b[fin]):.3f} | HMM {np.median(hm[fin]):.3f} | "
              f"per fly: HMM > sigma_MC in {np.mean(hm[fin]>mc[fin]):.2f} (paired Wilcoxon p {pw:.2g}), HMM >= bound in {np.mean(ok):.2f}, median HMM/bound {np.median(rat):.2f} [IQR {np.percentile(rat,25):.2f}-{np.percentile(rat,75):.2f}]")
