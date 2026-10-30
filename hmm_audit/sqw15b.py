# PRESPEC 20p amendment A: DD genotypes (per0DD, CS DD) analysed in 12-h clock halves; nulls A2 (two 12-h matrices), B (6-h).
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw
from sqw13 import mc_fit
from sqw15 import seqs, stats, gidx, NA, NB, I3, sim_rows
def fly(a):
    g,i=a; gi=gidx(g); y,_,_=seqs(g,i); miss=y==-1; ph2=R.phases(len(y),'12h'); ph6=R.phases(len(y),'6h')
    _,pa,Pa=mc_fit(y,ph2,2); _,pb,Pb=mc_fit(y,ph6,4); Pa,Pb=sim_rows(Pa),sim_rows(Pb); rng=np.random.default_rng(15500+1000*gi+i); sims={}
    for name,(pi,P,pp,ns) in (('A2',(pa,Pa,ph2,NA)),('B',(pb,Pb,ph6,NB))):
        L=[]
        for _ in range(ns):
            _,ys=R.simulate(pi,P,I3,pp,rng); ys[miss]=-1; L.append(stats(ys,ph2,2))
        sims[name]=L
    return dict(g=g,fly=i,obs=stats(y,ph2,2),sims=sims)
if __name__=='__main__':
    jobs=[(g,i) for g in ('per0DD','CSDD') for i in range(len(sqw.files(g)))]
    with Pool(4) as pl: out=pl.map(fly,jobs,chunksize=1)
    json.dump(out,open(sqw.OUTD+'p_blockkl_ddhalves.json','w')); print('done',len(out))
