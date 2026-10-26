# PRESPEC 20k: route vs bout duration, observed (permutation null) and hidden (Markov-model simulation null).
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
REL=list(R.RELATIONS); NSIM=400; NP_REAL=200; NP_SIM=20
def fly(d):
    g,i=d['g'],d['fly']; ld=g!='per0DD'; nph=2 if ld else 1
    x=R.read_counts(sqw.files(g)[i]); y,_=R.sqw_symbols(x,0.5,light_dark=ld); ph=R.phases(len(y),'12h' if ld else 'none')
    pi,P,B=(np.array(d[k]) for k in ('pi3','P3','B3')); rng=np.random.default_rng(10001+1000*sqw.GROUPS.index(g)+i)
    obs=R.route_duration_observed(x,0.5,ld)
    Z=R.ffbs(pi,P,B,y,ph,NP_REAL,int(rng.integers(1<<30)))
    def hidD(paths):
        out={}
        for rel in REL:
            for p in range(nph):
                v=[R.route_duration_D(R.route_duration_hidden(z,ph),rel,p) for z in paths]; v=np.array(v,float)
                out[f'{rel}|{p}']=float(np.nanmean(v)) if np.isfinite(v).mean()>=0.5 else None
        return out
    hid=hidD(Z); hidrec=[R.route_duration_hidden(z,ph) for z in Z[:20]]
    sims=[]
    for r in range(NSIM):
        _,ys=R.simulate(pi,P,B,ph,rng); Zs=R.ffbs(pi,P,B,ys,ph,NP_SIM,int(rng.integers(1<<30))); sims.append(hidD(Zs))
    return dict(g=g,fly=i,nph=nph,obs=obs,hid=hid,sims=sims,hidrec=hidrec)
if __name__=='__main__':
    F=sqw2.load(0.5)
    with Pool(4) as pl: res=pl.map(fly,[d for d in F.values() if not d['deg']],chunksize=1)
    json.dump(res,open(sqw.OUTD+'k_route_duration.json','w'),default=lambda o:bool(o) if isinstance(o,np.bool_) else float(o))
    print('done',len(res))
