# PRESPEC 20e section 1: no-change null for across-day sign flips.
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw, sqw2
from sqw2_common import *
def null_one(d):
    g=d['g']; nph=sqw2.nP(g); T=d['T']; ph=np.zeros(T,int) if g=='per0DD' else phases(T); h=((T//1440)//2)*1440
    rng=np.random.default_rng(11000+1000*sqw.GROUPS.index(g)+d['fly']); P=np.array(d['P3']); out=[]
    for r in range(8):
        _,y=simulate2(np.array(d['pi3']),P[ph],np.array(d['B3']),rng); halves=[]
        for yy,pp in ((y[:h],ph[:h]),(y[h:],ph[h:])):
            ll,pi,Pf,Bf,_=em(yy,pp,nph,3,10,r+100*d['fly']); Pf,Bf,_=label3b(Pf,Bf,pi)
            halves.append([float((yn(Pf[p])-xn(Pf[p]))/2) for p in range(nph)])
        out.append(halves)
    return f"{g}|{d['fly']}",out
if __name__=='__main__':
    R=sqw2.load(0.5)
    with Pool(4) as pl: res=dict(pl.map(null_one,[d for d in R.values() if not d['deg']],chunksize=1))
    json.dump(res,open(sqw.OUTD+'e_null.json','w'))
