# PRESPEC 20s-A (Plan A): observable-level quantities per fly x phase, with the HMM bracket. Writes s_observable.json.
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
from sqw13 import mc_fit
from sqw15 import seqs, DD
G=['CS','per0LD','per0DD','fumin','sss','CSDD']
def routes_obs(y,ph,nph):
    """Route counts (ROUTES order) on the observed symbols, each run of non-missing minutes separately."""
    y=np.asarray(y); c=np.zeros((nph,4)); ok=np.r_[False,y>=0,False]; d=np.diff(ok.astype(int)); st=np.flatnonzero(d==1); en=np.flatnonzero(d==-1)
    for s,e in zip(st,en):
        for r,t in R.route_events(y[s:e]): c[ph[s+t],R.ROUTES.index(r)]+=1
    return c
def block_stats(y,ph,nph,Pm,refok):
    T=len(y); N=R.observable_pair_counts(y,ph,nph); c=R.c_obs_minutes(y,ph,nph); rc=routes_obs(y,ph,nph); X,Y=R.xy_counts(rc); out=[]
    for p in range(nph):
        h=np.sum(ph==p)/60; e=R.model_energetics(Pm[p]) if Pm is not None else None; ls=R.landscape(Pm[p]) if Pm is not None else None
        out.append(dict(ok=bool(refok[p]),hours=float(h),c_obs=float(c[p]),J=float(R.j_cycle_count(N[p],h)),EPR=float(R.epr_count_bound(N[p],h)[0]),
                        r3=bool(R.epr_count_bound(N[p],h)[1]),N=N[p].tolist(),
                        J_chain=None if e is None else float(e['J']),A_chain=None if e is None else float(e['A']),EPR_chain=None if e is None else float(e['EPR']),
                        E21=None if ls is None else float(ls['E2_minus_E1']),E31=None if ls is None else float(ls['E3_minus_E1']),
                        F_SW=None if ls is None else float(ls['F_SW']),F_WS=None if ls is None else float(ls['F_WS']),Xn_cnt=float(X[p]),Yn_cnt=float(Y[p])))
    return out
def fly(a):
    g,i,hid=a; y,ph,ref=seqs(g,i); nP=1 if g in DD else 2; keys=['whole'] if g in DD else ['light','dark']
    refok=[ref[k] is not None and np.isfinite(ref[k]) for k in keys]
    _,_,Pm=mc_fit(y,ph,nP); main=block_stats(y,ph,nP,Pm,refok)
    six=None
    if g=='CS':
        ph6=R.phases(len(y),'6h'); _,_,P6=mc_fit(y,ph6,4); six=block_stats(y,ph6,4,P6,[refok[0],refok[0],refok[1],refok[1]])
    for p,m in enumerate(main): m.update(hid[p] if hid else {})
    return dict(g=g,fly=i,main=main,six=six)
if __name__=='__main__':
    F=sqw2.load(0.5); RT=json.load(open(sqw.OUTD+'route2_0.5.json')); jobs=[]
    for g in G:
        for i in range(len(sqw.files(g))):
            hid=None
            if g!='CSDD':
                d=F[(g,i)]; rt=RT[f'{g}|{i}']; P=np.array(d['P3']); hid=[]
                for p in range(len(rt['hours'])):
                    pp=np.array(rt['pairs_path'])[:,p]; h=rt['hours'][p]; e=R.model_energetics(P[p])
                    hid.append(dict(hmm_deg=bool(d['deg']),J_hmm=float(np.mean([R.j_cycle_count(n,h) for n in pp])),
                                    EPR_hmm=float(np.mean([R.epr_count_bound(n,h)[0] for n in pp])),EPR_hmm_model=float(e['EPR']),J_hmm_model=float(e['J'])))
            jobs.append((g,i,hid))
    with Pool(4) as pl: out=pl.map(fly,jobs,chunksize=2)
    json.dump(out,open(sqw.OUTD+'s_observable.json','w')); print('done',len(out))
