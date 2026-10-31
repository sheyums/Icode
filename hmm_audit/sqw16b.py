# PRESPEC 20q amendment E: rule of three for fly-phases whose observed counts have a one-way edge (cycle present).
# Writes q_tur_r3.json = q_tur.json with those fly-phases' chain quantities and chain nulls recomputed from P_r3.
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw
from sqw13 import mc_fit
from sqw15 import seqs, gidx, sim_rows, DD
from sqw16 import estimate, exact_B, exact_B_ct, DSGN, NMC
def counts(y,ph,p):
    C=np.zeros((3,3)); t=np.flatnonzero((ph[1:]==p)&(y[1:]>=0)&(y[:-1]>=0))+1
    np.add.at(C,(y[t-1],y[t]),1); return C
def kind(C,y,ph,p):
    if len({int(v) for v in y[(ph==p)&(y>=0)]})<3: return 'iii'
    off=[(a,b) for a in range(3) for b in range(a+1,3)]
    if any(C[a,b]==0 and C[b,a]==0 for a,b in off): return 'iii'
    if any((C[a,b]>0)!=(C[b,a]>0) for a,b in off): return 'i'
    return 'ii'
def r3(C):
    C=C.copy()
    for a in range(3):
        for b in range(a+1,3):
            if C[a,b]==0 and C[b,a]>0: C[a,b]=min(3.0,C[b,a])
            elif C[b,a]==0 and C[a,b]>0: C[b,a]=min(3.0,C[a,b])
    return C/C.sum(1,keepdims=True)
def redo(a):
    d,p=a; g,i=d['g'],d['fly']; gi=gidx(g); y,ph,_=seqs(g,i); nP=1 if g in DD else 2; miss=y==-1
    _,pm,Pm=mc_fit(y,ph,nP); Pm=sim_rows(Pm); Pr=r3(counts(y,ph,p)); Pm[p]=Pr
    rng=np.random.default_rng(16500+1000*gi+10*i+p); sims=[]
    for _ in range(NMC):
        _,ys=R.simulate(pm,Pm,np.eye(3),ph,rng); ys[miss]=-1; sims.append(estimate(ys,ph,nP)[p])
    Jm,Dm,Bm=exact_B(Pr,DSGN)
    return dict(g=g,fly=i,p=p,sig_mc=float(R.model_energetics(Pr)['EPR']),J_mc=Jm,D_mc=Dm,B_mc=Bm,B_mc_ct=exact_B_ct(Pr,DSGN)[2],sims=sims)
if __name__=='__main__':
    O=json.load(open(sqw.OUTD+'q_tur.json')); jobs=[]; K={}
    for d in O:
        g=d['g']; nP=1 if g in DD else 2; y,ph,_=seqs(g,d['fly'])
        for p in range(nP):
            key='whole' if g in DD else ['light','dark'][p]; e=d['ex'][p]
            if d['ref'][key] is None or (np.isfinite(e['sig_mc']) and e['sig_mc']>1e-9 and np.isfinite(e['B_mc_ct']) and e['B_mc_ct']>1e-12): continue
            k=kind(counts(y,ph,p),y,ph,p); K.setdefault((g,key),{}); K[(g,key)][k]=K[(g,key)].get(k,0)+1
            if k=='i': jobs.append((d,p))
    print('kinds of excluded fly-phases:',{f'{a} {b}':v for (a,b),v in K.items()})
    with Pool(4) as pl: res=pl.map(redo,jobs,chunksize=1)
    for r in res:
        d=[x for x in O if x['g']==r['g'] and x['fly']==r['fly']][0]; p=r['p']
        d['ex'][p].update({k:r[k] for k in ('sig_mc','J_mc','D_mc','B_mc','B_mc_ct')}); d['ex'][p]['r3']=True
        for s,new in zip(d['simMC'],r['sims']): s[p]=new
    json.dump(dict(data=O,kinds={f'{a}|{b}':v for (a,b),v in K.items()}),open(sqw.OUTD+'q_tur_r3.json','w')); print('done',len(res))
