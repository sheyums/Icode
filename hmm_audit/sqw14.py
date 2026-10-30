# PRESPEC 20o: memory (2nd-order chain) and drift (6-h rates) checks against the 20b HMM, CS and fumin.
# Usage: python3 sqw14.py check | run
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
from sqw13 import mc_fit, ppc_stats, NSIM
GS=('CS','fumin'); I3=np.eye(3)
B2=np.zeros((9,3)); B2[np.arange(9),np.arange(9)%3]=1          # pair state k = 3*y_{t-1} + y_t emits y_t
ALLOW=np.array([[k%3==l//3 for l in range(9)] for k in range(9)])   # (a,b) -> (b,c) only
def mc2_fit(y,ph,nP,tol=1e-10,maxit=1000):
    """2nd-order chain as a 9-state chain on pairs with fixed deterministic emissions; missing minutes marginalised."""
    y=np.asarray(y); P=np.where(ALLOW,1.0,0.0)[None].repeat(nP,0)
    for t in range(2,len(y)):
        if min(y[t-2],y[t-1],y[t])>=0: P[ph[t],3*y[t-2]+y[t-1],3*y[t-1]+y[t]]+=1
    P/=P.sum(2,keepdims=True); pi=np.full(9,1/9); old=-np.inf
    for it in range(maxit):
        pi,P,_,ll=R.em_step(pi,P,B2,y,ph); assert np.all(P[:,~ALLOW]==0)
        if ll-old<tol*abs(ll): break
        old=ll
    return R.loglik(pi,P,B2,y,ph),pi,P
K={'MC1_12':14,'HMM_12':20,'MC2_12':44,'MC1_6':26,'HMM_6':32,'MC2_6':80}
def fly(a):
    g,i,d=a; gi=sqw.GROUPS.index(g); y,ph,_=sqw2.seqs(g,i,0.5); ph6=R.phases(len(y),'6h'); miss=y==-1
    L={}; L['HMM_12']=d['F']['3']['ll']
    L['MC1_12'],pim,Pm=mc_fit(y,ph,2); L['MC1_6'],pi6,P6=mc_fit(y,ph6,4)
    L['MC2_12'],pi2,P2=mc2_fit(y,ph,2); L['MC2_6']=mc2_fit(y,ph6,4)[0]
    L['HMM_6']=R.em_fit(y,ph6,4,3,20,20+1000*gi+i)[0]
    rng=np.random.default_rng(14000+1000*gi+i); sims={'MC2_12':[],'MC1_6':[]}
    for _ in range(NSIM):
        for name,(pi,P,B,pp) in (('MC2_12',(pi2,P2,B2,ph)),('MC1_6',(pi6,P6,I3,ph6))):
            _,ys=R.simulate(pi,P,B,pp,rng); ys[miss]=-1; sims[name].append(ppc_stats(ys,ph,2))
    T=len(y); return dict(g=g,fly=i,T=T,LL=L,BIC={m:-2*L[m]+K[m]*np.log(T) for m in K},AIC={m:-2*L[m]+2*K[m] for m in K},sims=sims)
if __name__=='__main__':
    if sys.argv[1]=='check':
        y,ph,_=sqw2.seqs('CS',0,0.5); ll,pi,P=mc2_fit(y,ph,2)
        C=np.zeros((2,9,9))
        for t in range(2,len(y)): C[ph[t],3*y[t-2]+y[t-1],3*y[t-1]+y[t]]+=1
        Pc=C/np.maximum(C.sum(2,keepdims=True),1e-300); cond=np.sum(C[C>0]*np.log(Pc[C>0]))
        first=max(Pc[ph[1],3*a+y[0],3*y[0]+y[1]] for a in range(3))
        print('EM ll',ll,'| closed form (t>=2) + best first step',cond+np.log(first),'| diff',ll-cond-np.log(first))
        ys=y.copy(); ys[np.random.default_rng(1).random(len(y))<0.05]=-1; print('with 5% missing',mc2_fit(ys,ph,2)[0])
    else:
        F=sqw2.load(0.5); jobs=[(g,i,F[(g,i)]) for g in GS for i in range(len(sqw.files(g)))]
        with Pool(4) as pl: out=pl.map(fly,jobs,chunksize=1)
        json.dump(out,open(sqw.OUTD+'o_memory_drift.json','w')); print('done',len(out))
