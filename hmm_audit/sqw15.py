# PRESPEC 20p: block KL irreversibility of the S/Q/W symbols vs matched first-order-chain nulls; all genotypes + CS DD.
# Usage: python3 sqw15.py check | run
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
from sqw13 import mc_fit
from sqw14 import mc2_fit, B2, ALLOW
G=['CS','fumin','per0LD','per0DD','sss','CSDD']; DD={'per0DD','CSDD'}; I3=np.eye(3)
NA=200; NB=200; NM=50; NMAX=4
def gidx(g): return 5 if g=='CSDD' else sqw.GROUPS.index(g)
def seqs(g,i):
    """Symbols and phases as 20b; CS DD through the per0DD path (whole-record reference, one matrix)."""
    if g!='CSDD': return sqw2.seqs(g,i,0.5)
    x=np.loadtxt(sqw.files(g)[i]).astype(int); y,ref=sqw.symbols(x,'per0DD',0.5); y=y[:(len(y)//1440)*1440]
    return y,np.zeros(len(y),int),ref
def block_D(y,ph,nph,qw=False):
    """D_n (n = 1..NMAX) per phase: KL divergence between n-blocks and their reverses, pseudocount 0.5 on every block
    of the alphabet ({Q,W,S}, or {Q,W} if qw). Blocks lie in one phase and contain no missing minute."""
    y=np.asarray(y); ph=np.asarray(ph); T=len(y); out=np.zeros((nph,NMAX))
    ok=(y>=0)&((y<2) if qw else True)
    for n in range(2,NMAX+1):
        m=T-n+1; good=np.ones(m,bool); code=np.zeros(m,int); rc=np.zeros(m,int)
        for k in range(n):
            good&=ok[k:k+m]&(ph[k:k+m]==ph[:m]); yk=np.clip(y[k:k+m],0,2)
            code=code*3+yk; rc+=yk*3**k
        allc=np.arange(3**n); digits=np.array([[(c//3**k)%3 for k in range(n)] for c in allc])
        alph=(digits<2).all(1) if qw else np.ones(3**n,bool)
        rev=np.array([sum(digits[c,n-1-k]*3**k for k in range(n)) for c in allc])
        for p in range(nph):
            c=np.bincount(code[good&(ph[:m]==p)],minlength=3**n).astype(float)
            pr=np.where(alph,c+0.5,0.0); pr/=pr.sum(); s=alph
            out[p,n-1]=float(np.sum(pr[s]*np.log(pr[s]/pr[rev[s]])))
    return out
def stats(y,ph,nph):
    D=block_D(y,ph,nph); Dq=block_D(y,ph,nph,True); res=[]
    for p in range(nph):
        d=np.diff(D[p])*60; dq=np.diff(Dq[p])*60                 # d_2, d_3, d_4 in nats/h
        res.append(dict(d2=d[0],d3=d[1],d4=d[2],D3=d[1]-d[0],D4=d[2]-d[1],D3QW=dq[1]-dq[0]))
    return res
def sim_rows(P,allow=None):
    """For simulation only (amendment B): a row with zero expected count is replaced by its count-start row,
    uniform over the allowed transitions."""
    P=np.array(P,float); A=np.ones(P.shape[1:],bool) if allow is None else allow
    for p in range(P.shape[0]):
        for k in range(P.shape[1]):
            if P[p,k].sum()<0.5: P[p,k]=A[k]/A[k].sum()
    return P
def fly(a):
    g,i,d=a; gi=gidx(g); y,ph,_ref=seqs(g,i); nP=1 if g in DD else 2; miss=y==-1; ph6=R.phases(len(y),'6h')
    if d is None:
        ll,pi,P,B,_=R.em_fit(y,ph,nP,3,20,20+1000*gi+i); P,B,pi=R.label_states(P,B,pi); hmm=(pi,P,B); deg=bool(R.is_degenerate(B))
    else: hmm=tuple(np.array(d[k]) for k in ('pi3','P3','B3')); deg=d['deg']
    _,pa,Pa=mc_fit(y,ph,nP); _,pb,Pb=mc_fit(y,ph6,4); _,p2,P2=mc2_fit(y,ph,nP)
    Pa,Pb,P2=sim_rows(Pa),sim_rows(Pb),sim_rows(P2,ALLOW)
    rng=np.random.default_rng(15000+1000*gi+i); sims={}
    for name,(pi,P,B,pp,ns) in (('A',(pa,Pa,I3,ph,NA)),('B',(pb,Pb,I3,ph6,NB)),('HMM',hmm+(ph,NM)),('MC2',(p2,P2,B2,ph,NM))):
        L=[]
        for _ in range(ns):
            _,ys=R.simulate(pi,P,B,pp,rng); ys[miss]=-1; L.append(stats(ys,ph,nP))
        sims[name]=L
    return dict(g=g,fly=i,deg=deg,ref={k:(None if not np.isfinite(v) else float(v)) for k,v in _ref.items()},obs=stats(y,ph,nP),sims=sims,hmm=None if d is not None else dict(pi=hmm[0].tolist(),P=hmm[1].tolist(),B=hmm[2].tolist()))
if __name__=='__main__':
    if sys.argv[1]=='check':
        rng=np.random.default_rng(3)
        P=np.array([[[0.90,0.08,0.02],[0.02,0.95,0.03],[0.05,0.01,0.94]]]); ph=np.zeros(200000,int)
        _,y=R.simulate(np.ones(3)/3,P,I3,ph,rng); s=stats(y,ph,1)[0]
        print('first-order chain, T=2e5:',{k:round(v,4) for k,v in s.items()})
        B=np.array([[0.8,0.2,0],[0.1,0.9,0],[0,0,1]]); _,y=R.simulate(np.ones(3)/3,P,B,ph,rng); s=stats(y,ph,1)[0]
        print('same hidden chain + Q/W misclassification:',{k:round(v,4) for k,v in s.items()})
        y,ph,_=seqs('CSDD',0); print('CS DD fly 0: T',len(y),'missing',(y==-1).mean(),'Q/W/S',[round(float(np.mean(y==v)),3) for v in (0,1,2)])
    else:
        F=sqw2.load(0.5); jobs=[(g,i,None if g=='CSDD' else F[(g,i)]) for g in G for i in range(len(sqw.files(g)))]
        jobs.sort(key=lambda j:j[0]!='CSDD')                        # CS DD (HMM fits) first
        with Pool(4) as pl: out=pl.map(fly,jobs,chunksize=1)
        json.dump(out,open(sqw.OUTD+'p_blockkl.json','w')); print('done',len(out))
