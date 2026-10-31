# PRESPEC 20q: thermodynamic uncertainty relation per fly. Usage: python3 sqw16.py check | run
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
from sqw13 import mc_fit
from sqw15 import seqs, gidx, sim_rows, DD, G
DSGN=np.zeros((3,3)); Q_,W_,S_=0,1,2
for a,b in ((W_,Q_),(Q_,S_),(S_,W_)): DSGN[a,b]=1; DSGN[b,a]=-1   # forward cycle A->Q->S->A
TAUS=(180,120,360); NMC=200; NHMM=100
def estimate(y,ph,nph,taus=TAUS):
    """Per phase: J (/h) and, per tau, D (/h) and B = J^2 / D (k_B/h) from block variances of the counted current."""
    y=np.asarray(y); ph=np.asarray(ph); T=len(y)
    val=np.zeros(T,bool); val[1:]=(y[1:]>=0)&(y[:-1]>=0); inc=np.zeros(T); inc[1:]=np.where(val[1:],DSGN[np.clip(y[:-1],0,2),np.clip(y[1:],0,2)],0)
    seg=np.flatnonzero(np.diff(ph)!=0)+1; starts=np.r_[0,seg]; ends=np.r_[seg,T]; out=[]
    for p in range(nph):
        m=(ph==p)&val; hours=m.sum()/60; J=inc[m].sum()/hours if hours>0 else np.nan; r=dict(J=J)
        for tau in taus:
            Nb=[]
            for s,e in zip(starts,ends):
                if ph[s]!=p: continue
                for b0 in range(s,e-tau+1,tau):
                    if b0==0 or not val[b0:b0+tau].all(): continue
                    Nb.append(inc[b0:b0+tau].sum())
            D=np.var(Nb,ddof=1)/(2*tau/60) if len(Nb)>=4 else np.nan
            r[f'D{tau}']=D; r[f'B{tau}']=J**2/D if D and np.isfinite(D) and D>0 else np.nan; r[f'n{tau}']=len(Nb)
        out.append(r)
    return out
def _lam(M):
    w=np.linalg.eigvals(M); return float(np.log(np.max(w.real)))
def exact_B(T,d,h=1e-3):
    """Exact J, D, B (per hour) of the current with increments d on the chain T (per-minute steps)."""
    lp=_lam(T*np.exp(h*d)); lm=_lam(T*np.exp(-h*d)); l0=_lam(T)
    J=(lp-lm)/(2*h)*60; D=(lp-2*l0+lm)/h**2/2*60; return J,D,(J**2/D if D>0 else np.nan)      # amendment C: B = J^2/D
def exact_B_ct(T,d,h=1e-3):
    """Same for the continuous-time generator T - I (off-diagonal tilt), for the discrete-time inflation kappa."""
    def lam(s):
        M=T*np.exp(s*d); np.fill_diagonal(M,np.diag(T)); return float(np.max(np.linalg.eigvals(M).real))-1.0
    lp,lm=lam(h),lam(-h); J=(lp-lm)/(2*h)*60; D=(lp+lm)/h**2/2*60; return J,D,(J**2/D if D>0 else np.nan)
def hmm_ext(P,B):
    """9-state chain on (hidden z, symbol y): T[(z,y),(z',y')] = P[z,z'] B[z',y']; current increments d[y,y']."""
    T=np.einsum('ab,bc->abc',P,B)[:,None,:,:].repeat(3,1).reshape(9,9); d=np.tile(DSGN,(3,3)); return T,d
def fly(a):
    g,i,d,refs=a; gi=gidx(g); y,ph,_=seqs(g,i); nP=1 if g in DD else 2; miss=y==-1
    _,pm,Pm=mc_fit(y,ph,nP); Pm=sim_rows(Pm); hm=tuple(np.array(d[k]) for k in ('pi','P','B')); deg=d['deg']
    obs=estimate(y,ph,nP); rng=np.random.default_rng(16000+1000*gi+i); simMC=[];simH=[]
    for _ in range(NMC):
        _,ys=R.simulate(pm,Pm,np.eye(3),ph,rng); ys[miss]=-1; simMC.append(estimate(ys,ph,nP))
    for _ in range(NHMM):
        _,ys=R.simulate(*hm,ph,rng); ys[miss]=-1; simH.append(estimate(ys,ph,nP))
    ex=[]
    for p in range(nP):
        Jm,Dm,Bm=exact_B(Pm[p],DSGN); Th,dh=hmm_ext(hm[1][p],hm[2]); Jh,Dh,Bh=exact_B(Th,dh)
        Bm_ct=exact_B_ct(Pm[p],DSGN)[2]; Bh_ct=exact_B_ct(Th,dh)[2]
        ex.append(dict(sig_mc=float(R.model_energetics(Pm[p])['EPR']),J_mc=Jm,D_mc=Dm,B_mc=Bm,B_mc_ct=Bm_ct,B_hmm_ct=Bh_ct,
                       sig_hmm=float(R.model_energetics(hm[1][p])['EPR']),J_hmm=Jh,D_hmm=Dh,B_hmm=Bh))
    return dict(g=g,fly=i,deg=deg,ref=refs,obs=obs,ex=ex,simMC=simMC,simH=simH)
if __name__=='__main__':
    if sys.argv[1]=='check':
        rng=np.random.default_rng(5); worst=0
        for _ in range(1000):
            P=rng.dirichlet(np.ones(3)*rng.uniform(0.2,3),3); P=0.5*P+0.5*np.eye(3)*rng.uniform(0,1.9)
            P/=P.sum(1,keepdims=True); s=R.model_energetics(P)['EPR']; _,_,B=exact_B(P,DSGN)
            if np.isfinite(s) and s>1e-9: worst=max(worst,B/s); wct=max(locals().get('wct',0),exact_B_ct(P,DSGN)[2]/s)
        print('TUR check over 1000 random chains (discrete time; slack up to ~1.07 allowed): max B/sigma =',round(worst,4))
        print('continuous-time generator: max B/sigma =',round(wct,6),'(must be <= 1)')
        P=np.array([[[0.90,0.06,0.04],[0.03,0.95,0.02],[0.02,0.04,0.94]]]); ph=np.zeros(1_000_000,int)
        _,y=R.simulate(np.ones(3)/3,P,np.eye(3),ph,rng); e=estimate(y,ph,1)[0]; J,D,B=exact_B(P[0],DSGN)
        print('long record: J est',round(e['J'],4),'exact',round(J,4),'| D(3h) est',round(e['D180'],4),'exact',round(D,4),'| B est',round(e['B180'],4),'exact',round(B,4),'| sigma',round(R.model_energetics(P[0])['EPR'],4))
        Th,dh=hmm_ext(P[0],np.eye(3)); print('HMM extension with B = I reproduces the chain:',np.allclose(exact_B(Th,dh),exact_B(P[0],DSGN),rtol=1e-6))
    else:
        O={(d['g'],d['fly']):d for d in json.load(open(sqw.OUTD+'p_blockkl.json'))}; F=sqw2.load(0.5); jobs=[]
        for g in G:
            for i in range(len(sqw.files(g))):
                o=O[(g,i)]; h=o['hmm'] if o['hmm'] else dict(pi=F[(g,i)]['pi3'],P=F[(g,i)]['P3'],B=F[(g,i)]['B3'])
                jobs.append((g,i,dict(h,deg=o['deg']),o['ref']))
        with Pool(4) as pl: out=pl.map(fly,jobs,chunksize=1)
        json.dump(out,open(sqw.OUTD+'q_tur.json','w')); print('done',len(out))
