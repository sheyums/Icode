# PRESPEC 20n: Markov chain on S/Q/W symbols (HMM with B = I) vs the 20b 3-state HMM, CS and fumin.
# Usage: python3 sqw13.py main | boot
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
GS=('CS','fumin'); NSIM=20; NREP=10; NBOOT_FLIES=10; I3=np.eye(3)
def mc_fit(y,ph,nP,tol=1e-10,maxit=500):
    """MLE of a first-order Markov chain on the symbols = HMM with B fixed to I; missing minutes marginalised."""
    y=np.asarray(y); P=np.full((nP,3,3),1.0)
    for t in range(1,len(y)):
        if y[t-1]>=0 and y[t]>=0: P[ph[t],y[t-1],y[t]]+=1      # count start (+1 keeps every entry positive at start)
    P/=P.sum(2,keepdims=True); pi=np.full(3,1/3); old=-np.inf
    for it in range(maxit):
        pi,P,B,ll=R.em_step(pi,P,I3,y,ph); occ=[v for v in range(3) if np.any(y==v)]; assert np.allclose(B[occ],I3[occ],atol=1e-12)
        if ll-old<tol*abs(ll): break
        old=ll
    return R.loglik(pi,P,I3,y,ph),pi,P
def ppc_stats(y,ph,nph):
    """Per phase: changes/h, mean Q run, mean W run, mean S run, 1-min fraction of Q runs and of W runs, c_obs."""
    y=np.asarray(y); t=np.arange(1,len(y)); ch=t[(y[1:]>=0)&(y[:-1]>=0)&(y[1:]!=y[:-1])]
    rq,rw,rs=(R.runs_censored(y,ph,{v},nph) for v in (0,1,2)); c=R.c_obs_minutes(y,ph,nph); out=[]
    m=lambda d,e:float(np.mean(d[e==1])) if np.any(e==1) else np.nan
    f1=lambda d,e:float(np.mean(d[e==1]==1)) if np.any(e==1) else np.nan
    for p in range(nph):
        h=np.sum(ph==p)/60
        out.append(dict(chg_h=float(np.sum(ph[ch]==p)/h),q_run=m(*rq[p]),w_run=m(*rw[p]),s_run=m(*rs[p]),
                        q_1min=f1(*rq[p]),w_1min=f1(*rw[p]),c_obs=float(c[p])))
    return out
def thermo(P):
    e=R.model_energetics(P); return {k:float(v) for k,v in e.items()}
def fly(a):
    g,i,d=a; gi=sqw.GROUPS.index(g); y,ph,_=sqw2.seqs(g,i,0.5); nP=sqw2.nP(g); miss=y==-1
    llm,pim,Pm=mc_fit(y,ph,nP); llh=d['F']['3']['ll']; pih,Ph,Bh=(np.array(d[k]) for k in ('pi3','P3','B3'))
    T=len(y); nch=R.n_eff(y); kh=R.n_params(3,nP); km=kh-6
    rng=np.random.default_rng(13000+1000*gi+i); sims={'MC':[],'HMM':[]}
    for _ in range(NSIM):
        for name,(pi,P,B) in (('MC',(pim,Pm,I3)),('HMM',(pih,Ph,Bh))):
            _,ys=R.simulate(pi,P,B,ph,rng); ys[miss]=-1; sims[name].append(ppc_stats(ys,ph,nP))
    return dict(g=g,fly=i,T=T,nchg=nch,ll_mc=llm,ll_hmm=llh,deg=d['deg'],
                dBIC_T=float((-2*llm+km*np.log(T))-(-2*llh+kh*np.log(T))),
                dBIC_neff=float((-2*llm+km*np.log(nch))-(-2*llh+kh*np.log(nch))),
                pi_mc=pim.tolist(),P_mc=Pm.tolist(),real=ppc_stats(y,ph,nP),sims=sims,
                th_mc=[thermo(Pm[p]) for p in range(nP)],th_hmm=[thermo(Ph[p]) for p in range(nP)])
def boot(a):
    g,i,r,pim,Pm=a; gi=sqw.GROUPS.index(g); y,ph,_=sqw2.seqs(g,i,0.5); nP=sqw2.nP(g); miss=y==-1
    rng=np.random.default_rng(130000+10000*gi+100*i+r); _,ys=R.simulate(np.array(pim),np.array(Pm),I3,ph,rng); ys[miss]=-1
    llm=mc_fit(ys,ph,nP)[0]; llh=R.em_fit(ys,ph,nP,3,20,20+1000*gi+i+100000*(r+1))[0]
    return dict(g=g,fly=i,rep=r,ll_mc=llm,ll_hmm=llh,LR=float(2*(llh-llm)))
if __name__=='__main__':
    F=sqw2.load(0.5); fn=sqw.OUTD+'n_main.json'
    if sys.argv[1]=='main':
        jobs=[(g,i,F[(g,i)]) for g in GS for i in range(len(sqw.files(g)))]
        with Pool(4) as pl: out=pl.map(fly,jobs,chunksize=1)
        json.dump(out,open(fn,'w')); print('main done',len(out))
    else:
        M={(d['g'],d['fly']):d for d in json.load(open(fn))}; rng=np.random.default_rng(2026)
        sel={g:sorted(rng.choice(len(sqw.files(g)),NBOOT_FLIES,replace=False).tolist()) for g in GS}; print('boot flies',sel)
        jobs=[(g,i,r,M[(g,i)]['pi_mc'],M[(g,i)]['P_mc']) for g in GS for i in sel[g] for r in range(NREP)]
        bf=sqw.OUTD+'n_boot.jsonl'; done=set()
        if os.path.exists(bf): done={(d['g'],d['fly'],d['rep']) for d in map(json.loads,open(bf))}
        jobs=[j for j in jobs if j[:3] not in done]
        with Pool(4) as pl:
            for r in pl.imap_unordered(boot,jobs):
                open(bf,'a').write(json.dumps(r)+'\n'); print(r['g'],r['fly'],r['rep'],'LR*',round(r['LR'],2),flush=True)
