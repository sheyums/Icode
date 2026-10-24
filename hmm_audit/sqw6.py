# PRESPEC 20f: hidden vs observable EPR (difference, ratio) with calibration on synthetic flies.
import sys,os,json,numpy as np
from multiprocessing import Pool
from scipy.stats import wilcoxon
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw2
DRIVEN={('CS',1),('fumin',0),('fumin',1),('per0LD',0),('per0LD',1),('per0DD',0),('sss',0)}
BN=lambda g,p:'whole' if g=='per0DD' else ['light','dark'][p]
def obs_pairs(y,ph,nph):
    y=np.asarray(y); t=np.flatnonzero(y>=0); ys=y[t]; ps=np.asarray(ph)[t]
    ch=np.flatnonzero(np.diff(ys))+1; ok=(t[ch]-t[ch-1])==1
    N=np.zeros((nph,3,3)); np.add.at(N,(ps[ch[ok]],ys[ch-1][ok],ys[ch][ok]),1); return N
def real_one(d):
    g=d['g']; nph=sqw2.nP(g); y,ph,_=sqw2.seqs(g,d['fly'],0.5); out={}
    rt=RT[f"{g}|{d['fly']}"]; No=obs_pairs(y,ph,nph)
    for p in range(nph):
        h=np.sum(ph==p)/60; Eh=np.mean([R.epr_count_bound(n,h)[0] for n in np.array(rt['pairs_path'])[:,p]])
        out[p]=(Eh,R.epr_count_bound(No[p],h)[0])
    return f"{g}|{d['fly']}",out
def sim_one(d):
    g=d['g']; nph=sqw2.nP(g); T=d['T']; ph=R.phases(T,'none' if g=='per0DD' else '12h')
    pi,P,B=(np.array(d[k]) for k in ('pi3','P3','B3')); rng=np.random.default_rng(2020+1000*['CS','fumin','per0LD','per0DD','sss'].index(g)+d['fly'])
    hours=[np.sum(ph==p)/60 for p in range(nph)]; Em=[R.model_energetics(P[p])['EPR'] for p in range(nph)]; reps=[]
    for r in range(20):
        z,y=R.simulate(pi,P,B,ph,rng); Np=R.pair_counts(z,ph,nph); No=obs_pairs(y,ph,nph)
        Z=R.ffbs(pi,P,B,y,ph,100,int(rng.integers(1<<30)))
        rec={}
        for p in range(nph):
            Ef=np.mean([R.epr_count_bound(R.pair_counts(zz,ph,nph)[p],hours[p])[0] for zz in Z])
            rec[p]=(Em[p],R.epr_count_bound(Np[p],hours[p])[0],Ef,R.epr_count_bound(No[p],hours[p])[0])
        reps.append(rec)
    return f"{g}|{d['fly']}",reps
if __name__=='__main__':
    O=sqw2.OUTD; F=sqw2.load(0.5); RT=json.load(open(O+'route2_0.5.json'))
    D=[d for d in F.values() if not d['deg']]
    with Pool(4) as pl: real=dict(pl.map(real_one,D,chunksize=2))
    with Pool(4) as pl: sim=dict(pl.map(sim_one,D,chunksize=1))
    json.dump(dict(real={k:{str(p):v for p,v in o.items()} for k,o in real.items()},sim={k:[{str(p):v for p,v in r.items()} for r in reps] for k,reps in sim.items()}),
              open(O+'f_epr.json','w'),default=float)
    rng=np.random.default_rng(20); L=[]; tests=[]
    for g in ['CS','fumin','per0LD','per0DD','sss']:
        for p in range(sqw2.nP(g)):
            keys=[k for k in real if k.startswith(g+'|')]; Eh=np.array([real[k][p][0] for k in keys]); Eo=np.array([real[k][p][1] for k in keys])
            Dd=Eh-Eo; m=Eh>0; r=Eo[m]/Eh[m]
            bD=np.percentile([np.median(rng.choice(Dd,len(Dd))) for _ in range(2000)],[2.5,97.5]); br=np.percentile([np.median(rng.choice(r,len(r))) for _ in range(2000)],[2.5,97.5])
            pw=wilcoxon(Dd).pvalue if np.any(Dd!=0) else 1.0
            S=np.array([[sim[k][i][p] for k in keys] for i in range(20)])     # (rep, fly, 4)
            Em,Ep,Ef,Es=S[...,0],S[...,1],S[...,2],S[...,3]
            with np.errstate(all='ignore'):
                bias_path=np.nanmedian(np.where(Ep>0,(Ef-Ep)/Ep,np.nan)); bias_model=np.nanmedian(np.where(np.isfinite(Em)&(Em>0),(Ef-Em)/Em,np.nan))
                simr=[np.nanmedian(np.where(Ef[i]>0,Es[i]/Ef[i],np.nan)) for i in range(20)]
            lo,hi=np.percentile(simr,[2.5,97.5]); simneg=np.mean(Ef-Es<0)
            driven=(g,p) in DRIVEN
            L.append(f"{g:7s} {BN(g,p):5s} ({'driven' if driven else 'not driven: descriptive'}) flies {len(keys)} | median EPR_h {np.median(Eh):.3f} EPR_o {np.median(Eo):.3f} | D median {np.median(Dd):+.3f} [{bD[0]:+.3f},{bD[1]:+.3f}] Wilcoxon p {pw:.2g} | D<0 in {np.mean(Dd<0):.2f} (sim {simneg:.2f})"
                     f" | r median {np.median(r):.2f} [{br[0]:.2f},{br[1]:.2f}], simulated {np.median(simr):.2f} [{lo:.2f},{hi:.2f}] -> {'consistent' if lo<=np.median(r)<=hi else 'INCONSISTENT'}"
                     f" | E_ffbs bias vs true path {bias_path:+.2f}, vs model {bias_model:+.2f}")
            if driven: tests.append((len(L)-1,pw,bias_path))
    ps=np.array([t[1] for t in tests]); o=np.argsort(ps); adj=np.empty(len(ps)); run=0
    for i,j in enumerate(o): run=max(run,(len(ps)-i)*ps[j]); adj[j]=min(1,run)
    for (li,pw,bp),a in zip(tests,adj): L[li]+=f" | Holm p {a:.2g} -> hidden dissipation {'REPORTED' if a<0.05 and abs(bp)<0.25 else 'descriptive'}"
    open('sqw6.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
