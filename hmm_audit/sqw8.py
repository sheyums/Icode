# PRESPEC 20h: dwell times, observed S/Q/W (with model-predictive band) and hidden H1/H2/H3. theta 0.5.
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
OBS={'S':{2},'wake':{0,1},'Q':{0},'W':{1}}; HID={'H1':{0},'H2':{1},'H3':{2}}; NREP=20; NP=200
def fly(d):
    g,i=d['g'],d['fly']; nph=sqw2.nP(g); y,ph,_=sqw2.seqs(g,i,0.5); pi,P,B=(np.array(d[k]) for k in ('pi3','P3','B3'))
    out=dict(g=g,fly=i,obs={},obs_naive={},sim={},hid={},pstay=[[float(P[p][k,k]) for k in range(3)] for p in range(nph)])
    for nm,v in OBS.items():
        out['obs'][nm]=[(a.tolist(),b.tolist()) for a,b in R.runs_censored(y,ph,v,nph)]
        out['obs_naive'][nm]=[(a.tolist(),b.tolist()) for a,b in R.runs_censored(y,ph,v,nph,censor_at_phase=False)]
    rng=np.random.default_rng(8000+1000*sqw.GROUPS.index(g)+i); out['sim']={nm:[] for nm in OBS}
    for r in range(NREP):
        _,ys=R.simulate(pi,P,B,ph,rng)
        for nm,v in OBS.items(): out['sim'][nm].append([(a.tolist(),b.tolist()) for a,b in R.runs_censored(ys,ph,v,nph)])
    Z=R.ffbs(pi,P,B,y,ph,NP,9000+1000*sqw.GROUPS.index(g)+i)
    for nm,v in HID.items():
        acc=[([],[]) for _ in range(nph)]
        for z in Z:
            for p,(a,b) in enumerate(R.runs_censored(z,ph,v,nph)): acc[p][0].extend(a.tolist()); acc[p][1].extend(b.tolist())
        out['hid'][nm]=acc
    return out
if __name__=='__main__':
    F=sqw2.load(0.5)
    with Pool(4) as pl: res=pl.map(fly,[d for d in F.values() if not d['deg']],chunksize=1)
    json.dump(res,open(sqw.OUTD+'h_dwell.json','w'))
    E=R.DWELL_EDGES; L=[]; SUM={}
    for g in sqw.GROUPS:
        for p in range(sqw2.nP(g)):
            bn='whole' if g=='per0DD' else ['light','dark'][p]; Rg=[r for r in res if r['g']==g]
            for nm in OBS:
                d=np.concatenate([r['obs'][nm][p][0] for r in Rg]); e=np.concatenate([r['obs'][nm][p][1] for r in Rg])
                dn=np.concatenate([r['obs_naive'][nm][p][0] for r in Rg])
                t,S=R.km_survival(d,e); ev,ex,nr=R.life_table_hazard(d,e)
                Hs=[];Ss=[]
                for k in range(NREP):
                    ds=np.concatenate([r['sim'][nm][k][p][0] for r in Rg]); es=np.concatenate([r['sim'][nm][k][p][1] for r in Rg])
                    a,b,c=R.life_table_hazard(ds,es); Hs.append(a/np.maximum(b,1e-9))
                    ts,Sv=R.km_survival(ds,es); Ss.append(np.interp(np.arange(1,721),ts,Sv,left=1.0,right=Sv[-1]))
                Hs=np.array(Hs); lo,hi=np.percentile(Hs,[2.5,97.5],0); h=ev/np.maximum(ex,1e-9); show=nr>=10
                outside=np.mean(((h<lo)|(h>hi))[show]) if show.any() else np.nan
                med=[]
                for r in Rg:
                    tt,SS=R.km_survival(r['obs'][nm][p][0],r['obs'][nm][p][1]) if len(r['obs'][nm][p][0]) else (np.array([]),np.array([]))
                    if len(SS) and SS[-1]<=0.5: med.append(tt[np.argmax(SS<=0.5)])
                SUM[f'obs|{g}|{bn}|{nm}']=dict(n=len(d),censored=float(1-e.mean()),t=t.tolist(),S=S.tolist(),h=h.tolist(),h_lo=lo.tolist(),h_hi=hi.tolist(),n_risk=nr.tolist(),
                    sim_S_lo=np.percentile(Ss,2.5,0).tolist(),sim_S_hi=np.percentile(Ss,97.5,0).tolist(),sim_S_med=np.median(Ss,0).tolist(),
                    naive_mean=float(dn.mean()) if len(dn) else None,km_median_fly_median=float(np.median(med)) if med else None)
                L.append(f"obs {g:7s} {bn:5s} {nm:4s} runs {len(d):5d} (censored {1-e.mean():.2f}) | per-fly KM median {np.median(med) if med else float('nan'):.0f} min | hazard bins outside model 95% band {outside:.2f} of {show.sum()} | h per bin: "+' '.join(f'{x:.3f}' if s else '-' for x,s in zip(h,show)))
            for k,nm in enumerate(HID):
                d=[];e=[];w=[]
                for r in Rg:
                    a,b=r['hid'][nm][p]; d+=a; e+=b; w+=[1/NP]*len(a)
                d=np.array(d);e=np.array(e);w=np.array(w); t,S=R.km_survival(d,e,w); ev,ex,nr=R.life_table_hazard(d,e,w=w)
                ps=np.array([r['pstay'][p][k] for r in Rg]); tg=np.arange(0,721); wf=np.array([len(r['hid'][nm][p][0]) for r in Rg],float)
                Sg=np.average([R.geometric_survival(q,tg) for q in ps],axis=0,weights=np.maximum(wf,1e-9))   # weighted by each fly's decoded run count, like the pooled KM
                mean_dwell=1/(1-ps)
                hist,_=np.histogram(d,bins=E,weights=w)
                SUM[f'hid|{g}|{bn}|{nm}']=dict(n_per_path=float(len(d)/NP),t=t.tolist(),S=S.tolist(),geo_S=Sg.tolist(),hist=hist.tolist(),
                    mean_dwell_fly=mean_dwell.tolist(),h=(ev/np.maximum(ex,1e-9)).tolist(),n_risk=nr.tolist())
                L.append(f"hid {g:7s} {bn:5s} {nm} runs per path summed over flies {len(d)/NP:7.1f} | model mean dwell 1/(1-P_ii): median over flies {np.median(mean_dwell):.1f} min [IQR {np.percentile(mean_dwell,25):.1f}-{np.percentile(mean_dwell,75):.1f}] | decoded KM median {t[np.argmax(S<=0.5)] if (S<=0.5).any() else float('nan'):.0f} min")
    json.dump(SUM,open(sqw.OUTD+'h_summary.json','w'))
    open('sqw8.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
