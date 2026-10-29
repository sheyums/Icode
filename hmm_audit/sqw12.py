# PRESPEC 20m: 30-s vs 1-min steps, CS and fumin.
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
NP=200
def metrics(pi,P,B,y,ph,bin_s,seed):
    """Per phase: H1 mean dwell (min), fraction of decoded H1 runs of one step, count J (/h), EPR bound (/h),
    c_pred, FFBS c, hidden R1 D (matched definition rescaled)."""
    nph=P.shape[0]; spm=60//bin_s; Z=R.ffbs(pi,P,B,y,ph,NP,seed); out=[]
    for p in range(nph):
        hours=np.sum(ph==p)*bin_s/3600
        runs=[R.runs_censored(z,ph,{0},nph)[p] for z in Z]
        one=np.mean([np.mean(d==1) for d,e in runs if len(d)])
        pp=np.stack([R.pair_counts(z,ph,nph)[p] for z in Z])
        J=np.mean([R.j_cycle_count(n,hours) for n in pp]); E=np.mean([R.epr_count_bound(n,hours)[0] for n in pp])
        rp=np.stack([R.route_counts(z,ph,nph)[p] for z in Z]); X,Y=R.xy_counts(rp); C=(Y-X)/2
        cf=float(np.nanmean(C)) if np.isfinite(C).mean()>=0.5 else None
        Ds=np.array([R.route_duration_D(R.route_duration_hidden(z,ph,min_sleep=5*spm,min_gap=10*spm),'R1_exit_vs_sleep',p) for z in Z],float)
        D=float(np.nanmean(Ds)) if np.isfinite(Ds).mean()>=0.5 else None
        out.append(dict(h1_dwell_min=float(bin_s/60/(1-P[p][0,0])),h1_one_step=float(one),J=float(J),EPR=float(E),
                        c_pred=float(R.c_pred(P[p])),c_ffbs=cf,R1_D=D,noncycle=bool(abs(R.model_energetics(P[p],per_hour=3600/bin_s)['J'])<1e-12)))
    return out
def fly(a):
    g,i,d=a; gi=sqw.GROUPS.index(g); x=R.read_counts(sqw.files(g)[i]); res=dict(g=g,fly=i)
    for bin_s in (60,30):
        y,refs=R.sqw_symbols_bin(x,bin_s,0.5,True); ph=R.phases_bin(len(y),bin_s)
        if bin_s==60:
            pi,P,B=(np.array(d[k]) for k in ('pi3','P3','B3')); deg=d['deg']
        else:
            ll,pi,P,B,_=R.em_fit(y,ph,2,3,20,20+1000*gi+i); P,B,pi=R.label_states(P,B,pi); deg=bool(R.is_degenerate(B))
        res[str(bin_s)]=dict(deg=deg,refs={k:(None if not np.isfinite(v) else v) for k,v in refs.items()},
                             B=B.tolist(),m=None if deg else metrics(pi,P,B,y,ph,bin_s,7000+1000*gi+i))
    return res
if __name__=='__main__':
    F=sqw2.load(0.5); jobs=[(g,i,F[(g,i)]) for g in ('CS','fumin') for i in range(len(sqw.files(g)))]
    with Pool(4) as pl: out=pl.map(fly,jobs,chunksize=1)
    json.dump(out,open(sqw.OUTD+'m_resolution.json','w')); print('done',len(out))
