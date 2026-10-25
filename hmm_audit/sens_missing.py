# PRESPEC 20i sensitivity analysis.
import json,csv,numpy as np,sys
from scipy.stats import wilcoxon
sys.path.insert(0,'.'); import sqw_ref as R, sqw
F={(d['g'],d['fly']):d for d in map(json.loads,open('sqw_out/fit2_0.5.jsonl'))}; RT=json.load(open('sqw_out/route2_0.5.json'))
RTb=list(csv.DictReader(open('roundtrip_per_fly_theta0.5.csv'))); rng=np.random.default_rng(20); L=[]
def ci(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]; b=[rng.choice(v,len(v)).mean() for _ in range(2000)]; return v.mean(),*np.percentile(b,[2.5,97.5])
def noref(g,i,bn): return not np.isfinite(F[(g,i)]['ref'].get(bn,np.nan) if F[(g,i)]['ref'].get(bn) is not None else np.nan)
for g in ['fumin','sss']:
    for p,bn in enumerate(['light','dark']):
        out={}
        for excl in (False,True):
            C=[];J=[];E=[]
            for (gg,i),d in F.items():
                if gg!=g or d['deg'] or (excl and noref(g,i,bn)): continue
                rt=RT[f'{g}|{i}']; rp=np.array(rt['routes_path'])[:,p]; X,Y=R.xy_counts(rp); c=(Y-X)/2
                if np.isfinite(c).mean()>=0.5: C.append(np.nanmean(c))
                pp=np.array(rt['pairs_path'])[:,p]; h=rt['hours'][p]
                J.append(np.mean([R.j_cycle_count(n,h) for n in pp])); E.append(np.mean([R.epr_count_bound(n,h)[0] for n in pp]))
            rows=[r for r in RTb if r['genotype']==g and r['phase']==bn and not (excl and noref(g,int(r['fly_index']),bn))]
            rtv={lev:np.array([float(r['frac_forward'])-float(r['frac_reverse']) for r in rows if r['level']==lev and float(r['n_round_trips'])>=3]) for lev in ('observed','hidden_matched')}
            out[excl]=(len(C),ci(C),ci(J),np.median(E),{k:(v.mean(),wilcoxon(v).pvalue,len(v)) for k,v in rtv.items()})
        nex=sum(1 for (gg,i),d in F.items() if gg==g and not d['deg'] and noref(g,i,bn))
        for excl,lab in ((False,'all'),(True,f'excluding {nex} fly-phases')):
            n,c,j,e,rt=out[excl]
            L.append(f"{g:5s} {bn:5s} {lab:26s} | c_model {c[0]:+.3f} [{c[1]:+.3f},{c[2]:+.3f}] (n {n}) | J {j[0]:+.3f} [{j[1]:+.3f},{j[2]:+.3f}] | EPR median {e:.3f} | round trips fwd-rev: observed {rt['observed'][0]:+.3f} p {rt['observed'][1]:.2g} (n {rt['observed'][2]}), hidden {rt['hidden_matched'][0]:+.3f} p {rt['hidden_matched'][1]:.2g} (n {rt['hidden_matched'][2]})")
open('sens_missing.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
