# Presentation rule for Figs 13 / 13c (2026-10-29, user-approved after seeing the pooled curves; logged):
# each fly weighted equally within a duration bin (mean of per-fly proportions), a bin shown only if >= 10 flies
# contribute, 95% interval by bootstrap over flies (2000). Hidden level: a fly's records pooled over its 20 paths.
import numpy as np
MIN_FLIES=10
def fly_equal_bins(res, g, phase, lab, dur, edges, level, seed=20):
    rng=np.random.default_rng(seed); xs=[];ys=[];lo=[];hi=[];nf=[]
    for a,b in zip(edges[:-1],edges[1:]):
        props=[]
        for r in res:
            if r['g']!=g: continue
            recs=r['obs'] if level=='observed' else [x for path in r['hidrec'] for x in path]
            sel=[x for x in recs if x['phase_on']==phase and x[lab] is not None and a<=x[dur]<b]
            if sel: props.append(np.mean([bool(x[lab]) for x in sel]))
        if len(props)<MIN_FLIES: continue
        props=np.array(props); bs=[rng.choice(props,len(props)).mean() for _ in range(2000)]
        xs.append(np.sqrt(a*min(b,a*2))); ys.append(props.mean()); l,h=np.percentile(bs,[2.5,97.5]); lo.append(l); hi.append(h); nf.append(len(props))
    return np.array(xs),np.array(ys),np.array(lo),np.array(hi),nf
