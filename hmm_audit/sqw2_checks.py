# PRESPEC 20b section 5 summaries: day split (descriptive) and within-phase adequacy.
import json,numpy as np
from scipy.stats import spearmanr
O='sqw_out/'; SP=json.load(open(O+'split2.json')); AD=json.load(open(O+'adeq2.json')); rng=np.random.default_rng(20); L=[]
for g in ['CS','fumin','per0LD','per0DD','sss']:
    for p,bn in enumerate(['whole'] if g=='per0DD' else ['light','dark']):
        c1=[];c2=[]
        for k,v in SP.items():
            if not k.startswith(g+'|') or v['first']['deg'] or v['second']['deg']: continue
            c1.append((v['first']['yn'][p]-v['first']['xn'][p])/2); c2.append((v['second']['yn'][p]-v['second']['xn'][p])/2)
        c1=np.array(c1); c2=np.array(c2); m=np.isfinite(c1)&np.isfinite(c2); d=c2[m]-c1[m]
        b=np.percentile([rng.choice(d,len(d)).mean() for _ in range(2000)],[2.5,97.5]); rho=spearmanr(c1[m],c2[m])
        L.append(f"day split {g:7s} {bn:5s}: flies {m.sum()} | c_pred first {np.mean(c1[m]):+.3f} second {np.mean(c2[m]):+.3f} | change {d.mean():+.3f} [{b[0]:+.3f},{b[1]:+.3f}] | Spearman {rho.statistic:+.2f} (p {rho.pvalue:.2g}) | sign flips {np.mean(np.sign(c1[m])!=np.sign(c2[m])):.2f}")
L.append('within-phase adequacy: fraction of 4 h windows outside the fitted model\'s 95% predictive band, by ZT window (0-4..20-24); columns Q, W, S fraction, symbol changes; ~0.05 expected')
for g in ['CS','fumin','per0LD','per0DD','sss']:
    A=np.array([np.array(v) for k,v in AD.items() if k.startswith(g+'|')],dtype=object)
    rows=np.concatenate([np.array(a).reshape(-1,6,4) if len(a)%6==0 else np.array(a)[:len(a)//6*6].reshape(-1,6,4) for a in A])
    fr=rows.mean(0)
    L.append(f"  {g:7s} "+' | '.join(f"ZT{4*z}-{4*z+4}: "+' '.join(f'{v:.2f}' for v in fr[z]) for z in range(6)))
open('sqw2_checks.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
