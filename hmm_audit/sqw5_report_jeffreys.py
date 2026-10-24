# Sensitivity (post hoc, user question): null flip probability per fly smoothed as (flips + 0.5)/(n + 1) instead of flips/n, because with 8 replicates p_f = 0 is common and makes the test anti-conservative.
# PRESPEC 20e sections 1-2 (across days; light vs dark), and section 3 if e_within6_null.json exists.
import json,os,numpy as np, sqw_ref as R
from scipy.stats import spearmanr, norm
O='sqw_out/'; F={(d['g'],d['fly']):d for d in map(json.loads,open(O+'fit2_0.5.jsonl'))}
SP=json.load(open(O+'split2.json')); NU=json.load(open(O+'e_null.json')); L=[]; rng=np.random.default_rng(20)
def pbinom_tail(p,k):
    """P(X >= k) for X ~ Poisson-binomial(p), exact by convolution."""
    d=np.zeros(len(p)+1); d[0]=1
    for q in p: d[1:]=d[1:]*(1-q)+d[:-1]*q; d[0]*=(1-q)
    return d[k:].sum(), (d*np.arange(len(d))).sum()
def nocyc(c): return not np.isfinite(c) or abs(c)<1e-12
L.append('Section 1: across days (first vs second half), hidden-level c_pred; no-change null 8 replicates per fly')
for g in ['CS','fumin','per0LD','per0DD','sss']:
    for p,bn in enumerate(['whole'] if g=='per0DD' else ['light','dark']):
        obs=[];pf=[];dnull=[]; excl=0
        for k,v in SP.items():
            gg,f=k.split('|'); f=int(f)
            if gg!=g or v['first']['deg'] or v['second']['deg'] or F[(g,f)]['deg'] or k not in NU: continue
            c1=(v['first']['yn'][p]-v['first']['xn'][p])/2; c2=(v['second']['yn'][p]-v['second']['xn'][p])/2
            if nocyc(c1) or nocyc(c2): excl+=1; continue
            reps=[(r[0][p],r[1][p]) for r in NU[k] if not (nocyc(r[0][p]) or nocyc(r[1][p]))]
            if len(reps)<4: excl+=1; continue
            obs.append(np.sign(c1)!=np.sign(c2)); pf.append((np.sum([np.sign(a)!=np.sign(b) for a,b in reps])+0.5)/(len(reps)+1)); dnull.append([b-a for a,b in reps])
        k=int(np.sum(obs)); pval,exp=pbinom_tail(pf,k)
        L.append(f"  {g:7s} {bn:5s}: flies {len(obs)} (excluded {excl}) | observed flips {k} ({k/len(obs):.2f}) | expected under no change {exp:.1f} ({exp/len(obs):.2f}) | one-sided p {pval:.3g}")
# Section 2: light vs dark, same fly (12 h fits, theta 0.5), opposite-sign fraction vs bootstrap-SD expectation
L.append('Section 2: light vs dark within fly (12 h fits), descriptive')
B=json.load(open(O+'boot2.json'))
for g in ['CS','fumin','per0LD','sss']:
    sdc=np.nanmedian([np.sqrt(np.array(b['sdx'])**2+np.array(b['sdy'])**2)/2 for b in B if b['g']==g],0)  # rough c SD per phase (independence assumed)
    C=np.array([[R.c_pred(np.array(d['P3'])[p]) for p in (0,1)] for k,d in F.items() if k[0]==g and not d['deg']])
    C=C[np.all(np.isfinite(C)&(np.abs(C)>1e-12),1)]
    opp=np.mean(np.sign(C[:,0])!=np.sign(C[:,1]))
    pr=[(norm.cdf(0,c[0],sdc[0])*(1-norm.cdf(0,c[1],sdc[1]))+(1-norm.cdf(0,c[0],sdc[0]))*norm.cdf(0,c[1],sdc[1])) for c in C]
    rho=spearmanr(C[:,0],C[:,1])
    L.append(f"  {g:7s} flies {len(C)} | opposite signs light vs dark {opp:.2f} | expected from estimation noise if each fly's light and dark values were its true ones ~{np.mean(pr):.2f} (rough: c SD {sdc[0]:.2f}/{sdc[1]:.2f}, normal approx.) | Spearman light-dark {rho.statistic:+.2f} p {rho.pvalue:.2g}")
if os.path.exists(O+'e_within6_null.json'):
    from scipy.stats import wilcoxon
    W6=json.load(open(O+'e_within6_null.json')); F6={(d['g'],d['fly']):d for d in map(json.loads,open(O+'b6/fit2_0.5.jsonl'))}
    usable={l.split(':')[0]:('USABLE=True' in l) for l in open('sqw3_realistic.log') if 'USABLE=' in l}
    L.append('Section 3: within phase, 6 h blocks (block pairs 1-2 = ZT0-6 vs 6-12, 3-4 = ZT12-18 vs 18-24)')
    for g in ['CS','per0LD','fumin']:
        for (a,b) in ((0,1),(2,3)):
            obs=[];pf=[]
            for k,reps in W6.items():
                gg,f=k.split('|'); f=int(f)
                if gg!=g or (g,f) not in F6 or F6[(g,f)]['deg']: continue
                P=np.array(F6[(g,f)]['P3']); ca,cb=R.c_pred(P[a]),R.c_pred(P[b])
                if nocyc(ca) or nocyc(cb): continue
                rr=[r for r in reps if not (nocyc(r[a]) or nocyc(r[b]))]
                if len(rr)<4: continue
                obs.append(np.sign(ca)!=np.sign(cb)); pf.append((np.sum([np.sign(r[a])!=np.sign(r[b]) for r in rr])+0.5)/(len(rr)+1))
            k=int(np.sum(obs)); pval,exp=pbinom_tail(pf,k)
            # bias check: within-phase Wilcoxon on each null replicate dataset
            ps=[]
            for r in range(8):
                d=[W6[k][r][b]-W6[k][r][a] for k in W6 if k.startswith(g+'|') and len(W6[k])>r and not (nocyc(W6[k][r][a]) or nocyc(W6[k][r][b]))]
                if len(d)>5: ps.append(wilcoxon(d).pvalue)
            L.append(f"  {g:7s} blocks {a+1}-{b+1} ({'test' if usable[g] else 'descriptive'}): flies {len(obs)} | observed flips {k} ({k/max(1,len(obs)):.2f}) | expected {exp:.1f} | one-sided p {pval:.3g} | bias check: null-replicate Wilcoxon p<0.05 in {np.mean(np.array(ps)<0.05):.2f} of {len(ps)}")
open('sqw5_report_jeffreys.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
