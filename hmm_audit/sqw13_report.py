# PRESPEC 20n report: likelihood, predictive checks, thermodynamics (n_main.json) and bootstrap (n_boot.jsonl).
import os,sys,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
from scipy.stats import wilcoxon, binom
M=json.load(open(sqw.OUTD+'n_main.json')); PH=['light','dark']
STATS=['chg_h','q_run','w_run','s_run','q_1min','w_1min','c_obs']
print('## 1. likelihood')
for g in ('CS','fumin'):
    D=[d for d in M if d['g']==g]; dll=np.array([d['ll_hmm']-d['ll_mc'] for d in D]); per=1000*dll/np.array([d['T'] for d in D])
    bT=np.array([d['dBIC_T'] for d in D]); bn=np.array([d['dBIC_neff'] for d in D])
    print(f'{g:6s} n {len(D)} | dLL median {np.median(dll):.1f} [{dll.min():.1f}, {dll.max():.1f}] nats | per 1000 min {np.median(per):.1f}'
          f' | dBIC(n=T) median {np.median(bT):.1f} min {bT.min():.1f}, >10 in {np.mean(bT>10):.2f} | dBIC(n_eff) median {np.median(bn):.1f}, >10 in {np.mean(bn>10):.2f}')
print('## 3. posterior predictive: fraction of fly-phases with the real value outside the central 95% of 20 simulations')
for g in ('CS','fumin'):
    D=[d for d in M if d['g']==g]
    for p in range(2):
        row=[]
        for s in STATS:
            fr={}
            for mdl in ('MC','HMM'):
                o=[]
                for d in D:
                    r=d['real'][p][s]; v=np.array([x[p][s] for x in d['sims'][mdl]],float); v=v[np.isfinite(v)]
                    if np.isfinite(r) and len(v)>=10: lo,hi=np.percentile(v,[2.5,97.5]); o.append(r<lo or r>hi)
                fr[mdl]=np.mean(o) if o else np.nan
            row.append(f'{s} MC {fr["MC"]:.2f} HMM {fr["HMM"]:.2f}')
        print(f'{g:6s} {PH[p]:5s} | '+' | '.join(row))
    for p in range(2):
        for s in ('q_1min','w_1min','chg_h'):
            r=np.array([d['real'][p][s] for d in D]); mc=np.array([np.nanmean([x[p][s] for x in d['sims']['MC']]) for d in D])
            hm=np.array([np.nanmean([x[p][s] for x in d['sims']['HMM']]) for d in D])
            print(f'   {g} {PH[p]} {s}: real median {np.nanmedian(r):.3f} | MC sim median {np.nanmedian(mc):.3f} | HMM sim median {np.nanmedian(hm):.3f}')
print('## 4. thermodynamics, exact stationary, per fly (HMM non-degenerate flies)')
for g in ('CS','fumin'):
    D=[d for d in M if d['g']==g and not d['deg']]
    for p in range(2):
        for q in ('J','EPR'):
            a=np.array([d['th_mc'][p][q] for d in D]); b=np.array([d['th_hmm'][p][q] for d in D]); ok=np.isfinite(a)&np.isfinite(b)
            pw=wilcoxon(a[ok],b[ok]).pvalue if ok.sum()>5 else np.nan
            extra=f' | infinite MC {np.mean(~np.isfinite(a)):.2f} HMM {np.mean(~np.isfinite(b)):.2f}' if q=='EPR' else \
                  f' | same sign {np.mean(np.sign(a[ok])==np.sign(b[ok])):.2f}'
            print(f'{g:6s} {PH[p]:5s} {q:3s} n {ok.sum()} | MC median {np.median(a[ok]):+.3f} | HMM median {np.median(b[ok]):+.3f} | ratio of medians {np.median(a[ok])/np.median(b[ok]):.2f} | paired Wilcoxon p {pw:.2g}{extra}')
bf=sqw.OUTD+'n_boot.jsonl'
if os.path.exists(bf):
    B=[json.loads(l) for l in open(bf)]; Mi={(d['g'],d['fly']):d for d in M}
    print('## 2. bootstrap LRT (sign of J in the model per J/h convention: positive = forward)')
    n_top=0; nf=0
    for k in sorted({(b['g'],b['fly']) for b in B}):
        L=np.array([b['LR'] for b in B if (b['g'],b['fly'])==k]); obs=2*(Mi[k]['ll_hmm']-Mi[k]['ll_mc'])
        p=(1+np.sum(L>=obs))/(1+len(L)); nf+=1; n_top+=int(np.all(L<obs)) and len(L)==10
        print(f'  {k[0]:6s} fly {k[1]:2d} | reps {len(L)} | LR_obs {obs:8.1f} | LR* median {np.median(L):6.2f} max {L.max():6.2f} | p {p:.3f}')
    print(f'  flies with LR_obs above all 10 replicates: {n_top} of {nf}; P(X >= {n_top}) under Binomial({nf}, 1/11) = {binom.sf(n_top-1,nf,1/11):.2g}')
