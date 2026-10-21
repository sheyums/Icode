# UNPLANNED, DESCRIPTIVE (requested by user 2026-10-21): drift across recording days.
# (a) flat-HMM held-out score by held-out day (PRESPEC 11 stage 1, at K_g), within-fly centred.
# (b) raw daily behaviour from counts, all flies: movement fraction, sleep min/day (5-min rule),
#     sleep bouts/day, mean sleep bout (min). Within-fly slope per day; Wilcoxon of slopes vs 0.
import sys,json,glob; sys.path.insert(0,'.')
import numpy as np
from scipy.stats import wilcoxon
Kg=json.load(open('flatK_Kg.json')); R=[json.loads(l) for l in open('flatK_stage1.jsonl')]
print("(a) held-out log-lik per second at K_g, within-fly-bin centred, by held-out day (median over fly-bins)")
for g in ['CS','fumin']:
    rows=[]
    for r in [r for r in R if r['g']==g]:
        v=np.array(r['cv'][str(Kg[f"{g}_{r['bin']}"])],float)
        if np.all(np.isfinite(v)): rows.append(v-v.mean())
    A=np.array(rows); sl=[np.polyfit(np.arange(len(v)),v,1)[0] for v in A]
    print(f"  {g:6s} n={len(A)} fly-bins | by day {np.round(np.median(A,0)*3600,1).tolist()} nats/h | slope per day median {np.median(sl)*3600:+.1f} nats/h, Wilcoxon p {wilcoxon(sl).pvalue:.3g}")
CD='/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/counts/'
def daily(x):
    out=[]
    for d in range(int(round(len(x)/86400))):
        s=x[d*86400:(d+1)*86400]; z=(s==0).astype(int); ch=np.r_[0,np.flatnonzero(np.diff(z))+1,len(z)]
        runs=[b-a for a,b in zip(ch[:-1],ch[1:]) if z[a]==1 and b-a>=300]
        out.append(((s>0).mean(), sum(runs)/60, len(runs), np.mean(runs)/60 if runs else np.nan))
    return np.array(out)
print("\n(b) raw daily behaviour, median over flies by day | within-fly slope per day (median) and Wilcoxon p")
names=['movement fraction','sleep (min/day)','sleep bouts/day','mean sleep bout (min)']
for g in ['CS','fumin']:
    D=np.array([daily(np.loadtxt(f)) for f in sorted(glob.glob(CD+g+'/**/*.txt',recursive=True))])
    for k,n in enumerate(names):
        M=D[:,:,k]; sl=np.array([np.polyfit(np.arange(M.shape[1])[np.isfinite(m)],m[np.isfinite(m)],1)[0] for m in M if np.isfinite(m).sum()>=2])
        print(f"  {g:6s} {n:22s} by day {np.round(np.nanmedian(M,0),3).tolist()} | slope {np.median(sl):+.3f}/day, Wilcoxon p {wilcoxon(sl).pvalue:.3g} (n={len(sl)})")
    # excluding day 1 (acclimation) for the two main measures
    for k in (0,1):
        M=D[:,1:,k]; sl=np.array([np.polyfit(np.arange(M.shape[1]),m,1)[0] for m in M])
        print(f"  {g:6s} {names[k]:22s} days 2+ only: slope {np.median(sl):+.3f}/day, Wilcoxon p {wilcoxon(sl).pvalue:.3g}")
