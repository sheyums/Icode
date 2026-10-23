# Count-based cycle current and EPR from posterior hidden paths (analogue of compute_HMM_energetics.m GLOBAL:
# J_cycle = N32 - N23 on the run-compressed sequence, per hour; EPR_obs_full = sum_{i<j} (N_ij - N_ji) ln(N_ij/N_ji)
# per hour), beside the model's stationary flux. Sign: positive = forward A->Q->S->A (H2->H1->H3->H2).
# Added after the 20b results were seen. Usage: python3 sqw2_countJ.py THETA
import sys,os,json,numpy as np
from scipy.stats import spearmanr
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw, sqw2
from sqw2_common import stationary
theta=float(sys.argv[1]); R=sqw2.load(theta); RT=json.load(open(sqw.OUTD+f'route2_{theta}.json')); rng=np.random.default_rng(20); L=[]
def ci(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]
    return v.mean(),*np.percentile([rng.choice(v,len(v)).mean() for _ in range(2000)],[2.5,97.5])
def eprc(N,h):
    e=0.0
    for i in range(3):
        for j in range(i+1,3):
            a,b=N[i][j],N[j][i]
            if a>0 and b>0: e+=(a-b)*np.log(a/b)
            elif a!=b: return np.inf
    return e/h
for g in sqw.GROUPS:
    for p,bn in enumerate(['whole'] if g=='per0DD' else ['light','dark']):
        rows=[]
        for d in [R[k] for k in sorted(R) if k[0]==g and not R[k]['deg']]:
            rt=RT[f"{g}|{d['fly']}"]; h=rt['hours'][p]; P=np.array(d['P3'])[p]; pi=stationary(P)
            Jm=60*(pi[1]*P[1,0]-pi[0]*P[0,1])
            Nf=np.array(rt['pairs_ffbs'][p]); Nv=np.array(rt['pairs_vit'][p])
            rows.append((Jm,(Nf[2,1]-Nf[1,2])/h,(Nv[2,1]-Nv[1,2])/h,eprc(Nf,h),eprc(Nv,h)))
        r=np.array(rows); m=np.isfinite(r[:,0])&np.isfinite(r[:,1]); rho=spearmanr(r[m,0],r[m,1])
        fin=lambda x:x[np.isfinite(x)]
        L.append(f"{g:7s} {bn:5s} n {len(r)} | model J {ci(r[:,0])[0]:+.3f} [{ci(r[:,0])[1]:+.3f},{ci(r[:,0])[2]:+.3f}] | count J_cycle (FFBS) {ci(r[:,1])[0]:+.3f} [{ci(r[:,1])[1]:+.3f},{ci(r[:,1])[2]:+.3f}] | (Viterbi) {ci(r[:,2])[0]:+.3f} cycles/h"
                 f" | Spearman(model, count) {rho.statistic:+.2f} | same sign {np.mean(np.sign(r[m,0])==np.sign(r[m,1])):.2f} | count EPR (FFBS) median {np.median(fin(r[:,3])):.3f} /h, infinite {np.mean(~np.isfinite(r[:,3])):.2f}")
print('\n'.join(L)); open(f'sqw2_countJ_{theta}.log','w').write('\n'.join(L)+'\n')
# Per 4 h clock window: count J from the SAME whole-record posterior paths (no window edges), attributed to the
# window of the arrival minute; window counts sum to the global counts. Profile over ZT (mean over flies of each
# fly's mean over days), cycles/h.
L2=[]
for g in sqw.GROUPS:
    prof=[]
    for d in [R[k] for k in sorted(R) if k[0]==g and not R[k]['deg']]:
        W=np.array(RT[f"{g}|{d['fly']}"]['pairs_win']); jw=(W[:,2,1]-W[:,1,2])/4.0
        prof.append([jw[z::6].mean() for z in range(6)])
        assert np.allclose(W.sum(0),np.array(RT[f"{g}|{d['fly']}"]['pairs_ffbs']).sum(0))
    prof=np.array(prof); L2.append(f"{g:7s} ZT0-4..20-24 count J (cycles/h), mean over flies [SD across flies]: "+'  '.join(f"{m:+.3f} [{s:.2f}]" for m,s in zip(prof.mean(0),prof.std(0))))
    L2.append(f"{'':8s}fraction of flies with J > 0 per window: "+'  '.join(f"{v:.2f}" for v in (prof>0).mean(0)))
print('\n'.join(L2)); open(f'sqw2_countJ_{theta}.log','a').write('\n'.join(L2)+'\n')
# Jensen: nonlinear statistics as posterior means over paths, E[f(N)], beside the plug-in on mean counts, f(E[N]).
# Xn, Yn, c: per path (>= 3 routes per denominator, else undefined for that path); EPR per path per hour.
# J_cycle is linear in counts, so E[J] = J(E[N]) exactly.
def xy(c):
    X=(c[...,0]-c[...,1])/(c[...,0]+c[...,1]); Y=(c[...,2]-c[...,3])/(c[...,2]+c[...,3])
    X=np.where(c[...,0]+c[...,1]>=3,X,np.nan); Y=np.where(c[...,2]+c[...,3]>=3,Y,np.nan); return X,Y
L3=[]
for g in sqw.GROUPS:
    for p,bn in enumerate(['whole'] if g=='per0DD' else ['light','dark']):
        rows=[]
        for d in [R[k] for k in sorted(R) if k[0]==g and not R[k]['deg']]:
            rt=RT[f"{g}|{d['fly']}"]; h=rt['hours'][p]
            with np.errstate(all='ignore'):
                rp=np.array(rt['routes_path'])[:,p]; X,Y=xy(rp); C=(Y-X)/2
                Xm,Ym=xy(rp.mean(0)); Cm=(Ym-Xm)/2
                pp=np.array(rt['pairs_path'])[:,p]; E=np.array([eprc(n,h) for n in pp]); Em=eprc(pp.mean(0),h)
            rows.append((Cm,np.nanmean(C) if np.isfinite(C).mean()>=0.5 else np.nan,np.isfinite(C).mean(),
                         Xm,np.nanmean(X) if np.isfinite(X).mean()>=0.5 else np.nan,Ym,np.nanmean(Y) if np.isfinite(Y).mean()>=0.5 else np.nan,
                         Em,np.median(E),np.mean(~np.isfinite(E))))
        r=np.array(rows,float)
        f=lambda k:np.nanmean(r[:,k])
        L3.append(f"{g:7s} {bn:5s} | c plug-in {f(0):+.3f} vs per-path mean {f(1):+.3f} (paths with c defined {f(2):.2f}) | Xn {f(3):+.3f} vs {f(4):+.3f} | Yn {f(5):+.3f} vs {f(6):+.3f}"
                  f" | EPR plug-in median {np.nanmedian(np.where(np.isfinite(r[:,7]),r[:,7],np.nan)):.3f} vs per-path median {np.nanmedian(np.where(np.isfinite(r[:,8]),r[:,8],np.nan)):.3f} /h, per-path infinite {f(9):.2f}")
print('\n'.join(L3)); open(f'sqw2_countJ_{theta}.log','a').write('\n'.join(L3)+'\n')
# Rule of 3 (user suggestion): a pair (a, b) with b = 0 means the reverse rate is <= 3 counts at 95% (Poisson
# zero-count upper bound). The term f(a,b) = (a-b) ln(a/b) is minimised over b in [0, 3] at b* = min(3, a)
# (f falls to 0 at b = a), so substituting b* gives a 95% LOWER BOUND on that term -- a bound, not an estimate.
def eprc_lb(N,h):
    e=0.0; sub=False
    for i in range(3):
        for j in range(i+1,3):
            a,b=float(N[i][j]),float(N[j][i])
            if a==0 and b==0: continue
            if b==0: b=min(3.0,a); sub=True
            elif a==0: a=min(3.0,b); sub=True
            e+=(a-b)*np.log(a/b)
    return e/h,sub
L4=[]
for g in sqw.GROUPS:
    for p,bn in enumerate(['whole'] if g=='per0DD' else ['light','dark']):
        rows=[]
        for d in [R[k] for k in sorted(R) if k[0]==g and not R[k]['deg']]:
            rt=RT[f"{g}|{d['fly']}"]; h=rt['hours'][p]; pp=np.array(rt['pairs_path'])[:,p]
            lb=[eprc_lb(n,h) for n in pp]; plug=eprc_lb(pp.mean(0),h)[0]
            rows.append((np.mean([x for x,_ in lb]),np.mean([s for _,s in lb]),plug))
        r=np.array(rows)
        L4.append(f"{g:7s} {bn:5s} | EPR 95% lower bound (rule of 3), per-path posterior mean: median over flies {np.median(r[:,0]):.3f} /h, mean {r[:,0].mean():.3f}"
                  f" | paths needing substitution {r[:,1].mean():.2f} | plug-in on mean counts, median {np.median(r[:,2]):.3f} /h | flies {len(r)}")
print('\n'.join(L4)); open(f'sqw2_countJ_{theta}.log','a').write('\n'.join(L4)+'\n')
