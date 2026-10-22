# PRESPEC 19 (571f90b): signed cycle current S/Q/A from counts, +-5 min window, bridging (b).
import glob,json,sys,numpy as np
from multiprocessing import Pool
from scipy.stats import binomtest, wilcoxon, mannwhitneyu
sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from v18 import sleep_runs, CD, DAY, HALF
GROUPS=['CS','fumin','per0LD','per0DD','sss']
def blocks(g): return [('whole',None)] if g=='per0DD' else [('light',0),('dark',1)]
def fly(x):
    mv=(x>0).astype(float); sl=sleep_runs(x,True); ev=[]
    for (a1,n1),(a2,_) in zip(sl[:-1],sl[1:]):
        w0=a1+n1; w1=a2
        if w1-w0<600: continue
        mid=(w0+w1)//2; ref=mv[mid-150:mid+150].mean()
        ev.append(dict(ref=ref, wake=(int((w0%DAY)>=HALF), mv[w0+60:w0+300].mean()), onset=(int((w1%DAY)>=HALF), mv[w1-300:w1-60].mean())))
    return ev, len(x)
def stats(ev,ph,theta,hours):
    bb=[e for e in ev if (ph is None) or e['wake'][0]==ph or e['onset'][0]==ph]
    if len(bb)<3: return None
    r=np.median([e['ref'] for e in bb]); thr=theta*r
    on=[e['onset'][1]<thr for e in ev if ph is None or e['onset'][0]==ph]; wk=[e['wake'][1]<thr for e in ev if ph is None or e['wake'][0]==ph]
    if not on or not wk: return None
    nQS=sum(on); nAS=len(on)-nQS; nSQ=sum(wk); nSA=len(wk)-nSQ
    c=nQS/len(on)-nSQ/len(wk); J=((nQS-nSQ)+(nSA-nAS))/2/hours
    return dict(c=c,J=J,nQS=nQS,nAS=nAS,nSQ=nSQ,nSA=nSA)
if __name__=='__main__':
    rng=np.random.default_rng(19); OUT={}
    data={g:[np.loadtxt(p).astype(int) for p in sorted(glob.glob(CD+g+'/**/*.txt',recursive=True))] for g in GROUPS}
    EV={}
    for g in GROUPS:
        with Pool(4) as pl: EV[g]=pl.map(fly,data[g])
    for theta in (0.5,0.33,0.75):
        print(f"\n===== theta = {theta}{'  (PRIMARY)' if theta==0.5 else ''} =====")
        for g in GROUPS:
            for bn,ph in blocks(g):
                rows=[]
                for ev,L in EV[g]:
                    hours=(L/3600)/(1 if ph is None else 2)
                    s=stats(ev,ph,theta,hours)
                    if s: rows.append(s)
                c=np.array([r['c'] for r in rows]); J=np.array([r['J'] for r in rows])
                bs=np.array([rng.choice(c,len(c)).mean() for _ in range(2000)]); bJ=np.array([rng.choice(J,len(J)).mean() for _ in range(2000)])
                tot={k:sum(r[k] for r in rows) for k in ('nQS','nAS','nSQ','nSA')}
                OUT[f'{theta}|{g}|{bn}']=dict(c=c.tolist(),J=J.tolist(),c_ci=np.percentile(bs,[2.5,97.5]).tolist(),J_ci=np.percentile(bJ,[2.5,97.5]).tolist(),tot=tot,
                    sign_p=binomtest(int((c>0).sum()),int((c!=0).sum())).pvalue if (c!=0).any() else 1.0, wil_p=wilcoxon(c).pvalue if np.any(c!=0) else 1.0)
                o=OUT[f'{theta}|{g}|{bn}']
                print(f"{g:7s} {bn:5s} flies {len(c):2d} | c mean {c.mean():+.3f} [{o['c_ci'][0]:+.3f},{o['c_ci'][1]:+.3f}] median {np.median(c):+.3f} c>0 {np.mean(c>0):.2f} sign p {o['sign_p']:.2g} Wilcoxon p {o['wil_p']:.2g} | J {J.mean():+.3f}/h [{o['J_ci'][0]:+.3f},{o['J_ci'][1]:+.3f}] | onsets Q:A {tot['nQS']}:{tot['nAS']} wakings Q:A {tot['nSQ']}:{tot['nSA']}",flush=True)
        print("  sleepless vs others (Mann-Whitney on c):")
        for bn in ('light','dark'):
            s=OUT[f'{theta}|sss|{bn}']['c']
            for g in ['CS','fumin','per0LD']:
                print(f"    {bn}: sss vs {g}: p {mannwhitneyu(s,OUT[f'{theta}|{g}|{bn}']['c']).pvalue:.3g}")
            print(f"    {bn}: sss vs per0DD(whole): p {mannwhitneyu(s,OUT[f'{theta}|per0DD|whole']['c']).pvalue:.3g}")
    json.dump(OUT,open('cycle.json','w'),default=lambda o: o.item() if hasattr(o,'item') else str(o))
