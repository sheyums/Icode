# PRESPEC 20s-A report: A1 observable table (+ CS 6-h), A2 EPR vs J and Onsager L, A3 landscape.
import os,sys,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
from scipy.stats import wilcoxon, spearmanr, t as tdist, norm
O=json.load(open(sqw.OUTD+'s_observable.json')); rng=np.random.default_rng(20)
BL=[('CS',0,'light'),('CS',1,'dark'),('per0LD',0,'light'),('per0LD',1,'dark'),('per0DD',0,'whole'),('fumin',0,'light'),('fumin',1,'dark'),
    ('sss',0,'light'),('sss',1,'dark'),('CSDD',0,'whole')]
def ci(v,f=np.median,B=2000):
    v=np.asarray(v,float); v=v[np.isfinite(v)]
    if len(v)<3: return (np.nan,np.nan,np.nan)
    b=[f(rng.choice(v,len(v))) for _ in range(B)]; return (f(v),*np.percentile(b,[2.5,97.5]))
def rows(g,p,src='main'):
    out=[]
    for d in O:
        if d['g']!=g: continue
        m=d[src][p]
        if not m['ok']: continue
        out.append(dict(m,fly=d['fly']))
    return out
fmt=lambda t:f"{t[0]:+.3f} [{t[1]:+.3f},{t[2]:+.3f}]"
print('A1  observable quantities per fly x phase (median over flies [95% fly-bootstrap CI]); J cycles/h, EPR nats/h, A nats/cycle')
for g,p,lab in BL:
    R_=rows(g,p); cyc=[r for r in R_ if r['J_chain'] is not None and abs(r['J_chain'])>=1e-12]
    c=[r['c_obs'] for r in R_]; J=[r['J'] for r in R_]; E=[r['EPR'] for r in R_]
    print(f"{g:6s} {lab:5s} n {len(R_):2d} (chain with cycle {len(cyc)}) | c_obs {fmt(ci(c))} (>0 in {np.mean(np.array(c)>0):.2f}) | count J {fmt(ci(J))} (>0 in {np.mean(np.array(J)>0):.2f}) | count EPR {fmt(ci(E))} (rule of 3 used in {np.mean([r['r3'] for r in R_]):.2f})"
          f" | chain J {np.median([r['J_chain'] for r in cyc]):+.3f} A {np.nanmedian([r['A_chain'] for r in cyc if np.isfinite(r['A_chain'])]):+.2f} EPR {np.nanmedian([r['EPR_chain'] for r in cyc if np.isfinite(r['EPR_chain'])]):.3f} (inf {np.mean([not np.isfinite(r['EPR_chain']) for r in cyc]):.2f})")
    hb=[r for r in R_ if 'J_hmm' in r and not r['hmm_deg']]
    if hb:
        ratJ=[r['J_hmm']/r['J'] for r in hb if abs(r['J'])>1e-9]; ratE=[r['EPR_hmm']/r['EPR'] for r in hb if r['EPR']>1e-9]
        print(f"        HMM bracket (n {len(hb)}): count J {np.median([r['J_hmm'] for r in hb]):+.3f}, count EPR {np.median([r['EPR_hmm'] for r in hb]):.3f}, model EPR {np.nanmedian([r['EPR_hmm_model'] for r in hb if np.isfinite(r['EPR_hmm_model'])]):.3f} | per fly HMM/obs: J {np.median(ratJ):.2f} [IQR {np.percentile(ratJ,25):.2f}-{np.percentile(ratJ,75):.2f}], EPR {np.median(ratE):.2f} [IQR {np.percentile(ratE,25):.2f}-{np.percentile(ratE,75):.2f}] | sign(J) agree {np.mean([np.sign(r['J_hmm'])==np.sign(r['J']) for r in hb if r['J']!=0]):.2f}")
print('\nA1/C3  CS in 6-h ZT blocks (observable)')
names=['ZT0-6','ZT6-12','ZT12-18','ZT18-24']
for b in range(4):
    R_=rows('CS',b,'six'); print(f"  {names[b]:8s} n {len(R_)} | c_obs {fmt(ci([r['c_obs'] for r in R_]))} | count J {fmt(ci([r['J'] for r in R_]))} | count EPR {fmt(ci([r['EPR'] for r in R_]))}")
for k in ('c_obs','J','EPR'):
    dif=[]
    for d in O:
        if d['g']!='CS': continue
        s=d['six']
        if not all(s[b]['ok'] for b in range(4)): continue
        v=[s[b][k] for b in range(4)]
        if not np.all(np.isfinite(v)): continue
        dif.append((v[2]+v[3])/2-(v[0]+v[1])/2)
    print(f"  night - day per fly ({k}): median {np.median(dif):+.3f}, >0 in {np.mean(np.array(dif)>0):.2f}, paired Wilcoxon p {wilcoxon(dif).pvalue:.2g} (n {len(dif)})")
# ---------- A2 ----------
def lsq(x,y):
    A=np.c_[np.ones_like(x),x]; coef,_,_,_=np.linalg.lstsq(A,y,rcond=None); return coef[1],coef[0],float(((A@coef-y)**2).sum())
def models(rr,kJ,kE):
    J=np.array([r[kJ] for r in rr]); E=np.array([r[kE] for r in rr]); a,b1,r1=lsq(np.abs(J),E); c,b2,r2=lsq(J**2,E)
    m=(np.abs(J)>=0.05)&(E>=0.05); n=np.polyfit(np.log(np.abs(J[m])),np.log(E[m]-0.02),1)[0] if m.sum()>=5 else np.nan
    return dict(a=a,c=c,rss1=r1,rss2=r2,n=n)
def ols0(rr,kJ,kE):
    x=np.array([r[kJ] for r in rr])**2; y=np.array([r[kE] for r in rr]); s=(x@y)/(x@x); res=y-s*x; return s,np.sqrt((res@res)/(len(x)-1)/(x@x))
GA=['CS','per0LD','per0DD','fumin','sss','CSDD']
for kJ,kE,lab in (('J','EPR','count J and count EPR (rule of 3)'),('J_chain','EPR_chain','exact observed-chain J and EPR')):
    print(f'\nA2  EPR vs J at the observable level: {lab}. Points = fly x phase (both phases pooled), chains with a cycle and finite values.')
    S={}
    for g in GA:
        pts=[]
        for p in ([0] if g in ('per0DD','CSDD') else [0,1]):
            pts+= [dict(r,ph=p) for r in rows(g,p) if r['J_chain'] is not None and abs(r['J_chain'])>=1e-12 and np.isfinite(r[kE]) and np.isfinite(r[kJ])]
        flies=sorted({r['fly'] for r in pts}); est=models(pts,kJ,kE); s,se=ols0(pts,kJ,kE); Bm=[];Bl=[]
        for _ in range(2000):
            pick=rng.choice(flies,len(flies)); rr=[r for f in pick for r in pts if r['fly']==f]; Bm.append(models(rr,kJ,kE)); Bl.append(1/ols0(rr,kJ,kE)[0])
        S[g]=(s,se)
        print(f"  {g:6s} n {len(pts):3d} | M1 a {est['a']:.2f} RSS {est['rss1']:.3f} | M2 c {est['c']:.2f} RSS {est['rss2']:.3f} | better {'M1 linear' if est['rss1']<est['rss2'] else 'M2 quadratic'} (bootstrap share M1 {np.mean([b['rss1']<b['rss2'] for b in Bm]):.2f})"
              f" | exponent {est['n']:.2f} [{np.nanpercentile([b['n'] for b in Bm],2.5):.2f},{np.nanpercentile([b['n'] for b in Bm],97.5):.2f}] | Onsager L = {1/s:.3f} (fly bootstrap {np.percentile(Bl,2.5):.3f}-{np.percentile(Bl,97.5):.3f})")
    ps=[]
    for g in ('per0LD','per0DD','fumin','sss'):
        z=(S[g][0]-S['CS'][0])/np.hypot(S[g][1],S['CS'][1]); ps.append((g,z,2*norm.sf(abs(z))))
    p=np.array([x[2] for x in ps]); o=np.argsort(p); m=len(p); adj=np.empty(m); run=1
    for k in range(m-1,-1,-1): run=min(run,p[o[k]]*m/(k+1)); adj[o[k]]=run
    for (g,z,pp),aa in zip(ps,adj): print(f"    {g} vs CS (1/L) Wald z {z:+.2f} p {pp:.2g} BH {aa:.2g}")
# ---------- A3 ----------
print('\nA3  energy landscape from the observed chains (nats; E1 = quiet wake = 0)')
for g,p,lab in BL:
    R_=[r for r in rows(g,p) if r['J_chain'] is not None and abs(r['J_chain'])>=1e-12]
    e21=[r['E21'] for r in R_]; e31=[r['E31'] for r in R_]
    F=[(r['F_SW'],r['Xn_cnt']) for r in R_ if np.isfinite(r['Xn_cnt']) and np.isfinite(r['F_SW'])]; Fw=[(r['F_WS'],r['Yn_cnt']) for r in R_ if np.isfinite(r['Yn_cnt']) and np.isfinite(r['F_WS'])]
    sx=spearmanr(*zip(*F))[0] if len(F)>4 else np.nan; sy=spearmanr(*zip(*Fw))[0] if len(Fw)>4 else np.nan
    print(f"{g:6s} {lab:5s} n {len(R_):2d} | E2-E1 {fmt(ci(e21))} | E3-E1 {fmt(ci(e31))} | counted - predicted: Xn {np.median([b-a for a,b in F]):+.3f} (rho {sx:.2f}, n {len(F)}), Yn {np.median([b-a for a,b in Fw]):+.3f} (rho {sy:.2f}, n {len(Fw)})")
