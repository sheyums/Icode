# PRESPEC 20d: (a) reversible null, (b) per-fly attenuation / posterior-predictive check.
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw, sqw2
from sqw2_common import *
from sqw2_synth import PL1,PD1,B1
NM=['N312','N32','N213','N23']; T=5760; rngb=np.random.default_rng(20)
def from_flows(pi,f01,f02,f12):
    F=np.array([[0,f01,f02],[f01,0,f12],[f02,f12,0]]); P=F/np.array(pi)[:,None]; np.fill_diagonal(P,1-P.sum(1)); return P
RL=from_flows((0.1,0.4,0.5),0.0165,0.0165,0.0115); RD=from_flows((0.1,0.3,0.6),0.0165,0.0165,0.0045)
def ener(P):
    P=np.array(P); Z=P<1e-10; pi=stationary(P); J=60*(pi[0]*P[0,1]-pi[1]*P[1,0])
    if np.any(Z!=Z.T): return np.inf,J
    return 60*sum((pi[i]*P[i,j]-pi[j]*P[j,i])*np.log(pi[i]*P[i,j]/(pi[j]*P[j,i])) for i in range(3) for j in range(i+1,3) if not Z[i,j]),J
def cXY(c):
    X=(c[0]-c[1])/(c[0]+c[1]) if c[0]+c[1]>=3 else np.nan; Y=(c[2]-c[3])/(c[2]+c[3]) if c[2]+c[3]>=3 else np.nan; return (Y-X)/2
def route_counts(z,ph,nph):
    c=np.zeros((nph,4))
    for rt,t in sqw.runs_routes(z): c[ph[t],NM.index(rt)]+=1
    return c
def c_obs_min(y,ph,nph):
    t=np.arange(1,len(y)); on=t[(y[:-1]>=0)&(y[:-1]<2)&(y[1:]==2)]; wk=t[(y[:-1]==2)&(y[1:]>=0)&(y[1:]<2)]   # missing (-1) is not wake; out=[]
    for p in range(nph):
        o=on[ph[on]==p]; w=wk[ph[wk]==p]
        out.append(np.mean(y[o-1]==0)-np.mean(y[w]==0) if len(o) and len(w) else np.nan)
    return out
# ---- (a) ----
def a_one(args):
    case,f=args; rng=np.random.default_rng(8000+100*(case=='V1b')+f); ph=phases(T)
    P2=np.stack([RL,RD]) if case=='rev' else np.stack([PL1,PD1])
    _,y=simulate2(stationary(P2[0]),P2[ph],B1,rng); ll,pi,P,B,_=em(y,ph,2,3,20,50+f); P,B,pi=label3b(P,B,pi)
    cf,_=sqw2.decode(pi,P,B,y,ph,f)
    return dict(case=case,f=f,cpred=[(yn(P[p])-xn(P[p]))/2 for p in range(2)],cffbs=[cXY(np.array(cf[p])) for p in range(2)],
                E=[ener(P[p]) for p in range(2)],deg=bool(degenerate(B)))
def stage_a():
    with Pool(4) as pl: res=pl.map(a_one,[(c,f) for c in ('rev','V1b') for f in range(20)],chunksize=1)
    json.dump(res,open(sqw.OUTD+'d_a.json','w'),default=float); L=[]
    def ci(v):
        v=np.array(v,float); v=v[np.isfinite(v)]; b=[rngb.choice(v,len(v)).mean() for _ in range(2000)]; return v.mean(),*np.percentile(b,[2.5,97.5])
    ok=True; null95={}
    for case in ('rev','V1b'):
        R=[r for r in res if r['case']==case]
        for p,bn in enumerate(['light','dark']):
            tP=(RL,RD)[p] if case=='rev' else (PL1,PD1)[p]; tc=(yn(tP)-xn(tP))/2; tE=ener(tP)
            cp=ci([r['cpred'][p] for r in R]); cf=ci([r['cffbs'][p] for r in R]); E=np.array([r['E'][p][0] for r in R])
            L.append(f"{case} {bn}: true c {tc:+.3f} EPR {tE[0]:.3f}/h | c_pred {cp[0]:+.3f} [{cp[1]:+.3f},{cp[2]:+.3f}] | c_ffbs {cf[0]:+.3f} [{cf[1]:+.3f},{cf[2]:+.3f}] | fitted EPR median {np.median(E):.3f} 95th pct {np.percentile(E,95):.3f} /h (inf {np.mean(~np.isfinite(E)):.2f})")
            if case=='rev':
                null95[bn]=float(np.percentile(E[np.isfinite(E)],95))
                for c in (cp,cf): ok&= abs(c[0])<0.05 and c[1]<=0<=c[2]
    L.append(f"(a) PASS={ok}"); json.dump(null95,open(sqw.OUTD+'d_null95.json','w'))
    # fly EPRs against the null's 95th percentile
    R20=sqw2.load(0.5)
    for g in sqw.GROUPS:
        for p,bn in enumerate(['whole'] if g=='per0DD' else ['light','dark']):
            E=[ener(np.array(d['P3'])[p])[0] for d in R20.values() if d['g']==g and not d['deg']]
            thr=null95['light'] if bn in ('light','whole') else null95['dark']
            L.append(f"   fly EPR {g} {bn}: fraction above reversible-null 95th pct ({thr:.3f}/h): {np.mean(np.array(E)>thr):.2f} (n {len(E)})")
    open('sqw4_a.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
# ---- (b) ----
def b_one(d):
    g=d['g']; nph=sqw2.nP(g); y,ph,_=sqw2.seqs(g,d['fly'],0.5); rng=np.random.default_rng(9000+1000*sqw.GROUPS.index(g)+d['fly'])
    P=np.array(d['P3']); real=c_obs_min(y,ph,nph); so=[];st=[]
    for r in range(20):
        z,ys=simulate2(np.array(d['pi3']),P[ph],np.array(d['B3']),rng); so.append(c_obs_min(ys,ph,nph))
        c=route_counts(z,ph,nph); st.append([cXY(c[p]) for p in range(nph)])
    return dict(g=g,fly=d['fly'],real=real,sim_obs=so,sim_true=st)
def stage_b():
    R20=sqw2.load(0.5); RT=json.load(open(sqw.OUTD+'route2_0.5.json'))
    with Pool(4) as pl: res=pl.map(b_one,[d for d in R20.values() if not d['deg']],chunksize=2)
    json.dump(res,open(sqw.OUTD+'d_b.json','w'),default=float); L=[]
    for g in sqw.GROUPS:
        Rg=[r for r in res if r['g']==g]
        for p,bn in enumerate(['whole'] if g=='per0DD' else ['light','dark']):
            real=np.array([r['real'][p] for r in Rg]); so=np.array([[s[p] for s in r['sim_obs']] for r in Rg]); st=np.array([[s[p] for s in r['sim_true']] for r in Rg])
            gm=np.nanmean(so,0); lo,hi=np.percentile(gm,[2.5,97.5]); rm=np.nanmean(real)
            band=np.nanpercentile(so,[2.5,97.5],axis=1); out=np.mean([(x<a or x>b) for x,a,b in zip(real,band[0],band[1]) if np.isfinite(x)])
            cm=np.nanmean([cXY(np.array(RT[f"{g}|{r['fly']}"]['ffbs'][p])) for r in Rg])
            opp=np.nanmean(np.where(np.isfinite(st)&np.isfinite(so),np.sign(st)!=np.sign(so),np.nan))
            L.append(f"{'':13s} simulated records with opposite signs of c_true and c_obs_min: {opp:.2f}")
            L.append(f"{g:7s} {bn:5s}: real c_obs_min {rm:+.3f}; simulated group mean {np.mean(gm):+.3f} [{lo:+.3f},{hi:+.3f}] -> {'consistent' if lo<=rm<=hi else 'INCONSISTENT'}; "
                     f"flies outside own 95% band {out:.2f}; sim c_true {np.nanmean(st):+.3f} / sim c_obs {np.nanmean(so):+.3f} = {np.nanmean(st)/np.nanmean(so):.2f}; real c_model {cm:+.3f} / c_obs {rm:+.3f} = {cm/rm:.2f}")
    open('sqw4_b.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
if __name__=='__main__':
    {'a':stage_a,'b':stage_b}[sys.argv[1]]()
