# PRESPEC 17 (88c6ea7): per-fly E2 effect sizes, taper depth, relaxation times in/out of sleep.
import glob,json,numpy as np
from multiprocessing import Pool
from scipy.stats import wilcoxon, binomtest
from scipy.optimize import curve_fit
CD='/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/counts/'
GROUPS=['CS','fumin','per0LD','per0DD','sss']; DAY=86400; HALF=43200
def blocks(g): return [('whole',None)] if g=='per0DD' else [('light',0),('dark',1)]
def zruns(x):
    z=(x==0).astype(int); ch=np.r_[0,np.flatnonzero(np.diff(z))+1,len(z)]; return [(a,b-a) for a,b in zip(ch[:-1],ch[1:]) if z[a]==1]
def fly_bouts(x):
    mv=(x>0).astype(float); sl=[(a,n) for a,n in zruns(x) if n>=300]; out=[]
    for (a1,n1),(a2,_) in zip(sl[:-1],sl[1:]):
        w0=a1+n1; w1=a2; L=w1-w0
        if L<1200: continue
        mid=(w0+w1)//2; ref=mv[mid-150:mid+150].mean()
        ex10=[mv[w0+60*k:w0+60*(k+1)].mean() for k in range(1,10)]; ap10=[mv[w1-60*(k+1):w1-60*k].mean() for k in range(1,10)]
        ex20=[mv[w0+60*k:w0+60*(k+1)].mean() for k in range(1,20)] if L>=2400 else None
        ap20=[mv[w1-60*(k+1):w1-60*k].mean() for k in range(1,20)] if L>=2400 else None
        out.append(dict(ph=int((w0%DAY)>=HALF),ref=ref,ex10=ex10,ap10=ap10,ex20=ex20,ap20=ap20,depth=ref-np.mean(mv[w1-180:w1-60])))
    return out
f=lambda t,A,tau: 1-A*np.exp(-(t-1)/tau)
def fit(prof):
    t=np.arange(2,21); 
    try: p,_=curve_fit(f,t,prof,p0=[0.3,3],bounds=([-5,0.1],[5,100])); return p
    except Exception: return np.array([np.nan,np.nan])
if __name__=='__main__':
    rng=np.random.default_rng(17); OUT={}
    for g in GROUPS:
        X=[np.loadtxt(p).astype(int) for p in sorted(glob.glob(CD+g+'/**/*.txt',recursive=True))]
        with Pool(4) as pl: B=pl.map(fly_bouts,X)
        for bn,ph in blocks(g):
            fb=[[b for b in bs if ph is None or b['ph']==ph] for bs in B]
            # A: per-fly D
            Ds=[];ps=[];deps=[];rdeps=[]
            for bs in fb:
                if len(bs)<3: continue
                d=np.array([np.mean(np.array(b['ex10'])-np.array(b['ap10'])) for b in bs]); D=d.mean()
                nul=np.abs((d[None,:]*np.where(rng.random((2000,len(d)))<.5,-1,1)).mean(1)); ps.append((1+np.sum(nul>=abs(D)))/2001); Ds.append(D)
                dp=np.mean([b['depth'] for b in bs]); rf=np.mean([b['ref'] for b in bs]); deps.append(dp); rdeps.append(dp/rf if rf>0 else np.nan)
            Ds=np.array(Ds); ps=np.array(ps); deps=np.array(deps); rdeps=np.array(rdeps)
            st=binomtest(int(np.sum(Ds>0)),len(Ds)).pvalue; wd=wilcoxon(deps).pvalue
            # C: normalised 2-20 min profiles, per fly then mean; fly bootstrap
            P=[]
            for bs in fb:
                bb=[b for b in bs if b['ex20'] is not None and b['ref']>0]
                if bb: P.append((np.mean([np.array(b['ap20'])/b['ref'] for b in bb],0),np.mean([np.array(b['ex20'])/b['ref'] for b in bb],0)))
            res=dict(n_flies=len(Ds),D=Ds.tolist(),p_fly=ps.tolist(),depth=deps.tolist(),rel_depth=rdeps.tolist(),sign_p=st,depth_wilcoxon_p=wd,n_fit=len(P))
            if len(P)>=8:
                AP=np.array([p[0] for p in P]); EX=np.array([p[1] for p in P]); fi=fit(AP.mean(0)); fo=fit(EX.mean(0))
                bs_=[]
                for _ in range(1000):
                    i=rng.integers(0,len(P),len(P)); a=fit(AP[i].mean(0)); o=fit(EX[i].mean(0)); bs_.append([a[0],a[1],o[0],o[1]])
                bs_=np.array(bs_); q=lambda v: np.nanpercentile(v,[2.5,97.5]).tolist()
                res.update(ap_mean=AP.mean(0).tolist(),ex_mean=EX.mean(0).tolist(),ap_sem=(AP.std(0,ddof=1)/np.sqrt(len(P))).tolist(),ex_sem=(EX.std(0,ddof=1)/np.sqrt(len(P))).tolist(),
                           A_in=fi[0],tau_in=fi[1],A_out=fo[0],tau_out=fo[1],tau_in_ci=q(bs_[:,1]),tau_out_ci=q(bs_[:,3]),A_in_ci=q(bs_[:,0]),A_out_ci=q(bs_[:,2]),
                           dtau_ci=q(bs_[:,1]-bs_[:,3]),dA_ci=q(bs_[:,0]-bs_[:,2]))
            OUT[f'{g}|{bn}']=res
            r=res; print(f"{g:7s} {bn:5s} flies {r['n_flies']:2d} | D median {np.median(Ds):+.4f} IQR [{np.percentile(Ds,25):+.4f},{np.percentile(Ds,75):+.4f}] D>0 {np.mean(Ds>0):.2f} sign p {st:.3g} | fly perm p<.05 {np.mean(ps<0.05):.2f} | depth median {np.median(deps):+.4f} rel {np.nanmedian(rdeps):+.2f} Wilcoxon p {wd:.3g}",flush=True)
            if 'tau_in' in r: print(f"         fit n={r['n_fit']}: tau_in {r['tau_in']:.2f} min {np.round(r['tau_in_ci'],2).tolist()} A_in {r['A_in']:+.2f} | tau_out {r['tau_out']:.2f} {np.round(r['tau_out_ci'],2).tolist()} A_out {r['A_out']:+.2f} | dtau CI {np.round(r['dtau_ci'],2).tolist()} dA CI {np.round(r['dA_ci'],2).tolist()}",flush=True)
            else: print(f"         fit skipped (flies with >=40 min bouts: {r['n_fit']})")
    json.dump(OUT,open('taper.json','w'),default=float)
