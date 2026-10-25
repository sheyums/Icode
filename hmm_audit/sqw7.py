# PRESPEC 20g: round-trip types per sleep bout, observed (1-s counts) and hidden (posterior paths).
import sys,os,json,numpy as np
from multiprocessing import Pool
from scipy.stats import wilcoxon
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
TYPES=['forward','reverse','direct','quiet_both']     # (in, out) = (Q,A), (A,Q), (A,A), (Q,Q)
def tidx(qin,qout): return {(1,0):0,(0,1):1,(0,0):2,(1,1):3}[(int(qin),int(qout))]
def observed(x,g,theta):
    """Per phase: counts of the 4 types (observed)."""
    nph=1 if g=='per0DD' else 2; bouts=R.sleep_bouts(x,True); mv=(np.asarray(x)>0).astype(float)
    ev=R.wake_bout_events(x,bouts); ref=[R.wake_reference(ev,None)] if nph==1 else [R.wake_reference(ev,0),R.wake_reference(ev,1)]
    C=np.zeros((nph,4))
    for k in range(1,len(bouts)-1):
        s,n=bouts[k]; e=s+n; w_before=s-(bouts[k-1][0]+bouts[k-1][1]); w_after=bouts[k+1][0]-e
        if w_before<600 or w_after<600: continue
        p=0 if nph==1 else int((s%R.DAY)>=R.HALF); r=ref[p]
        if not np.isfinite(r): continue
        qin=mv[s-300:s-60].mean()<theta*r; qout=mv[e+60:e+300].mean()<theta*r
        C[p,tidx(qin,qout)]+=1
    return C
def hidden_path(z,ph,nph):
    C=np.zeros((nph,4)); ch=np.r_[0,np.flatnonzero(np.diff(z))+1]; st=z[ch]
    for k in range(1,len(st)-1):
        if st[k]!=2 or st[k-1]==2 or st[k+1]==2: continue
        C[ph[ch[k]],tidx(st[k-1]==0,st[k+1]==0)]+=1
    return C
def one(a):
    g,i,theta,fitfile=a; x=np.loadtxt(sqw.files(g)[i]).astype(int); nph=sqw2.nP(g)
    out=dict(g=g,fly=i,obs=observed(x,g,theta).tolist(),obs_hours=[float(np.sum(sqw2.seqs(g,i,theta)[1]==q)/60) for q in range(nph)])
    d=FITS[fitfile].get((g,i))
    if d is not None and 'P3' in d and not d['deg']:
        y,ph,_=sqw2.seqs(g,i,theta); pi,P,B=(np.array(d[k]) for k in ('pi3','P3','B3'))
        Z=R.ffbs(pi,P,B,y,ph,200,7000+1000*sqw.GROUPS.index(g)+i)
        H=np.stack([hidden_path(z,ph,nph) for z in Z])              # paths x phase x 4
        out['hid_frac']=[(np.nanmean(np.where(H[:,p].sum(1,keepdims=True)>0,H[:,p]/np.maximum(H[:,p].sum(1,keepdims=True),1),np.nan),0)).tolist() for p in range(nph)]
        out['hid_n']=[float(H[:,p].sum(1).mean()) for p in range(nph)]
        out['hid_counts']=[H[:,p].mean(0).tolist() for p in range(nph)]
        out['hours']=[float(np.sum(ph==q)/60) for q in range(nph)]
        def dep(Hp):   # per path: joint - product of marginals, per type; then mean over paths
            n=Hp.sum(1,keepdims=True); ok=n[:,0]>0; J=Hp[ok]/n[ok]; qin=J[:,0]+J[:,3]; qout=J[:,1]+J[:,3]
            prod=np.stack([qin*(1-qout),(1-qin)*qout,(1-qin)*(1-qout),qin*qout],1)
            return (J-prod).mean(0).tolist(), float(qin.mean()), float(qout.mean())
        out['hid_dep']=[dep(H[:,p]) for p in range(nph)]
    return out
FITS={}
if __name__=='__main__':
    res={}
    for theta in (0.5,0.33,0.75):
        fn=sqw.OUTD+f'fit2_{theta}.jsonl'; FITS[fn]={(d['g'],d['fly']):d for d in map(json.loads,open(fn))}
        jobs=[(g,i,theta,fn) for g in sqw.GROUPS for i in range(len(sqw.files(g)))]
        with Pool(4) as pl: res[str(theta)]=pl.map(one,jobs,chunksize=2)
        print('theta',theta,'done',flush=True)
    json.dump(res,open(sqw.OUTD+'g_roundtrips.json','w'))
    rng=np.random.default_rng(20); L=[]; SUM={}
    for theta in ('0.5','0.33','0.75'):
        L.append(f'===== theta {theta}'+(' (PRIMARY)' if theta=='0.5' else ''))
        for level in ('observed','hidden'):
            ps=[];lines=[]
            for g in sqw.GROUPS:
                for p in range(sqw2.nP(g)):
                    bn='whole' if g=='per0DD' else ['light','dark'][p]; Fr=[]
                    for r in res[theta]:
                        if r['g']!=g: continue
                        if level=='observed':
                            c=np.array(r['obs'][p]);
                            if c.sum()>=3: Fr.append(c/c.sum())
                        elif 'hid_frac' in r and r['hid_n'][p]>=3 and np.all(np.isfinite(r['hid_frac'][p])): Fr.append(np.array(r['hid_frac'][p]))
                    # counts per hour, marginals, dependence
                    CH=[];QI=[];QO=[];DP=[]
                    for r in res[theta]:
                        if r['g']!=g: continue
                        if level=='observed':
                            c=np.array(r['obs'][p])
                            if c.sum()<3: continue
                            f=c/c.sum(); qi=f[0]+f[3]; qo=f[1]+f[3]
                            CH.append(c/r['obs_hours'][p]); QI.append(qi); QO.append(qo)
                            DP.append(f-np.array([qi*(1-qo),(1-qi)*qo,(1-qi)*(1-qo),qi*qo]))
                        elif 'hid_frac' in r and r['hid_n'][p]>=3 and np.all(np.isfinite(r['hid_frac'][p])):
                            CH.append(np.array(r['hid_counts'][p])/r['hours'][p]); DP.append(np.array(r['hid_dep'][p][0])); QI.append(r['hid_dep'][p][1]); QO.append(r['hid_dep'][p][2])
                    CH=np.array(CH); DP=np.array(DP)
                    Fr=np.array(Fr); m=Fr.mean(0); ci=np.percentile([Fr[rng.integers(0,len(Fr),len(Fr))].mean(0) for _ in range(2000)],[2.5,97.5],axis=0)
                    d=Fr[:,0]-Fr[:,1]; pw=wilcoxon(d).pvalue if np.any(d!=0) else 1.0; ps.append(pw)
                    top=TYPES[int(np.argmax(m))]; fwdtop=np.mean(np.argmax(Fr,1)==0)
                    pdep=[wilcoxon(DP[:,j]).pvalue if np.any(DP[:,j]!=0) else 1.0 for j in range(4)]
                    SUM[f'{theta}|{level}|{g}|{bn}']=dict(n=len(Fr),mean=m.tolist(),lo=ci[0].tolist(),hi=ci[1].tolist(),fwd_minus_rev=float(d.mean()),p=pw,
                        rate_per_h=CH.mean(0).tolist(),p_q_entry=float(np.mean(QI)),p_q_exit=float(np.mean(QO)),dependence=DP.mean(0).tolist(),dependence_p=pdep)
                    lines.append([f"  {level:8s} {g:7s} {bn:5s} n {len(Fr):2d} | "+' '.join(f"{t} {m[j]:.2f} [{ci[0][j]:.2f},{ci[1][j]:.2f}]" for j,t in enumerate(TYPES))
                                  +f" | most common: {top} | flies with forward most common {fwdtop:.2f} | fwd-rev {d.mean():+.3f} Wilcoxon p {pw:.2g}"
                                  +f"\n      rates/h "+' '.join(f"{t} {CH.mean(0)[j]:.3f}" for j,t in enumerate(TYPES))
                                  +f" | P(Q entry) {np.mean(QI):.2f} P(Q exit) {np.mean(QO):.2f} | joint - product: "+' '.join(f"{t} {DP.mean(0)[j]:+.3f} (p {pdep[j]:.2g})" for j,t in enumerate(TYPES)),f'{theta}|{level}|{g}|{bn}'])
            ps=np.array(ps); o=np.argsort(ps); adj=np.empty(len(ps)); run=0
            for k,j in enumerate(o): run=max(run,(len(ps)-k)*ps[j]); adj[j]=min(1,run)
            for (ln,key),a in zip(lines,adj): L.append(ln+f" Holm {a:.2g}"); SUM[key]['p_holm']=float(a)
    json.dump(SUM,open(sqw.OUTD+'g_summary.json','w'))
    open('sqw7.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
