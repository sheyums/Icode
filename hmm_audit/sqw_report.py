# PRESPEC 20 summaries. Usage: python3 sqw_report.py THETA  (needs fit_THETA.jsonl, route_THETA.json)
import sys,os,json,numpy as np
from scipy.stats import spearmanr
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
from sqw_common import *
WIN=sqw.WIN
from cycle import fly as cyc_fly, stats as cyc_stats
theta=float(sys.argv[1]); R=sqw.load_fit(theta); OUTD=sqw.OUTD
RT=json.load(open(OUTD+f'route_{theta}.json')) if os.path.exists(OUTD+f'route_{theta}.json') else {}
rng=np.random.default_rng(20); L=[]; OUT={}
def ci(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]
    if len(v)<2: return (np.nan,np.nan,np.nan)
    b=[rng.choice(v,len(v)).mean() for _ in range(2000)]; return (v.mean(),*np.percentile(b,[2.5,97.5]))
def ener(P):
    P=np.array(P); Z=P<1e-10; pi=stationary(P)
    if np.any(Z!=Z.T) or (Z[0,1]|Z[1,2]|Z[2,0]): return np.inf,pi[0]*P[0,1]-pi[1]*P[1,0],np.inf
    epr=sum((pi[i]*P[i,j]-pi[j]*P[j,i])*np.log(pi[i]*P[i,j]/(pi[j]*P[j,i])) for i in range(3) for j in range(i+1,3) if not Z[i,j])
    return 60*epr,60*(pi[0]*P[0,1]-pi[1]*P[1,0]),np.log(P[0,1]*P[1,2]*P[2,0]/(P[0,2]*P[2,1]*P[1,0]))
NM=['N312','N32','N213','N23']
for g in sqw.GROUPS:
    for bn in (['whole'] if g=='per0DD' else ['light','dark']):
        W=[d for k,ws in R.items() if k[0]==g for d in ws if d['block']==bn]; Fw=[d for d in W if d['fit']]
        o={}
        if Fw and len(Fw[0]['F'])==3:
            best=[];beste=[];flag=0;db2=[];db4=[]
            for d in Fw:
                F={int(k):v['ll'] for k,v in d['F'].items()}; b={K:bic(F[K],K,WIN) for K in F}
                be={K:bic(F[K],K,max(2,d['nchg'])) for K in F}
                best.append(min(b,key=b.get)); beste.append(min(be,key=be.get)); db2.append(b[2]-b[3]); db4.append(b[4]-b[3])
                flag+=min(b[2],b[4])<b[3]-10
            fr={K:np.mean([x==K for x in best]) for K in (2,3,4)}; fre={K:np.mean([x==K for x in beste]) for K in (2,3,4)}
            o.update(bestK=fr,bestK_eff=fre,flag=flag/len(Fw),Kdiff=np.mean([a!=b for a,b in zip(best,beste)]),
                     dbic2_med=float(np.median(db2)),dbic4_med=float(np.median(db4)),nchg_med=float(np.median([d['nchg'] for d in Fw])))
            L.append(f"{g:7s} {bn:5s} windows {len(W)} fitted {len(Fw)} | best K (n=240) 2:{fr[2]:.2f} 3:{fr[3]:.2f} 4:{fr[4]:.2f} | (n_eff) 2:{fre[2]:.2f} 3:{fre[3]:.2f} 4:{fre[4]:.2f} | K differs {o['Kdiff']:.2f} | flagged {o['flag']:.2f} | median dBIC(2-3) {o['dbic2_med']:+.1f} (4-3) {o['dbic4_med']:+.1f} | median symbol changes {o['nchg_med']:.0f}")
        nd=[d for d in Fw if 'P3' in d and not d['deg']]; o['deg']=1-len(nd)/max(1,len(Fw))
        L.append(f"{'':13s} H1 degenerate in {o['deg']:.2f} of fitted windows")
        # per fly
        flies=sorted({d['fly'] for d in W}); rows=[]
        for f in flies:
            fw=[d for d in nd if d['fly']==f]; r=dict(fly=f)
            r['xp']=np.nanmean([d['xn'] for d in fw]) if fw else np.nan; r['yp']=np.nanmean([d['yn'] for d in fw]) if fw else np.nan
            rt=RT.get(f'{g}|{f}')
            for kind in ('ffbs','vit'):
                if rt is None: r[kind]=None; continue
                c=np.sum([rt[kind][d['w']] for d in fw],0) if fw else np.zeros(4)
                X=(c[0]-c[1])/(c[0]+c[1]) if c[0]+c[1]>=3 else np.nan; Y=(c[2]-c[3])/(c[2]+c[3]) if c[2]+c[3]>=3 else np.nan
                r[kind]=(X,Y,c.tolist())
            x=np.loadtxt(sqw.files(g)[f]).astype(int); ev,Ln=cyc_fly(x); ph=None if bn=='whole' else (0 if bn=='light' else 1)
            s=cyc_stats(ev,ph,theta,(Ln/3600)/(1 if ph is None else 2))
            r['obs']=None if s is None else (2*s['nSQ']/(s['nSQ']+s['nSA'])-1,2*s['nQS']/(s['nQS']+s['nAS'])-1,s['c'])
            E=[ener(d['P3']) for d in fw]; r['E']=E; rows.append(r)
        def col(fn): return [fn(r) for r in rows]
        xp,yp=col(lambda r:r['xp']),col(lambda r:r['yp'])
        L.append(f"{'':13s} Xn_pred {ci(xp)[0]:+.3f} [{ci(xp)[1]:+.3f},{ci(xp)[2]:+.3f}]  Yn_pred {ci(yp)[0]:+.3f} [{ci(yp)[1]:+.3f},{ci(yp)[2]:+.3f}]")
        for kind in ('ffbs','vit'):
            X=col(lambda r:r[kind][0] if r[kind] else np.nan); Y=col(lambda r:r[kind][1] if r[kind] else np.nan)
            cm=[(y-x)/2 for x,y in zip(X,Y)]
            L.append(f"{'':13s} {kind:4s} Xn {ci(X)[0]:+.3f} [{ci(X)[1]:+.3f},{ci(X)[2]:+.3f}] Yn {ci(Y)[0]:+.3f} [{ci(Y)[1]:+.3f},{ci(Y)[2]:+.3f}] c_model {ci(cm)[0]:+.3f} [{ci(cm)[1]:+.3f},{ci(cm)[2]:+.3f}] (flies {np.isfinite(cm).sum()})")
            o[kind]=dict(X=X,Y=Y,c=cm)
        cobs=col(lambda r:r['obs'][2] if r['obs'] else np.nan); cff=o['ffbs']['c']; m=np.isfinite(cobs)&np.isfinite(cff)
        rho=spearmanr(np.array(cobs)[m],np.array(cff)[m]) if m.sum()>4 else None
        L.append(f"{'':13s} model-free c {ci(cobs)[0]:+.3f} [{ci(cobs)[1]:+.3f},{ci(cobs)[2]:+.3f}]; Spearman(c_obs, c_model ffbs) "+(f"{rho.statistic:+.2f} p {rho.pvalue:.2g} (n {m.sum()})" if rho else 'n/a'))
        Ew=[e for r in rows for e in r['E']]; inf=np.mean([not np.isfinite(e[0]) for e in Ew]) if Ew else np.nan
        fin=[e for e in Ew if np.isfinite(e[0])]
        if fin:
            L.append(f"{'':13s} energetics over {len(Ew)} windows: infinite EPR/affinity in {inf:.2f}; finite: median EPR {np.median([e[0] for e in fin]):.3f} nats/h, median J {np.median([e[1] for e in fin]):+.4f} cycles/h, median affinity {np.median([e[2] for e in fin]):+.2f}")
        o['obs']=col(lambda r:r['obs']); o['xp']=xp; o['yp']=yp; OUT[f'{g}|{bn}']=o
print('\n'.join(L)); open(f'sqw_{theta}.log','w').write('\n'.join(L)+'\n')
json.dump(OUT,open(OUTD+f'report_{theta}.json','w'),default=lambda v:None if isinstance(v,float) and not np.isfinite(v) else (v.item() if hasattr(v,'item') else str(v)))
