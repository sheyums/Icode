# PRESPEC 20b summaries. Usage: python3 sqw2_report.py THETA
import sys,os,json,numpy as np
from scipy.stats import spearmanr
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw, sqw2
from sqw2_common import kpar2, stationary
from cycle import fly as cyc_fly, stats as cyc_stats
theta=float(sys.argv[1]); R=sqw2.load(theta); OUTD=sqw.OUTD
RT=json.load(open(OUTD+f'route2_{theta}.json')) if os.path.exists(OUTD+f'route2_{theta}.json') else {}
rng=np.random.default_rng(20); L=[]; OUT={}
def ci(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]
    if len(v)<2: return (np.nan,np.nan,np.nan)
    return (v.mean(),*np.percentile([rng.choice(v,len(v)).mean() for _ in range(2000)],[2.5,97.5]))
def fmt(v): m,a,b=ci(v); return f"{m:+.3f} [{a:+.3f},{b:+.3f}]"
def ener(P):
    P=np.array(P); Z=P<1e-10; pi=stationary(P)
    J=60*(pi[0]*P[0,1]-pi[1]*P[1,0])
    if np.any(Z!=Z.T) or Z[0,1] or Z[1,2] or Z[2,0]: return np.inf,J,np.inf
    epr=sum((pi[i]*P[i,j]-pi[j]*P[j,i])*np.log(pi[i]*P[i,j]/(pi[j]*P[j,i])) for i in range(3) for j in range(i+1,3))
    return 60*epr,J,np.log(P[0,1]*P[1,2]*P[2,0]/(P[0,2]*P[2,1]*P[1,0]))
for g in sqw.GROUPS:
    D=[R[k] for k in sorted(R) if k[0]==g]; nP=sqw2.nP(g)
    if len(D[0]['F'])==3:
        best=[];beste=[];flag=0
        for d in D:
            F={int(k):v['ll'] for k,v in d['F'].items()}; b={K:-2*F[K]+kpar2(K,nP)*np.log(d['T']) for K in F}
            be={K:-2*F[K]+kpar2(K,nP)*np.log(max(2,d['nchg'])) for K in F}
            best.append(min(b,key=b.get)); beste.append(min(be,key=be.get)); flag+=min(b[2],b[4])<b[3]-10
        L.append(f"{g}: flies {len(D)} | best K (n=T) "+' '.join(f"{K}:{np.mean([x==K for x in best]):.2f}" for K in (2,3,4))
                 +" | (n_eff) "+' '.join(f"{K}:{np.mean([x==K for x in beste]):.2f}" for K in (2,3,4))
                 +f" | rules disagree {np.mean([a!=b for a,b in zip(best,beste)]):.2f} | flagged {flag/len(D):.2f} | H1 degenerate {np.mean([d['deg'] for d in D]):.2f}")
        OUT[g]=dict(best=best,best_eff=beste,flag=flag/len(D))
    nd=[d for d in D if not d['deg']]
    for p,bn in enumerate(['whole'] if nP==1 else ['light','dark']):
        xp=[d['xn'][p] for d in nd]; yp=[d['yn'][p] for d in nd]; o=dict(xp=xp,yp=yp)
        L.append(f"  {bn:5s} (non-degenerate flies {len(nd)}) Xn_pred {fmt(xp)}  Yn_pred {fmt(yp)}  c_pred {fmt([(b-a)/2 for a,b in zip(xp,yp)])}")
        for kind in ('ffbs','vit'):
            X=[];Y=[]
            for d in nd:
                c=np.array(RT[f"{g}|{d['fly']}"][kind][p]) if f"{g}|{d['fly']}" in RT else np.zeros(4)
                X.append((c[0]-c[1])/(c[0]+c[1]) if c[0]+c[1]>=3 else np.nan); Y.append((c[2]-c[3])/(c[2]+c[3]) if c[2]+c[3]>=3 else np.nan)
            cm=[(b-a)/2 for a,b in zip(X,Y)]; o[kind]=dict(X=X,Y=Y,c=cm)
            L.append(f"        {kind:4s} Xn {fmt(X)}  Yn {fmt(Y)}  c_model {fmt(cm)}")
        cobs=[];obs=[]
        for d in nd:
            x=np.loadtxt(sqw.files(g)[d['fly']]).astype(int); ev,Ln=cyc_fly(x); ph=None if bn=='whole' else p
            s=cyc_stats(ev,ph,theta,(Ln/3600)/(1 if ph is None else 2))
            obs.append(None if s is None else (2*s['nSQ']/(s['nSQ']+s['nSA'])-1,2*s['nQS']/(s['nQS']+s['nAS'])-1)); cobs.append(np.nan if s is None else s['c'])
        o['obs']=obs; o['cobs']=cobs
        cm=np.array(o['ffbs']['c']); cb=np.array(cobs); m=np.isfinite(cm)&np.isfinite(cb)
        rho=spearmanr(cb[m],cm[m]) if m.sum()>4 else None
        L.append(f"        model-free c {fmt(cobs)}; Spearman(c_obs, c_model ffbs) "+(f"{rho.statistic:+.2f} p {rho.pvalue:.2g} n {m.sum()}" if rho else 'n/a'))
        E=[ener(np.array(d['P3'])[p]) for d in nd]; inf=np.mean([not np.isfinite(e[0]) for e in E]); fin=[e for e in E if np.isfinite(e[0])]
        if fin: L.append(f"        energetics: infinite EPR in {inf:.2f} of flies; finite median EPR {np.median([e[0] for e in fin]):.3f} nats/h, J {np.median([e[1] for e in fin]):+.4f} cycles/h, affinity {np.median([e[2] for e in fin]):+.2f}")
        o['E']=E; OUT[f'{g}|{bn}']=o
print('\n'.join(L)); open(f'sqw2_{theta}.log','w').write('\n'.join(L)+'\n')
json.dump(OUT,open(OUTD+f'report2_{theta}.json','w'),default=lambda v:(v.item() if hasattr(v,'item') else str(v)))
