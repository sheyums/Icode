# PRESPEC 20g amendment B (post hoc): hidden round trips matched to the observed definition.
import sys,os,json,numpy as np
from multiprocessing import Pool
from scipy.stats import wilcoxon
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
TYPES=['forward','reverse','direct','quiet_both']
def matched(z,ph,nph,min_s=5,gap=10):
    C=np.zeros((nph,4)); ch=np.r_[0,np.flatnonzero(np.diff(z))+1,len(z)]; st=z[ch[:-1]]; ln=np.diff(ch)
    s_idx=np.flatnonzero(st==2)
    for a,k in enumerate(s_idx):
        if ln[k]<min_s or k==0 or k==len(st)-1: continue
        prev_end=ch[k]; prev_s=ch[s_idx[a-1]+1] if a>0 else 0          # end of previous H3 run
        next_s=ch[s_idx[a+1]] if a+1<len(s_idx) else len(z)
        if prev_end-prev_s<gap or next_s-ch[k+1]<gap: continue
        if a==0 or a+1==len(s_idx): continue
        C[ph[ch[k]],{(1,0):0,(0,1):1,(0,0):2,(1,1):3}[(int(st[k-1]==0),int(st[k+1]==0))]]+=1
    return C
def one(d):
    g,i=d['g'],d['fly']; nph=sqw2.nP(g); y,ph,_=sqw2.seqs(g,i,0.5); pi,P,B=(np.array(d[k]) for k in ('pi3','P3','B3'))
    Z=R.ffbs(pi,P,B,y,ph,200,7000+1000*sqw.GROUPS.index(g)+i); H=np.stack([matched(z,ph,nph) for z in Z])
    out=[]
    for p in range(nph):
        Hp=H[:,p]; n=Hp.sum(1); ok=n>0
        if n.mean()<3 or ok.mean()<0.5: out.append(None); continue
        J=Hp[ok]/n[ok,None]; qi=J[:,0]+J[:,3]; qo=J[:,1]+J[:,3]
        prod=np.stack([qi*(1-qo),(1-qi)*qo,(1-qi)*(1-qo),qi*qo],1)
        out.append(dict(frac=J.mean(0).tolist(),rate=(Hp.mean(0)/(np.sum(ph==p)/60)).tolist(),qin=float(qi.mean()),qout=float(qo.mean()),dep=(J-prod).mean(0).tolist(),n=float(n.mean())))
    return f'{g}|{i}',out
if __name__=='__main__':
    F=sqw2.load(0.5)
    with Pool(4) as pl: res=dict(pl.map(one,[d for d in F.values() if not d['deg']],chunksize=2))
    json.dump(res,open(sqw.OUTD+'g_matched.json','w'))
    rng=np.random.default_rng(20); L=[]; ps=[]; lines=[]
    for g in sqw.GROUPS:
        for p in range(sqw2.nP(g)):
            bn='whole' if g=='per0DD' else ['light','dark'][p]; rows=[v[p] for k,v in res.items() if k.startswith(g+'|') and v[p] is not None]
            Fr=np.array([r['frac'] for r in rows]); m=Fr.mean(0); ci=np.percentile([Fr[rng.integers(0,len(Fr),len(Fr))].mean(0) for _ in range(2000)],[2.5,97.5],axis=0)
            d=Fr[:,0]-Fr[:,1]; pw=wilcoxon(d).pvalue; ps.append(pw); DP=np.array([r['dep'] for r in rows])
            lines.append(f"  hidden(matched) {g:7s} {bn:5s} n {len(Fr):2d} (round trips/fly {np.mean([r['n'] for r in rows]):.1f}) | "+' '.join(f"{t} {m[j]:.2f} [{ci[0][j]:.2f},{ci[1][j]:.2f}]" for j,t in enumerate(TYPES))
                         +f" | most common: {TYPES[int(np.argmax(m))]} | fwd-rev {d.mean():+.3f} p {pw:.2g}\n      rates/h "+' '.join(f"{t} {np.mean([r['rate'][j] for r in rows]):.3f}" for j,t in enumerate(TYPES))
                         +f" | P(Q entry) {np.mean([r['qin'] for r in rows]):.2f} P(Q exit) {np.mean([r['qout'] for r in rows]):.2f} | joint-product forward {DP[:,0].mean():+.3f} (p {wilcoxon(DP[:,0]).pvalue:.2g})")
    ps=np.array(ps); o=np.argsort(ps); adj=np.empty(len(ps)); run=0
    for k,j in enumerate(o): run=max(run,(len(ps)-k)*ps[j]); adj[j]=min(1,run)
    L=[ln.replace('\n',f' Holm {a:.2g}\n',1) for ln,a in zip(lines,adj)]
    open('sqw7b.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
