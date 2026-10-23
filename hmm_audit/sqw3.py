# PRESPEC 20c: 6 h transition blocks (secondary design). Reuses sqw2's stages with 4 phases.
# Usage: python3 sqw3.py synth | realistic | fit 0.5 | route 0.5
import sys,os,json,time,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw, sqw2
from sqw2_common import *
GROUPS=['CS','per0LD','fumin']; NB=4; T=5760
def phases6(T): return (np.arange(T)%1440)//360
def seqs6(g,i,theta):
    x=np.loadtxt(sqw.files(g)[i]).astype(int); y,ref=sqw.symbols(x,g,theta); y=y[:(len(y)//1440)*1440]; return y,phases6(len(y)),ref
sqw2.GROUPS=GROUPS; sqw2.seqs=seqs6; sqw2.nP=lambda g:NB
sqw2.OUTD=sqw.OUTD+'b6/'; os.makedirs(sqw2.OUTD,exist_ok=True)
def kpar6(K): return (K-1)+NB*K*(K-1)+2*K
def fit_eval(y,ph,seed):
    F={K:em(y,ph,NB,K,20,seed) for K in (2,3,4)}; b={K:-2*F[K][0]+kpar6(K)*np.log(len(y)) for K in F}
    P,B,pi=label3b(F[3][2],F[3][3],F[3][1])
    return min(b,key=b.get),[float(xn(P[p])) for p in range(NB)],[float(yn(P[p])) for p in range(NB)],bool(degenerate(B))
# ---- step 1: V1c, V2c ----
def synth_one(a):
    from sqw2_synth import PL1,PD1,B1,PL2,PD2,B2
    case,f=a; rng=np.random.default_rng(5000+100*(case=='V2c')+f); ph=phases6(T)
    Pb=np.stack([PL1,PL1,PD1,PD1]) if case=='V1c' else np.stack([PL2,PL2,PD2,PD2]); B=B1 if case=='V1c' else B2
    _,y=simulate2(stationary(Pb[0]),Pb[ph],B,rng); k,x,yy,dg=fit_eval(y,ph,30+f); return dict(case=case,best=k,xn=x,yn=yy,deg=dg)
def stage_synth():
    from sqw2_synth import TRUTH
    t0=time.time()
    with Pool(4) as p: res=p.map(synth_one,[(c,f) for c in ('V1c','V2c') for f in range(20)],chunksize=1)
    json.dump(res,open('sqw3_synth.json','w')); L=[]
    tx=[TRUTH['light'][0]]*2+[TRUTH['dark'][0]]*2; ty=[TRUTH['light'][1]]*2+[TRUTH['dark'][1]]*2
    for c in ('V1c','V2c'):
        R=[r for r in res if r['case']==c]; fr={K:np.mean([r['best']==K for r in R]) for K in (2,3,4)}
        X=np.array([r['xn'] for r in R]); Y=np.array([r['yn'] for r in R])
        L.append(f"{c}: best K {fr}; H1 degenerate {np.mean([r['deg'] for r in R]):.2f}")
        L.append('   Xn mean '+' '.join(f'{v:+.3f}' for v in np.nanmean(X,0))+' SD '+' '.join(f'{v:.3f}' for v in np.nanstd(X,0)))
        L.append('   Yn mean '+' '.join(f'{v:+.3f}' for v in np.nanmean(Y,0))+' SD '+' '.join(f'{v:.3f}' for v in np.nanstd(Y,0)))
        if c=='V1c':
            ok=fr[3]>=0.7 and np.all(np.abs(np.nanmean(X,0)-tx)<0.1) and np.all(np.abs(np.nanmean(Y,0)-ty)<0.1) and np.all(np.nanstd(X,0)<0.3) and np.all(np.nanstd(Y,0)<0.3)
            L.append(f'   V1c PASS={ok}')
        else: L.append(f"   V2c PASS={fr[2]>=0.8}")
    L.append(f'elapsed {time.time()-t0:.0f} s'); open('sqw3_synth.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
# ---- step 2: realistic truths from the 20b fits ----
def real_one(a):
    d,r=a; rng=np.random.default_rng(6000+r+100*d['fly']); P2=np.array(d['P3']); Pb=np.stack([P2[0],P2[0],P2[1],P2[1]])
    ph=phases6(d['T']); _,y=simulate2(np.array(d['pi3']),Pb[ph],np.array(d['B3']),rng); k,x,yy,dg=fit_eval(y,ph,40+r)
    return dict(g=d['g'],fly=d['fly'],r=r,best=k,xn=x,yn=yy,deg=dg,tx=[xn(P2[0])]*2+[xn(P2[1])]*2,ty=[yn(P2[0])]*2+[yn(P2[1])]*2)
def stage_realistic():
    R={(d['g'],d['fly']):d for d in map(json.loads,open(sqw.OUTD+'fit2_0.5.jsonl'))}; rng=np.random.default_rng(21); jobs=[]
    for g in GROUPS:
        pool=[d for k,d in sorted(R.items()) if k[0]==g and 'P3' in d and not d['deg']]
        for k in rng.choice(len(pool),min(5,len(pool)),replace=False): jobs+=[(pool[k],r) for r in range(20)]
    with Pool(4) as p: res=p.map(real_one,jobs,chunksize=1)
    json.dump(res,open(sqw2.OUTD+'realistic.json','w'),default=float); L=[]
    for g in GROUPS:
        sds=[]
        for f in sorted({r['fly'] for r in res if r['g']==g}):
            Rf=[r for r in res if r['g']==g and r['fly']==f]; X=np.array([r['xn'] for r in Rf]); Y=np.array([r['yn'] for r in Rf])
            sx,sy=np.nanstd(X,0),np.nanstd(Y,0); sds.append(np.r_[sx,sy])
            L.append(f"{g} fly {f}: K=3 chosen {np.mean([r['best']==3 for r in Rf]):.2f}; Xn SD "+' '.join(f'{v:.2f}' for v in sx)+'; Yn SD '+' '.join(f'{v:.2f}' for v in sy)
                     +'; Xn bias '+' '.join(f'{v:+.2f}' for v in np.nanmean(X,0)-Rf[0]['tx'])+'; Yn bias '+' '.join(f'{v:+.2f}' for v in np.nanmean(Y,0)-Rf[0]['ty']))
        med=np.median(np.array(sds),0); L.append(f"{g}: median SD over flies (Xn blocks 1-4, Yn blocks 1-4) "+' '.join(f'{v:.2f}' for v in med)+f" -> USABLE={bool(np.all(med<0.3))}")
    open('sqw3_realistic.log','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
if __name__=='__main__':
    st=sys.argv[1]
    if st=='synth': stage_synth()
    elif st=='realistic': stage_realistic()
    elif st=='fit': sqw2.stage_fit(float(sys.argv[2]),(2,3,4))
    elif st=='route': sqw2.stage_route(float(sys.argv[2]))
