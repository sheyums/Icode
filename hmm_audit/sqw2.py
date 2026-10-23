# PRESPEC 20b: S/Q/W HMM per fly, shared emissions, light/dark transition matrices.
# Usage: python3 sqw2.py fit THETA [K3only] | route THETA | boot | split | adequacy | mat THETA
import sys,os,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
from sqw2_common import *
from sqw2_common import _E
from sqw_common import bic
GROUPS=sqw.GROUPS; ST=20; OUTD=sqw.OUTD; NM=['N312','N32','N213','N23']
def nP(g): return 1 if g=='per0DD' else 2
def seqs(g,i,theta):
    x=np.loadtxt(sqw.files(g)[i]).astype(int); y,ref=sqw.symbols(x,g,theta); y=y[:(len(y)//1440)*1440]
    ph=np.zeros(len(y),int) if g=='per0DD' else phases(len(y)); return y,ph,ref
def fit_fly(a):
    g,i,theta,Ks=a; y,ph,ref=seqs(g,i,theta); r=dict(g=g,fly=i,T=len(y),ref=ref,nchg=int(np.sum(np.diff(y)!=0)),
        frac={s:float(np.mean(y==v)) for s,v in (('Q',0),('W',1),('S',2))},F={})
    for K in Ks:
        ll,pi,P,B,it=em(y,ph,nP(g),K,ST,20+1000*GROUPS.index(g)+i); r['F'][K]=dict(ll=ll,pi=pi.tolist(),P=P.tolist(),B=B.tolist(),it=it)
    if 3 in r['F']:
        F=r['F'][3]; P,B,pi=label3b(np.array(F['P']),np.array(F['B']),np.array(F['pi']))
        r.update(P3=P.tolist(),B3=B.tolist(),pi3=pi.tolist(),deg=bool(degenerate(B)),
                 xn=[float(xn(P[p])) for p in range(nP(g))],yn=[float(yn(P[p])) for p in range(nP(g))])
    return r
def stage_fit(theta,Ks):
    fn=OUTD+f'fit2_{theta}.jsonl'; done=set()
    if os.path.exists(fn): done={(d['g'],d['fly']) for d in map(json.loads,open(fn))}
    jobs=[(g,i,theta,Ks) for g in GROUPS for i in range(len(sqw.files(g))) if (g,i) not in done]
    with Pool(4) as pl:
        for r in pl.imap_unordered(fit_fly,jobs):
            open(fn,'a').write(json.dumps(r)+'\n'); print(r['g'],r['fly'],'done',flush=True)
def load(theta):
    return {(d['g'],d['fly']):d for d in map(json.loads,open(OUTD+f'fit2_{theta}.jsonl'))}
# ---- route counting: FFBS (200 paths) and Viterbi, attributed to the phase of the completion minute ----
def pair_counts(z,ph,nph):
    """Single-step counts N_ij on the run-compressed path, attributed to the phase of the arrival minute."""
    ch=np.flatnonzero(np.diff(z))+1; c=np.zeros((nph,3,3))
    np.add.at(c,(ph[ch],z[ch-1],z[ch]),1); return c
def decode(pi,P,B,y,ph,seed,NP=200,pairs=False):
    T=len(y); E=_E(B[None],y)[:,0]; A=np.empty((T,3)); a=pi*E[0]; A[0]=a/a.sum()
    for t in range(1,T): a=(A[t-1]@P[ph[t]])*E[t]; A[t]=a/a.sum()
    rng=np.random.default_rng(seed); Z=np.empty((NP,T),int); c=np.cumsum(A[-1]); Z[:,-1]=np.searchsorted(c/c[-1],rng.random(NP))
    for t in range(T-2,-1,-1):
        pr=np.cumsum(A[t][:,None]*P[ph[t+1]][:,Z[:,t+1]],0); pr/=pr[-1]; Z[:,t]=(rng.random(NP)[None,:]>pr).sum(0)
    with np.errstate(divide='ignore'):
        lE=np.log(E); d=np.log(pi)+lE[0]; bp=np.empty((T,3),int); lP=np.log(P)
        for t in range(1,T): m=d[:,None]+lP[ph[t]]; bp[t]=m.argmax(0); d=m.max(0)+lE[t]
    zv=np.empty(T,int); zv[-1]=d.argmax()
    for t in range(T-1,0,-1): zv[t-1]=bp[t,zv[t]]
    nph=ph.max()+1; cf=np.zeros((nph,4)); cv=np.zeros((nph,4))
    rp=np.zeros((NP,nph,4))                   # per-path route counts (kept so nonlinear statistics are
    for k,z in enumerate(Z):                  # averaged per path, E[f(N)], not computed on mean counts, f(E[N]))
        for rt,t in sqw.runs_routes(z): rp[k,ph[t],NM.index(rt)]+=1
    cf=rp.sum(0)
    for rt,t in sqw.runs_routes(zv): cv[ph[t],NM.index(rt)]+=1
    if pairs:
        pp=np.stack([pair_counts(z,ph,nph) for z in Z]); pf=pp.mean(0); pv=pair_counts(zv,ph,nph)
        win=np.arange(T)//240; nw=T//240     # same counts, attributed to the 4 h clock window of the arrival minute
        pw=sum(pair_counts(z,win,nw) for z in Z)/NP
        return (cf/NP).tolist(),cv.tolist(),pf.tolist(),pv.tolist(),pw.tolist(),rp.astype(int).tolist(),pp.astype(int).tolist()
    return (cf/NP).tolist(),cv.tolist()
def route_one(a):
    d,theta=a; g,i=d['g'],d['fly']; y,ph,_=seqs(g,i,theta)
    f,v,pf,pv,pw,rp,pp=decode(np.array(d['pi3']),np.array(d['P3']),np.array(d['B3']),y,ph,20+1000*GROUPS.index(g)+i,pairs=True)
    hours=[float(np.sum(ph==p))/60 for p in range(ph.max()+1)]
    return f'{g}|{i}',dict(ffbs=f,vit=v,pairs_ffbs=pf,pairs_vit=pv,pairs_win=pw,hours=hours,routes_path=rp,pairs_path=pp)
def stage_route(theta):
    R=load(theta)
    with Pool(4) as pl: out=dict(pl.map(route_one,[(d,theta) for d in R.values() if 'P3' in d]))
    json.dump(out,open(OUTD+f'route2_{theta}.json','w'))
# ---- parametric bootstrap (section 3) ----
def boot_one(a):
    d,seed=a; g=d['g']; rng=np.random.default_rng(seed); T=d['T']
    ph=np.zeros(T,int) if g=='per0DD' else phases(T); P=np.array(d['P3']); Pt=P[ph]; xs=[];ys=[]
    for b in range(30):
        _,yb=simulate2(np.array(d['pi3']),Pt,np.array(d['B3']),rng); ll,pi,Pb,Bb,_=em(yb,ph,nP(g),3,10,seed*100+b)
        Pb,Bb,_=label3b(Pb,Bb,pi); xs.append([xn(Pb[p]) for p in range(nP(g))]); ys.append([yn(Pb[p]) for p in range(nP(g))])
    xs=np.array(xs); ys=np.array(ys)
    return dict(g=g,fly=d['fly'],xn=d['xn'],yn=d['yn'],biasx=(np.nanmean(xs,0)-d['xn']).tolist(),biasy=(np.nanmean(ys,0)-d['yn']).tolist(),
                sdx=np.nanstd(xs,0).tolist(),sdy=np.nanstd(ys,0).tolist())
def stage_boot():
    R=load(0.5); rng=np.random.default_rng(20); jobs=[]
    for g in GROUPS:
        pool=[d for k,d in sorted(R.items()) if k[0]==g and 'P3' in d and not d['deg']]
        for k in rng.choice(len(pool),min(5,len(pool)),replace=False): jobs.append((pool[k],20+len(jobs)))
    with Pool(4) as pl: json.dump(pl.map(boot_one,jobs,chunksize=1),open(OUTD+'boot2.json','w'))
# ---- day split (section 5) ----
def split_one(a):
    g,i=a; y,ph,_=seqs(g,i,0.5); days=len(y)//1440; h=(days//2)*1440; out={}
    for nm,(yy,pp) in (('first',(y[:h],ph[:h])),('second',(y[h:],ph[h:]))):
        ll,pi,P,B,_=em(yy,pp,nP(g),3,ST,20+i); P,B,pi=label3b(P,B,pi)
        out[nm]=dict(deg=bool(degenerate(B)),xn=[float(xn(P[p])) for p in range(nP(g))],yn=[float(yn(P[p])) for p in range(nP(g))])
    return f'{g}|{i}',out
def stage_split():
    with Pool(4) as pl: json.dump(dict(pl.map(split_one,[(g,i) for g in GROUPS for i in range(len(sqw.files(g)))])),open(OUTD+'split2.json','w'))
# ---- within-phase adequacy (section 5) ----
def stat4h(y):
    W4=y[:(len(y)//240)*240].reshape(-1,240)
    return np.stack([(W4==0).mean(1),(W4==1).mean(1),(W4==2).mean(1),(np.diff(W4,axis=1)!=0).sum(1)],1)
def adeq_one(a):
    d=a; g,i=d['g'],d['fly']; y,ph,_=seqs(g,i,0.5); rng=np.random.default_rng(20+i); P=np.array(d['P3'])[ph]
    sims=np.stack([stat4h(simulate2(np.array(d['pi3']),P,np.array(d['B3']),rng)[1]) for _ in range(100)])
    lo,hi=np.percentile(sims,[2.5,97.5],axis=0); obs=stat4h(y); out=(obs<lo)|(obs>hi)
    return f'{g}|{i}',out.astype(int).tolist()
def stage_adequacy():
    R=load(0.5)
    with Pool(4) as pl: json.dump(dict(pl.map(adeq_one,[d for d in R.values() if 'P3' in d])),open(OUTD+'adeq2.json','w'))
# ---- export ----
def stage_mat(theta):
    import scipy.io as sio
    R=load(theta)
    for g in GROUPS:
        keys=sorted(k for k in R if k[0]==g); nw=max(R[k]['T'] for k in keys)//240
        Rc=np.empty((len(keys),nw),object); Cc=np.empty((3*len(keys),nw),object)
        for a in (Rc,Cc):
            for idx in np.ndindex(a.shape): a[idx]=np.zeros((0,0))
        for j,k in enumerate(keys):
            d=R[k]; y,ph,_=seqs(g,k[1],theta)
            for w in range(len(y)//240):
                Rc[j,w]=(y[w*240:(w+1)*240]+1).astype(float)[None,:]
                if 'P3' in d: Cc[3*j,w]=np.array([[3.0]]); Cc[3*j+1,w]=np.array(d['P3'][ph[w*240]]); Cc[3*j+2,w]=np.array(d['B3'])
        sio.savemat(OUTD+f'SQW2_{g}_theta{theta}.mat',{'R':Rc,'C':Cc,'symbols':'1=Q 2=W 3=S; states H1,H2,H3; P per window = its phase'})
if __name__=='__main__':
    st=sys.argv[1]
    if st=='fit': stage_fit(float(sys.argv[2]),(3,) if 'K3only' in sys.argv else (2,3,4))
    elif st=='route': stage_route(float(sys.argv[2]))
    elif st=='boot': stage_boot()
    elif st=='split': stage_split()
    elif st=='adequacy': stage_adequacy()
    elif st=='mat': stage_mat(float(sys.argv[2]))
