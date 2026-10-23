# PRESPEC 20 (+ amendment A): S/Q/W per-minute HMM on count data.
# Usage: python3 sqw.py fit THETA [K3only] | route THETA | boot | lrt | mat THETA | report
import glob,json,sys,os,numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from v18 import sleep_runs, CD, DAY, HALF
from cycle import fly as cyc_fly
from sqw_common import *
GROUPS=['CS','fumin','per0LD','per0DD','sss']; WIN=240; STARTS=20
OUTD=os.path.dirname(os.path.abspath(__file__))+'/sqw_out/'; os.makedirs(OUTD,exist_ok=True)
def files(g): return sorted(glob.glob(CD+g+'/**/*.txt',recursive=True))
def wblock(g,w):   # block of window index w (4 h windows from ZT0)
    return 'whole' if g=='per0DD' else ('light' if (w%6)<3 else 'dark')
def symbols(x,g,theta):
    """1-min S/Q/W sequence (0=Q,1=W,2=S), per PRESPEC 20 section 1."""
    nm=len(x)//60; x=x[:nm*60]; sl=np.zeros(len(x),bool)
    for s,n in sleep_runs(x,True): sl[s:s+n]=True
    ssec=sl.reshape(nm,60).sum(1); act=(x>0).reshape(nm,60).mean(1)
    ev,_=cyc_fly(x); ref={}
    for bn,ph in ([('whole',None)] if g=='per0DD' else [('light',0),('dark',1)]):
        bb=[e['ref'] for e in ev if ph is None or e['wake'][0]==ph or e['onset'][0]==ph]
        ref[bn]=np.median(bb) if len(bb)>=3 else np.nan
    y=np.full(nm,W); y[ssec>=30]=S
    for m in np.flatnonzero(ssec<30):
        bn='whole' if g=='per0DD' else ('light' if (m*60)%DAY<HALF else 'dark')
        r=ref[bn]
        if np.isnan(r): y[m]=-1          # no reference: minute unclassifiable
        elif act[m]<theta*r: y[m]=Q
    return y,ref
# ---------------- stage: fit ----------------
def fit_fly(args):
    g,i,path,theta,Ks=args
    x=np.loadtxt(path).astype(int); y,ref=symbols(x,g,theta); nw=len(y)//WIN; out=[]
    for w in range(nw):
        yy=y[w*WIN:(w+1)*WIN]; r=dict(g=g,fly=i,w=w,block=wblock(g,w),ref=ref.get(wblock(g,w)),
            nS=int((yy==S).sum()),nQ=int((yy==Q).sum()),nW=int((yy==W).sum()),nchg=int(np.sum(np.diff(yy)!=0)))
        if (yy<0).any() or r['nS']==0 or r['nS']==WIN: r['fit']=False; out.append(r); continue
        r['fit']=True; seed=20+10000*GROUPS.index(g)+100*i+w; F={}
        for K in Ks:
            f=fit(yy,K,STARTS,seed)
            if f is None: continue
            F[K]=dict(ll=f[0],pi=f[1].tolist(),P=f[2].tolist(),B=f[3].tolist())
        r['F']=F
        if 3 in F:
            P,B,pi=label3(np.array(F[3]['P']),np.array(F[3]['B']),np.array(F[3]['pi']))
            r.update(P3=P.tolist(),B3=B.tolist(),pi3=pi.tolist(),deg=bool(degenerate(B)),xn=float(xn(P)),yn=float(yn(P)))
        out.append(r)
    return out
def stage_fit(theta,Ks):
    fn=OUTD+f'fit_{theta}.jsonl'; done=set()
    if os.path.exists(fn):
        for l in open(fn): d=json.loads(l); done.add((d['g'],d['fly']))
    jobs=[(g,i,p,theta,Ks) for g in GROUPS for i,p in enumerate(files(g)) if (g,i) not in done]
    with Pool(4) as pl:
        for res in pl.imap_unordered(fit_fly,jobs):
            with open(fn,'a') as fh:
                for r in res: fh.write(json.dumps(r)+'\n')
            print(res[0]['g'],res[0]['fly'],'done',flush=True)
def load_fit(theta):
    R={}
    for l in open(OUTD+f'fit_{theta}.jsonl'):
        d=json.loads(l); R.setdefault((d['g'],d['fly']),[]).append(d)
    for k in R: R[k].sort(key=lambda d:d['w'])
    return R

# ---------------- stage: route (whole-record decoding, time-varying parameters) ----------------
def fly_params(ws):
    """Per-window (pi,P,B) at K=3, label-aligned; unfit windows borrow the nearest fitted one (previous first)."""
    have=[d['w'] for d in ws if 'P3' in d]
    if not have: return None
    par=[]
    for d in ws:
        if 'P3' in d: src=d
        else:
            prev=[w for w in have if w<d['w']]; src=ws[max(prev)] if prev else ws[min(have)]
        par.append((np.array(src['pi3']),np.array(src['P3']),np.array(src['B3'])))
    return par
def emis(B,y): return np.where(y[:,None]>=0,B[:,np.clip(y,0,2)].T,1.0)
def runs_routes(z):
    """Routes in a state path: list of (route, completion step). Sleep->wake via H1: 2,0,1 / direct 2,1;
    wake->sleep via H1: 1,0,2 / direct 1,2 (indices H1=0,H2=1,H3=2)."""
    ch=np.r_[0,np.flatnonzero(np.diff(z))+1]; st=z[ch]; out=[]
    for k in range(1,len(st)):
        a,b=st[k-1],st[k]
        if (a,b)==(2,1): out.append(('N32',ch[k]))
        elif (a,b)==(1,2): out.append(('N23',ch[k]))
        elif k>=2 and st[k-1]==0:
            p=st[k-2]
            if (p,b)==(2,1): out.append(('N312',ch[k]))
            elif (p,b)==(1,2): out.append(('N213',ch[k]))
    return out
def route_fly(args):
    key,ws,seed=args; par=fly_params(ws)
    if par is None: return key,None
    y=np.concatenate([np.array(d['y']) for d in ws]); T=len(y); wi=np.arange(T)//WIN
    E=np.vstack([emis(par[w][2],y[w*WIN:(w+1)*WIN]) for w in range(len(ws))])
    # forward (scaled)
    A=np.empty((T,3)); a=par[0][0]*E[0]; A[0]=a/a.sum()
    for t in range(1,T): a=(A[t-1]@par[wi[t]][1])*E[t]; A[t]=a/a.sum()
    rng=np.random.default_rng(seed); NP=200; Z=np.empty((NP,T),int)
    c=np.cumsum(A[-1]); Z[:,-1]=np.searchsorted(c/c[-1],rng.random(NP))
    for t in range(T-2,-1,-1):
        pr=A[t][:,None]*par[wi[t+1]][1][:,Z[:,t+1]]; pr=np.cumsum(pr,0); pr/=pr[-1]
        Z[:,t]=(rng.random(NP)[None,:]>pr).sum(0)
    # Viterbi
    with np.errstate(divide='ignore'):
        lE=np.log(E); d=np.log(par[0][0])+lE[0]; bp=np.empty((T,3),int)
        for t in range(1,T):
            m=d[:,None]+np.log(par[wi[t]][1]); bp[t]=m.argmax(0); d=m.max(0)+lE[t]
    zv=np.empty(T,int); zv[-1]=d.argmax()
    for t in range(T-1,0,-1): zv[t-1]=bp[t,zv[t]]
    names=['N312','N32','N213','N23']; cnt=np.zeros((len(ws),4)); cv=np.zeros((len(ws),4))
    for z in Z:
        for r,t in runs_routes(z): cnt[t//WIN,names.index(r)]+=1
    cnt/=NP
    for r,t in runs_routes(zv): cv[t//WIN,names.index(r)]+=1
    return key,dict(ffbs=cnt.tolist(),vit=cv.tolist())
def stage_route(theta):
    R=load_fit(theta)
    # sequences are needed again; rebuild them
    jobs=[]
    for (g,i),ws in R.items():
        x=np.loadtxt(files(g)[i]).astype(int); y,_=symbols(x,g,theta)
        for d in ws: d['y']=y[d['w']*WIN:(d['w']+1)*WIN].tolist()
        jobs.append(((g,i),ws,20+1000*GROUPS.index(g)+i))
    out={}
    with Pool(4) as pl:
        for key,res in pl.imap_unordered(route_fly,jobs): out[f'{key[0]}|{key[1]}']=res
    json.dump(out,open(OUTD+f'route_{theta}.json','w'))
# ---------------- stage: boot (section 5) ----------------
def boot_one(args):
    d,seed=args; rng=np.random.default_rng(seed); pi,P,B=(np.array(d[k]) for k in ('pi3','P3','B3')); xs=[];ys=[]
    for b in range(50):
        _,yb=simulate(pi,P,B,WIN,rng); f=fit(yb,3,STARTS,seed*100+b)
        if f is None: continue
        Pb,Bb,_=label3(f[2],f[3]); xs.append(xn(Pb)); ys.append(yn(Pb))
    return dict(g=d['g'],fly=d['fly'],w=d['w'],block=d['block'],xn=d['xn'],yn=d['yn'],
                biasx=float(np.nanmean(xs)-d['xn']),biasy=float(np.nanmean(ys)-d['yn']),sdx=float(np.nanstd(xs)),sdy=float(np.nanstd(ys)))
def stage_boot():
    R=load_fit(0.5); rng=np.random.default_rng(20); jobs=[]
    for g in GROUPS:
        for bn in (['whole'] if g=='per0DD' else ['light','dark']):
            pool=[d for (gg,_),ws in sorted(R.items()) if gg==g for d in ws if 'P3' in d and d['block']==bn and not d['deg']]
            for k in rng.choice(len(pool),min(6,len(pool)),replace=False): jobs.append((pool[k],20+len(jobs)))
    with Pool(4) as pl: res=pl.map(boot_one,jobs,chunksize=1)
    json.dump(res,open(OUTD+'boot.json','w'))
# ---------------- stage: lrt (amendment A; null = the smaller model) ----------------
def lrt_one(args):
    d,seed=args; F={int(k):v for k,v in d['F'].items()}; b2,b3,b4=(bic(F[K]['ll'],K,WIN) for K in (2,3,4))
    Kc=2 if b2<b4 else 4; K0,K1=min(3,Kc),max(3,Kc); lr0=2*(F[K1]['ll']-F[K0]['ll'])
    rng=np.random.default_rng(seed); n0=F[K0]; ge=0; lrs=[]
    for b in range(100):
        _,yb=simulate(np.array(n0['pi']),np.array(n0['P']),np.array(n0['B']),WIN,rng)
        f0=fit(yb,K0,10,seed*1000+b); f1=fit(yb,K1,10,seed*1000+b)
        if f0 is None or f1 is None: continue
        lr=2*max(0.0,f1[0]-f0[0]); lrs.append(lr); ge+=lr>=lr0
    return dict(g=d['g'],fly=d['fly'],w=d['w'],block=d['block'],Kc=Kc,dbic=float(min(b2,b4)-b3),LR=float(lr0),
                p=(1+ge)/(1+len(lrs)),nrep=len(lrs))
def stage_lrt():
    R=load_fit(0.5); rng=np.random.default_rng(20); jobs=[]
    for g in GROUPS:
        for bn in (['whole'] if g=='per0DD' else ['light','dark']):
            pool=[]
            for (gg,_),ws in sorted(R.items()):
                if gg!=g: continue
                for d in ws:
                    if d['block']!=bn or not d['fit'] or len(d['F'])<3: continue
                    F={int(k):v for k,v in d['F'].items()}; b=[bic(F[K]['ll'],K,WIN) for K in (2,3,4)]
                    if min(b[0],b[2])<b[1]-10: pool.append(d)
            for k in rng.choice(len(pool),min(20,len(pool)),replace=False) if pool else []: jobs.append((pool[k],20+len(jobs)))
    print(len(jobs),'flagged windows tested',flush=True)
    with Pool(4) as pl: res=pl.map(lrt_one,jobs,chunksize=1)
    json.dump(res,open(OUTD+'lrt.json','w'))
# ---------------- stage: mat export ----------------
def stage_mat(theta):
    import scipy.io as sio
    R=load_fit(theta)
    for g in GROUPS:
        keys=sorted(k for k in R if k[0]==g); nw=max(len(R[k]) for k in keys)
        Rc=np.empty((len(keys),nw),object); Cc=np.empty((3*len(keys),nw),object)
        Rc[:]=[[np.zeros((0,0))]*nw]; Cc[:]=[[np.zeros((0,0))]*nw]
        for i,k in enumerate(keys):
            x=np.loadtxt(files(g)[k[1]]).astype(int); y,_=symbols(x,g,theta)
            for d in R[k]:
                w=d['w']; Rc[i,w]=(y[w*WIN:(w+1)*WIN]+1).astype(float)[None,:]
                if 'P3' in d: Cc[3*i,w]=np.array([[3.0]]); Cc[3*i+1,w]=np.array(d['P3']); Cc[3*i+2,w]=np.array(d['B3'])
        sio.savemat(OUTD+f'SQW_{g}_theta{theta}.mat',{'R':Rc,'C':Cc,'symbols':'1=Q 2=W 3=S; states H1,H2,H3'})
if __name__=='__main__':
    st=sys.argv[1]
    if st=='fit': stage_fit(float(sys.argv[2]),(3,) if 'K3only' in sys.argv else (2,3,4))
    elif st=='route': stage_route(float(sys.argv[2]))
    elif st=='boot': stage_boot()
    elif st=='lrt': stage_lrt()
    elif st=='mat': stage_mat(float(sys.argv[2]))
