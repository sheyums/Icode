# PRESPEC 11: flat count HMM K=2..6 per fly x 4 h ZT bin, days as sequences.
import os; os.environ['OMP_NUM_THREADS']='1'
import sys, json, glob; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
import numpy as np
from multiprocessing import Pool
from hmmlearn.hmm import CategoricalHMM
CD='/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/counts/'
GEN=['CS','fumin']; DAY=86400; BIN=14400
def load(g): return [np.loadtxt(f).astype(int) for f in sorted(glob.glob(CD+g+'/**/*.txt',recursive=True))]
def segs(x,b):
    out=[]
    for d in range(int(np.ceil(len(x)/DAY))):
        s=x[d*DAY+b*BIN:d*DAY+(b+1)*BIN]
        if len(s)>=0.9*BIN: out.append(np.minimum(s,3))
    return out
def hfit(S,K,seed):
    m=CategoricalHMM(n_components=K,n_features=4,n_iter=500,tol=0.01,random_state=seed,init_params='ste',params='ste')
    m.fit(np.concatenate(S).reshape(-1,1),[len(s) for s in S]); return m
def score(m,S): return m.score(np.concatenate(S).reshape(-1,1),[len(s) for s in S])
def stage1(a):
    g,fi,b,x=a; S=segs(x,b); out={'g':g,'fly':fi,'bin':b,'ndays':len(S),'cv':{}}
    for K in range(2,7):
        sc=[]
        for d in range(len(S)):
            tr=[s for k,s in enumerate(S) if k!=d]
            try:
                m=max((hfit(tr,K,r) for r in range(3)),key=lambda m:score(m,tr)); sc.append(score(m,[S[d]])/len(S[d]))
            except Exception as e: sc.append(float('nan')); out.setdefault('errors',[]).append(f'K{K} d{d}: {e}')
        out['cv'][K]=sc
    return out
def pcca(P,m):
    w,V=np.linalg.eig(P); o=np.argsort(-np.abs(w)); w=w[o]; V=V[:,o]
    cplx=bool(np.any(np.abs(w[:m+1].imag)>1e-9)); X=V[:,:m+1].real
    X[:,0]=1.0
    rep=[int(np.argmax(np.linalg.norm(X,axis=1)))]
    for _ in range(m):
        Q,_=np.linalg.qr(X[rep].T); R=X-(X@Q)@Q.T; rep.append(int(np.argmax(np.linalg.norm(R,axis=1))))
    A=np.linalg.inv(X[rep]); chi=X@A; return chi,cplx
def stationary(P):
    w,V=np.linalg.eig(P.T); v=np.abs(V[:,np.argmin(np.abs(w-1))].real); return v/v.sum()
def committor(P,A,B,avoid=()):
    n=len(P); q=np.zeros(n); q[list(B)]=1; F=[i for i in range(n) if i not in A and i not in B and i not in avoid]
    if F:
        M=np.eye(len(F))-P[np.ix_(F,F)]; rhs=P[np.ix_(F,list(B))].sum(1); q[F]=np.linalg.solve(M,rhs)
    return q
def via_T(P,pi,A,B,T):
    q=committor(P,A,B); h=committor(P,A,B,avoid=T)
    tot=sum(pi[a]*P[a,j]*q[j] for a in A for j in range(len(P)) if j not in A)
    dirc=sum(pi[a]*P[a,j]*h[j] for a in A for j in range(len(P)) if j not in A and j not in T)
    return 1-dirc/tot if tot>0 else float('nan')
def describe(m):
    P=m.transmat_; B=m.emissionprob_; pi=stationary(P); K=len(P)
    w=np.linalg.eigvals(P); w=w[np.argsort(-np.abs(w))][1:]; tau=[(-1/np.log(abs(x)) if 0<abs(x)<1 else float('inf')) for x in w]
    d={'tau':tau,'m':{th:int(sum(t>=th for t in tau)) for th in (30,60,300)},'B':B.round(4).tolist(),'pi':pi.round(5).tolist(),'dwell':(1/(1-np.diag(P))).round(1).tolist()}
    mm=d['m'][60]; mf=1-B[:,0]
    if 1<=mm<=K-1:
        chi,cplx=pcca(P,mm); D=np.diag(pi); Pc=np.linalg.solve(chi.T@D@chi, chi.T@D@P@chi)
        macro_mf=[(pi*chi[:,I]*mf).sum()/(pi*chi[:,I]).sum() for I in range(mm+1)]; occ=[(pi*chi[:,I]).sum() for I in range(mm+1)]
        order=np.argsort(macro_mf); lab={int(order[0]):'S',int(order[-1]):'W'}
        d['macro']={'mf':[macro_mf[i] for i in order],'occ':[occ[i] for i in order],'dwell':[1/(1-Pc[i,i]) if Pc[i,i]<1 else float('inf') for i in order],'complex_flag':cplx}
        if mm==2:
            sets={L:[k for k in range(K) if chi[k,I]>0.5] for I,L in lab.items()}; Tm=int(order[1]); sets['T']=[k for k in range(K) if chi[k,Tm]>0.5]
            if sets['S'] and sets['W']:
                d['fracT_SW']=via_T(P,pi,set(sets['S']),set(sets['W']),set(sets['T'])); d['fracT_WS']=via_T(P,pi,set(sets['W']),set(sets['S']),set(sets['T']))
    return d
def stage2(a):
    g,fi,b,x,K=a; S=segs(x,b); out={'g':g,'fly':fi,'bin':b,'K':K}
    fits=[]
    for r in range(10):
        try: m=hfit(S,K,r); fits.append((score(m,S),m))
        except Exception as e: out.setdefault('errors',[]).append(str(e))
    fits.sort(key=lambda t:-t[0]); best=fits[0][0]; near=[f for f in fits if f[0]>=best-0.5]
    out['logL']=best; out['n_near']=len(near); out['best']=describe(near[0][1])
    fr=[describe(f[1]) for f in near[1:]]
    for key in ('fracT_SW','fracT_WS'):
        vals=[out['best'].get(key)]+[d.get(key) for d in fr]; vals=[v for v in vals if v is not None and np.isfinite(v)]
        if vals: out[key+'_range']=max(vals)-min(vals)
    return out
def run(fn,jobs,path,key):
    done=set()
    if os.path.exists(path): done={tuple(json.loads(l)[k] for k in key) for l in open(path)}
    todo=[j for j in jobs if tuple(j[i] for i in range(len(key))) not in done]
    print(path,'remaining',len(todo),flush=True)
    with Pool(4) as pl:
        for r in pl.imap_unordered(fn,todo):
            with open(path,'a') as f: f.write(json.dumps(r,default=lambda o: float(o) if not isinstance(o,(list,dict)) else o)+'\n')
    return [json.loads(l) for l in open(path)]
if __name__=='__main__':
    data={g:load(g) for g in GEN}
    sub={g:sorted(np.random.default_rng(11).choice(len(data[g]),10,replace=False).tolist()) for g in GEN}
    print('CV subset (0-based file index):',sub,flush=True)
    R1=run(stage1,[(g,fi,b,data[g][fi]) for g in GEN for fi in sub[g] for b in range(6)],'flatK_stage1.jsonl',('g','fly','bin'))
    Kg={}
    for g in GEN:
        for b in range(6):
            ks=[]
            for r in [r for r in R1 if r['g']==g and r['bin']==b]:
                sc={int(K):np.nanmean(v) for K,v in r['cv'].items()}; ks.append(max(sc,key=sc.get))
            vals,cnt=np.unique(ks,return_counts=True); Kg[f'{g}_{b}']=int(vals[np.flatnonzero(cnt==cnt.max())[0]])
    json.dump(Kg,open('flatK_Kg.json','w')); print('K_g:',Kg,flush=True)
    R2=run(stage2,[(g,fi,b,data[g][fi],Kg[f'{g}_{b}']) for g in GEN for fi in range(len(data[g])) for b in range(6)],'flatK_stage2.jsonl',('g','fly','bin'))
    print('stage 2 done',len(R2),flush=True)
