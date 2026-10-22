# PRESPEC 18: bridged scoring and +-5 min window; versions V0, V1, V2; figures for V2.
import glob,json,sys,numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.stats import binomtest
CD='/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/counts/'
GROUPS=['CS','fumin','per0LD','per0DD','sss']; DAY=86400; HALF=43200
VERS={'V0':(False,60,600,60),'V1b':(True,60,600,60),'V2b':(True,60,300,30)}   # bridge, start s, end s, bin s
def blocks(g): return [('whole',None)] if g=='per0DD' else [('light',0),('dark',1)]
def sleep_runs(x,bridge):
    a=(x>0).astype(int)
    if bridge:   # option (b): bridge a 1 s active run only if BOTH flanking inactive runs are >= 300 s
        z0=1-a; ch=np.r_[0,np.flatnonzero(np.diff(z0))+1,len(z0)]; runs=[(s,e-s,z0[s]) for s,e in zip(ch[:-1],ch[1:])]
        for i in range(1,len(runs)-1):
            s_,n_,zz=runs[i]
            if zz==0 and n_==1 and runs[i-1][2]==1 and runs[i+1][2]==1 and runs[i-1][1]>=300 and runs[i+1][1]>=300: a[s_]=0
    z=1-a; ch=np.r_[0,np.flatnonzero(np.diff(z))+1,len(z)]
    return [(s,b-s) for s,b in zip(ch[:-1],ch[1:]) if z[s]==1 and b-s>=300]
def fly_bouts(args):
    x,bridge,s0,s1,bw=args; mv=(x>0).astype(float); sl=sleep_runs(x,bridge); out=[]
    nb=(s1-s0)//bw
    for (a1,n1),(a2,_) in zip(sl[:-1],sl[1:]):
        w0=a1+n1; w1=a2
        if w1-w0<1200: continue
        mid=(w0+w1)//2; ref=mv[mid-150:mid+150].mean()
        ex=np.array([mv[w0+s0+bw*k:w0+s0+bw*(k+1)].mean() for k in range(nb)])
        ap=np.array([mv[w1-s0-bw*(k+1):w1-s0-bw*k].mean() for k in range(nb)])
        out.append((int((w0%DAY)>=HALF),ref,ex,ap))
    return out
if __name__=='__main__':
    rng=np.random.default_rng(18); R={}
    data={g:[np.loadtxt(p).astype(int) for p in sorted(glob.glob(CD+g+'/**/*.txt',recursive=True))] for g in GROUPS}
    for v,(br,s0,s1,bw) in VERS.items():
        for g in GROUPS:
            with Pool(4) as pl: B=pl.map(fly_bouts,[(x,br,s0,s1,bw) for x in data[g]])
            for bn,ph in blocks(g):
                fl=[]
                for bs in B:
                    bb=[b for b in bs if ph is None or b[0]==ph]
                    if len(bb)>=3: fl.append(bb)
                D=np.array([np.mean([np.mean(b[2]-b[3]) for b in bb]) for bb in fl])
                obs=np.mean(D); null=[]
                for _ in range(2000):
                    null.append(np.mean([np.mean([np.mean(b[2]-b[3])*(1 if rng.random()<.5 else -1) for b in bb]) for bb in fl]))
                pp=(1+np.sum(np.abs(null)>=abs(obs)))/2001
                xin=np.array([np.mean([b[1]-b[3].mean() for b in bb]) for bb in fl]); yout=np.array([np.mean([b[1]-b[2].mean() for b in bb]) for bb in fl])
                EX=np.array([np.mean([b[2] for b in bb],0) for bb in fl]); AP=np.array([np.mean([b[3] for b in bb],0) for bb in fl])
                R[f'{v}|{g}|{bn}']=dict(n=len(fl),bouts=sum(map(len,fl)),D=D.tolist(),sign_p=binomtest(int((D>0).sum()),len(D)).pvalue,perm_p=pp,
                    x=xin.tolist(),y=yout.tolist(),EX=EX.tolist(),AP=AP.tolist())
                r=R[f'{v}|{g}|{bn}']
                print(f"{v} {g:7s} {bn:5s} flies {r['n']:2d} bouts {r['bouts']:4d} | median D {100*np.median(D):+.2f} pp, D>0 {np.mean(D>0):.2f}, sign p {r['sign_p']:.2g}, perm p {pp:.4f} | x {100*xin.mean():+.2f} y {100*yout.mean():+.2f}",flush=True)
    json.dump(R,open('v18b.json','w'))
