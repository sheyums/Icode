# PRESPEC 15 E2/E3 on count data.
import glob,numpy as np,itertools
from multiprocessing import Pool
CD='/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/counts/'
GEN=['CS','fumin']; DAY=86400; HALF=43200
def load(g): return [np.loadtxt(f).astype(int) for f in sorted(glob.glob(CD+g+'/**/*.txt',recursive=True))]
def zruns(x):
    z=(x==0).astype(int); ch=np.r_[0,np.flatnonzero(np.diff(z))+1,len(z)]; return [(a,b-a) for a,b in zip(ch[:-1],ch[1:]) if z[a]==1]
# ---- E2
def e2_fly(x):
    mv=(x>0).astype(float); sl=[(a,n) for a,n in zruns(x) if n>=300]; out=[]
    for (a1,n1),(a2,_) in zip(sl[:-1],sl[1:]):
        w0=a1+n1; w1=a2
        if w1-w0<1200: continue
        ex=np.array([mv[w0+60*k:w0+60*(k+1)].mean() for k in range(1,10)])
        ap=np.array([mv[w1-60*(k+1):w1-60*k].mean() for k in range(1,10)])
        out.append((int((w0%DAY)>=HALF),ex,ap))
    return out
# ---- E3
def runsym(seg):
    m=(seg>0).astype(int); ch=np.r_[0,np.flatnonzero(np.diff(m))+1,len(m)]; s=[]
    for a,b in zip(ch[:-1],ch[1:]):
        n=b-a
        if m[a]==0: s.append(0 if n<10 else (1 if n<300 else 2))
        else: s.append(3 if n<=2 else 4)
    return s
TRI={t:i for i,t in enumerate(itertools.product(range(5),repeat=3))}
def tri_counts(s):
    f=np.zeros(125); r=np.zeros(125)
    for i in range(len(s)-2):
        t=(s[i],s[i+1],s[i+2]); f[TRI[t]]+=1; r[TRI[t[::-1]]]+=1
    return f,r
def kl(f,r): a=f+.5; b=r+.5; a/=a.sum(); b/=b.sum(); return float((a*np.log(a/b)).sum())
if __name__=='__main__':
    rng=np.random.default_rng(15)
    for g in GEN:
        X=load(g)
        print(f"\n===== {g} ({len(X)} flies) =====",flush=True)
        with Pool(4) as p: E=p.map(e2_fly,X)
        for ph,pn in [(0,'light'),(1,'dark ')]:
            bouts=[[(ex,ap) for q,ex,ap in E[i] if q==ph] for i in range(len(X))]; keep=[b for b in bouts if len(b)>=3]
            def stat(sw):
                Ds=[]
                for b,s in zip(keep,sw):
                    d=[(ap-ex) if f else (ex-ap) for (ex,ap),f in zip(b,s)]; Ds.append(np.mean(d))
                return np.mean(Ds)
            obs=stat([np.zeros(len(b),bool) for b in keep])
            null=np.array([stat([rng.random(len(b))<.5 for b in keep]) for _ in range(2000)])
            pv=(1+np.sum(np.abs(null)>=abs(obs)))/2001
            EX=np.mean([np.mean([ex for ex,ap in b],0) for b in keep],0); AP=np.mean([np.mean([ap for ex,ap in b],0) for b in keep],0)
            print(f"E2 {pn} flies {len(keep)} bouts {sum(map(len,keep))} | mean D (exit - approach) {obs:+.4f}  perm p {pv:.4f}")
            print(f"      exit profile min2-10     {np.round(EX,3).tolist()}\n      approach profile min2-10 {np.round(AP,3).tolist()}",flush=True)
        for ph,pn in [(0,'light'),(1,'dark ')]:
            segs=[]
            for x in X:
                for d in range(int(np.ceil(len(x)/DAY))):
                    for b in range(3):
                        s=x[d*DAY+ph*HALF+b*14400:d*DAY+ph*HALF+(b+1)*14400]
                        if len(s)>=0.9*14400: segs.append(tri_counts(runsym(s)))
            F=np.array([s[0] for s in segs]); R=np.array([s[1] for s in segs])
            obs=kl(F.sum(0),R.sum(0)); null=[]
            for _ in range(1000):
                fl=rng.random(len(segs))<.5; null.append(kl(np.where(fl[:,None],R,F).sum(0),np.where(fl[:,None],F,R).sum(0)))
            null=np.array(null); pv=(1+np.sum(null>=obs))/1001
            print(f"E3 {pn} segments {len(segs)} triples {int(F.sum())} | D {obs:.5f} null mean {null.mean():.5f} D* {max(obs-null.mean(),0):.5f} nats/triple | p {pv:.4f}",flush=True)
