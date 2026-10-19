import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from common import *
from analysis_C import ffbs
from analysis_D import SS, ll
def stat(P):
    n=P.shape[0]; A=np.transpose(P,(0,2,1))-np.eye(3); A[:,2,:]=1
    b=np.zeros((n,3)); b[:,2]=1; return np.linalg.solve(A,b)
def rates(P):
    pi=stat(P); den=P[:,0,1]+P[:,0,2]
    return pi[:,2]*P[:,2,0]*P[:,0,1]/den, pi[:,1]*P[:,1,0]*P[:,0,2]/den
rng=np.random.default_rng(1)
for g in ['CS','fumin','sss','per0']:
    W=windows(g)
    for blk,name in [((0,1),'H1==H2'),((0,2),'H1==H3')]:
        S,Si=SS[blk]; R=[]; spot=[]
        for w in W:
            if np.abs(w['B'][blk[0]]-w['B'][blk[1]]).sum()>=1e-6: continue
            Pp=S@w['P']@Si; ok=np.flatnonzero(Pp.min(axis=(1,2))>=0); Pq=Pp[ok]
            a,b=rates(Pq); K=(a-b)/(a+b)
            f=np.isfinite(K)
            if not f.any(): R.append([np.nan]*6); continue
            R.append((K[f].min(),K[f].max(),a[f].min(),a[f].max(),b[f].min(),b[f].max()))
            if name=='H1==H2' and len(spot)<3 and w['win']==min(v['win'] for v in W if v['fly']==w['fly'] and np.abs(v['B'][0]-v['B'][1]).sum()<1e-6):
                Kf=np.where(f,K,np.nan); res=[]
                for lab,k in [('fit',None),('maxK',np.nanargmax(Kf)),('minK',np.nanargmin(Kf))]:
                    P=w['P'] if k is None else Pq[k]; pi0=np.ones(3)/3 if k is None else np.ones(3)/3@Si[ok[k]]
                    c=ffbs(w['y'],P,w['B'],rng)
                    res.append(f"{lab}: N312 {c[1]:.2f} N213 {c[3]:.2f} logL {ll(P,w['B'],w['y'],pi0):.4f}")
                spot.append(f"fly{w['fly']} win{w['win']}: "+' | '.join(res))
        R=np.array(R)
        if not len(R): continue
        fin=np.isfinite(R[:,0]); lo,hi=R[fin,0],R[fin,1]
        with np.errstate(all='ignore'): ra=np.nanmedian(R[fin,3]/R[fin,2]); rb=np.nanmedian(R[fin,5]/R[fin,4])
        print(f"\n{g:6s} {name}: {len(R)} windows ({(~fin).sum()} undefined)")
        print(f"   K=(r312-r213)/(r312+r213): invariant {np.mean(hi-lo<0.05):.3f}  sign determined {np.mean((lo>0)|(hi<0)):.3f}  median width {np.median(hi-lo):.2f}  window-mean set [{lo.mean():+.3f}, {hi.mean():+.3f}]")
        print(f"   median max/min over set: r312 {ra:.3g}   r213 {rb:.3g}")
        for s in spot: print('   spot',s)
