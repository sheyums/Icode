import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from common import *
from analysis_C import ffbs, counts
g_=np.linspace(-3,4,701); A,Bg=np.meshgrid(g_,g_,indexing='ij'); A=np.r_[A.ravel(),1.]; Bg=np.r_[Bg.ravel(),0.]
keep=np.abs(A*(1-Bg)-(1-A)*Bg)>1e-9; A,Bg=A[keep],Bg[keep]
def Smats(blk):
    i,j=blk; o=3-i-j; S=np.zeros((A.size,3,3)); S[:,i,i]=A; S[:,i,j]=1-A; S[:,j,i]=Bg; S[:,j,j]=1-Bg; S[:,o,o]=1
    return S,np.linalg.inv(S)
SS={(0,1):Smats((0,1)),(0,2):Smats((0,2))}
def ll(P,B,y,pi0):
    f=pi0*B[:,y[0]]; s=f.sum(); L=np.log(s); f/=s
    for t in range(1,len(y)): f=(f@P)*B[:,y[t]]; s=f.sum(); L+=np.log(s); f/=s
    return L
rng=np.random.default_rng(1); out={}
for g in ['CS','fumin','sss','per0']:
    W=windows(g)
    for blk,name in [((0,1),'H1==H2'),((0,2),'H1==H3')]:
        S,Si=SS[blk]; R=[]; spot=[]
        for w in W:
            if np.abs(w['B'][blk[0]]-w['B'][blk[1]]).sum()>=1e-6: continue
            Pp=S@w['P']@Si; ok=np.flatnonzero(Pp.min(axis=(1,2))>=0)
            x=xn(Pp[ok]); y=yn(Pp[ok]); d=x-y; s=x+y
            rr=lambda v:(np.nanmin(v),np.nanmax(v)) if np.isfinite(v).any() else (np.nan,np.nan)
            R.append((*rr(x),*rr(y),*rr(d),*rr(s)))
            if name=='H1==H2' and len(spot)<3 and w['win']==min(v['win'] for v in W if v['fly']==w['fly'] and np.abs(v['B'][0]-v['B'][1]).sum()<1e-6) and np.isfinite(d).any():
                res=[]
                for lab,P,pi0 in [('fit',w['P'],np.ones(3)/3),('maxD',Pp[ok[np.nanargmax(d)]],np.ones(3)/3@Si[ok[np.nanargmax(d)]]),('minD',Pp[ok[np.nanargmin(d)]],np.ones(3)/3@Si[ok[np.nanargmin(d)]])]:
                    c=ffbs(w['y'],P,w['B'],rng); X=(c[1]-c[0])/(c[1]+c[0]); Y=(c[3]-c[2])/(c[3]+c[2])
                    res.append(f"{lab}: Xn {X:+.2f} Yn {Y:+.2f} D {X-Y:+.2f} logL {ll(P,w['B'],w['y'],pi0):.4f}")
                spot.append(f"fly{w['fly']} win{w['win']}: "+' | '.join(res))
        R=np.array(R); 
        if not len(R): print(f"{g:6s} {name}: 0 windows"); continue
        print(f"\n{g:6s} {name}: {len(R)} windows")
        for k,lab in enumerate(['Xn','Yn','Delta=Xn-Yn','Sigma=Xn+Yn']):
            lo,hi=R[:,2*k],R[:,2*k+1]; fin=np.isfinite(lo)
            print(f"   {lab:12s} invariant {np.mean((hi-lo)[fin]<0.05):.3f}  sign determined {np.mean(((lo>0)|(hi<0))[fin]):.3f}  median width {np.median((hi-lo)[fin]):.2f}  undefined {np.sum(~fin)}  | window-mean identified set [{lo[fin].mean():+.3f}, {hi[fin].mean():+.3f}]")
        for sp in spot: print('   spot',sp)
        out[(g,name)]=R
np.save('D_results.npy',out,allow_pickle=True)
