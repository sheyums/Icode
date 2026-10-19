# PRESPEC 8: model-free time-irreversibility of sleep/wake pairs, reversal null, fly bootstrap.
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from pooled import GEN, np
from pooled_LD import bouts_t
def kl(a,b):
    a=a+0.5; b=b+0.5; a=a/a.sum(-1,keepdims=True); b=b/b.sum(-1,keepdims=True); return (a*np.log(a/b)).sum(-1)
def build(B,nb,phase_sel,qS,qW,edgefree=False):
    # per window: SW and WS count matrices (flattened nb*nb), fly id, cycle durations
    ph=((B[:,6]/3600)%24>=12).astype(int); W={}; cyc=[]
    edge=np.zeros(len(B),bool)
    if edgefree:
        for key in set(zip(B[:,0],B[:,1])):
            m=np.flatnonzero((B[:,0]==key[0])&(B[:,1]==key[1])); edge[m[0]]=edge[m[-1]]=True
    for k in range(len(B)-1):
        if B[k,2]!=B[k+1,2] or B[k+1,5]!=B[k,5]+1 or B[k,1]!=B[k+1,1]: continue   # same window, consecutive
        if phase_sel is not None and ph[k]!=phase_sel: continue
        if edge[k] or edge[k+1]: continue
        key=(B[k,0],B[k,1]); W.setdefault(key,[np.zeros(nb*nb),np.zeros(nb*nb)])
        if B[k,3]==0: s,w,t=B[k,4],B[k+1,4],0; cyc.append(s+w)
        else: s,w,t=B[k+1,4],B[k,4],1
        i=np.searchsorted(qS,s,side='left'); j=np.searchsorted(qW,w,side='left')
        W[key][t][i*nb+j]+=1
    keys=list(W); SW=np.array([W[k][0] for k in keys]); WS=np.array([W[k][1] for k in keys]); fl=np.array([k[0] for k in keys])
    return SW,WS,fl,np.mean(cyc)
def null_D(SW,WS,flips):        # flips: (R, nW) 0/1 ; 1 = reversed (swap)
    d=WS-SW; cSW=SW.sum(0)+flips@d; cWS=WS.sum(0)-flips@d; return kl(cSW,cWS)
EDGEFREE=len(sys.argv)>1 and sys.argv[1]=='edgefree'
if __name__=='__main__':
    rng=np.random.default_rng(8)
    for g in GEN:
        B=bouts_t(g)
        for blk,sel in [('whole',None),('light',0),('dark',1)]:
            line=[]; nullsets=[]; obsD=[]
            flips=None
            for nb in [3,4,5]:
                qS=np.quantile(B[B[:,3]==0,4],np.linspace(0,1,nb+1)[1:-1]); qW=np.quantile(B[B[:,3]==1,4],np.linspace(0,1,nb+1)[1:-1])
                SW,WS,fl,mc=build(B,nb,sel,qS,qW,EDGEFREE)
                if flips is None: flips=(rng.random((2000,len(SW)))<0.5).astype(float)
                Dobs=kl(SW.sum(0),WS.sum(0)); Dn=null_D(SW,WS,flips); p=(1+np.sum(Dn>=Dobs))/2001
                Dstar=max(Dobs-Dn.mean(),0.)
                flies=np.unique(fl); bs=[]
                for _ in range(2000):
                    pick=rng.choice(flies,len(flies)); idx=np.concatenate([np.flatnonzero(fl==f) for f in pick])
                    s,w=SW[idx],WS[idx]; fr=(rng.random((500,len(idx)))<0.5).astype(float)
                    bs.append(kl(s.sum(0),w.sum(0))-null_D(s,w,fr).mean())
                ci=np.percentile(bs,[2.5,97.5]); obsD.append(Dobs); nullsets.append(Dn)
                line.append(f"{nb}x{nb}: D {Dobs:.4f} null {Dn.mean():.4f} D* {Dstar:.4f} CI [{ci[0]:+.4f},{ci[1]:+.4f}] p {p:.4f} | D*/h {Dstar/(mc/3600):.4f}")
            N=np.array(nullsets); thr=np.percentile(N,95,axis=1)[:,None]; glob=np.mean(np.all(N>thr,axis=0))
            passed=all((1+np.sum(N[i]>=obsD[i]))/2001<0.05 for i in range(3))
            print(f"{g:6s} {blk:5s} pairs SW {int(SW.sum())} WS {int(WS.sum())} windows {len(SW)} | PASS(all 3 bin sets) {passed} | null pass rate {glob:.4f}\n    "+"\n    ".join(line),flush=True)
