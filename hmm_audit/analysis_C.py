import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from common import *
from multiprocessing import Pool
M=200
def counts(z):
    c=z[np.r_[0,np.flatnonzero(np.diff(z))+1]]
    n=lambda *p: int(np.sum(np.all([c[k:len(c)-len(p)+1+k]==p[k] for k in range(len(p))],axis=0))) if len(c)>=len(p) else 0
    return [n(2,1),n(2,0,1),n(1,2),n(1,0,2)]
def viterbi(y,P,B):
    lP=np.log(np.maximum(P,1e-300)); lB=np.log(np.maximum(B,1e-300)); T=len(y)
    d=np.log(1/3)+lB[:,y[0]]; psi=np.zeros((T,3),int)
    for t in range(1,T):
        s=d[:,None]+lP; psi[t]=s.argmax(0); d=s.max(0)+lB[:,y[t]]
    z=np.empty(T,int); z[-1]=d.argmax()
    for t in range(T-1,0,-1): z[t-1]=psi[t,z[t]]
    return z
def ffbs(y,P,B,rng):
    T=len(y); F=np.empty((T,3)); f=B[:,y[0]]/3
    if f.sum()==0: return None
    F[0]=f/f.sum()
    for t in range(1,T):
        f=(F[t-1]@P)*B[:,y[t]]; s=f.sum()
        if s==0: return None
        F[t]=f/s
    Z=np.empty((M,T),int); Z[:,-1]=np.searchsorted(F[-1].cumsum(),rng.random(M)*F[-1].sum())
    for t in range(T-2,-1,-1):
        Q=F[t][:,None]*P[:,Z[:,t+1]]; Q/=Q.sum(0)
        Z[:,t]=np.minimum((rng.random(M)[None,:]>Q.cumsum(0)).sum(0),2)
    return np.mean([counts(z) for z in Z],0)
def job(a):
    k,w=a; rng=np.random.default_rng([1,k])
    v=counts(viterbi(w['y'],w['P'],w['B'])); f=ffbs(w['y'],w['P'],w['B'],rng)
    L12=np.abs(w['B'][0]-w['B'][1]).sum()
    return (w['fly'],w['win'],L12,*v,*(f if f is not None else [np.nan]*4))
if __name__=='__main__':
    for g in (sys.argv[1:] or ['CS','fumin']):
        W=windows(g)
        with Pool(4) as p: r=np.array(p.map(job,list(enumerate(W)),chunksize=4))
        np.save(f'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit/C_{g}.npy',r)
        print(f'{g}: done {len(r)} windows, FFBS failures {np.isnan(r[:,7]).sum()}',flush=True)
