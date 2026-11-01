import sys, time, numpy as np
sys.path.insert(0, '/home/user/Icode/hmm_audit')
import planB_20u as M
rng = np.random.default_rng(7)
# known reversible model: Ka=1, Ki=2 (inactive fast + slow)
pi = np.array([0.08, 0.32, 0.60]); F = np.zeros((3, 3))
F[0,1]=F[1,0]=0.08*0.05; F[1,2]=F[2,1]=0.32*0.002; F[0,2]=F[2,0]=0.08*0.0005
for i in range(3): F[i,i]=pi[i]-F[i].sum()
X0 = F/F.sum(); P = X0/X0.sum(1)[:,None]
# brute-force check on a short sequence
M._seed(1); a = M._sim(P, X0.sum(1)/X0.sum(), 1, 3000).astype(int)
sym, L = M.runs_of(a); ll, C = M.loglik_counts(X0, 1, sym, L)
E = np.array([[1,0],[0,1],[0,1]],float)   # state -> symbol (col0 inactive? use emission prob of a)
em = lambda z: np.array([1.0 if (k<1)==(z==1) else 0.0 for k in range(3)])
al = np.zeros((len(a),3)); c = np.zeros(len(a)); p0 = X0.sum(1)
al[0] = p0*em(a[0]); c[0]=al[0].sum(); al[0]/=c[0]
for t in range(1,len(a)):
    al[t] = (al[t-1]@P)*em(a[t]); c[t]=al[t].sum(); al[t]/=c[t]
be = np.ones((len(a),3))
for t in range(len(a)-2,-1,-1):
    be[t] = P@(em(a[t+1])*be[t+1]); be[t]/=c[t+1]
Cb = np.zeros((3,3))
for t in range(len(a)-1):
    Cb += np.outer(al[t], em(a[t+1])*be[t+1])*P/c[t+1]
print('loglik kernel', ll, 'brute', np.log(c).sum(), '| max |C-Cb|', np.abs(C-Cb).max(), 'sum C', C.sum())
# recovery + timing on 300k s
M._seed(2); a = M._sim(P, X0.sum(1)/X0.sum(), 1, 300000).astype(int); sym, L = M.runs_of(a)
t=time.time(); M.loglik_counts(X0,1,sym,L); print('one E-step s', time.time()-t, 'runs', len(L))
t=time.time(); ll, X, it = M.fit_order(sym, L, 1, 2, rng, starts=3)
print('fit s', time.time()-t, 'iters', it)
Ph = X/X.sum(1)[:,None]; print('true P offdiag', np.round(P,5)); print('fit  P', np.round(Ph,5))
# monotonicity
X = M.random_start(3, rng); lls=[]
for k in range(30):
    l, Cc = M.loglik_counts(X,1,sym,L); lls.append(l); X = M.reversible_mle(Cc, X)
print('min EM step', np.min(np.diff(lls)))
