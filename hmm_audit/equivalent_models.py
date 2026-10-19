# Degenerate H1/H2 emissions: construct an observationally EQUIVALENT model
# P' = S P S^-1 (S block-diagonal on {H1,H2}, rows sum to 1) and compare routing.
import numpy as np
exec(open(__file__.replace('equivalent_models.py','viterbi_vs_ffbs.py')).read().split('for label,B')[0].replace('REPS = 200','REPS = 1'))
B = np.array([[.4,.6],[.4,.6],[.98,.02]])
def xnpred(P): k=P[0,1]/(P[0,1]+P[0,2]); return (P[2,0]*k-P[2,1])/(P[2,0]*k+P[2,1])
def loglik(P,y,pi0):
    f = pi0*B[:,y[0]]; ll=np.log(f.sum()); f/=f.sum()
    for t in range(1,len(y)):
        f=(f@P)*B[:,y[t]]; s=f.sum(); ll+=np.log(s); f/=s
    return ll
best=None
for a in np.linspace(0.5,1.5,201):
  for b in np.linspace(-0.5,0.5,201):
    S=np.array([[a,1-a,0],[b,1-b,0],[0,0,1.]])
    if abs(np.linalg.det(S))<1e-6: continue
    Pp=S@P@np.linalg.inv(S)
    if Pp.min()<0: continue
    if best is None or xnpred(Pp)>xnpred(best[0]): best=(Pp,S)
Pp,S=best
print("Original P:\n",np.array2string(P,precision=6)); print("Xn_pred =",round(xnpred(P),3))
print("Equivalent P':\n",np.array2string(Pp,precision=6)); print("Xn_pred' =",round(xnpred(Pp),3))
z,y=simulate(B); pi0=np.ones(3)/3; pi0p=pi0@np.linalg.inv(S)
print("pi0' =",pi0p, " (valid distribution:", bool(pi0p.min()>=0),")")
print("logL under P  :",loglik(P,y,pi0)); print("logL under P' :",loglik(Pp,y,pi0p))
print("true counts (N32,N312):",counts(z))
for nm,PP in [("P ",P),("P'",Pp)]:
    globals()['P']=PP; D=ffbs(y,B); print(f"FFBS under {nm}: mean N32={D[:,0].mean():.2f} N312={D[:,1].mean():.2f}")
