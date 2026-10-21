import sys,time; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
import numpy as np
from hhmm_fly import *
from synth_test import P3, truth_model
from validate13 import perturb
KS=np.exp([-.4,-.13,.13,.4]); rs=np.random.RandomState(7); S=[];fl=[];dy=[]
for i in range(8):
    base=perturb(P3,rs); eps,A,G,Bin=base._struct
    for d in range(4):
        e=eps.copy(); e[:2]*=KS[d]; m=truth_model(3,e,A,G,Bin,base.emissionprob_); S.append(m.sample(14400,random_state=rs)[0].ravel()); fl.append(i); dy.append(d)
t=time.time(); fits=[fit_fly(S,fl,dy,8,4,3,r) for r in range(3)]; best=max(fits,key=lambda m:m.monitor_.history[-1]); cs=summary_fly(best)
h=np.array(best.monitor_.history); print(f"{time.time()-t:.0f}s for 3 fits | iters {len(h)} | largest logL decrease {max(0,-np.diff(h).min()):.4f} nats")
print("true kappa_S",np.round(KS,2).tolist(),"| fitted kappa (S,T,W)",np.round(cs['kappa'],2).tolist())
print("Xn_c %+.3f (truth +0.333) | Yn_c %+.3f (truth +0.077) | mf %s"%(cs['Xn_c'],cs['Yn_c'],np.round(cs['mf'],3).tolist()))
print("restart final logLs",[round(m.monitor_.history[-1],1) for m in fits])
