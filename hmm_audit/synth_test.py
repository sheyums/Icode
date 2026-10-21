# PRESPEC 12 synthetic validation (criteria in PRESPEC_12.md, commit b3b43df).
import sys,time; sys.path.insert(0,'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit')
from hhmm import *
def truth_model(Nc,eps,A,G,Bin,E):
    m=HHMM(Nc=Nc); m.startprob_=np.ones(Nc*2)/(Nc*2); m.transmat_=m.build(eps,A,G,Bin); m.emissionprob_=E; m._struct=(eps,A,G,Bin); return m
ES=[[.999,.001,0,0],[.70,.0,.10,.20]]; ET=[[.90,.10,0,0],[.70,.28,.02,0]]; EW=[[.85,.13,.02,0],[.25,.60,.12,.03]]
AS=np.array([[.995,.005],[.2,.8]]); AT=np.array([[.9,.1],[.15,.85]]); AW=np.array([[.92,.08],[.08,.92]])
P3=truth_model(3,np.array([1/900,1/900,1/120,1/120,1/300,1/300]),[AS,AT,AW],np.array([[0,.8,.2],[.5,0,.5],[.3,.7,0]]),
               [np.array([.9,.1]),np.array([.5,.5]),np.array([.5,.5])],np.array(ES+ET+EW))
P2=truth_model(2,np.array([1/900,1/900,1/300,1/300]),[AS,AW],np.array([[0,1.],[1.,0]]),[np.array([.9,.1]),np.array([.5,.5])],np.array(ES+EW))
def sample(m,seed): rs=np.random.RandomState(seed); return [m.sample(14400,random_state=rs)[0].ravel() for _ in range(4)]
if __name__=="__main__":
  ok=True
  for name,tm,seed in [('P3',P3,1),('P2',P2,2)]:
      S=sample(tm,seed); tsum=coarse_summary(tm); t=time.time()
      fits=[fit(S,3,r) for r in range(6)]; best=max(fits,key=lambda m:score(m,S)); cs=coarse_summary(best)
      mono=all(np.all(np.diff(f.monitor_.history)>=-1e-6) for f in fits)
      tr,te=S[:3],S[3:]; f2=max([fit(tr,2,r) for r in range(6)],key=lambda m:score(m,tr)); f3=max([fit(tr,3,r) for r in range(6)],key=lambda m:score(m,tr))
      gain=(score(f3,te)-score(f2,te))/len(te[0])*3600
      print(f"== {name} ({time.time()-t:.0f}s) truth mf {np.round(tsum['mf'],3).tolist()} dwell {np.round(tsum['dwell'],0).tolist()}"+(f" Xn_c {tsum['Xn_c']:+.3f} Yn_c {tsum['Yn_c']:+.3f}" if 'Xn_c' in tsum else ''))
      print(f"   fitted Nc=3: mf {np.round(cs['mf'],3).tolist()} dwell {np.round(cs['dwell'],0).tolist()} Xn_c {cs['Xn_c']:+.3f} Yn_c {cs['Yn_c']:+.3f} | restart logLs {sorted([round(score(m,S),1) for m in fits])}")
      print(f"   monotone EM: {mono} | held-out gain Nc3-Nc2: {gain:+.1f} nats/h")
      if name=='P3':
          c1=all(abs(cs[k]-tsum[k])<=0.15 and np.sign(cs[k])==np.sign(tsum[k]) for k in ('Xn_c','Yn_c'))
          c2=all(1/1.5<=a/b<=1.5 for a,b in zip(cs['dwell'],tsum['dwell'])); c3=gain>0
          print(f"   criterion 1 {c1} | 2 {c2} | 3 {c3} | 5 {mono}"); ok&=c1 and c2 and c3 and mono
      else:
          c4=gain<=5; print(f"   criterion 4 (negative control) {c4} | 5 {mono}"); ok&=c4 and mono
  print('SYNTHETIC VALIDATION', 'PASS' if ok else 'FAIL')
