# Checks that sqw_ref.py section 15 reproduces the stored Plan A results (20n, 20o, 20p, 20q, 20s-A). Run from hmm_audit/.
import sys,os,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, sqw, sqw2
from sqw15 import seqs, DD
ok=True
def check(name,cond):
    global ok; ok&=bool(cond); print(('PASS ' if cond else 'FAIL ')+name,flush=True)
N={(d['g'],d['fly']):d for d in json.load(open(sqw.OUTD+'n_main.json'))}
Om={(d['g'],d['fly']):d for d in json.load(open(sqw.OUTD+'o_memory_drift.json'))}
Pk={(d['g'],d['fly']):d for d in json.load(open(sqw.OUTD+'p_blockkl.json'))}
Q={(d['g'],d['fly']):d for d in json.load(open(sqw.OUTD+'q_tur.json'))}
S={(d['g'],d['fly']):d for d in json.load(open(sqw.OUTD+'s_observable.json'))}
picks=[('CS',0),('CS',7),('fumin',12),('fumin',30),('per0LD',3),('per0DD',4),('sss',1),('sss',34),('CSDD',2)]
for g,i in picks:
    y,ph,_=seqs(g,i); nP=1 if g in DD else 2; tag=f'{g} {i}'
    ll,pi,P=R.mc_fit(y,ph,nP)
    if (g,i) in N: check(f'{tag}: mc_fit loglik = 20n ({ll:.6f})',abs(ll-N[(g,i)]['ll_mc'])<1e-6)
    if (g,i) in Om:
        check(f'{tag}: mc_fit 6h loglik = 20o',abs(R.mc_fit(y,R.phases(len(y),'6h'),4)[0]-Om[(g,i)]['LL']['MC1_6'])<1e-6)
        check(f'{tag}: mc2_fit loglik = 20o',abs(R.mc2_fit(y,ph,nP)[0]-Om[(g,i)]['LL']['MC2_12'])<1e-6)
    st=R.block_kl_stats(y,ph,nP)
    check(f'{tag}: block_kl_stats = 20p',all(abs(st[p][k]-Pk[(g,i)]['obs'][p][k])<1e-12 for p in range(nP) for k in st[p]))
    te=R.tur_estimate(y,ph,nP); qo=Q[(g,i)]['obs']
    check(f'{tag}: tur_estimate = 20q',all((np.isnan(te[p][k]) and np.isnan(qo[p][k])) or abs(te[p][k]-qo[p][k])<1e-12 for p in range(nP) for k in te[p]))
    Ps=R.sim_rows(P); qe=Q[(g,i)]['ex']
    check(f'{tag}: tur_exact / tur_exact_ct = 20q',all((not np.isfinite(qe[p]['B_mc']) and not np.isfinite(R.tur_exact(Ps[p])[2])) or abs(R.tur_exact(Ps[p])[2]-qe[p]['B_mc'])<1e-9 for p in range(nP)))
    check(f'{tag}: routes_observed = route_counts on complete records',(y>=0).all()==False or np.allclose(R.routes_observed(y,ph,nP),R.route_counts(y,ph,nP)))
    oq=R.observable_quantities(y,ph,nP,P); so=S[(g,i)]['main']
    keys=['c_obs','J','EPR','J_chain','A_chain','EPR_chain','E21','E31','F_SW','F_WS','Xn_cnt','Yn_cnt']
    same=lambda a,b:(a is None and b is None) or (isinstance(a,float) and isinstance(b,float) and ((not np.isfinite(a) and not np.isfinite(b) and str(a)==str(b)) or abs(a-b)<1e-9))
    check(f'{tag}: observable_quantities = 20s-A',all(same(oq[p][k],so[p][k]) for p in range(nP) for k in keys))
T,d=R.hmm_observed_chain(np.array([[0.9,0.06,0.04],[0.03,0.95,0.02],[0.02,0.04,0.94]]),np.eye(3))
check('hmm_observed_chain with B = I reproduces the chain',np.allclose(R.tur_exact(T,d),R.tur_exact(np.array([[0.9,0.06,0.04],[0.03,0.95,0.02],[0.02,0.04,0.94]]))))
print('ALL PASS' if ok else 'SOME FAILED')
