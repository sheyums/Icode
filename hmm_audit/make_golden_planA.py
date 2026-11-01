# Golden test vectors G12 for sqw_ref.py section 15 (Plan A). Synthetic input, no fly data. 0-BASED indices (see PORTING_GUIDE).
import os,json,numpy as np
import sqw_ref as R
from make_golden import clean, dump
rng=np.random.default_rng(1212)
P=np.array([[[0.80,0.15,0.05],[0.06,0.90,0.04],[0.02,0.03,0.95]],      # light: reverse-biased cycle (J < 0; label corrected 2026-11-01, values unchanged)
            [[0.70,0.20,0.10],[0.10,0.85,0.05],[0.03,0.05,0.92]]])     # dark
T=2*1440; ph=R.phases(T); z,_=R.simulate(np.ones(3)/3,P,np.eye(3),ph,rng); y=z.copy(); y[500:530]=-1; y[2000:2003]=-1
ll,pi,Pm=R.mc_fit(y,ph,2); ll2,_,_=R.mc2_fit(y,ph,2)
out=dict(note='y: Q=0, W=1, S=2, missing=-1; ph: 0 = ZT0-12, 1 = ZT12-24 (phase of the ARRIVAL minute). Tolerances: logliks 1e-6, others 1e-9.',
         y=y,ph=ph,mc_fit=dict(loglik=ll,P=Pm),mc2_fit=dict(loglik=ll2),
         block_kl=R.block_kl(y,ph,2),block_kl_qw=R.block_kl(y,ph,2,qw=True),block_kl_stats=R.block_kl_stats(y,ph,2),
         tur_estimate=R.tur_estimate(y,ph,2),
         tur_exact_true_light=R.tur_exact(P[0]),tur_exact_ct_true_light=R.tur_exact_ct(P[0]),sigma_true_light=R.model_energetics(P[0])['EPR'],
         hmm_observed_chain_check=R.tur_exact(*R.hmm_observed_chain(P[0],np.array([[0.9,0.1,0],[0.05,0.95,0],[0,0,1]]))),
         routes_observed=R.routes_observed(y,ph,2),observable_quantities=R.observable_quantities(y,ph,2,Pm),
         sim_rows_example=dict(input=[[[0.5,0.5,0],[0,0,0],[0.2,0.3,0.5]]],output=R.sim_rows([[[0.5,0.5,0],[0,0,0],[0.2,0.3,0.5]]])))
dump('G12_planA.json',out); print('written golden/G12_planA.json; mc_fit loglik',round(ll,4),'| B/sigma light (discrete, continuous):',round(R.tur_exact(P[0])[2]/R.model_energetics(P[0])['EPR'],4),round(R.tur_exact_ct(P[0])[2]/R.model_energetics(P[0])['EPR'],4))
