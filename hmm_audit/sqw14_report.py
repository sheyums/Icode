# PRESPEC 20o report (o_memory_drift.json; MC1_12 / HMM_12 predictive results from n_main.json, 20n).
import os,sys,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
O=json.load(open(sqw.OUTD+'o_memory_drift.json')); N={(d['g'],d['fly']):d for d in json.load(open(sqw.OUTD+'n_main.json'))}
MS=['MC1_12','HMM_12','MC2_12','MC1_6','HMM_6','MC2_6']; PH=['light','dark']
q=lambda a:f'{np.median(a):+.1f} [{np.percentile(a,25):+.1f}, {np.percentile(a,75):+.1f}]'
for g in ('CS','fumin'):
    D=[d for d in O if d['g']==g]; L=lambda m:np.array([d['LL'][m] for d in D]); T=np.array([d['T'] for d in D])
    print(f'## {g} (n {len(D)})')
    for m in MS[1:]: print(f'  dLL {m:7s} - MC1_12: median [IQR] {q(L(m)-L("MC1_12"))} nats | per 1000 min {np.median(1000*(L(m)-L("MC1_12"))/T):+.2f}')
    gain=L('HMM_12')-L('MC1_12'); ok=gain>3
    Rm=(L('MC2_12')-L('MC1_12'))[ok]/gain[ok]; Rd=1-(L('HMM_6')-L('MC1_6'))[ok]/gain[ok]
    print(f'  flies with HMM gain > 3 nats: {ok.sum()} | R_mem median {np.median(Rm):.2f} IQR [{np.percentile(Rm,25):.2f}, {np.percentile(Rm,75):.2f}], >= 1 in {np.mean(Rm>=1):.2f}'
          f' | R_drift median {np.median(Rd):.2f} IQR [{np.percentile(Rd,25):.2f}, {np.percentile(Rd,75):.2f}], >= 0.5 in {np.mean(Rd>=0.5):.2f}')
    print(f'  HMM_6 - MC1_6 (gain left after drift): {q(L("HMM_6")-L("MC1_6"))} | MC2_12 - HMM_12: {q(L("MC2_12")-L("HMM_12"))} | MC2_6 - HMM_6: {q(L("MC2_6")-L("HMM_6"))}')
    for crit in ('BIC','AIC'):
        best=[min(MS,key=lambda m:d[crit][m]) for d in D]; print(f'  best by {crit}: '+', '.join(f'{m} {best.count(m)/len(D):.2f}' for m in MS))
    print('  predictive: fraction of fly-phases outside the central 95% of 20 simulations')
    for p in range(2):
        row=[]
        for s in ('chg_h','q_run','w_run','s_run','q_1min','w_1min','c_obs'):
            fr=[]
            for name,src in (('MC1_12','n'),('HMM_12','n'),('MC2_12','o'),('MC1_6','o')):
                o=[]
                for d in D:
                    dn=N[(g,d['fly'])]; r=dn['real'][p][s]
                    sims=dn['sims']['MC' if name=='MC1_12' else 'HMM'] if src=='n' else d['sims'][name]
                    v=np.array([x[p][s] for x in sims],float); v=v[np.isfinite(v)]
                    if np.isfinite(r) and len(v)>=10: lo,hi=np.percentile(v,[2.5,97.5]); o.append(r<lo or r>hi)
                fr.append(f'{np.mean(o):.2f}')
            row.append(f'{s} '+'/'.join(fr))
        print(f'   {PH[p]:5s} (MC1_12/HMM_12/MC2_12/MC1_6) | '+' | '.join(row))
