# PRESPEC 20q report (q_tur.json), with amendment C. All quantities per fly x phase; group statements by paired tests.
import os,sys,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
from scipy.stats import wilcoxon, binomtest
R3=len(sys.argv)>1 and sys.argv[1]=='r3'
O=json.load(open(sqw.OUTD+'q_tur_r3.json'))['data'] if R3 else json.load(open(sqw.OUTD+'q_tur.json'))   # amendment E
BL=[('per0DD',0,'whole',True),('CSDD',0,'whole',True),('CS',0,'light',False),('CS',1,'dark',False),('per0LD',0,'light',False),
    ('per0LD',1,'dark',False),('fumin',0,'light',False),('fumin',1,'dark',False),('sss',0,'light',False),('sss',1,'dark',False)]
def holm(ps):
    ps=np.array(ps); o=np.argsort(ps); m=len(ps); adj=np.empty(m); run=0
    for k,i in enumerate(o): run=max(run,(m-k)*ps[i]); adj[i]=min(1,run)
    return adj
NOCYC={}; NOEST={}
def block(g,p,key,tau=180):
    rows=[]
    for d in O:
        if d['g']!=g or d['ref'][key] is None: continue
        Bo=d['obs'][p][f'B{tau}']; Bs=np.array([s[p][f'B{tau}'] for s in d['simMC']],float); Bs=Bs[np.isfinite(Bs)]
        Bh=np.array([s[p][f'B{tau}'] for s in d['simH']],float); Bh=Bh[np.isfinite(Bh)]; e=d['ex'][p]
        if not (np.isfinite(e['sig_mc']) and e['sig_mc']>1e-9 and np.isfinite(e['B_mc_ct']) and e['B_mc_ct']>1e-12):   # amendment D: no cycle in the chain
            NOCYC[(g,key)]=NOCYC.get((g,key),0)+1; continue
        if not (np.isfinite(Bo) and len(Bs)>=100): NOEST[(g,key)]=NOEST.get((g,key),0)+1; continue
        kap=e['B_mc']/e['B_mc_ct']; Bc=Bo*e['B_mc']/np.median(Bs)/kap
        rows.append(dict(fly=d['fly'],deg=d['deg'],Bo=Bo,medB=np.median(Bs),p1=(1+np.sum(Bs>=Bo))/(1+len(Bs)),Bc=Bc,kap=kap,
                         sig_mc=e['sig_mc'],B_mc=e['B_mc'],sig_h=e['sig_hmm'],
                         inH=(np.percentile(Bh,2.5)<=Bo<=np.percentile(Bh,97.5)) if len(Bh)>=50 and not d['deg'] else np.nan,
                         J=d['obs'][p]['J'],D=d['obs'][p][f'D{tau}'],nb=d['obs'][p][f'n{tau}']))
    return rows
out=[]
for g,p,key,prim in BL:
    R_=block(g,p,key); r=dict(g=g,ph=key,prim=prim,n=len(R_))
    a=np.array([np.log(x['Bo']/x['medB']) for x in R_]); r['S1_med']=np.median(np.exp(a)); r['S1_p']=wilcoxon(a).pvalue
    p1=np.array([x['p1'] for x in R_]); r['S1_frac']=np.mean(p1<0.05); r['S1_bin']=binomtest(int(np.sum(p1<0.05)),len(p1),0.05,alternative='greater').pvalue
    b=np.array([np.log(x['Bc']/x['sig_mc']) for x in R_]); r['S2_med']=np.median(np.exp(b)); r['S2_p']=wilcoxon(b).pvalue; r['S2_frac']=np.mean(b>0)
    nd=[x for x in R_ if not x['deg']]; sh=np.array([x['sig_h'] for x in nd]); bc=np.array([x['Bc'] for x in nd])
    r['HMM_enough']=np.mean(sh>=bc); r['HMM_inside']=np.nanmean([x['inH'] for x in nd])
    r['Q_mc']=np.median([x['B_mc']/x['sig_mc'] for x in R_]); r['Q_obs']=np.median([x['Bc']/x['sig_mc'] for x in R_])
    qh=np.array([x['Bc']/x['sig_h'] if np.isfinite(x['sig_h']) and x['sig_h']>0 else 0.0 for x in nd]); r['Q_hmm']=np.median(qh)
    r['kap_max']=max(x['kap'] for x in R_); r['J']=np.median([x['J'] for x in R_]); r['D']=np.median([x['D'] for x in R_]); r['nb']=int(np.median([x['nb'] for x in R_]))
    r['sig_mc']=np.median([x['sig_mc'] for x in R_]); r['Bo']=np.median([x['Bo'] for x in R_])
    for tau in (120,360):
        Rt=block(g,p,key,tau); r[f'S1_{tau}']=np.median([x['Bo']/x['medB'] for x in Rt]); r[f'S2_{tau}']=np.median([x['Bc']/x['sig_mc'] for x in Rt])
    out.append(r)
prim=[r for r in out if r['prim']]
for k in ('S1_p','S2_p'):
    h=holm([r[k] for r in prim])
    for r,x in zip(prim,h): r[k+'_holm']=x
print(('AMENDMENT E (rule of three for one-way edges). ' if R3 else '')+'20q. All per fly x phase; tau = 3 h primary. B = J^2/D (k_B/h). Ratios are medians over flies of per-fly ratios.')
for r in out:
    tag='PRIMARY' if r['prim'] else 'descriptive (LD)'
    print(f"   [{r['g']} {r['ph']}: excluded, observed chain without a cycle: {NOCYC.get((r['g'],r['ph']),0)//3}; no estimate (blocks < 4 or D <= 0): {NOEST.get((r['g'],r['ph']),0)//3}]")
    s1=f"S1 B_obs/median(B* chain) {r['S1_med']:.2f} (Wilcoxon p {r['S1_p']:.2g}{', Holm '+format(r['S1_p_holm'],'.2g') if r['prim'] else ''}; flies p<.05 {r['S1_frac']:.2f}, binom {r['S1_bin']:.2g})"
    s2=f"S2 B_c/sigma_MC {r['S2_med']:.2f} (p {r['S2_p']:.2g}{', Holm '+format(r['S2_p_holm'],'.2g') if r['prim'] else ''}; >1 in {r['S2_frac']:.2f})"
    print(f"{r['g']:6s} {r['ph']:5s} n {r['n']:2d} [{tag}] | J {r['J']:+.3f}/h D {r['D']:.3f}/h blocks {r['nb']} | B_obs {r['Bo']:.3f} sigma_MC {r['sig_mc']:.3f} | {s1} | {s2} | "
          f"sigma_HMM >= B_c in {r['HMM_enough']:.2f}, HMM sims bracket B_obs {r['HMM_inside']:.2f} | Q_MC {r['Q_mc']:.2f} Q_obs {r['Q_obs']:.2f} Q_HMM {r['Q_hmm']:.2f} | kappa max {r['kap_max']:.4f} | "
          f"tau 2h/6h: S1 {r['S1_120']:.2f}/{r['S1_360']:.2f} S2 {r['S2_120']:.2f}/{r['S2_360']:.2f}")
    if r['prim']:
        v1=r['S1_p_holm']<0.05 and r['S1_med']>1 and r['S1_bin']<0.05; v2=r['S2_p_holm']<0.05 and r['S2_med']>1
        print(f"        verdict: S1 precision beyond observed chain {'SUPPORTED' if v1 else 'not supported'}; S2 TUR demands hidden dissipation {'SUPPORTED' if v2 else 'not supported'}")
