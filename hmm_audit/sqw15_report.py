# PRESPEC 20p report (p_blockkl.json).
import os,sys,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sqw
from scipy.stats import wilcoxon, binomtest
O=json.load(open(sqw.OUTD+'p_blockkl.json'))
BL=[(g,p) for g in ('CS','per0LD','fumin','sss') for p in (0,1)]+[('per0DD',0),('CSDD',0)]   # genotype order rule
PN=lambda g,p:'whole' if g in ('per0DD','CSDD') else ['light','dark'][p]
rows=[]
for g,p in BL:
    D=[d for d in O if d['g']==g and d['ref'][PN(g,p)] is not None]; r=dict(g=g,ph=PN(g,p),n=len(D))   # amendment B: no Q/W reference -> excluded
    for k in ('d2','d3','d4'): r[k]=np.median([d['obs'][p][k] for d in D])
    for nul in ('A','B'):
        E=[];pv=[]
        for d in D:
            o=d['obs'][p]['D3']; s=np.array([x[p]['D3'] for x in d['sims'][nul]])
            E.append(o-np.median(s)); pv.append((1+np.sum(s>=o))/(1+len(s)))
        E=np.array(E); pv=np.array(pv); r[f'E_{nul}']=np.median(E); r[f'pW_{nul}']=wilcoxon(E).pvalue
        r[f'frac_{nul}']=np.mean(pv<0.05); r[f'pbin_{nul}']=binomtest(int(np.sum(pv<0.05)),len(pv),0.05,alternative='greater').pvalue
    for k in ('D4','D3QW'):
        E=np.array([d['obs'][p][k]-np.median([x[p][k] for x in d['sims']['A']]) for d in D]); r[f'E_{k}']=np.median(E); r[f'pW_{k}']=wilcoxon(E).pvalue
    for m in ('HMM','MC2'):
        inside=[];pred=[]
        for d in D:
            s=np.array([x[p]['D3'] for x in d['sims'][m]]); lo,hi=np.percentile(s,[2.5,97.5]); o=d['obs'][p]['D3']
            inside.append(lo<=o<=hi); pred.append(np.median(s))
        r[f'in_{m}']=np.mean(inside); r[f'pred_{m}']=np.median(pred)
    r['D3obs']=np.median([d['obs'][p]['D3'] for d in D])
    r['ratio']=np.median([d['obs'][p]['d3']/d['obs'][p]['d2'] for d in D if d['obs'][p]['d2']>0])
    rows.append(r)
def holm(ps):
    ps=np.array(ps); o=np.argsort(ps); m=len(ps); adj=np.empty(m); run=0
    for k,i in enumerate(o): run=max(run,(m-k)*ps[i]); adj[i]=min(1,run)
    return adj
for nul in ('A','B'):
    h=holm([r[f'pW_{nul}'] for r in rows])
    for r,x in zip(rows,h): r[f'holm_{nul}']=x
print('d_n in nats/h, medians over flies. E = Delta3_obs - median Delta3 under the null. Delta3 = d3 - d2.')
for r in rows:
    q=all(r[f'holm_{n}']<0.05 and r[f'E_{n}']>0 and r[f'pbin_{n}']<0.05 for n in ('A','B'))
    print(f"{r['g']:6s} {r['ph']:5s} n {r['n']:2d} | d2 {r['d2']:.3f} d3 {r['d3']:.3f} d4 {r['d4']:.3f} | Delta3 obs {r['D3obs']:+.3f} | "
          f"null A: E {r['E_A']:+.3f} Holm {r['holm_A']:.2g} p<.05 in {r['frac_A']:.2f} (binom {r['pbin_A']:.2g}) | "
          f"null B: E {r['E_B']:+.3f} Holm {r['holm_B']:.2g} p<.05 in {r['frac_B']:.2f} (binom {r['pbin_B']:.2g}) | "
          f"Delta4 E {r['E_D4']:+.3f} (p {r['pW_D4']:.2g}) | Delta3_QW E {r['E_D3QW']:+.3f} (p {r['pW_D3QW']:.2g}) | "
          f"HMM pred {r['pred_HMM']:+.3f} inside {r['in_HMM']:.2f} | MC2 pred {r['pred_MC2']:+.3f} inside {r['in_MC2']:.2f} | d3/d2 {r['ratio']:.2f} | QUALIFIES {q}")

# amendment A: DD genotypes in 12-h clock halves, nulls A2 (two 12-h matrices) and B (6-h)
fn=sqw.OUTD+'p_blockkl_ddhalves.json'
if os.path.exists(fn):
    H=json.load(open(fn)); REF={(d['g'],d['fly']):d['ref']['whole'] for d in O if d['g'] in ('per0DD','CSDD')}; hr=[]
    for g in ('per0DD','CSDD'):
        for p,lab in ((0,'ZT0-12'),(1,'ZT12-24')):
            D=[d for d in H if d['g']==g and REF[(g,d['fly'])] is not None]; r=dict(g=g,ph=lab,n=len(D))
            for k in ('d2','d3','d4','D3'): r[k]=np.median([d['obs'][p][k] for d in D])
            for nul in ('A2','B'):
                E=[];pv=[]
                for d in D:
                    o=d['obs'][p]['D3']; s_=np.array([x[p]['D3'] for x in d['sims'][nul]]); E.append(o-np.median(s_)); pv.append((1+np.sum(s_>=o))/(1+len(s_)))
                E=np.array(E); pv=np.array(pv); r[f'E_{nul}']=np.median(E); r[f'pW_{nul}']=wilcoxon(E).pvalue; r[f'frac_{nul}']=np.mean(pv<0.05)
                r[f'pbin_{nul}']=binomtest(int(np.sum(pv<0.05)),len(pv),0.05,alternative='greater').pvalue
            hr.append(r)
    for nul in ('A2','B'):
        h=holm([r[f'pW_{nul}'] for r in hr])
        for r,x in zip(hr,h): r[f'holm_{nul}']=x
    print('AMENDMENT A: DD genotypes in 12-h clock halves')
    for r in hr:
        q=all(r[f'holm_{n}']<0.05 and r[f'E_{n}']>0 and r[f'pbin_{n}']<0.05 for n in ('A2','B'))
        print(f"{r['g']:6s} {r['ph']:7s} n {r['n']:2d} | d2 {r['d2']:.3f} d3 {r['d3']:.3f} d4 {r['d4']:.3f} | Delta3 obs {r['D3']:+.3f} | "
              f"null A2: E {r['E_A2']:+.3f} Holm {r['holm_A2']:.2g} p<.05 in {r['frac_A2']:.2f} (binom {r['pbin_A2']:.2g}) | "
              f"null B: E {r['E_B']:+.3f} Holm {r['holm_B']:.2g} p<.05 in {r['frac_B']:.2f} (binom {r['pbin_B']:.2g}) | QUALIFIES {q}")
