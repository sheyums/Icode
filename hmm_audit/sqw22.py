# PRESPEC 22 (482898d): within-genotype Spearman of sleep (min/h) against c_obs (primary) and chain J (co-primary,
# fly-phases with a cycle), permutation null (10000, seed 22), Holm over 14. Secondary: no-cycle sleep, stratified
# Fisher z, fumin/sss overlap, pooled figure. Writes sqw_out/s22_sleep_vs_J.json, figures/fig22_J_vs_sleep.{png,pdf}.
import json, numpy as np, sqw
from sqw15 import seqs, DD
from scipy.stats import rankdata
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = json.load(open(sqw.OUTD + 's_observable.json')); rng = np.random.default_rng(22); NP = 10000
PRIM = [('CS', 0), ('CS', 1), ('per0DD', 0), ('fumin', 0), ('fumin', 1), ('sss', 0), ('sss', 1)]
SEC = [('CSDD', 0), ('per0LD', 0), ('per0LD', 1)]
rows = []
for d in O:
    g, i = d['g'], d['fly']; y, ph, _ = seqs(g, i); y = np.asarray(y); ph = np.asarray(ph)
    for p, r in enumerate(d['main']):
        if not r['ok']: continue
        m = (ph == p) & (y >= 0); sl = 60.0 * np.mean(y[m] == 2)
        Jc = r['J_chain']; cyc = Jc is not None and abs(Jc) >= 1e-12
        rows.append(dict(g=g, fly=i, p=p, sleep=float(sl), c_obs=r['c_obs'], J=Jc if cyc else 0.0, cycle=bool(cyc)))
def rk(v): return rankdata(np.round(np.asarray(v, float), 10))   # ties that differ only by round-off rank as ties (2026-11-04, lClaude)
def rho(x, y): return np.corrcoef(rk(x), rk(y))[0, 1]
def perm(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); r0 = rho(x, y); rx = rk(x); ry = rk(y)
    null = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(NP)])
    return r0, (1 + np.sum(np.abs(null) >= abs(r0) - 1e-12)) / (NP + 1)
def block(g, p): return [r for r in rows if r['g'] == g and r['p'] == p]
res = {}
for g, p in PRIM + SEC:
    B = block(g, p); bc = [r for r in B if r['c_obs'] is not None and np.isfinite(r['c_obs'])]; bj = [r for r in B if r['cycle']]
    rc, pc = perm([r['sleep'] for r in bc], [r['c_obs'] for r in bc]); rj, pj = perm([r['sleep'] for r in bj], [r['J'] for r in bj])
    nc = [r['sleep'] for r in B if not r['cycle']]
    res[f'{g}|{p}'] = dict(n=len(B), n_c=len(bc), rho_c=rc, p_c=pc, n_J=len(bj), rho_J=rj, p_J=pj, primary=(g, p) in PRIM,
                           sleep_med=float(np.median([r['sleep'] for r in B])), n_nocycle=len(nc),
                           sleep_nocycle=float(np.median(nc)) if nc else None, sleep_cycle=float(np.median([r['sleep'] for r in bj])) if bj else None)
# Holm over the 14 primary tests
fam = [(k, t) for k, v in res.items() if v['primary'] for t in ('c', 'J')]; ps = np.array([res[k]['p_' + t] for k, t in fam])
o = np.argsort(ps); adj = np.empty(len(ps)); run = 0
for j, idx in enumerate(o): run = max(run, min(1, (len(ps) - j) * ps[idx])); adj[idx] = run
for (k, t), a in zip(fam, adj): res[k]['holm_' + t] = float(a)
# stratified Fisher z for c_obs, permutation within blocks
def strat(perm_it):
    zs = []; ws = []
    for g, p in PRIM:
        bc = [r for r in block(g, p) if r['c_obs'] is not None and np.isfinite(r['c_obs'])]
        x = np.array([r['sleep'] for r in bc]); yv = np.array([r['c_obs'] for r in bc])
        if perm_it: x = rng.permutation(x)
        zs.append(np.arctanh(np.clip(rho(x, yv), -0.999999, 0.999999))); ws.append(len(bc) - 3)
    return np.average(zs, weights=ws)
z0 = strat(False); zn = np.array([strat(True) for _ in range(2000)]); pz = (1 + np.sum(np.abs(zn) >= abs(z0))) / 2001
# fumin/sss overlap in sleep (per phase): fly-phases whose sleep lies inside both genotypes' ranges
ov = {}
for p in (0, 1):
    F = [r for r in block('fumin', p) if r['c_obs'] is not None and np.isfinite(r['c_obs'])]; S = [r for r in block('sss', p) if r['c_obs'] is not None and np.isfinite(r['c_obs'])]
    lo = max(min(r['sleep'] for r in F), min(r['sleep'] for r in S)); hi = min(max(r['sleep'] for r in F), max(r['sleep'] for r in S))
    f = [r for r in F if lo <= r['sleep'] <= hi]; s = [r for r in S if lo <= r['sleep'] <= hi]
    ov[p] = dict(range=[lo, hi], n_fumin=len(f), n_sss=len(s), sleep_med_fumin=float(np.median([r['sleep'] for r in f])) if f else None,
                 sleep_med_sss=float(np.median([r['sleep'] for r in s])) if s else None,
                 c_med_fumin=float(np.median([r['c_obs'] for r in f])) if f else None, c_med_sss=float(np.median([r['c_obs'] for r in s])) if s else None,
                 J_med_fumin=float(np.median([r['J'] for r in f if r['cycle']])) if f else None, J_med_sss=float(np.median([r['J'] for r in s if r['cycle']])) if s else None)
json.dump(dict(rows=rows, blocks=res, strat=dict(z=float(z0), rho_equiv=float(np.tanh(z0)), p=float(pz)), overlap=ov),
          open(sqw.OUTD + 's22_sleep_vs_J.json', 'w'))
print('block        n  | sleep med | rho(sleep,c_obs) p Holm (n) | rho(sleep,J | cycle) p Holm (n) | no-cycle n, sleep med vs cycle')
for k, v in res.items():
    print(f"{k:9s} {v['n']:3d} | {v['sleep_med']:5.1f} | {v['rho_c']:+.2f} p {v['p_c']:.4f} Holm {v.get('holm_c', float('nan')):.3f} ({v['n_c']}) | "
          f"{v['rho_J']:+.2f} p {v['p_J']:.4f} Holm {v.get('holm_J', float('nan')):.3f} ({v['n_J']}) | {v['n_nocycle']}, "
          f"{v['sleep_nocycle'] if v['sleep_nocycle'] is None else round(v['sleep_nocycle'], 1)} vs {v['sleep_cycle'] if v['sleep_cycle'] is None else round(v['sleep_cycle'], 1)}" + ('' if v['primary'] else '  [secondary]'))
print(f"stratified c_obs: mean z {z0:+.3f} (rho ~ {np.tanh(z0):+.2f}), within-block permutation p {pz:.4f}")
for p in (0, 1): print('overlap', ['light', 'dark'][p], ov[p])
# figure
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
COL = {'CS': '#2a78d6', 'per0DD': '#eda100', 'fumin': '#d43d51', 'sss': '#6a3d9a'}; NAME = {'CS': 'CS', 'per0DD': '$per^0$ DD', 'fumin': '$fumin$', 'sss': '$sss$'}
fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
for ax, key, yl in ((axs[0], 'J', 'cycle flux $J$ (cycles h$^{-1}$)\n+ forward W→Q→S→W'), (axs[1], 'c_obs', 'direction per transition $c_{obs}$\n+ forward')):
    for g, p in PRIM:
        B = [r for r in block(g, p) if (r['cycle'] if key == 'J' else (r['c_obs'] is not None and np.isfinite(r['c_obs'])))]
        fill = p == 1 or g == 'per0DD'; x = np.array([r['sleep'] for r in B]); yv = np.array([r[key] for r in B])
        ax.plot(x, yv, 'o', ms=3, mfc=COL[g] if fill else 'white', mec=COL[g], mew=0.6, alpha=0.85)
        if key == 'J':
            Z = [r for r in block(g, p) if not r['cycle']]
            ax.plot([r['sleep'] for r in Z], [0] * len(Z), 'x', ms=3, color=COL[g], mew=0.6, alpha=0.7)
    ax.axhline(0, color='0.6', lw=0.5); ax.set_xlabel('sleep (min per hour of the phase)'); ax.set_ylabel(yl)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
for g in COL: axs[0].plot([], [], 's', color=COL[g], label=NAME[g])
axs[0].plot([], [], 'o', mfc='white', mec='k', mew=0.6, label='light'); axs[0].plot([], [], 'o', color='k', label='dark / DD'); axs[0].plot([], [], 'x', color='k', label='no cycle')
fig.legend(loc='outside upper center', ncol=7, frameon=False, fontsize=7)
for a, l in zip(axs, 'ab'): a.text(-0.14, 1.02, l, transform=a.transAxes, fontsize=10, fontweight='bold')
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig22_J_vs_sleep.{ext}', dpi=300)
print('figure saved')
