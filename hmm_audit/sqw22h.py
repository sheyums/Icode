# POST HOC to PRESPEC 22 (Sheyum's request, 2026-11-03, after seeing the observed-level result): the same question at the
# hidden level. x = hidden sleep, 60 * stationary pi_S of the fly's phase matrix P3 (min/h). Outcomes: J_hmm (count J
# averaged over 200 posterior paths, cycles/h, as in fig_sss_reversal), and a rate-free direction, net forward turns per
# sleep onset, d = (N_SW - N_WS) / (N_QS + N_WS), computed on each path and then averaged (Jensen). Same blocks
# (minus CS DD, not fitted), same permutation null and Holm structure as 22. Labelled post hoc throughout.
import json, numpy as np, sqw, sqw2
from scipy.stats import rankdata
import sqw_ref as R
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = json.load(open(sqw.OUTD + 's_observable.json')); F = sqw2.load(0.5); RT = json.load(open(sqw.OUTD + 'route2_0.5.json'))
rng = np.random.default_rng(221); NP = 10000
PRIM = [('CS', 0), ('CS', 1), ('per0DD', 0), ('fumin', 0), ('fumin', 1), ('sss', 0), ('sss', 1)]
rows = []
for d in O:
    g, i = d['g'], d['fly']; h = F.get((g, i)); rt = RT.get(f'{g}|{i}')
    if g == 'CSDD' or h is None or h['deg'] or rt is None: continue
    for p, r in enumerate(d['main']):
        if not r['ok'] or r.get('J_hmm') is None: continue
        P = np.array(h['P3'])[p]; sl = 60.0 * R.stationary(P)[2]
        N = np.array(rt['pairs_path'])[:, p]                    # 200 x 3 x 3, indices Q=0, W=1, S=2
        on = N[:, 0, 2] + N[:, 1, 2]; net = N[:, 2, 1] - N[:, 1, 2]
        dd = np.where(on > 0, net / np.maximum(on, 1), np.nan)
        dval = float(np.nanmean(dd)) if np.mean(on > 0) >= 0.5 else np.nan
        rows.append(dict(g=g, fly=i, p=p, sleep=float(sl), J=float(r['J_hmm']), d=dval, cycle=abs(r['J_hmm']) >= 1e-12))
def rho(x, y): return np.corrcoef(rankdata(x), rankdata(y))[0, 1]
def perm(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); r0 = rho(x, y); rx = rankdata(x); ry = rankdata(y)
    null = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(NP)])
    return r0, (1 + np.sum(np.abs(null) >= abs(r0) - 1e-12)) / (NP + 1)
def block(g, p): return [r for r in rows if r['g'] == g and r['p'] == p]
res = {}
for g, p in PRIM:
    B = block(g, p); bd = [r for r in B if np.isfinite(r['d'])]; bj = [r for r in B if r['cycle']]
    rd, pd = perm([r['sleep'] for r in bd], [r['d'] for r in bd]); rj, pj = perm([r['sleep'] for r in bj], [r['J'] for r in bj])
    res[f'{g}|{p}'] = dict(n=len(B), n_d=len(bd), rho_d=rd, p_d=pd, n_J=len(bj), rho_J=rj, p_J=pj, sleep_med=float(np.median([r['sleep'] for r in B])))
fam = [(k, t) for k in res for t in ('d', 'J')]; ps = np.array([res[k]['p_' + t] for k, t in fam]); o = np.argsort(ps); adj = np.empty(len(ps)); run = 0
for j, idx in enumerate(o): run = max(run, min(1, (len(ps) - j) * ps[idx])); adj[idx] = run
for (k, t), a in zip(fam, adj): res[k]['holm_' + t] = float(a)
def strat(perm_it):
    zs, ws = [], []
    for g, p in PRIM:
        bd = [r for r in block(g, p) if np.isfinite(r['d'])]; x = np.array([r['sleep'] for r in bd]); yv = np.array([r['d'] for r in bd])
        if perm_it: x = rng.permutation(x)
        zs.append(np.arctanh(np.clip(rho(x, yv), -0.999999, 0.999999))); ws.append(len(bd) - 3)
    return np.average(zs, weights=ws)
z0 = strat(False); zn = np.array([strat(True) for _ in range(2000)]); pz = (1 + np.sum(np.abs(zn) >= abs(z0))) / 2001
ov = {}
for p in (0, 1):
    f = [r for r in block('fumin', p) if np.isfinite(r['d'])]; s = [r for r in block('sss', p) if np.isfinite(r['d'])]
    lo = max(min(r['sleep'] for r in f), min(r['sleep'] for r in s)); hi = min(max(r['sleep'] for r in f), max(r['sleep'] for r in s))
    f = [r for r in f if lo <= r['sleep'] <= hi]; s = [r for r in s if lo <= r['sleep'] <= hi]
    ov[p] = {k: round(float(v), 3) for k, v in dict(lo=lo, hi=hi, n_f=len(f), n_s=len(s), sleep_f=np.median([r['sleep'] for r in f]), sleep_s=np.median([r['sleep'] for r in s]),
             d_f=np.median([r['d'] for r in f]), d_s=np.median([r['d'] for r in s]), J_f=np.median([r['J'] for r in f]), J_s=np.median([r['J'] for r in s])).items()}
json.dump(dict(rows=rows, blocks=res, strat=dict(z=float(z0), p=float(pz)), overlap=ov), open(sqw.OUTD + 's22h_hidden.json', 'w'))
print('POST HOC hidden | block n | hidden sleep med | rho(sleep, d) p Holm (n) | rho(sleep, J_hmm | cycle) p Holm (n)')
for k, v in res.items():
    print(f"{k:9s} {v['n']:3d} | {v['sleep_med']:5.1f} | {v['rho_d']:+.2f} p {v['p_d']:.4f} Holm {v['holm_d']:.3f} ({v['n_d']}) | {v['rho_J']:+.2f} p {v['p_J']:.4f} Holm {v['holm_J']:.3f} ({v['n_J']})")
print(f'stratified d: z {z0:+.3f} (rho ~ {np.tanh(z0):+.2f}) p {pz:.4f}'); print('overlap', ov)
# figure: observed (top, from s22) and hidden (bottom)
S22 = json.load(open(sqw.OUTD + 's22_sleep_vs_J.json'))['rows']
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in'})
COL = {'CS': '#2a78d6', 'per0DD': '#eda100', 'fumin': '#d43d51', 'sss': '#6a3d9a'}; NAME = {'CS': 'CS', 'per0DD': '$per^0$ DD', 'fumin': '$fumin$', 'sss': '$sss$'}
fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.6), constrained_layout=True, sharex=True)
def panel(ax, R_, key, cyc_key, yl):
    for g, p in PRIM:
        B = [r for r in R_ if r['g'] == g and r['p'] == p]; fill = p == 1 or g == 'per0DD'
        V = [r for r in B if (r[cyc_key] if key == 'J' else (r[key] is not None and np.isfinite(r[key])))]
        ax.plot([r['sleep'] for r in V], [r[key] for r in V], 'o', ms=3, mfc=COL[g] if fill else 'white', mec=COL[g], mew=0.6, alpha=0.85)
        if key == 'J':
            Z = [r for r in B if not r[cyc_key]]; ax.plot([r['sleep'] for r in Z], [0] * len(Z), 'x', ms=3, color=COL[g], mew=0.6, alpha=0.7)
    ax.axhline(0, color='0.6', lw=0.5); ax.set_ylabel(yl)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
panel(axs[0, 0], S22, 'J', 'cycle', 'observed cycle flux $J$ (cycles h$^{-1}$)'); panel(axs[0, 1], S22, 'c_obs', 'cycle', 'observed direction per transition $c_{obs}$')
panel(axs[1, 0], rows, 'J', 'cycle', 'hidden cycle flux $J$ (cycles h$^{-1}$)'); panel(axs[1, 1], rows, 'd', 'cycle', 'hidden net forward turns\nper sleep onset')
for a in axs[1]: a.set_xlabel('sleep (min per hour of the phase)')
for a, l in zip(axs.flat, 'abcd'): a.text(-0.14, 1.02, l, transform=a.transAxes, fontsize=10, fontweight='bold')
axs[0, 0].set_title('observed chain', fontsize=8.5, loc='left'); axs[1, 0].set_title('hidden model (post hoc; x = model sleep)', fontsize=8.5, loc='left')
for g in COL: axs[0, 0].plot([], [], 's', color=COL[g], label=NAME[g])
axs[0, 0].plot([], [], 'o', mfc='white', mec='k', mew=0.6, label='light'); axs[0, 0].plot([], [], 'o', color='k', label='dark / DD'); axs[0, 0].plot([], [], 'x', color='k', label='no cycle')
fig.legend(loc='outside upper center', ncol=7, frameon=False, fontsize=7)
for ext in ('png', 'pdf'): fig.savefig(f'figures/fig22_J_vs_sleep_obs_hidden.{ext}', dpi=300)
print('figure saved')
