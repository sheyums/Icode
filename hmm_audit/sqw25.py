# PRESPEC 25 (f939ba6): direction of MEASURABLE cycles in 4-h / 6-h single-day windows within flies, against a
# clock-matched (primary) and a 12-h (secondary) first-order chain null, 500 simulated records per fly.
import json, numpy as np, sqw
from sqw15 import seqs
from scipy.stats import binom, wilcoxon
NR = 500; ALPHA = 0.05; rng = np.random.default_rng(25)
PRIM = ['per0LD', 'per0DD']; REF = ['CS', 'CSDD']; DESC = ['sss']; G = PRIM + REF + DESC
SIZES = [360, 240]   # minutes: 6 h, 4 h
def nfly(g): return len(sqw.files(g))
def chain(y, idx, nm, fallback=None):
    """Row-normalised transition counts per matrix index idx[t] (arrival minute). Empty rows -> fallback[m] row."""
    C = np.zeros((nm, 3, 3)); ok = (y[:-1] >= 0) & (y[1:] >= 0); np.add.at(C, (idx[1:][ok], y[:-1][ok], y[1:][ok]), 1); P = np.empty_like(C)
    for m in range(nm):
        for s in range(3):
            r = C[m, s].sum()
            P[m, s] = C[m, s] / r if r > 0 else (fallback[m][s] if fallback is not None else np.full(3, 1 / 3))
    return P
def simulate(P, idx, y0, T, R):
    Cc = np.cumsum(P, 2); Cc[..., -1] = 1.0; S = np.empty((R, T), np.int8); cur = np.full(R, y0); S[:, 0] = cur
    for t in range(1, T):
        cur = (rng.random(R)[:, None] > Cc[idx[t], cur]).sum(1); S[:, t] = cur
    return S
def window_counts(Y, L):
    """Y (R x T). Per window: N_SW, N_WS over consecutive minute pairs inside the window."""
    Y = np.atleast_2d(Y); R, T = Y.shape; sw = (Y[:, :-1] == 2) & (Y[:, 1:] == 1); ws = (Y[:, :-1] == 1) & (Y[:, 1:] == 2)
    t = np.arange(1, T); inside = (t % L) != 0          # pair (t-1, t) inside one window
    nw = T // L; a = np.zeros((R, nw), int); b = np.zeros((R, nw), int); w = t // L
    np.add.at(a, (slice(None), w[inside]), sw[:, inside]) if False else None
    for k in range(nw):
        m = inside & (w == k); a[:, k] = sw[:, m].sum(1); b[:, k] = ws[:, m].sum(1)
    return a, b
def classify(a, b):
    n = a + b; p = np.minimum(1.0, 2 * binom.cdf(np.minimum(a, b), n, 0.5)); meas = (p < ALPHA) & (n > 0)
    return np.where(meas & (a > b), 1, np.where(meas & (b > a), -1, 0))   # +1 forward, -1 backward, 0 not measurable
out = {}
for g in G:
    for L in SIZES:
        key = f'{g}|{L // 60}h'; flies = []; Btot_n = {'clock': np.zeros(NR), 'phase': np.zeros(NR)}; Ftot_n = {'clock': np.zeros(NR), 'phase': np.zeros(NR)}
        for i in range(nfly(g)):
            y, ph, _ = seqs(g, i); y = np.asarray(y).astype(int); ph = np.asarray(ph); T = len(y)
            miss = y < 0   # Amendment A (sss only has any): windows with a missing minute are excluded in data and simulations
            nP = int(ph.max()) + 1; Pph = chain(y, ph, nP)
            wpos = (np.arange(T) % 1440) // L; nwp = 1440 // L
            fb = np.array([Pph[ph[np.argmax(wpos == m)]] for m in range(nwp)])     # 12-h row for empty clock-window rows
            Pclk = chain(y, wpos, nwp, fb)
            a, b = window_counts(y, L); cls = classify(a, b)[0]; wbad = miss.reshape(-1, L).any(1); cls = np.where(wbad, 9, cls)
            rec = dict(fly=i, n_excluded=int(wbad.sum()), cls=cls.tolist(), nSW=a[0].tolist(), nWS=b[0].tolist(), B=int(np.sum(cls == -1)), F=int(np.sum(cls == 1)))
            for nm, P, idx in (('clock', Pclk, wpos), ('phase', Pph, ph)):
                S = simulate(P, idx, y[np.argmax(y >= 0)] if np.any(y >= 0) else 0, T, NR); sa, sb = window_counts(S, L); sc = classify(sa, sb); sc[:, wbad] = 9
                Bs = (sc == -1).sum(1); Fs = (sc == 1).sum(1); Btot_n[nm] += Bs; Ftot_n[nm] += Fs
                rec[f'B_null_{nm}'] = float(Bs.mean()); rec[f'F_null_{nm}'] = float(Fs.mean()); rec[f'anyB_null_{nm}'] = float(np.mean(Bs > 0))
            flies.append(rec)
        B = sum(f['B'] for f in flies); F = sum(f['F'] for f in flies); res = dict(n_flies=len(flies), windows_per_fly=len(flies[0]['cls']), B=B, F=F,
                 frac_B=B / (B + F) if B + F else None, flies_anyB=sum(f['B'] > 0 for f in flies))
        for nm in ('clock', 'phase'):
            Bn, Fn = Btot_n[nm], Ftot_n[nm]; fr = np.where(Bn + Fn > 0, Bn / np.maximum(Bn + Fn, 1), np.nan)
            d = np.array([f['B'] - f[f'B_null_{nm}'] for f in flies])
            try: pw = float(wilcoxon(d).pvalue) if np.any(d != 0) else 1.0
            except ValueError: pw = 1.0
            res[nm] = dict(B_null_mean=float(Bn.mean()), B_null_95=[float(x) for x in np.percentile(Bn, [2.5, 97.5])], p_B=float((1 + np.sum(Bn >= B)) / (NR + 1)),
                           F_null_mean=float(Fn.mean()), F_null_95=[float(x) for x in np.percentile(Fn, [2.5, 97.5])], p_F=float((1 + np.sum(Fn >= F)) / (NR + 1)),
                           frac_B_null_median=float(np.nanmedian(fr)), frac_B_null_95=[float(x) for x in np.nanpercentile(fr, [2.5, 97.5])],
                           flies_anyB_null=float(sum(f[f'anyB_null_{nm}'] for f in flies)), wilcoxon_p=pw,
                           B_null_dist=Bn.tolist(), F_null_dist=Fn.tolist())
        out[key] = dict(res=res, flies=flies); print(key, 'done', flush=True)
def holm(keys, nm='clock'):
    ps = np.array([out[k]['res'][nm]['p_B'] for k in keys]); o = np.argsort(ps); adj = np.empty(len(ps)); run = 0
    for j, i in enumerate(o): run = max(run, min(1, (len(ps) - j) * ps[i])); adj[i] = run
    for k, a in zip(keys, adj): out[k]['res'][nm]['holm_B'] = float(a)
for fam in (PRIM, REF):
    for nm in ('clock', 'phase'): holm([f'{g}|{L // 60}h' for g in fam for L in SIZES], nm)
json.dump(out, open(sqw.OUTD + 's25_windows.json', 'w'))
print('\nblock        flies win/fly | measurable F  B  (B frac) flies>=1B | CLOCK null: B mean [95%] p_B Holm | F mean p_F | Bfrac null | PHASE null: B mean p_B Holm, Bfrac null')
for k, v in out.items():
    r = v['res']; c = r['clock']; p = r['phase']
    print(f"{k:10s} {r['n_flies']:3d} {r['windows_per_fly']:3d} | {r['F']:4d} {r['B']:3d} ({(r['frac_B'] or 0):.3f}) {r['flies_anyB']:3d} | "
          f"{c['B_null_mean']:5.1f} [{c['B_null_95'][0]:.0f},{c['B_null_95'][1]:.0f}] p {c['p_B']:.3f} Holm {c.get('holm_B', float('nan')):.3f} | {c['F_null_mean']:5.1f} p {c['p_F']:.3f} | {c['frac_B_null_median']:.3f} | "
          f"{p['B_null_mean']:5.1f} p {p['p_B']:.3f} Holm {p.get('holm_B', float('nan')):.3f}, {p['frac_B_null_median']:.3f} | flies>=1B null {c['flies_anyB_null']:.1f}")
