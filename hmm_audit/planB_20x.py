"""PRESPEC 20x: sleep deprivation, within-fly contrasts at ZT0-6 (day 4 vs days 0-2; placebo day 2 vs days 0-1).
Usage: python3 -I planB_20x.py <deprivation folder>"""
import sys, os, glob, json, numpy as np
from scipy.stats import wilcoxon, spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R
Q, W, S = 0, 1, 2; DAY = 86400; HALF = 43200

def files(d): return sorted(f for f in glob.glob(os.path.join(d, '**', '*'), recursive=True) if os.path.isfile(f))

def symbols_baseline_ref(x, theta=0.5):
    """sqw_symbols with light/dark references from days 0-2 only, applied to every day (PRESPEC 20x)."""
    x = np.asarray(x); nm = len(x) // 60; x = x[:nm * 60]
    sl = np.zeros(len(x), bool)
    for s, n in R.sleep_bouts(x, True): sl[s:s + n] = True
    xb = x[:3 * DAY]; ev = R.wake_bout_events(xb, R.sleep_bouts(xb, True))
    refs = {'light': R.wake_reference(ev, 0), 'dark': R.wake_reference(ev, 1)}
    ssec = sl.reshape(nm, 60).sum(1); act = (x > 0).reshape(nm, 60).mean(1)
    y = np.full(nm, W); y[ssec >= 30] = S
    for m in np.flatnonzero(ssec < 30):
        r = refs['light'] if (m * 60) % DAY < HALF else refs['dark']
        if np.isnan(r): y[m] = -1
        elif act[m] < theta * r: y[m] = Q
    nd = len(y) // 1440
    return y[:nd * 1440], refs, sl[:nd * DAY]

def mask(T, days, lo, hi):
    t = np.arange(T); d = t // 1440; m = t % 1440
    return np.isin(d, days) & (m >= lo) & (m < hi)

def readouts(y, msk):
    ph = np.where(msk, 0, 1); hours = msk.sum() / 60
    c = R.c_obs_minutes(y, ph, 2)[0]; N = R.observable_pair_counts(y, ph, 2)[0]
    J = R.j_cycle_count(N, hours); E, sub = R.epr_count_bound(N, hours)
    kind = None
    pairs = [(0, 1), (0, 2), (1, 2)]
    if any(N[i][j] == 0 and N[j][i] == 0 for i, j in pairs): kind = 'iii'
    elif all(N[i][j] == N[j][i] for i, j in pairs): kind = 'ii'
    elif sub: kind = 'i'
    return dict(c=float(c), J=float(J), EPR=float(E), kind=kind)

def pair_counts_all(y, msk):
    """All consecutive non-missing pairs (incl. self) arriving in the window."""
    t = np.arange(1, len(y)); ok = (y[:-1] >= 0) & (y[1:] >= 0) & msk[1:]
    C = np.zeros((3, 3)); np.add.at(C, (y[:-1][ok], y[1:][ok]), 1); return C

def sleep_min(sl, day, lo_min, hi_min):
    return sl[day * DAY + lo_min * 60: day * DAY + hi_min * 60].sum() / 60

def boot_ci_median(v, rng, B=2000):
    v = v[np.isfinite(v)]; bs = [np.median(v[rng.integers(len(v), size=len(v))]) for _ in range(B)]
    return np.percentile(bs, [2.5, 97.5])

def wil(v):
    v = v[np.isfinite(v)]
    return wilcoxon(v).pvalue if np.any(v != 0) else 1.0

def main(d):
    fl = files(d); rng = np.random.default_rng(2026); rows = []; noref = {'light': 0, 'dark': 0}; Tmin = 10**9
    CB = []; CR = []
    for f in fl:
        x = R.read_counts(f); Tmin = min(Tmin, len(x))
        y, refs, sl = symbols_baseline_ref(x); T = len(y)
        for k in noref: noref[k] += int(np.isnan(refs[k]))
        r = dict(file=os.path.basename(f))
        r['sleepB'] = np.mean([sleep_min(sl, dd, 0, 360) for dd in (0, 1, 2)]); r['sleepR'] = sleep_min(sl, 4, 0, 360)
        r['sleep2'] = sleep_min(sl, 2, 0, 360); r['sleep01'] = np.mean([sleep_min(sl, dd, 0, 360) for dd in (0, 1)])
        r['lost'] = np.mean([sleep_min(sl, dd, 900, 1440) for dd in (0, 1, 2)]) - sleep_min(sl, 3, 900, 1440)
        r['depr_sleep'] = sleep_min(sl, 3, 900, 1440)
        r['B'] = readouts(y, mask(T, [0, 1, 2], 0, 360)); r['R'] = readouts(y, mask(T, [4], 0, 360))
        r['P2'] = readouts(y, mask(T, [2], 0, 360)); r['P01'] = readouts(y, mask(T, [0, 1], 0, 360))
        r['desc'] = {}
        for nm, days_r, lo, hi in (('d4_ZT6-12', [4], 360, 720), ('d4_ZT12-24', [4], 720, 1440)) + tuple(
                (f'd{k}_ZT0-6', [k], 0, 360) for k in (5, 6, 7, 8)):
            rb = readouts(y, mask(T, [0, 1, 2], lo, hi)); rr = readouts(y, mask(T, days_r, lo, hi))
            sb = np.mean([sleep_min(sl, dd, lo, hi) for dd in (0, 1, 2)]); sr = sleep_min(sl, days_r[0], lo, hi)
            r['desc'][nm] = {k: rr[k] - rb[k] for k in ('c', 'J', 'EPR')} | {'sleep': sr - sb}
        CB.append(pair_counts_all(y, mask(T, [0, 1, 2], 0, 360))); CR.append(pair_counts_all(y, mask(T, [4], 0, 360)))
        rows.append(r)
    n = len(rows)
    print(f"20x: {n} flies | shortest trace {Tmin/86400:.4f} d | flies without a baseline reference: light {noref['light']}, dark {noref['dark']}")
    # ---- stage 0
    dS = np.array([r['sleepR'] - r['sleepB'] for r in rows]); p0 = wil(dS); ci = boot_ci_median(dS, rng)
    gate = np.median(dS) > 0 and p0 < 0.05
    print(f"STAGE 0: sleep in ZT0-6, day 4 - days 0-2: median {np.median(dS):+.1f} min [95% CI {ci[0]:+.1f}, {ci[1]:+.1f}] "
          f"(baseline median {np.median([r['sleepB'] for r in rows]):.0f} min) | Wilcoxon p {p0:.3g} -> {'PASS' if gate else 'FAIL'}")
    out = dict(rows=rows, stage0=dict(dS=dS.tolist(), p=float(p0), PASS=bool(gate)))
    if not gate:
        print('no measurable rebound: stop (PRESPEC 20x)'); json.dump(out, open(os.path.join(HERE, 'planB_20x.json'), 'w'), default=float); return
    # ---- primary
    names = {'c': 'c_obs', 'J': 'count J (cycles/h)', 'EPR': 'count EPR (nats/h)'}; res = {}
    for k in names:
        D = np.array([r['R'][k] - r['B'][k] for r in rows]); Dp = np.array([r['P2'][k] - r['P01'][k] for r in rows])
        g = np.isfinite(D) & np.isfinite(Dp)
        res[k] = dict(D=D, Dp=Dp, p=wil(D), ppl=wil(Dp), pvs=wil((D - Dp)[g]), medD=float(np.nanmedian(D)),
                      medDp=float(np.nanmedian(Dp)), medDiff=float(np.median((D - Dp)[g])), ci=boot_ci_median(D, rng),
                      n=int(np.isfinite(D).sum()), base=float(np.nanmedian([r['B'][k] for r in rows])))
    ps = sorted((res[k]['p'], k) for k in names); run = 0
    for i, (pv, k) in enumerate(ps): run = max(run, min(1, (3 - i) * pv)); res[k]['holm'] = run
    for k in names:
        a = res[k]; sup = a['holm'] < 0.05 and a['pvs'] < 0.05 and np.sign(a['medDiff']) == np.sign(a['medD']) and a['medD'] != 0
        a['verdict'] = 'SUPPORTED' if sup else 'NOT SUPPORTED'
        print(f"{names[k]:20s} n {a['n']} | baseline median {a['base']:+.3f} | Delta median {a['medD']:+.3f} [95% CI {a['ci'][0]:+.3f}, {a['ci'][1]:+.3f}] "
              f"p {a['p']:.3g} Holm {a['holm']:.3g} | vs placebo: Delta_pl median {a['medDp']:+.3f}, paired p {a['pvs']:.3g} | "
              f"time-trend check Delta_pl vs 0 p {a['ppl']:.3g}{' FLAG' if a['ppl'] < 0.05 else ''} -> {a['verdict']}")
    kinds = {}
    for r in rows:
        for c in ('B', 'R', 'P2', 'P01'):
            kd = r[c]['kind']
            if kd: kinds[f'{c}:{kd}'] = kinds.get(f'{c}:{kd}', 0) + 1
    print(f"  no-cycle / rule-of-three kinds (condition:kind -> flies): {kinds or 'none'}")
    # ---- transitions
    CB = np.array(CB); CR = np.array(CR); off = [(i, j) for i in range(3) for j in range(3) if i != j]
    def ratios(idx):
        b = CB[idx].sum(0); r_ = CR[idx].sum(0); Pb = b / b.sum(1, keepdims=True); Pr = r_ / r_.sum(1, keepdims=True)
        return np.array([Pr[i, j] / Pb[i, j] for i, j in off])
    pt = ratios(np.arange(n)); bs = np.array([ratios(rng.integers(n, size=n)) for _ in range(2000)])
    pb = [min(1.0, 2 * min(np.mean(bs[:, k] <= 1), np.mean(bs[:, k] >= 1))) for k in range(6)]
    order = np.argsort(pb); holm = np.zeros(6); run = 0
    for i, k in enumerate(order): run = max(run, min(1, (6 - i) * pb[k])); holm[k] = run
    lab = 'QWS'
    for k, (i, j) in enumerate(off):
        lo, hi = np.percentile(bs[:, k], [2.5, 97.5])
        print(f"  P({lab[i]}->{lab[j]}) recovery/baseline {pt[k]:.2f} [95% CI {lo:.2f}, {hi:.2f}] bootstrap p {pb[k]:.3f} Holm {holm[k]:.3f}"
              f"{'  -> sleep pressure changes this transition' if holm[k] < 0.05 else ''}")
    # ---- dose association
    lost = np.array([r['lost'] for r in rows]); eff = np.array([r['depr_sleep'] <= 108 for r in rows])
    for k in names:
        D = res[k]['D']; g = np.isfinite(D)
        ra = spearmanr(lost[g], D[g]); re_ = spearmanr(lost[g & eff], D[g & eff])
        print(f"  dose association {names[k]:20s}: all flies rho {ra.statistic:+.2f} (p {ra.pvalue:.2f}, n {g.sum()}) | "
              f"efficacious rho {re_.statistic:+.2f} (p {re_.pvalue:.2f}, n {(g & eff).sum()})")
    ra = spearmanr(lost, dS); print(f"  dose association sleep rebound: all flies rho {ra.statistic:+.2f} (p {ra.pvalue:.2f})")
    # ---- descriptive windows
    print("  descriptive (median Delta vs the same window of days 0-2): window | sleep min | c_obs | J | EPR")
    for nm in rows[0]['desc']:
        v = {k: np.nanmedian([r['desc'][nm][k] for r in rows]) for k in ('sleep', 'c', 'J', 'EPR')}
        print(f"    {nm:11s} | {v['sleep']:+6.1f} | {v['c']:+.3f} | {v['J']:+.3f} | {v['EPR']:+.3f}")
    out.update(primary={k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv) for kk, vv in res[k].items()} for k in res},
               transitions=dict(ratio=pt.tolist(), p=pb, holm=holm.tolist()))
    json.dump(out, open(os.path.join(HERE, 'planB_20x.json'), 'w'), default=float)

if __name__ == '__main__':
    main(sys.argv[1])
