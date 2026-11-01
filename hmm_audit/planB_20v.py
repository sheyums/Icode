"""PRESPEC 20v: waiting-time irreversibility beyond the observed first-order chain.
Usage: python3 -I planB_20v.py <data root with genotype subfolders> <CS DD folder>
Null A exactly as 20p (sqw15.fly): mc_fit per phase, sim_rows, sqw_ref.simulate with B = I, missing minutes reinserted.
LD descriptive: also the 6-h-rate null (mc_fit on phases(T, '6h'))."""
import sys, os, glob, json, numpy as np
from multiprocessing import Pool
from scipy.stats import wilcoxon, binom
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R
Q, W, S = 0, 1, 2
EDGES = np.array([1, 2, 3, 4, 6, 10, 10**9])           # bins 1, 2, 3, 4-5, 6-9, >= 10 minutes
PAIRS = {Q: [(W, S), (S, W)], W: [(Q, S), (S, Q)]}
NNULL = 200
SETS = {'per0DD': ('per0DD_activity_files_1sec_res', False), 'CSDD': (None, False),
        'CS': ('CS_activity_files_1sec_res', True), 'per0LD': ('per0LD_activity_1sec_res', True),
        'fumin': ('Fmn_activity_files_1sec', True), 'sss': ('sleepless_activity_1sec_res', True)}
PRIMARY = ('per0DD', 'CSDD')
I3 = np.eye(3)

def files(d): return sorted(f for f in glob.glob(os.path.join(d, '**', '*'), recursive=True) if os.path.isfile(f))

def sim_rows(P):
    P = np.array(P, float)
    for p in range(P.shape[0]):
        for k in range(3):
            if P[p, k].sum() < 0.5: P[p, k] = 1 / 3
    return P

def dwell_hists(y, ph, nph):
    """Per phase: dict (s, e, x) -> histogram of dwell lengths over the 6 bins; runs + neighbours inside one phase."""
    y = np.asarray(y); ph = np.asarray(ph); T = len(y)
    ch = np.flatnonzero(np.diff(y)) + 1; st = np.r_[0, ch]; en = np.r_[ch, T]
    out = [{(s, e, x): np.zeros(6) for s in PAIRS for (e, x) in PAIRS[s]} for _ in range(nph)]
    for a, b in zip(st, en):
        s = y[a]
        if s not in PAIRS or a == 0 or b == T: continue
        e, x = y[a - 1], y[b]
        if (e, x) not in PAIRS[s]: continue
        p = ph[a - 1]
        if ph[b] != p or np.any(ph[a:b] != p): continue
        out[p][(s, e, x)][np.searchsorted(EDGES, b - a, side='right') - 1] += 1
    return out

def sigma_wt(hists, hours):
    tot = 0.0
    for s in PAIRS:
        (a, b) = PAIRS[s]
        ha, hb = hists[(s,) + a] + 0.5, hists[(s,) + b] + 0.5
        pa, pb = ha / ha.sum(), hb / hb.sum()
        na, nb = hists[(s,) + a].sum(), hists[(s,) + b].sum()
        tot += na / hours * np.sum(pa * np.log(pa / pb)) + nb / hours * np.sum(pb * np.log(pb / pa))
    return tot

def stats(y, ph, nph):
    H = dwell_hists(y, ph, nph); ph = np.asarray(ph)
    return [sigma_wt(H[p], np.sum(ph == p) / 60) for p in range(nph)]

def fly(args):
    g, i, path, ld = args
    y, ref = R.sqw_symbols(R.read_counts(path), theta=0.5, light_dark=ld)
    T = len(y); nP = 2 if ld else 1; ph = R.phases(T, '12h') if ld else np.zeros(T, int); miss = y == -1
    okph = ([np.isfinite(ref['light']), np.isfinite(ref['dark'])] if ld else [np.isfinite(ref['whole'])])
    obs = stats(y, ph, nP)
    nulls = {}
    specs = [('A', ph, nP)] + ([('B6', R.phases(T, '6h'), 4)] if ld else [])
    rng = np.random.default_rng(20_000 + 1000 * list(SETS).index(g) + i)
    for name, pp, n in specs:
        _, pi, P = R.mc_fit(y, pp, n); P = sim_rows(P); L = []
        for _ in range(NNULL):
            _, ys = R.simulate(pi, P, I3, pp, rng); ys[miss] = -1; L.append(stats(ys, ph, nP))
        nulls[name] = np.array(L)
    res = dict(g=g, fly=i, file=os.path.basename(path), ok=[bool(v) for v in okph], obs=obs)
    for name, L in nulls.items():
        res['E_' + name] = [float(obs[p] - np.median(L[:, p])) for p in range(nP)]
        res['p_' + name] = [float((1 + np.sum(L[:, p] >= obs[p])) / (NNULL + 1)) for p in range(nP)]
    return res

def group(rows, ph, key):
    r = [x for x in rows if x['ok'][ph]]
    E = np.array([x['E_' + key][ph] for x in r]); p = np.array([x['p_' + key][ph] for x in r])
    w = wilcoxon(E).pvalue if len(E) > 0 and np.any(E != 0) else np.nan
    k = int(np.sum(p < 0.05)); bp = binom.sf(k - 1, len(p), 0.05)
    return dict(n=len(E), E_med=float(np.median(E)), wil=float(w), frac=float(k / len(p)), binom=float(bp),
                obs_med=float(np.median([x['obs'][ph] for x in r])))

def main(root, csdd):
    jobs = []
    for g, (sub, ld) in SETS.items():
        d = csdd if g == 'CSDD' else os.path.join(root, sub)
        for i, f in enumerate(files(d)): jobs.append((g, i, f, ld))
    with Pool(4) as pool:
        rows = pool.map(fly, jobs, chunksize=1)
    json.dump(rows, open(os.path.join(HERE, 'planB_20v.json'), 'w'))
    prim = {g: group([r for r in rows if r['g'] == g], 0, 'A') for g in PRIMARY}
    ps = sorted((prim[g]['wil'], g) for g in PRIMARY); m = len(ps); holm = {}; run = 0
    for k, (pv, g) in enumerate(ps):
        run = max(run, min(1, (m - k) * pv)); holm[g] = run
    for g in PRIMARY:
        s = prim[g]; sup = holm[g] < 0.05 and s['E_med'] > 0 and s['binom'] < 0.05
        print(f"{g:7s} whole n {s['n']} [PRIMARY] | sigma_WT obs median {s['obs_med']:.4f} nats/h | E_WT median {s['E_med']:+.4f} | "
              f"Wilcoxon p {s['wil']:.3g} Holm {holm[g]:.3g} | flies p<.05 {s['frac']:.2f} (binom {s['binom']:.3g}) -> "
              f"{'SUPPORTED' if sup else 'NOT SUPPORTED'}")
    for g in ('CS', 'per0LD', 'fumin', 'sss'):
        for ph, nm in ((0, 'light'), (1, 'dark')):
            a = group([r for r in rows if r['g'] == g], ph, 'A'); b = group([r for r in rows if r['g'] == g], ph, 'B6')
            print(f"{g:7s} {nm:5s} n {a['n']} [descriptive] | sigma_WT obs {a['obs_med']:.4f} | 12-h null: E {a['E_med']:+.4f} "
                  f"p {a['wil']:.3g} frac {a['frac']:.2f} | 6-h null: E {b['E_med']:+.4f} p {b['wil']:.3g} frac {b['frac']:.2f}")

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
