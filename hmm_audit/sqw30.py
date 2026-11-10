# PRESPEC 30 (a439a0c): long-bout (>= 25 / 30 min) sleep after 9-h deprivation; within-fly contrasts with the
# day-to-day placebo (day 2 vs days 0-1). Symbols for routes as 20x (baseline-referenced).
import sys, glob, json, numpy as np
sys.path.insert(0, '/home/user/Icode/hmm_audit')
import planB_20x as P20
import sqw_ref as R
from scipy.stats import wilcoxon
D = '/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/data_depr'; DAY = 86400; Q, W, S = 0, 1, 2
def overlap(bouts, lo, hi): return sum(max(0, min(s + n, hi) - max(s, lo)) for s, n in bouts) / 60.0
def fly(f):
    x = np.loadtxt(f).astype(int); bouts = R.sleep_bouts(x, True); y, refs, _ = P20.symbols_baseline_ref(x); nd = len(x) // DAY; out = {}
    for d in range(min(nd, 9)):
        lo, hi = d * DAY, (d + 1) * DAY; L25 = [b for b in bouts if b[1] >= 1500]; L30 = [b for b in bouts if b[1] >= 1800]
        start = [b for b in bouts if lo <= b[0] < hi]; lstart = [b for b in L25 if lo <= b[0] < hi]
        def route_exit(b):
            m = (b[0] + b[1]) // 60
            while m < len(y) and (y[m] == S or y[m] < 0): m += 1
            return None if m >= len(y) else int(y[m] == W)
        def route_entry(b):
            m = b[0] // 60 - 1
            while m >= 0 and (y[m] == S or y[m] < 0): m -= 1
            return None if m < 0 else int(y[m] == Q)
        ex = [r for r in map(route_exit, lstart) if r is not None]; en = [r for r in map(route_entry, lstart) if r is not None]
        out[d] = dict(L25=overlap(L25, lo, hi), L30=overlap(L30, lo, hi), all=overlap(bouts, lo, hi), L25_zt06=overlap(L25, lo, lo + 21600),
                      dur=float(np.mean([b[1] for b in start]) / 60) if start else np.nan, nlong=len(lstart),
                      exitW=float(np.mean(ex)) if ex else np.nan, entryQ=float(np.mean(en)) if en else np.nan)
        out[d]['share'] = out[d]['L25'] / out[d]['all'] if out[d]['all'] > 0 else np.nan
    return out
def contrast(F, key, rec):
    D_ = np.array([f[rec][key] - np.nanmean([f[d][key] for d in (0, 1, 2)]) for f in F]); Dp = np.array([f[2][key] - np.nanmean([f[d][key] for d in (0, 1)]) for f in F])
    g = np.isfinite(D_) & np.isfinite(Dp); D_, Dp = D_[g], Dp[g]
    p1 = wilcoxon(np.round(D_, 10)).pvalue if np.any(D_ != 0) else 1.0; p2 = wilcoxon(np.round(D_ - Dp, 10)).pvalue if np.any(D_ != Dp) else 1.0
    return dict(key=key, day=rec, n=int(g.sum()), Delta_med=float(np.median(D_)), Delta_pl_med=float(np.median(Dp)), excess_med=float(np.median(D_ - Dp)), p_Delta=float(p1), p_excess=float(p2))
def holm(ps):
    ps = np.asarray(ps); o = np.argsort(ps); a = np.empty(len(ps)); run = 0
    for j, i in enumerate(o): run = max(run, min(1, (len(ps) - j) * ps[i])); a[i] = run
    return a
if __name__ == '__main__':
    from multiprocessing import Pool
    fl = sorted(glob.glob(D + '/**/*', recursive=True)); fl = [f for f in fl if not f.endswith('/') and '.' in f.rsplit('/', 1)[1] or True]
    import os; fl = [f for f in fl if os.path.isfile(f)]
    with Pool(4) as pool: F = pool.map(fly, fl)
    print('flies', len(F))
    prim = contrast(F, 'L25', 4); prim['rebound'] = bool(prim['p_Delta'] < 0.05 and prim['p_excess'] < 0.05 and prim['Delta_med'] > 0 and prim['excess_med'] > 0)
    sec = [contrast(F, 'L25', 5), contrast(F, 'L25', 6), contrast(F, 'L30', 4), contrast(F, 'share', 4), contrast(F, 'dur', 4), contrast(F, 'L25_zt06', 4), contrast(F, 'exitW', 4), contrast(F, 'entryQ', 4)]
    for s, a in zip(sec, holm([s['p_excess'] for s in sec])): s['holm_excess'] = float(a)
    course = {k: [[float(np.nanmedian([f[d][k] for f in F])) for d in range(9)]] for k in ('L25', 'all', 'dur', 'nlong', 'share')}
    json.dump(dict(primary=prim, secondary=sec, course=course, flies=[{str(d): v for d, v in f.items()} for f in F]), open('sqw_out/s30_longbouts.json', 'w'), default=float)
    print('PRIMARY', prim)
    for s in sec: print('SEC', s)
    for k, v in course.items(): print('course', k, [round(x, 2) for x in v[0]])
