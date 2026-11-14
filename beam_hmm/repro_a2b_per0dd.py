"""Reproduce hmm_audit PRESPEC 21 (planB_21.py) on the per0 DD 1-s files, with the ORIGINAL code and with the A2b port.
Usage: python -I repro_a2b_per0dd.py <dir with hmm_audit/{sqw_ref,planB_21,planB_20u}.py and analysis_20r/*.mat> <per0DD data dir>
The reference files come from origin/claude/trusting-allen-7nxrsf (git show). Expected (planB_21.log, per0DD):
L300 z med +0.58, Wilcoxon 0.0136, flies p<.05 0.32, day-block z med +0.93; L120 z med +0.95, Wilcoxon 0.000912, 0.44."""
import sys, os, glob, numpy as np
ref, data = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(ref, 'hmm_audit')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R, planB_21 as Bo, planB_20u as U, pipeline as pl
_, _, names, tl = U.load_ref(); up = sorted(glob.glob(os.path.join(data, '*')))
X = []
for i, nm in enumerate(names):
    m = [f for f in up if os.path.basename(f).startswith(nm)]; assert len(m) == 1, (nm, m)
    x = R.read_counts(m[0])[:tl[i]]; X.append(x[:(len(x) // 86400) * 86400])

def port_windows(x, L, b):
    a = (np.asarray(x) > 0).astype(np.int8); c = np.asarray(x, float); sb = pl.sleep_bouts(a)
    D, Dc, blk, lev = [], [], [], []
    for k in range(len(sb) - 1):
        w0 = sb[k][0] + sb[k][1]; w1 = sb[k + 1][0]
        if w1 - w0 < 2 * L: continue
        for arr, lst in ((a.astype(float), D), (c, Dc)):
            lst.append(arr[w0:w0 + L].reshape(-1, b).mean(1) - arr[w1 - L:w1][::-1].reshape(-1, b).mean(1))
        lev.append(a[w0:w0 + 60].mean() - a[w1 - 60:w1].mean()); blk.append(((w0 + w1) // 2) // 86400)
    return np.array(D), np.array(Dc), np.array(blk), np.array(lev)

for label, fn in (('ORIGINAL', Bo.windows), ('A2B PORT', port_windows)):
    for st, (L, b) in Bo.SETTINGS.items():
        rows = []
        for i, x in enumerate(X):
            D, Dc, blk, lev = fn(x, L, b)
            if len(D) >= 5:
                r = Bo.fly_stats(D, blk, np.random.default_rng(27_000 + 1000 * Bo.GIDX['per0DD'] + i)); r['md'] = D.mean(0); rows.append(r)
        g = Bo.group(rows, 'z', 'p'); gb = Bo.group(rows, 'zb', 'pb'); md = np.median([r['md'] for r in rows], 0)
        print(f"{label} {st}: flies {len(rows)} | z med {g['zmed']:+.2f} Wilcoxon {g['wil']:.3g} flies p<.05 {g['frac']:.2f} | "
              f"day-block z med {gb['zmed']:+.2f} p {gb['wil']:.3g} | d(tau) " + ' '.join(f'{v:+.3f}' for v in md[:6]))
