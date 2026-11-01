"""PRESPEC 20x post hoc P1 (placebo arm for the six transitions) and P2 (active fraction within wake).
POST HOC, not confirmatory. Usage: python3 -I planB_20x_posthoc.py <deprivation folder>"""
import sys, os, json, numpy as np
from scipy.stats import wilcoxon
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R
from planB_20x import files, symbols_baseline_ref, mask, pair_counts_all

def wil(v):
    v = v[np.isfinite(v)]; return wilcoxon(v).pvalue if np.any(v != 0) else 1.0

def main(d):
    C = {k: [] for k in ('B', 'R', 'P2', 'P01')}; A = []
    for f in files(d):
        x = R.read_counts(f); y, refs, sl = symbols_baseline_ref(x); T = len(y)
        act = (x[:T * 60] > 0).reshape(T, 60).mean(1); nonS = y != R.S
        for k, days in (('B', [0, 1, 2]), ('R', [4]), ('P2', [2]), ('P01', [0, 1])):
            C[k].append(pair_counts_all(y, mask(T, days, 0, 360)))
        a = {}
        for w, (lo, hi) in (('ZT0-6', (0, 360)), ('ZT6-12', (360, 720)), ('ZT12-24', (720, 1440))):
            for k, days in (('B', [0, 1, 2]), ('R', [4]), ('P2', [2]), ('P01', [0, 1])):
                m = mask(T, days, lo, hi) & nonS; a[(w, k)] = float(np.median(act[m])) if m.any() else np.nan
        A.append(a)
    C = {k: np.array(v) for k, v in C.items()}; n = len(A); off = [(i, j) for i in range(3) for j in range(3) if i != j]
    def P(k, idx):
        c = C[k][idx].sum(0); return c / c.sum(1, keepdims=True)
    def all3(idx):
        rr = np.array([P('R', idx)[i, j] / P('B', idx)[i, j] for i, j in off])
        rp = np.array([P('P2', idx)[i, j] / P('P01', idx)[i, j] for i, j in off])
        return rr, rp, np.log(rr) - np.log(rp)
    rr, rp, ct = all3(np.arange(n)); rng = np.random.default_rng(2027); bs = []
    for _ in range(2000):
        idx = rng.integers(n, size=n); bs.append(np.array(all3(idx)))
    bs = np.array(bs)                    # (2000, 3, 6)
    pb = [min(1.0, 2 * min(np.mean(bs[:, 2, k] <= 0), np.mean(bs[:, 2, k] >= 0))) for k in range(6)]
    order = np.argsort(pb); holm = np.zeros(6); run = 0
    for i, k in enumerate(order): run = max(run, min(1, (6 - i) * pb[k])); holm[k] = run
    lab = 'QWS'; out = dict(P1=[], P2={})
    print("P1 (POST HOC): transition | recovery/baseline [CI] | placebo day2/days0-1 [CI] | log contrast [CI] p Holm -> reading")
    for k, (i, j) in enumerate(off):
        c1 = np.percentile(bs[:, 0, k], [2.5, 97.5]); c2 = np.percentile(bs[:, 1, k], [2.5, 97.5]); c3 = np.percentile(bs[:, 2, k], [2.5, 97.5])
        beyond = holm[k] < 0.05 and (c3[0] > 0 or c3[1] < 0)
        print(f"  {lab[i]}->{lab[j]} | {rr[k]:.2f} [{c1[0]:.2f}, {c1[1]:.2f}] | {rp[k]:.2f} [{c2[0]:.2f}, {c2[1]:.2f}] | "
              f"{ct[k]:+.2f} [{c3[0]:+.2f}, {c3[1]:+.2f}] p {pb[k]:.3f} Holm {holm[k]:.3f} -> "
              f"{'changes beyond the day-to-day trend (post hoc hint)' if beyond else 'not separable from the trend in morning sleep'}")
        out['P1'].append(dict(t=f'{lab[i]}->{lab[j]}', rec=rr[k], pla=rp[k], contrast=ct[k], p=pb[k], holm=holm[k], beyond=bool(beyond)))
    print("P2 (POST HOC): median active fraction of non-S minutes")
    for w in ('ZT0-6', 'ZT6-12', 'ZT12-24'):
        D = np.array([a[(w, 'R')] - a[(w, 'B')] for a in A]); Dp = np.array([a[(w, 'P2')] - a[(w, 'P01')] for a in A])
        g = np.isfinite(D) & np.isfinite(Dp); p1 = wil(D); p2 = wil((D - Dp)[g])
        ch = p1 < 0.05 and p2 < 0.05 and np.sign(np.median((D - Dp)[g])) == np.sign(np.nanmedian(D))
        tag = ('CHANGE beyond day-to-day variation (post hoc)' if ch else 'no change beyond day-to-day variation') if w == 'ZT0-6' else 'descriptive'
        print(f"  {w:7s} baseline median {np.nanmedian([a[(w,'B')] for a in A]):.3f} | Delta median {np.nanmedian(D):+.4f} p {p1:.3g} | "
              f"Delta_pl median {np.nanmedian(Dp):+.4f} | Delta vs Delta_pl paired p {p2:.3g} -> {tag}")
        out['P2'][w] = dict(D=D.tolist(), Dp=Dp.tolist(), p=float(p1), p_vs_pl=float(p2))
    json.dump(out, open(os.path.join(HERE, 'planB_20x_posthoc.json'), 'w'), default=float)

if __name__ == '__main__':
    main(sys.argv[1])
