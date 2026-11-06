# POST HOC (Sheyum, 2026-11-06, after PRESPEC 24's result): sleep per day and phase, 23E10>Kir vs NC; same death rule.
import glob, numpy as np, sqw_ref as R
from multiprocessing import Pool
from scipy.stats import mannwhitneyu
from sqw24 import trim, G
def one(a):
    g, f = a; x, nd, _ = trim(np.loadtxt(f).astype(int))
    if nd < 1: return g, None
    y, _ = R.sqw_symbols(x, 0.5, light_dark=True); out = np.full((8, 2), np.nan)
    for d in range(nd):
        for p in range(2):
            seg = y[d * 1440 + p * 720: d * 1440 + (p + 1) * 720]; out[d, p] = 60 * np.mean(seg == 2)
    return g, out
if __name__ == '__main__':
    jobs = [(g, f) for g, d in G.items() for f in sorted(glob.glob(d + '/**/*.txt', recursive=True))]
    with Pool(4) as pool: res = pool.map(one, jobs)
    A = {g: np.array([o for gg, o in res if gg == g and o is not None]) for g in G}
    print('day  phase | Kir median (n) | NC median (n) | diff | MW two-sided p')
    for d in range(8):
        for p, ph in enumerate(['light', 'dark']):
            a = A['kir'][:, d, p]; b = A['nc'][:, d, p]; a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
            print(f'd{d}  {ph:5s} | {np.median(a):5.1f} ({len(a)}) | {np.median(b):5.1f} ({len(b)}) | {np.median(a) - np.median(b):+5.1f} | {mannwhitneyu(a, b).pvalue:.3g}')
    for lab, ds in (('days 0-3', slice(0, 4)), ('days 4-7', slice(4, 8))):
        for p, ph in enumerate(['light', 'dark']):
            a = np.nanmean(A['kir'][:, ds, p], 1); b = np.nanmean(A['nc'][:, ds, p], 1); a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
            print(f'{lab} {ph:5s}: Kir {np.median(a):.1f} ({len(a)}) vs NC {np.median(b):.1f} ({len(b)}) diff {np.median(a) - np.median(b):+.1f} p {mannwhitneyu(a, b).pvalue:.3g}')
