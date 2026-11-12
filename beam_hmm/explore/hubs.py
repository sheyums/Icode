"""E0 follow-up: structure of the jump artifacts in the artifact-heavy channels (week 1, descriptive)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, pipeline as pl
X, t = pl.load_week('data/week1/week1.npz')
names = [c + 1 for c in pl.CH]
for j, c in enumerate(names):
    if c == 30: continue
    x, _ = pl.spike_filter(X[:, j]); d = np.diff(x.astype(int)); big = np.abs(d) >= 4
    if big.sum() < 50: continue
    into = x[1:][big]; v, n = np.unique(into, return_counts=True); o = np.argsort(-n)
    hub = v[o[0]]; share = n[o[0]] / big.sum()
    # stays at the hub after a big jump into it, and where it goes next
    idx = np.nonzero(big & (x[1:] == hub))[0] + 1
    stays, back = [], 0
    for i in idx[:5000]:
        k = i
        while k + 1 < len(x) and x[k + 1] == hub: k += 1
        stays.append((k - i + 1) * 0.4)
        if k + 1 < len(x) and x[k + 1] == x[i - 1]: back += 1
    nat = (x == hub).mean()
    print(f'ch{c}: big jumps {big.sum():6d} ({big.sum()/max((d!=0).sum(),1):.0%} of changes); top target {hub} takes {share:.0%}; '
          f'stay at hub after jump: median {np.median(stays):.1f}s, 90% {np.percentile(stays,90):.1f}s; returns to pre-jump value {back/len(stays):.0%}; '
          f'time at hub value overall {nat:.1%}; 2nd/3rd targets {list(v[o[1:3]])} ({n[o[1]]/big.sum():.0%}, {n[o[2]]/big.sum():.0%})')
