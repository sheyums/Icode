import numpy as np
d = np.load('data/week1/week1.npz'); X, T = d['X'], d['T']
neg = (X == -1).all(1)
print('all -1 lines', neg.sum(), 'time NaN', np.isnan(T).sum(), 'overlap', (neg & np.isnan(T)).sum())
print('lines with time but all -1:', (neg & ~np.isnan(T)).sum(), ' no time but data:', (~neg & np.isnan(T)).sum())
# pattern of line types
k = (~np.isnan(T)).astype(int); print('first 20 types', k[:20], 'strictly alternating:', (np.diff(k)!=0).all())
D = X[~np.isnan(T)]; t = T[~np.isnan(T)]
dt = np.diff(t); u, c = np.unique(np.round(dt, 4), return_counts=True); print('time steps', dict(zip(u, c)))
print('t first/last', t[:3], t[-1], 'n', len(t))
print('values present', np.unique(D))
partial = ((X == -1).any(1) & ~neg).sum(); print('lines with some -1:', partial)
for j in range(64):
    x = D[:, j]; v = np.unique(x)
    ch = np.count_nonzero(np.diff(x))
    jumps = np.abs(np.diff(x.astype(int))); big = (jumps > 1).sum()
    print(f'ch{j+1:2d} nvals {len(v):2d} range {x.min()}-{x.max()} changes {ch:7d} jumps>1 {big:6d} maxjump {jumps.max():2d} mean pos {x.mean():5.2f}')
