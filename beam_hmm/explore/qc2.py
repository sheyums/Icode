import numpy as np
d = np.load('data/week1/week1.npz'); X, T = d['X'], d['T']
ok = ~np.isnan(T); D = X[ok]; t = T[ok]
for j in list(range(0,16,3))+list(range(16,32,2)):
    x = D[:, j]; v, c = np.unique(x, return_counts=True)
    print(f'ch{j+1:2d}', ' '.join(f'{a}:{b/len(x):.3f}' for a, b in zip(v, c)))
# transition matrix of values for ch1 and ch17
for j in (0, 5, 20, 30):
    x = D[:, j].astype(int); a, b = x[:-1], x[1:]; m = a != b
    pairs, c = np.unique(np.c_[a[m], b[m]], axis=0, return_counts=True)
    top = np.argsort(-c)[:14]
    print(f'ch{j+1} top transitions', [(tuple(pairs[i]), c[i]) for i in top])
# irregular gaps: where
g = np.where(np.round(np.diff(t),4) > 0.4)[0]; print('gaps at t=', t[g], np.diff(t)[g])
# 33-64 changes: when
for j in (32, 48):
    x = D[:, j]; ix = np.nonzero(np.diff(x))[0]; print(f'ch{j+1} changes at t', t[ix], x[ix], x[ix+1])
