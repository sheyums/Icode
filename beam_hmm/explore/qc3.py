import numpy as np
d = np.load('data/week1/week1.npz'); X, T = d['X'], d['T']
ok = ~np.isnan(T); D = X[ok].astype(int)
print('|jump| distribution among changes (fraction), monitor 2 channels 17-32')
for j in range(16, 32):
    x = D[:, j]; dx = np.diff(x); m = dx != 0; a = np.abs(dx[m])
    h = np.bincount(a, minlength=17)[1:]/m.sum()
    # values involved in big jumps (>=4)
    big = np.where(np.abs(dx) >= 4)[0]
    vals = np.concatenate([x[big], x[big+1]]); vv, cc = np.unique(vals, return_counts=True); top = vv[np.argsort(-cc)[:3]]
    # A-B-A spikes: single-sample excursions
    spike = ((x[1:-1] != x[:-2]) & (x[2:] == x[:-2]) & (np.abs(x[1:-1]-x[:-2]) >= 4)).sum()
    print(f'ch{j+1}: 1:{h[0]:.2f} 2:{h[1]:.2f} 3:{h[2]:.2f} >=4:{h[3:].sum():.2f}  n_changes {m.sum():6d}  big-jump hub values {top}  single-sample spikes>=4: {spike}')
print('monitor 1: |jump| with values read as-is')
for j in range(16):
    x = D[:, j]; dx = np.diff(x); m = dx != 0; a = np.abs(dx[m]); h = np.bincount(a, minlength=4)[1:]/max(m.sum(),1)
    print(f'ch{j+1}: 1:{h[0]:.2f} 2:{h[1]:.2f} 3:{h[2]:.2f}  n {m.sum()}')
