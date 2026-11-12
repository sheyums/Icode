import numpy as np
d = np.load('data/week1/week1.npz'); X, T = d['X'], d['T']
ok = ~np.isnan(T); D = X[ok].astype(int); t = T[ok]
for j in range(32):
    x = D[:, j]; ix = np.nonzero(np.diff(x))[0]
    if len(ix) < 3: print(f'ch{j+1}: {len(ix)} changes in the week; value {np.unique(x)}'); continue
    tc = t[ix+1]; gaps = np.diff(np.r_[t[0], tc, t[-1]])
    k = np.argmax(gaps); start = np.r_[t[0], tc][k]
    print(f'ch{j+1}: longest silence {gaps[k]/3600:5.1f} h starting {start/3600:6.1f} h; last change {tc[-1]/3600:6.1f} h; silences>=12h: {(gaps>=43200).sum()}')
# what happens near hour 150-158 across channels: changes per 10 min
w = (t//600).astype(int); M = (np.diff(D[:, :32], axis=0) != 0)
cnt = np.array([np.bincount(w[1:], weights=M[:, j], minlength=w.max()+1) for j in range(32)])
seg = slice(150*6, 160*6)
print('10-min change counts, hours 150-160, channels with any change:')
for j in range(32): print(f'ch{j+1:2d}', ''.join('.' if c == 0 else ('o' if c < 50 else 'O') for c in cnt[j, seg]))
