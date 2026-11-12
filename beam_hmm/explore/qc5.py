import numpy as np
d = np.load('data/week1/week1.npz'); X, T = d['X'], d['T']
D = X[~np.isnan(T)].astype(int)
for j in range(16, 32):
    x = D[:, j]; z = x == 0
    if z.sum() == 0: print(f'ch{j+1}: no zeros'); continue
    st = np.nonzero(z[1:] & ~z[:-1])[0]; en = np.nonzero(~z[1:] & z[:-1])[0]
    before = x[st]; after = x[en+1] if len(en) else np.array([])
    runs = []
    s = np.r_[0, np.nonzero(np.diff(z.astype(int)))[0]+1, len(x)]
    runs = [(s[i+1]-s[i])*0.4 for i in range(len(s)-1) if z[s[i]]]
    vb, cb = np.unique(before, return_counts=True); va, ca = np.unique(after, return_counts=True)
    ones = np.nonzero(x[1:] == 1)[0]; prev1 = x[ones]; p1 = prev1[prev1 != 1]
    v1, c1 = np.unique(p1, return_counts=True)
    print(f'ch{j+1}: zero frac {z.mean():.3f} runs {len(runs)} median {np.median(runs):.1f}s max {max(runs)/60:.0f}min | before 0: {dict(zip(vb[np.argsort(-cb)][:3], np.sort(cb)[::-1][:3]))} after 0: {dict(zip(va[np.argsort(-ca)][:3], np.sort(ca)[::-1][:3]))} | entries to 1 from: {dict(zip(v1[np.argsort(-c1)][:3], np.sort(c1)[::-1][:3]))}')
