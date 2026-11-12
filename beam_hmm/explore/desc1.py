import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = np.load('data/week1/week1.npz'); X, T = d['X'], d['T']
ok = ~np.isnan(T); D = X[ok].astype(int); t = T[ok]
clean = [18, 19, 22, 24, 26, 28, 29]   # monitor-2 channels with >=95% of changes <= 2 beams (qc3)
keep = t < 150*3600                     # before the hour-153 event
zt = ((t/3600) % 24).astype(int)
fig, axs = plt.subplots(1, len(clean), figsize=(16, 4), sharey=True)
for a, c in zip(axs, clean):
    x = D[keep, c-1]; z = zt[keep]
    M = np.array([np.bincount(x[z == h], minlength=17)/max((z == h).sum(), 1) for h in range(24)]).T
    a.imshow(M, aspect='auto', origin='lower', cmap='magma', vmax=0.5); a.set_title(f'ch{c}'); a.set_xlabel('ZT h'); a.axvline(11.5, color='c', lw=0.8)
axs[0].set_ylabel('beam value (0-16)'); plt.suptitle('Occupancy by ZT hour, week 1 days 1-6 (exploration, descriptive)'); plt.tight_layout(); plt.savefig('explore/fig_occupancy_zt.png', dpi=110)
# pauses: durations of runs of identical position, pooled clean channels, by sub-5-min classes
dur = []
for c in clean:
    x = D[keep, c-1]; ch = np.r_[0, np.nonzero(np.diff(x))[0]+1, len(x)]; r = np.diff(ch)*0.4; dur.append(r)
dur = np.concatenate(dur)
edges = [0, 1, 5, 30, 60, 300, 1800, 1e9]
h = np.histogram(dur, edges)[0]; tw = np.histogram(dur, edges, weights=dur)[0]
for i in range(len(h)): print(f'still {edges[i]:>5}-{edges[i+1]:<6}s : runs {h[i]:7d}  share of time {tw[i]/tw.sum():.3f}')
