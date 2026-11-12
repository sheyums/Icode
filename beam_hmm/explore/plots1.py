import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = np.load('data/week1/week1.npz'); X, T = d['X'], d['T']
ok = ~np.isnan(T); D = X[ok].astype(int); t = T[ok]
zt = (t/3600) % 24; day = (t//86400).astype(int)
mv = np.vstack([np.zeros(32, bool), D[1:, :32] != D[:-1, :32]])
# activity: fraction of samples with a change, per hour of the week
hr = (t//3600).astype(int); H = hr.max()+1
act = np.array([np.bincount(hr, weights=mv[:, j], minlength=H)/np.bincount(hr, minlength=H) for j in range(32)])
fig, ax = plt.subplots(figsize=(12, 7)); im = ax.imshow(act, aspect='auto', cmap='viridis', vmax=np.percentile(act, 99))
for k in range(0, H, 12): ax.axvline(k-0.5, color='w' if (k//12) % 2 else 'k', lw=0.5)
ax.set_xlabel('hour from file start (lines every 12 h: lights on/off as stated)'); ax.set_ylabel('channel'); ax.set_yticks(range(32)); ax.set_yticklabels(range(1, 33), fontsize=6)
plt.colorbar(im, label='fraction of 0.4-s samples with a position change'); plt.tight_layout(); plt.savefig('explore/fig_activity_heatmap.png', dpi=110)
# raw traces day 2 for some channels
fig, axs = plt.subplots(8, 1, figsize=(13, 12), sharex=True)
for a, j in zip(axs, [0, 5, 9, 12, 18, 20, 28, 30]):
    m = day == 1; a.plot(zt[m], D[m, j], lw=0.3, color='k'); a.set_ylabel(f'ch{j+1}'); a.set_ylim(-0.5, 16.5); a.axvspan(12, 24, color='0.85', zorder=0)
axs[-1].set_xlabel('ZT (h), day 2 of week 1; grey = dark'); plt.tight_layout(); plt.savefig('explore/fig_traces_day2.png', dpi=110)
# zoomed 10 min traces
fig, axs = plt.subplots(4, 1, figsize=(13, 8), sharex=True)
for a, j in zip(axs, [9, 20, 28, 30]):
    m = (t > 86400+3*3600) & (t < 86400+3*3600+600); a.step(t[m]-t[m][0], D[m, j], lw=0.6, color='k', where='post'); a.set_ylabel(f'ch{j+1}'); a.set_ylim(-0.5, 16.5)
axs[-1].set_xlabel('s (day 2, ZT3 + 10 min)'); plt.tight_layout(); plt.savefig('explore/fig_traces_10min.png', dpi=110)
np.save('explore/_act.npy', act)
