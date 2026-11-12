"""Run PRESPEC P1 on one week file. Usage: python3 -I run_p1.py <week.npz> <week label> <outdir>"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl

npz, label, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)
X, t = pl.load_week(npz)
nch = X.shape[1]; names = [f'ch{c+1}' for c in pl.CH]

# preprocessing
P = np.empty_like(X); nfix = []
for j in range(nch):
    P[:, j], k = pl.spike_filter(X[:, j]); nfix.append(k)

# inclusion
first24 = t < 86400
nchg24 = (np.diff(P[first24], axis=0) != 0).sum(0)
empty = nchg24 < 100
bigfrac = []
for j in range(nch):
    d = np.diff(P[:, j]); c = d != 0
    bigfrac.append(float((np.abs(d[c]) >= 4).mean()) if c.any() else np.nan)
bigfrac = np.array(bigfrac); artifact = (bigfrac > 0.05) & ~empty
C = pl.change_counts_10min(P, t)
hw = pl.handling_start(C, ~empty)
t_handle = hw * 600.0 if hw is not None else np.inf
valid = (t >= 7200) & (t < t_handle)

rows = []
for j in range(nch):
    sil = pl.longest_silences(P[:, j], t, valid)
    long12 = [(s / 3600, L / 3600) for s, L in sil if L >= 43200]
    status = 'empty/dead at load' if empty[j] else ('artifact-heavy' if artifact[j] else 'included')
    rows.append(dict(ch=names[j], changes_first24h=int(nchg24[j]), spikes_fixed=nfix[j],
                     frac_jumps_ge4=round(float(bigfrac[j]), 4), silences_ge12h_h=long12, status=status))

summary = dict(week=label, t_handling_start_h=None if hw is None else t_handle / 3600,
               valid_window_h=[2.0, (t_handle if np.isfinite(t_handle) else t[-1]) / 3600], channels=rows)

# scoring + descriptors for non-empty channels
per = {}
for j in range(nch):
    if empty[j]: continue
    x = P[:, j]
    a = pl.activity_1s(x, t, valid)
    y, refs, bouts = pl.sqw_minutes(a)
    D = pl.minute_descriptors(x, t, valid)
    n = min(len(y), len(D['changes'])); y = y[:n]
    light = ((np.arange(n) * 60) % pl.DAY) < pl.HALF
    res = dict(refs=refs, n_bouts=len(bouts), n_bouts_uncensored=sum(1 for b in bouts if not b[2]),
               sleep_min_per_day=float((y == pl.S).sum() / ((y >= 0).sum() / 1440)) if (y >= 0).any() else None)
    for sym, nm in ((pl.S, 'S'), (pl.Q, 'Q'), (pl.W, 'W')):
        for ph, pn in ((True, 'L'), (False, 'D')):
            sel = (y == sym) & (light == ph)
            if sel.sum() < 10: continue
            r = {k: float(np.nanmedian(v[:n][sel])) for k, v in D.items()}
            r['mean_changes'] = float(np.nanmean(D['changes'][:n][sel]))
            act = np.nansum(D['act_s'][:n][sel]); fl = np.nansum(D['flickonly_s'][:n][sel])
            r['share_active_s_flicker_only'] = float(fl / act) if act > 0 else None
            r['n_min'] = int(sel.sum())
            res[f'{nm}_{pn}'] = r
    per[names[j]] = res
summary['per_fly'] = per
with open(os.path.join(outdir, f'p1_{label}.json'), 'w') as f:
    json.dump(summary, f, indent=1, default=float)
print(json.dumps({k: summary[k] for k in ('week', 't_handling_start_h', 'valid_window_h')}))
for r in rows: print(r)
