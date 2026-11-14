"""A1 — where flies rest (PRESPEC A1). Model-free.
Per fly x block (light/dark): E_S = share of SLEEP-stillness samples at a tube END, E_P = the same for PAUSE
stillness (30 s - < 5 min without a position change, outside sleep bouts), D = E_S - E_P.
Usage (descriptive, one week file): python -I a1.py <week.npz> <label> <valid_start_h> <valid_end_h> <outdir>
"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl

DT = 0.4
END = lambda x: (x <= 2) | (x >= 14)


def stillness_masks(x, t, valid):
    """Per sample: in a sleep bout (P1), in a pause run (30 s to < 5 min, no change, outside sleep)."""
    a = pl.activity_1s(x, t, valid)
    bouts = pl.sleep_bouts(a)
    sec = np.floor(t).astype(int); insleep_s = np.zeros(len(a), bool)
    for s, n, _ in bouts:
        insleep_s[s:s + n] = True
    sleep = insleep_s[np.minimum(sec, len(a) - 1)] & valid
    chg = np.r_[True, x[1:] != x[:-1]] | ~valid
    run_id = np.cumsum(chg)
    run_len = np.bincount(run_id)[run_id] * DT
    pause = (run_len >= 30) & (run_len < 300) & ~sleep & valid
    return sleep, pause


def measures(x, t, valid, min_minutes=30):
    sleep, pause = stillness_masks(x, t, valid)
    dark = (t % 86400) >= 43200; end = END(x); b1 = x == 0
    out = {}
    for name, blk in (('light', ~dark), ('dark', dark)):
        s = sleep & blk; p = pause & blk
        ms, mp = s.sum() * DT / 60, p.sum() * DT / 60
        r = dict(sleep_min=float(ms), pause_min=float(mp))
        if ms >= min_minutes and mp >= min_minutes:
            r.update(E_S=float(end[s].mean()), E_P=float(end[p].mean()), D=float(end[s].mean() - end[p].mean()),
                     B1_S=float(b1[s].mean()), B1_P=float(b1[p].mean()),
                     low_end_share_of_end_sleep=float((x[s & end] <= 2).mean()) if (s & end).any() else None)
        out[name] = r
    return out, sleep, pause


def day_bootstrap(x, t, sleep, pause, B=1000, seed=0):
    """95% interval of D per block, resampling whole days (fly-level, descriptive)."""
    rng = np.random.default_rng(seed); day = (t // 86400).astype(int); days = np.unique(day)
    dark = (t % 86400) >= 43200; end = END(x); res = {}
    for name, blk in (('light', ~dark), ('dark', dark)):
        s_end = np.bincount(day, weights=(sleep & blk & end), minlength=days.max() + 1)
        s_all = np.bincount(day, weights=(sleep & blk), minlength=days.max() + 1)
        p_end = np.bincount(day, weights=(pause & blk & end), minlength=days.max() + 1)
        p_all = np.bincount(day, weights=(pause & blk), minlength=days.max() + 1)
        D = []
        for _ in range(B):
            d = rng.choice(days, len(days), replace=True)
            if s_all[d].sum() > 0 and p_all[d].sum() > 0:
                D.append(s_end[d].sum() / s_all[d].sum() - p_end[d].sum() / p_all[d].sum())
        res[name] = [float(np.percentile(D, 2.5)), float(np.percentile(D, 97.5))] if D else None
    return res


if __name__ == '__main__':
    npz, label, v0, v1, outdir = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    os.makedirs(outdir, exist_ok=True)
    X, t = pl.load_week(npz)
    valid = (t >= v0 * 3600) & (t < v1 * 3600)
    rows = {}
    for j, c in enumerate(cc + 1 for cc in pl.CH):
        x, _ = pl.spike_filter(X[:, j]); x, info = pl.hub_repair(x)
        d = np.abs(np.diff(x.astype(int))); ch = d != 0
        if (np.diff(x[t < 86400]) != 0).sum() < 100:
            rows[c] = dict(set='empty/dead'); continue
        big = float((d[ch] >= 4).mean()) if ch.any() else 1.0
        m, sl, pa = measures(x, t, valid)
        rows[c] = dict(set='primary' if big <= 0.05 else 'artifact-heavy', frac_big_after_repair=big, hub=info,
                       measures=m, D_day_bootstrap95=day_bootstrap(x, t, sl, pa))
    json.dump(dict(week=label, valid_hours=[v0, v1], flies=rows), open(os.path.join(outdir, f'a1_{label}.json'), 'w'),
              indent=1, default=float)
    for c, r in rows.items():
        if r['set'] != 'primary':
            print(f'ch{c}: {r["set"]}'); continue
        L, Dk = r['measures']['light'], r['measures']['dark']
        f = lambda b: (f"E_S {b['E_S']:.2f} E_P {b['E_P']:.2f} D {b['D']:+.2f}" if 'D' in b else f"insufficient ({b['sleep_min']:.0f}/{b['pause_min']:.0f} min)")
        print(f'ch{c}: light {f(L)} | dark {f(Dk)} | boot95 {r["D_day_bootstrap95"]}')
