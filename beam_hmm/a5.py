"""A5 — ageing and the approach to death (PRESPEC A5 + P1-F). Usage: python -I a5.py <outdir>"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.stats import wilcoxon
import pipeline as pl, p1e, a1

DT = 0.4; WEEKS = list(range(1, 10)); TREND = [2, 3, 4, 5]; STAMP_OVERRIDE = {6: 23.83}
SHEET = {'died_in_after_wk5': [25], 'alive_in_after_wk6': [19, 20, 23, 24, 27, 30, 31, 32], 'alive_in_after_wk8': [19]}
rng = np.random.default_rng(5)


def status(x):
    n = int((np.diff(x) != 0).sum())
    return 'DEAD-LIKE' if n < 100 else None


def day_components(x, t, v, dark):
    """Per-day sums from which M1-M4 are rebuilt (for the day bootstrap)."""
    sleep, pause = a1.stillness_masks(x, t, v)
    chg = np.r_[False, (x[1:] != x[:-1]) & v[1:] & v[:-1]]
    end = a1.END(x); day = (t // 86400).astype(int); rows = []
    # walking bouts (A4 definition)
    cidx = np.nonzero(chg)[0]; speeds_by_day = {}
    if len(cidx):
        tc = t[cidx]; step = np.abs(x[cidx].astype(int) - x[cidx - 1]); brk = np.r_[True, np.diff(tc) > 2.0]
        bid = np.cumsum(brk) - 1; path = np.bincount(bid, weights=step)
        first = tc[brk]; last = np.r_[tc[np.nonzero(brk)[0][1:] - 1], tc[-1]]
        for p_, f_, l_ in zip(path, first, last):
            if p_ >= 3:
                speeds_by_day.setdefault(int(f_ // 86400), []).append(p_ / (l_ - f_ + DT))
    for d in np.unique(day[v]):
        m = (day == d) & v; ml = m & ~dark; mdk = m & dark
        rows.append(dict(day=int(d), valid_s=float(m.sum() * DT), sleep_s=float((sleep & m).sum() * DT),
                         sleep_light_s=float((sleep & ml).sum() * DT), sleep_dark_s=float((sleep & mdk).sum() * DT),
                         wake_s=float((~sleep & m).sum() * DT), changes_wake=int((chg & ~sleep & m).sum()),
                         sL=int((sleep & ml).sum()), sL_end=int((sleep & ml & end).sum()),
                         pL=int((pause & ml).sum()), pL_end=int((pause & ml & end).sum()), speeds=speeds_by_day.get(int(d), [])))
    return rows


def measures(rows):
    vs = sum(r['valid_s'] for r in rows)
    if vs < 3 * 86400:
        return None
    m = dict(valid_days=vs / 86400, M1=sum(r['sleep_s'] for r in rows) / 60 / (vs / 86400))
    ws = sum(r['wake_s'] for r in rows)
    m['M2'] = sum(r['changes_wake'] for r in rows) / (ws / 3600) if ws > 0 else None
    sL, pL = sum(r['sL'] for r in rows), sum(r['pL'] for r in rows)
    m['M3'] = (sum(r['sL_end'] for r in rows) / sL - sum(r['pL_end'] for r in rows) / pL) if (sL * DT >= 1800 and pL * DT >= 1800) else None
    sp = [s for r in rows for s in r['speeds']]
    m['M4'] = float(np.median(sp)) if len(sp) >= 50 else None
    sd = sum(r['sleep_dark_s'] for r in rows); sl = sum(r['sleep_light_s'] for r in rows)
    m['LD_sleep_ratio'] = sl / sd if sd > 0 else None
    return m


def slope(weeks, vals):
    w = np.array(weeks, float); y = np.array(vals, float)
    return float(np.polyfit(w, y, 1)[0])


if __name__ == '__main__':
    outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
    F = {}
    for w in WEEKS:
        X, t = pl.load_week(f'data/week{w}/week{w}.npz'); P, info = p1e.preprocess(X, t)
        st = p1e.phase_stamp(P, t, info)
        if w in STAMP_OVERRIDE:
            st = dict(st, lights_on_h=STAMP_OVERRIDE[w], flagged=False, override='hand check (PRESPEC A5)')
        assert not st['flagged'], f'week {w} phase flagged'
        F[w] = dict(P=P, t=t, info=info, stamp=st)
    out = dict(stamps={w: F[w]['stamp'] for w in WEEKS}, flies={})
    for j, c in enumerate(range(17, 33)):
        stat = []
        for w in WEEKS:
            x = F[w]['P'][:, j]; t = F[w]['t']; hub = F[w]['info'][j]['hub']
            n = int((np.diff(x) != 0).sum())
            if n < 100:
                s = 'DEAD-LIKE'
            elif hub['replaced'] / len(x) > 0.5 or (x.max() - x.min()) < 8:
                s = 'UNREADABLE'
            else:
                s = 'READABLE'
            stat.append(s)
        readable = [i for i, s in enumerate(stat) if s == 'READABLE']
        death = None; undetermined = None
        if readable:
            r = readable[-1]; w = WEEKS[r]; later = stat[r + 1:]
            x = F[w]['P'][:, j]; t = F[w]['t']; tc = t[1:][x[1:] != x[:-1]]
            edges = np.r_[t[0], tc, t[-1]]; gaps = np.diff(edges); big = np.nonzero(gaps >= 43200)[0]
            cut = float(edges[big[-1]]) if len(big) else float(t[-1])
            if r < len(WEEKS) - 1 and all(s == 'DEAD-LIKE' for s in later):
                death = dict(week=w, t_h=cut / 3600)
            elif any(s == 'UNREADABLE' for s in later):
                undetermined = dict(week=w, t_h=cut / 3600)
        weekly = {}; comps = {}
        for i, w in enumerate(WEEKS):
            if stat[i] != 'READABLE' or not F[w]['info'][j]['included']:
                continue
            x = F[w]['P'][:, j]; t = F[w]['t']
            v = t >= 86400
            if w == 1:
                v &= t < 150 * 3600
            stop = death or undetermined
            if stop is not None:
                if w > stop['week']:
                    continue
                if w == stop['week']:
                    v &= t < stop['t_h'] * 3600 - (86400 if death else 0)
            if v.sum() == 0 or (x[v].max() - x[v].min()) < 8:
                continue
            dark = ((t / 3600 - F[w]['stamp']['lights_on_h']) % 24) >= 12
            rows = day_components(x, t, v, dark); m = measures(rows)
            if m:
                weekly[w] = m; comps[w] = rows
        rec = dict(status=stat, death=death, undetermined=undetermined, weekly=weekly)
        # trend over weeks 2-5
        tw = [w for w in TREND if w in weekly]
        rec['slopes'] = {}; rec['slope_boot95'] = {}
        for k in ('M1', 'M2', 'M3', 'M4'):
            ww = [w for w in tw if weekly[w].get(k) is not None]
            if len(ww) >= 3:
                rec['slopes'][k] = slope(ww, [weekly[w][k] for w in ww])
                bs = []
                for _ in range(1000):
                    vals = []
                    for w in ww:
                        rr = comps[w]; pick = [rr[i] for i in rng.integers(0, len(rr), len(rr))]
                        mm = measures(pick); vals.append(mm.get(k) if mm else None)
                    if all(v_ is not None for v_ in vals):
                        bs.append(slope(ww, vals))
                rec['slope_boot95'][k] = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))] if len(bs) > 50 else None
        # death-aligned daily M1/M2 (descriptive)
        if death:
            daily = []
            for w in sorted(comps):
                base = sum(F[ww]['t'][-1] for ww in WEEKS if ww < w)
                for r_ in comps[w]:
                    daily.append(dict(abs_day=(base + r_['day'] * 86400) / 86400, sleep_min=r_['sleep_s'] / 60 * 86400 / r_['valid_s'],
                                      changes_per_wake_h=r_['changes_wake'] / (r_['wake_s'] / 3600) if r_['wake_s'] else None))
            dabs = sum(F[ww]['t'][-1] for ww in WEEKS if ww < death['week']) / 86400 + death['t_h'] / 24
            for d_ in daily:
                d_['days_before_death'] = dabs - d_['abs_day']
            rec['death_aligned'] = [d_ for d_ in daily if d_['days_before_death'] <= 21]
        out['flies'][c] = rec
    # tests
    tests = {}
    for k in ('M1', 'M2', 'M3', 'M4'):
        S = np.round([out['flies'][c]['slopes'][k] for c in out['flies'] if k in out['flies'][c]['slopes']], 10)
        if len(S) >= 5:
            ties = len(np.unique(np.abs(S))) < len(S) or (S == 0).any()
            p = float(wilcoxon(S, method='approx' if ties else 'exact').pvalue)
            tests[k] = dict(n=int(len(S)), median_slope=float(np.median(S)), n_positive=int((S > 0).sum()), sd=float(np.std(S, ddof=1)), p=p,
                            method='approx' if ties else 'exact')
    order = sorted(tests, key=lambda k: tests[k]['p']); prev = 0
    for i, k in enumerate(order):
        prev = min(1.0, max(prev, (len(order) - i) * tests[k]['p'])); tests[k]['p_holm'] = prev
    out['tests'] = tests
    # survival sheet comparison
    cmp = []
    for c, r in out['flies'].items():
        dw = r['death']['week'] if r['death'] else None
        if c in SHEET['alive_in_after_wk6'] and dw is not None and dw < 6:
            cmp.append(f'ch{c}: sheet alive in/after wk6, rule says dead in wk{dw}')
        if c in SHEET['alive_in_after_wk8'] and dw is not None and dw < 8:
            cmp.append(f'ch{c}: sheet alive in/after wk8, rule says dead in wk{dw}')
        if c in SHEET['died_in_after_wk5'] and (dw is None or dw < 5):
            cmp.append(f'ch{c}: sheet died in/after wk5, rule says {dw}')
    out['sheet_disagreements'] = cmp
    json.dump(out, open(os.path.join(outdir, 'a5.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    for c, r in out['flies'].items():
        print(f"ch{c}: {' '.join(s[:4] for s in r['status'])} | death {r['death']} undetermined {r['undetermined']} | weeks {sorted(r['weekly'])}")
    print('TESTS', json.dumps(tests, indent=1)); print('SHEET DISAGREEMENTS', cmp)
