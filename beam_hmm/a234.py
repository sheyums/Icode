"""A2 (pauses around sleep), A3 (edge flicker day/night), A4 (walking day/night) — PRESPEC A2-A4 + P1-E.
Usage: python -I a234.py <outdir> <analysis weeks, comma-separated> <all supplied weeks...>
  e.g. tests:        python -I a234.py results/a234 2,3 1 2 3 4 5 6 7 8 9
       descriptive:  python -I a234.py results/a234 1 1 2 3 4 5 6 7 8 9"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.stats import wilcoxon
import pipeline as pl, p1e, a1, d1d

DT = 0.4
rng0 = np.random.default_rng(0)


def per_fly_week(x, t, v, dark):
    """Per-unit records keyed by day, for one fly-week (x repaired; v valid; dark per sample)."""
    day = (t // 86400).astype(int)
    out = dict(A2=[], A3=[], A4=[], sleep=[])
    # ---- A2
    a = pl.activity_1s(x, t, v); bouts = pl.sleep_bouts(a)
    sleep, pause = a1.stillness_masks(x, t, v)
    for (s0, n0, c0), (s1, n1, c1) in zip(bouts[:-1], bouts[1:]):
        w0, w1 = s0 + n0, s1
        if w1 - w0 < 1800:
            continue
        m = (t >= w0) & (t < w1)
        if not v[m].all():
            continue
        def frac(lo, hi):
            mm = (t >= lo) & (t < hi); return float(pause[mm].mean())
        mid = (w0 + w1) / 2
        out['A2'].append(dict(day=int(w0 // 86400), pre_mid=frac(w1 - 600, w1) - frac(mid - 300, mid + 300),
                              post_mid=frac(w0, w0 + 600) - frac(mid - 300, mid + 300), onset_dark=bool(dark[np.searchsorted(t, w1) - 1])))
    # ---- A3
    i = np.nonzero((x[:-2] == x[2:]) & (np.abs(x[1:-1].astype(int) - x[:-2]) == 1) & v[:-2] & v[1:-1] & v[2:])[0] + 1
    flick = np.zeros(len(x), bool); flick[i] = True
    xf = d1d.remove_flicker(x)
    chg = np.r_[True, xf[1:] != xf[:-1]] | ~v
    rid = np.cumsum(chg); rlen = np.bincount(rid)[rid] * DT
    ftr = (rlen >= 300) & v
    mid_zone = (xf > 2) & (xf < 14)
    for d in np.unique(day[v]):
        md = day == d
        for blk, bm in (('light', ~dark), ('dark', dark)):
            m = md & bm & ftr
            out['A3'].append(dict(day=int(d), block=blk, events=int((flick & m).sum()), hours=float(m.sum() * DT / 3600),
                                  events_mid=int((flick & m & mid_zone).sum()), hours_mid=float((m & mid_zone).sum() * DT / 3600)))
    af = pl.activity_1s(xf, t, v); bf = pl.sleep_bouts(af)
    vd = v.sum() * DT / 86400
    out['sleep'].append(dict(valid_days=float(vd), sleep_min_std=float(sum(n for _, n, _ in bouts) / 60),
                             sleep_min_noflick=float(sum(n for _, n, _ in bf) / 60), n_bouts_std=len(bouts), n_bouts_noflick=len(bf)))
    # ---- A4
    ok = v[1:] & v[:-1]; cidx = np.nonzero((x[1:] != x[:-1]) & ok)[0] + 1
    if len(cidx):
        tc = t[cidx]; step = np.abs(x[cidx].astype(int) - x[cidx - 1])
        brk = np.r_[True, np.diff(tc) > 2.0]; bid = np.cumsum(brk) - 1
        path = np.bincount(bid, weights=step); first = tc[brk]; last = np.r_[tc[np.nonzero(brk)[0][1:] - 1], tc[-1]]
        dur = last - first + DT
        dk = dark[cidx[brk]]
        for p_, du, f_, dkk in zip(path, dur, first, dk):
            if p_ >= 3:
                out['A4'].append(dict(day=int(f_ // 86400), block='dark' if dkk else 'light', path=float(p_), speed=float(p_ / du)))
    wake_h = {blk: float(((~sleep) & v & (dark if blk == 'dark' else ~dark)).sum() * DT / 3600) for blk in ('light', 'dark')}
    out['wake_hours'] = wake_h
    return out


def stats(recs, days=None):
    """Per-fly primary statistics from pooled records (optionally restricted to a list of (week, day) keys)."""
    sel = (lambda r: True) if days is None else (lambda r: (r['wk'], r['day']) in days)
    a2 = [r for r in recs['A2'] if sel(r)]; a3 = [r for r in recs['A3'] if sel(r)]; a4 = [r for r in recs['A4'] if sel(r)]
    s = {}
    if len(a2) >= 5:
        s['A2_pre_mid'] = float(np.mean([r['pre_mid'] for r in a2])); s['A2_post_mid'] = float(np.mean([r['post_mid'] for r in a2]))
        s['A2_n_bouts'] = len(a2)
    ev = {b: sum(r['events'] for r in a3 if r['block'] == b) for b in ('light', 'dark')}
    hr = {b: sum(r['hours'] for r in a3 if r['block'] == b) for b in ('light', 'dark')}
    if min(hr.values()) >= 1:
        s['A3_rate_light'] = ev['light'] / hr['light']; s['A3_rate_dark'] = ev['dark'] / hr['dark']
        s['A3_diff'] = s['A3_rate_dark'] - s['A3_rate_light']
        evm = {b: sum(r['events_mid'] for r in a3 if r['block'] == b) for b in ('light', 'dark')}
        hrm = {b: sum(r['hours_mid'] for r in a3 if r['block'] == b) for b in ('light', 'dark')}
        if min(hrm.values()) >= 1:
            s['A3_diff_middle'] = evm['dark'] / hrm['dark'] - evm['light'] / hrm['light']
    nb = {b: [r for r in a4 if r['block'] == b] for b in ('light', 'dark')}
    if min(len(v) for v in nb.values()) >= 50:
        s['A4_speed_diff'] = float(np.median([r['speed'] for r in nb['dark']]) - np.median([r['speed'] for r in nb['light']]))
        s['A4_path_diff'] = float(np.median([r['path'] for r in nb['dark']]) - np.median([r['path'] for r in nb['light']]))
        s['A4_n_bouts'] = {b: len(v) for b, v in nb.items()}
    return s


def boot_all(recs, keys, B=1000):
    """Day-bootstrap 95% intervals for several per-fly statistics at once (resampling (week, day) units)."""
    units = sorted({(r['wk'], r['day']) for part in ('A2', 'A3', 'A4') for r in recs[part]})
    if not units:
        return {}
    by = {p: {} for p in ('A2', 'A3', 'A4')}
    for p in by:
        for r in recs[p]:
            by[p].setdefault((r['wk'], r['day']), []).append(r)
    vals = {k: [] for k in keys}
    for _ in range(B):
        pick = [units[i] for i in rng0.integers(0, len(units), len(units))]
        rr = {p: [r for u in pick for r in by[p].get(u, [])] for p in by}
        st = stats(rr)
        for k in keys:
            if k in st:
                vals[k].append(st[k])
    return {k: ([float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if len(v) > 50 else None) for k, v in vals.items()}


def wtest(vals):
    D = np.round(np.array(vals, float), 10)
    if len(D) < 5:
        return dict(n=int(len(D)), p=None)
    ties = len(np.unique(np.abs(D))) < len(D) or (D == 0).any()
    p = float(wilcoxon(D, alternative='two-sided', method='approx' if ties else 'exact').pvalue)
    return dict(n=int(len(D)), median=float(np.median(D)), n_positive=int((D > 0).sum()), sd=float(D.std(ddof=1)),
                p=p, method='approx (ties/zeros)' if ties else 'exact')


def holm(tests):
    order = sorted([k for k in tests if tests[k].get('p') is not None], key=lambda k: tests[k]['p']); prev = 0
    for i, k in enumerate(order):
        prev = min(1.0, max(prev, (len(order) - i) * tests[k]['p'])); tests[k]['p_holm'] = prev
    return tests


if __name__ == '__main__':
    outdir = sys.argv[1]; aw = [int(w) for w in sys.argv[2].split(',')]; weeks = [int(w) for w in sys.argv[3:]]
    os.makedirs(outdir, exist_ok=True)
    F = {}
    for w in weeks:
        X, t = pl.load_week(f'data/week{w}/week{w}.npz'); P, info = p1e.preprocess(X, t)
        F[w] = dict(P=P, t=t, info=info, stamp=p1e.phase_stamp(P, t, info) if w in aw else None)
    for w in aw:
        assert not F[w]['stamp']['flagged'], f'week {w} failed the phase guard'
    chans = [c + 1 for c in pl.CH]; res = {}
    for j, c in enumerate(chans):
        chg = [p1e.change_times(F[w]['P'][:, j], F[w]['t']) for w in weeks]
        death = p1e.death_times(chg, [F[w]['t'][-1] for w in weeks])
        recs = dict(A2=[], A3=[], A4=[]); sleep_rows = []; used = []; never_art = True; wake_h = {'light': 0, 'dark': 0}
        for fi, w in enumerate(weeks):
            if w not in aw:
                continue
            inf = F[w]['info'][j]; never_art &= inf['never_artifact']
            if not inf['included']:
                continue
            x = F[w]['P'][:, j]; t = F[w]['t']; v = p1e.valid_mask(t, w, death, fi)
            if v.sum() == 0 or (x[v].max() - x[v].min()) < 8:            # narrow-range exclusion (A2-A4 common)
                continue
            dark = ((t / 3600 - F[w]['stamp']['lights_on_h']) % 24) >= 12
            o = per_fly_week(x, t, v, dark)
            for p in ('A2', 'A3', 'A4'):
                for r in o[p]:
                    r['wk'] = w; recs[p].append(r)
            sleep_rows += o['sleep']; used.append(w)
            for b in wake_h: wake_h[b] += o['wake_hours'][b]
        vd = sum(r['valid_days'] for r in sleep_rows)
        st = stats(recs) if vd >= 3 else {}
        r = dict(weeks_used=used, valid_days=vd, never_artifact=bool(never_art), death=death, stats=st, sleep=sleep_rows,
                 secondary=dict(bout_rate_per_wake_h={b: (sum(1 for q in recs['A4'] if q['block'] == b) / wake_h[b]) if wake_h[b] else None for b in wake_h},
                                cross_share=float(np.mean([q['path'] >= 12 for q in recs['A4']])) if recs['A4'] else None))
        if st:
            r['boot95'] = boot_all(recs, [k for k in ('A2_pre_mid', 'A2_post_mid', 'A3_diff', 'A4_speed_diff', 'A4_path_diff') if k in st], B=1000)
        res[c] = r
    prim = [c for c in chans if res[c]['stats']]; sens = [c for c in prim if res[c]['never_artifact']]
    def run(sel):
        A2 = holm({k: wtest([res[c]['stats'][k] for c in sel if k in res[c]['stats']]) for k in ('A2_pre_mid', 'A2_post_mid')})
        A3 = {'A3_diff': wtest([res[c]['stats']['A3_diff'] for c in sel if 'A3_diff' in res[c]['stats']])}
        A4 = holm({k: wtest([res[c]['stats'][k] for c in sel if k in res[c]['stats']]) for k in ('A4_speed_diff', 'A4_path_diff')})
        sec = {'A3_diff_middle': wtest([res[c]['stats']['A3_diff_middle'] for c in sel if 'A3_diff_middle' in res[c]['stats']])}
        return dict(A2=A2, A3=A3, A4=A4, secondary=sec)
    summ = dict(analysis_weeks=aw, supplied_weeks=weeks, stamps={w: F[w]['stamp'] for w in aw}, primary=prim, sensitivity=sens,
                tests_primary=run(prim), tests_sensitivity=run(sens) if len(sens) >= 6 else f'{len(sens)} flies: not tested', flies=res)
    name = f"a234_weeks{'-'.join(map(str, aw))}.json"
    json.dump(summ, open(os.path.join(outdir, name), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print('primary flies', prim, '| sensitivity', sens)
    print(json.dumps(summ['tests_primary'], indent=1, default=float))
