"""A1 primary test (PRESPEC A1 + P1-E): weeks 2-3 pooled per fly; death rule over all supplied files.
Usage: python -I a1_test.py <outdir> <week numbers of all supplied files, e.g. 2 3 4 5>"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.stats import wilcoxon
import pipeline as pl, p1e, a1

outdir = sys.argv[1]; weeks = [int(w) for w in sys.argv[2:]]; TEST = [2, 3]
os.makedirs(outdir, exist_ok=True)
F = {}
for w in weeks:
    X, t = pl.load_week(f'data/week{w}/week{w}.npz'); P, info = p1e.preprocess(X, t)
    F[w] = dict(P=P, t=t, info=info, stamp=p1e.phase_stamp(P, t, info))
    print(f"week{w}: lights-on {F[w]['stamp']['lights_on_h']:.2f} h into file, flagged {F[w]['stamp']['flagged']}")
assert not any(F[w]['stamp']['flagged'] for w in TEST), 'a test file failed the phase guard: check by hand'
chans = [c + 1 for c in pl.CH]
res = {}
for j, c in enumerate(chans):
    chg = [p1e.change_times(F[w]['P'][:, j], F[w]['t']) for w in weeks]
    death = p1e.death_times(chg, [F[w]['t'][-1] for w in weeks])
    acc = {b: dict(s=0, s_end=0, p=0, p_end=0, s1=0, p1=0) for b in ('light', 'dark')}
    used = []; never_art = True
    for fi, w in enumerate(weeks):
        if w not in TEST:
            continue
        inf = F[w]['info'][j]; never_art &= inf['never_artifact']
        if not inf['included']:
            continue
        x = F[w]['P'][:, j]; t = F[w]['t']
        v = p1e.valid_mask(t, w, death, fi)
        if v.sum() == 0:
            continue
        used.append(dict(week=w, valid_days=float(v.sum() * 0.4 / 86400)))
        sleep, pause = a1.stillness_masks(x, t, v)
        zt = (t / 3600 - F[w]['stamp']['lights_on_h']) % 24; dark = zt >= 12
        end = a1.END(x); b1z = x == 0
        for b, blk in (('light', ~dark), ('dark', dark)):
            s = sleep & blk; p = pause & blk
            acc[b]['s'] += s.sum(); acc[b]['s_end'] += (s & end).sum(); acc[b]['p'] += p.sum(); acc[b]['p_end'] += (p & end).sum()
            acc[b]['s1'] += (s & b1z).sum(); acc[b]['p1'] += (p & b1z).sum()
    days = sum(u['valid_days'] for u in used)
    r = dict(death=death, fly_weeks_used=used, valid_days=days, never_artifact=bool(never_art), blocks={})
    for b, a in acc.items():
        ms, mp = a['s'] * 0.4 / 60, a['p'] * 0.4 / 60
        d = dict(sleep_min=float(ms), pause_min=float(mp))
        if days >= 3 and ms >= 30 and mp >= 30:
            d.update(E_S=a['s_end'] / a['s'], E_P=a['p_end'] / a['p'], D=a['s_end'] / a['s'] - a['p_end'] / a['p'],
                     B1_S=a['s1'] / a['s'], B1_P=a['p1'] / a['p'])
        r['blocks'][b] = d
    res[c] = r


def test(sel):
    out = {}
    for b in ('light', 'dark'):
        D = np.array([round(res[c]['blocks'][b]['D'], 10) for c in sel if 'D' in res[c]['blocks'][b]])
        ties = len(np.unique(np.abs(D))) < len(D) or (D == 0).any()
        p = float(wilcoxon(D, alternative='two-sided', method='approx' if ties else 'exact').pvalue) if len(D) >= 5 else None
        out[b] = dict(n=int(len(D)), median_D=float(np.median(D)) if len(D) else None, n_positive=int((D > 0).sum()),
                      sd=float(D.std(ddof=1)) if len(D) > 1 else None, p=p, method='approx (ties)' if ties else 'exact')
    ps = [(b, out[b]['p']) for b in out if out[b]['p'] is not None]
    order = sorted(ps, key=lambda z: z[1]); m = len(order); prev = 0
    for k, (b, p) in enumerate(order):
        h = min(1.0, max(prev, (m - k) * p)); out[b]['p_holm'] = h; prev = h
        out[b]['supported'] = bool(h < 0.05 and out[b]['median_D'] > 0)
        out[b]['opposite'] = bool(h < 0.05 and out[b]['median_D'] < 0)
    return out


primary = [c for c in chans if res[c]['fly_weeks_used']]
sensitivity = [c for c in primary if res[c]['never_artifact']]
summary = dict(stamps={w: F[w]['stamp'] for w in weeks}, primary=primary, sensitivity=sensitivity,
               test_primary=test(primary), test_sensitivity=test(sensitivity), flies=res)
json.dump(summary, open(os.path.join(outdir, 'a1_test_weeks2-3.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
for name in ('test_primary', 'test_sensitivity'):
    print(name, json.dumps(summary[name], default=float))
for c in chans:
    r = res[c]
    print(f"ch{c}: days {r['valid_days']:.1f}, death {r['death']}, " + ' | '.join(
        f"{b} D {r['blocks'][b].get('D', float('nan')):+.3f}" for b in ('light', 'dark')))
