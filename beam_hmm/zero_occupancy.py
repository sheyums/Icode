"""DESCRIPTIVE (post hoc, 2026-11-16): position-0 occupancy per fly-week, weeks 1-9, ch17-32, to see whether tube
orientation (food next to beam 1 vs beam 16; user: mixed in the pilot) leaves a mark that does NOT use where flies
sleep. Per fly-week, after spike filter + P1-C hub repair, on P1-F READABLE fly-weeks:
  share of samples at 0; entries into 0 per day; median visit length at 0 (s); the reading just before each entry
  into 0 (1 expected; anything else is reported); end-beam hub (1/2 or 15/16) if any.
Nothing here is a test, and nothing here may be used to assign a food end for the A1 food-end split.
Usage: python -I zero_occupancy.py <out.json>"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pipeline as pl

here = os.path.dirname(os.path.abspath(__file__))


def main(out):
    a5 = json.load(open(os.path.join(here, 'results/a5/a5.json')))
    rows = []
    for w in range(1, 10):
        X, T = pl.load_week(os.path.join(here, f'data/week{w}/week{w}.npz'))
        days = (T[-1] - T[0]) / 86400
        for j in range(X.shape[1]):
            ch = 17 + j
            st = a5['flies'].get(str(ch), {}).get('status', [None] * 9)[w - 1]
            xs, _ = pl.spike_filter(X[:, j]); xr, info = pl.hub_repair(xs)
            x = xr.astype(int); at0 = x == 0
            ent = np.nonzero(at0[1:] & ~at0[:-1])[0] + 1
            prev = np.bincount(x[ent - 1], minlength=17).tolist() if len(ent) else [0] * 17
            ex = np.nonzero(~at0[1:] & at0[:-1])[0] + 1
            lens = []
            for e in ent:
                k = ex[ex > e]
                if len(k): lens.append((k[0] - e) * 0.4)
            hub = info.get('hub') if info.get('replaced', 0) > 0 else None
            rows.append(dict(week=w, ch=ch, status=st, share0=float(at0.mean()), entries_per_day=float(len(ent) / days),
                             median_visit_s=float(np.median(lens)) if lens else None, prev_before_0=prev,
                             share16=float((x == 16).mean()), share1=float((x == 1).mean()), dead_hub=hub))
    R = [r for r in rows if r['status'] == 'READABLE']
    s0 = np.array([r['share0'] for r in R]); e0 = np.array([r['entries_per_day'] for r in R])
    print(f'READABLE fly-weeks {len(R)}; share at 0: quantiles 10/25/50/75/90 % =',
          np.round(np.percentile(s0, [10, 25, 50, 75, 90]) * 100, 2))
    print('entries into 0 per day: quantiles =', np.round(np.percentile(e0, [10, 25, 50, 75, 90]), 1))
    print('fly-weeks with < 1 entry/day:', int((e0 < 1).sum()), ' with 0 entries:', int((e0 == 0).sum()))
    prev = np.sum([r['prev_before_0'] for r in rows], 0)
    print('reading before entering 0, all fly-weeks (index = beam):', {i: int(v) for i, v in enumerate(prev) if v})
    for ch in range(17, 33):
        rr = [r for r in R if r['ch'] == ch]
        print(f'ch{ch}: ' + '  '.join(f"w{r['week']} {100*r['share0']:.1f}%/{r['entries_per_day']:.0f}" for r in rr))
    json.dump(rows, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
