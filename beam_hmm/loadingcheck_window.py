"""LoadingCheck first-screen window: 1 h vs 2 h (the recommended default) on the pilot, weeks 1-9, ch17-32.
Descriptive validation (thresholds unchanged, not tuned). References, as in the 2-h validation (ledger 2026-11-14):
  - full-week hub declared dead (spike filter + P1-C criteria) -> was BLOCKED BEAM flagged in the window?
  - the 14 fly-weeks with a live fly hidden behind a constant 16 -> was FROZEN AT END BEAM 16 flagged?
  - flags on READABLE fly-weeks (P1-F, results/a5/a5.json) with no full-week hub and not hidden: candidate false
    alarms (e.g. a resting fly read as FROZEN / LOW ACTIVITY). Not proof of a false alarm: a blockage can be
    transient and gone by the full-week view.
Usage: python -I loadingcheck_window.py <out.json>"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pipeline as pl, labconv as lc

HIDDEN = {(1, 30)} | {(2, c) for c in (19, 27, 29, 30, 31, 32)} | {(w, 27) for w in range(3, 9)} | {(5, 32)}
here = os.path.dirname(os.path.abspath(__file__))


def main(out):
    a5 = json.load(open(os.path.join(here, 'results/a5/a5.json')))
    rows = []
    for w in range(1, 10):
        X, T = pl.load_week(os.path.join(here, f'data/week{w}/week{w}.npz'))
        for j in range(X.shape[1]):
            ch = 17 + j; x = X[:, j]
            xs, _ = pl.spike_filter(x); _, info = pl.hub_repair(xs)
            st = a5['flies'].get(str(ch), {}).get('status', [None] * 9)[w - 1]
            r = dict(week=w, ch=ch, week_hub_dead=bool(info.get('replaced', 0) > 0), hidden=(w, ch) in HIDDEN, status=st)
            for h in (1.0, 2.0):
                r[f'flags_{int(h)}h'] = lc.loading_check(x, T, hours=h)['flags']
            rows.append(r)
    summ = {}
    for h in ('1h', '2h'):
        fl = lambda r: r[f'flags_{h}']
        blk = lambda r: any(f.startswith('BLOCKED') for f in fl(r))
        n_both = sum(r['week_hub_dead'] and blk(r) for r in rows)
        n_win_only = sum((not r['week_hub_dead']) and blk(r) for r in rows)
        n_week_only = sum(r['week_hub_dead'] and not blk(r) for r in rows)
        n_week_only_anyflag = sum(r['week_hub_dead'] and not blk(r) and bool(fl(r)) for r in rows)
        n_neither = sum((not r['week_hub_dead']) and not blk(r) for r in rows)
        hid = sum(r['hidden'] and 'FROZEN AT END BEAM 16' in fl(r) for r in rows)
        clean = [r for r in rows if r['status'] == 'READABLE' and not r['week_hub_dead'] and not r['hidden']]
        cand = [r for r in clean if fl(r)]
        kinds = {}
        for r in cand:
            for f in fl(r):
                k = f.split(' (')[0]; k = 'BLOCKED BEAM' if k.startswith('BLOCKED') else k
                k = 'FROZEN' if k.startswith('FROZEN') else k
                kinds[k] = kinds.get(k, 0) + 1
        summ[h] = dict(blocked_both=n_both, blocked_window_only=n_win_only, blocked_week_only=n_week_only,
                       week_only_but_other_flag=n_week_only_anyflag, neither=n_neither, hidden_caught=hid,
                       hidden_total=len(HIDDEN), clean_readable_fly_weeks=len(clean),
                       clean_flagged=len(cand), clean_flag_kinds=kinds)
        print(h, summ[h], flush=True)
    json.dump(dict(summary=summ, rows=rows), open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
