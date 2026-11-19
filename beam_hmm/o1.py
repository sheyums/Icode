"""PRESPEC O1 (2026-11-19). O1-A: extract fully observed sleep and wake bouts (weeks 2-3, primary set) for the MATLAB
hazard / model fits (o1a_run.m). O1-B: T1-C1 rotation test per week (descriptive).
Usage: python -I o1.py <outdir> [A|B|AB]"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl, p1e, t1c

here = os.path.dirname(os.path.abspath(__file__))
WEEKS = list(range(1, 10))
STAMP = {int(w): s['lights_on_h'] for w, s in json.load(open(os.path.join(here, 'results/a5/a5.json')))['stamps'].items()}


def load_all():
    F = {}
    for w in WEEKS:
        X, t = pl.load_week(os.path.join(here, f'data/week{w}/week{w}.npz')); P, info = p1e.preprocess(X, t)
        F[w] = dict(P=P, t=t, info=info)
    deaths = {}
    for j, c in enumerate(c + 1 for c in pl.CH):
        chg = [p1e.change_times(F[w]['P'][:, j], F[w]['t']) for w in WEEKS]
        deaths[c] = p1e.death_times(chg, [F[w]['t'][-1] for w in WEEKS])
    return F, deaths


def fly_weeks(F, deaths, weeks):
    """Yield (ch, week, activity a) for included, non-narrow fly-weeks with valid data (t1c_real.build rules)."""
    for j, c in enumerate(c + 1 for c in pl.CH):
        for w in weeks:
            if not F[w]['info'][j]['included']:
                continue
            x = F[w]['P'][:, j]; t = F[w]['t']; v = p1e.valid_mask(t, w, deaths[c], WEEKS.index(w))
            if v.sum() == 0 or (x[v].max() - x[v].min()) < 8:
                continue
            yield c, w, pl.activity_1s(x, t, v)


def bouts_of(a, lights_on_h):
    """Fully observed sleep and wake bouts. Returns rows (kind, start_s, length_s, start_zt) and drop counts."""
    sb = pl.sleep_bouts(a); rows = []; drop = dict(sleep=0, wake=0)
    zt = lambda s: float(((s / 3600.0 - lights_on_h) % 24))
    for s, n, cens in sb:
        if cens: drop['sleep'] += 1
        else: rows.append(('sleep', s, n, zt(s)))
    for (s1, n1, c1), (s2, n2, c2) in zip(sb[:-1], sb[1:]):
        w0, w1 = s1 + n1, s2
        if c1 or c2 or (a[w0:w1] == -1).any():
            drop['wake'] += 1
        else:
            rows.append(('wake', w0, w1 - w0, zt(w0)))
    return rows, drop


def run_A(F, deaths, outdir):
    per = {}; rows = []; drops = dict(sleep=0, wake=0)
    acc = {}
    for c, w, a in fly_weeks(F, deaths, [2, 3]):
        acc.setdefault(c, []).append((w, a))
    for c, lst in acc.items():
        if sum((a >= 0).sum() for _, a in lst) < 3 * 86400:
            continue
        for w, a in lst:
            r, d = bouts_of(a, STAMP[w]); drops['sleep'] += d['sleep']; drops['wake'] += d['wake']
            rows += [(c, w) + x for x in r]
        per[c] = [w for w, _ in lst]
    with open(os.path.join(outdir, 'bouts_weeks2-3.csv'), 'w') as f:
        f.write('ch,week,kind,start_s,length_s,start_zt\n')
        for c, w, k, s, n, z in rows:
            f.write(f'{c},{w},{k},{s},{n},{z:.4f}\n')
    summ = dict(flies=sorted(per), fly_weeks=per, dropped=drops)
    for k in ('sleep', 'wake'):
        L = np.array([n for _, _, kk, _, n, _ in rows if kk == k]); summ[k] = dict(n=int(len(L)), dropped=drops[k],
            dropped_share=float(drops[k] / (len(L) + drops[k])), min=int(L.min()), max=int(L.max()))
    json.dump(summ, open(os.path.join(outdir, 'bouts_summary.json'), 'w'), indent=1)
    print('O1-A extraction:', json.dumps({k: summ[k] for k in ('flies', 'sleep', 'wake')}))


def run_B(F, deaths, outdir):
    out = {}
    for w in WEEKS:
        B, used = [], []
        for c, ww, a in fly_weeks(F, deaths, [w]):
            if (a >= 0).sum() < 3 * 86400:
                continue
            B.append(t1c.fly_bins(t1c.labels_from_activity(a), STAMP[w])); used.append(c)
        res = dict(flies=used, n=len(B))
        if len(B) >= 4:
            r = t1c.rotation_test(B, np.random.SeedSequence([11301, w]), excl_h=2.0)
            r.update(peaks_h=t1c.within_phase_peak(r['h']), peaks_g=t1c.within_phase_peak(r['g']),
                     LR_over_null_median=r['LR'] / r['null_median'])
            res['E2'] = r
            print(f"week {w}: {len(B)} flies {used}  LR {r['LR']:.1f}  p {r['p']:.3f}  null median {r['null_median']:.1f}  "
                  f"ratio {r['LR'] / r['null_median']:.1f} | wake-ups peak L {r['peaks_h']['light_peak_zt']:.1f} "
                  f"D {r['peaks_h']['dark_peak_zt']:.1f} range {r['peaks_h']['dark_range']:.2f} | onsets peak L "
                  f"{r['peaks_g']['light_peak_zt']:.1f} D {r['peaks_g']['dark_peak_zt']:.1f} range {r['peaks_g']['dark_range']:.2f}",
                  flush=True)
        else:
            print(f'week {w}: {len(B)} qualifying flies -- not computed', flush=True)
        out[w] = res
    json.dump(out, open(os.path.join(outdir, 'o1b_weeks.json'), 'w'), indent=1,
              default=lambda o: o.item() if hasattr(o, 'item') else str(o))


if __name__ == '__main__':
    outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True); which = sys.argv[2] if len(sys.argv) > 2 else 'AB'
    F, deaths = load_all()
    if 'A' in which: run_A(F, deaths, outdir)
    if 'B' in which: run_B(F, deaths, outdir)
