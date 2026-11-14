"""T1-C1 real-data run (PRESPEC T1-C + T1-C1): weeks 2-3 pooled per fly, primary set, E = 2 h (decision), 1 and 3 h
(sensitivity), C1 vs C0 (secondary); week 1 descriptive. Usage: python -I t1c_real.py <outdir>"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl, p1e, t1c

WEEKS = list(range(1, 10))


def build(aw, F):
    B, used, deaths = [], [], {}
    for j, c in enumerate(c + 1 for c in pl.CH):
        chg = [p1e.change_times(F[w]['P'][:, j], F[w]['t']) for w in WEEKS]
        death = p1e.death_times(chg, [F[w]['t'][-1] for w in WEEKS]); deaths[c] = death
        acc = None; ok_weeks = []
        for fi, w in enumerate(WEEKS):
            if w not in aw or not F[w]['info'][j]['included']:
                continue
            x = F[w]['P'][:, j]; t = F[w]['t']; v = p1e.valid_mask(t, w, death, fi)
            if v.sum() == 0 or (x[v].max() - x[v].min()) < 8:
                continue
            b = t1c.fly_bins(t1c.labels_from_activity(pl.activity_1s(x, t, v)), F[w]['stamp']['lights_on_h'])
            acc = b if acc is None else {k: acc[k] + b[k] for k in acc}; ok_weeks.append(w)
        if acc is not None and (acc['es'].sum() + acc['ew'].sum()) >= 3 * 86400:
            B.append(acc); used.append(dict(ch=c, weeks=ok_weeks, death=death))
    return B, used


if __name__ == '__main__':
    outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
    F = {}
    for w in WEEKS:
        X, t = pl.load_week(f'data/week{w}/week{w}.npz'); P, info = p1e.preprocess(X, t)
        F[w] = dict(P=P, t=t, info=info, stamp=p1e.phase_stamp(P, t, info) if w in (1, 2, 3) else None)
    out = {}
    for label, aw in (('weeks2-3_TEST', [2, 3]), ('week1_descriptive', [1])):
        B, used = build(aw, F)
        assert all(u['death'] is None or u['death']['file'] > max(aw) - 1 for u in used), 'a used fly dies inside the window'
        res = dict(flies=used, n=len(B))
        for E in (2.0, 1.0, 3.0):
            r = t1c.rotation_test(B, np.random.SeedSequence([11001, 1]), excl_h=E)
            r.update(peaks_h=t1c.within_phase_peak(r['h']), peaks_g=t1c.within_phase_peak(r['g'])); res[f'E{E:g}'] = r
        r1 = t1c.rotation_test(B, np.random.SeedSequence([11001, 1]), excl_h=2.0, model_id=1); res['C1_vs_C0_E2'] = r1
        out[label] = res
        print(f"== {label}: {len(B)} flies {[u['ch'] for u in used]}")
        for k in ('E2', 'E1', 'E3', 'C1_vs_C0_E2'):
            r = res[k]; print(f"   {k}: LR {r['LR']:.1f}  rotation p {r['p']:.3f}  (null median LR {r['null_median']:.1f})")
        for nm, pk in (('wake-ups h', res['E2']['peaks_h']), ('onsets g', res['E2']['peaks_g'])):
            print(f"   {nm}: light-phase peak ZT {pk['light_peak_zt']:.1f}, dark-phase peak ZT {pk['dark_peak_zt']:.1f}, dark range {pk['dark_range']:.2f}")
    json.dump(out, open(os.path.join(outdir, 't1c1_real.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
