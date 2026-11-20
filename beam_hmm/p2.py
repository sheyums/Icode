"""PRESPEC P2 (2026-11-20): sleep-bout duration by place at sleep, END vs MIDDLE (orientation-free).
Usage: python -I p2.py <outdir> check|run"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl, p1e, o1, q1
from scipy.stats import kstest

END = {0, 1, 2, 14, 15, 16}


def bout_table(F, deaths):
    """Rows (ch, kind, start_zt, T, place) for the O1-A bouts; place = modal reading inside the bout (sleep) or of the
    preceding sleep bout (wake)."""
    rows = []
    for j, c in enumerate(c + 1 for c in pl.CH):
        lst = []
        for w in (2, 3):
            if not F[w]['info'][j]['included']:
                continue
            x = F[w]['P'][:, j]; t = F[w]['t']; v = p1e.valid_mask(t, w, deaths[c], o1.WEEKS.index(w))
            if v.sum() == 0 or (x[v].max() - x[v].min()) < 8:
                continue
            lst.append((w, x, t, pl.activity_1s(x, t, v)))
        if not lst or sum((a >= 0).sum() for *_, a in lst) < 3 * 86400:
            continue
        for w, x, t, a in lst:
            br, _ = o1.bouts_of(a, o1.STAMP[w]); place_end = {}
            for kind, s, n, zt in br:
                if kind != 'sleep':
                    continue
                i0, i1 = np.searchsorted(t, s), np.searchsorted(t, s + n)
                pl_ = int(np.bincount(x[i0:i1].astype(int), minlength=17).argmax())
                place_end[s + n] = pl_; rows.append((c, 'sleep', zt, n, pl_))
            for kind, s, n, zt in br:
                if kind == 'wake' and s in place_end:
                    rows.append((c, 'wake', zt, n, place_end[s]))
    return rows


def analyse(rows, kind, phase, nperm, seed, labels=None):
    R = [r for r in rows if r[1] == kind and (phase is None or (r[2] < 12) == (phase == 'light'))]
    fly = np.array([r[0] for r in R]); zt = np.array([r[2] for r in R]); T = np.array([r[3] for r in R], float)
    place = np.array([r[4] for r in R])
    g = np.isin(place, list(END)).astype(int) if labels is None else labels
    strata = fly * 10 + np.floor(zt / 4).astype(int)
    Z, Zp = q1.logrank_perm(T, g, strata, nperm, np.random.default_rng(np.random.SeedSequence(seed)))
    p = (1 + int((np.abs(Zp) >= abs(Z)).sum())) / (nperm + 1)
    per = {}
    for c in np.unique(fly):
        e = T[(fly == c) & (g == 1)]; m = T[(fly == c) & (g == 0)]
        per[int(c)] = dict(n_end=int(len(e)), n_mid=int(len(m)), share_end=float(len(e) / (len(e) + len(m))),
                           med_end=float(np.median(e)) if len(e) else None, med_mid=float(np.median(m)) if len(m) else None)
    cmp_ = {c: v for c, v in per.items() if v['n_end'] >= 5 and v['n_mid'] >= 5}
    return dict(kind=kind, phase=phase or 'both', n=int(len(T)), n_end=int(g.sum()), n_mid=int((g == 0).sum()),
                n_beyond1=int((place == 0).sum()), median_end=float(np.median(T[g == 1])), median_mid=float(np.median(T[g == 0])),
                Z=Z, p=p, nperm=nperm, null_absZ_95=float(np.percentile(np.abs(Zp), 95)), per_fly=per,
                flies_end_longer=int(sum(v['med_end'] > v['med_mid'] for v in cmp_.values())), flies_compared=len(cmp_))


def main(outdir, mode):
    os.makedirs(outdir, exist_ok=True); F, deaths = o1.load_all(); rows = bout_table(F, deaths)
    if mode == 'check':
        sl = [r for r in rows if r[1] == 'sleep']; share = float(np.mean([r[4] in END for r in sl]))
        rng = np.random.default_rng(np.random.SeedSequence([13005])); ps = []
        for i in range(200):
            lab = (rng.random(len(sl)) < share).astype(int)
            ps.append(analyse(rows, 'sleep', None, 999, [13006, i], labels=lab)['p'])
        out = dict(draws=200, nperm=999, end_share=share, ks_p=float(kstest(ps, 'uniform').pvalue),
                   share_p_below_005=float(np.mean(np.array(ps) < 0.05)))
        print('CHECK (random labels):', out); json.dump(dict(out, p=ps), open(os.path.join(outdir, 'p2_check.json'), 'w'))
        return
    res = dict(primary=analyse(rows, 'sleep', None, 9999, [13001]),
               sleep_light=analyse(rows, 'sleep', 'light', 9999, [13002]),
               sleep_dark=analyse(rows, 'sleep', 'dark', 9999, [13003]),
               wake_by_preceding_place=analyse(rows, 'wake', None, 9999, [13004]))
    for k, r in res.items():
        print(f"{k}: n {r['n']} (END {r['n_end']}, MIDDLE {r['n_mid']}, at 0 {r['n_beyond1']}) median END {r['median_end']:.0f} s "
              f"MIDDLE {r['median_mid']:.0f} s | Z {r['Z']:.2f} p {r['p']:.4f} (null |Z| 95% {r['null_absZ_95']:.2f}) | "
              f"flies END longer {r['flies_end_longer']}/{r['flies_compared']}")
    json.dump(res, open(os.path.join(outdir, 'p2_results.json'), 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
