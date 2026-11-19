"""PRESPEC Q1 (2026-11-19): Q at the start of a wake bout vs its remaining duration. Landmark L; stratified log-rank
(strata = fly x 4-h start-ZT block); within-stratum permutation null.
Usage: python -I q1.py <outdir> check|run"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl, o1
from scipy.stats import kstest

THETA = 0.5


def zt_of(s, lon):
    return ((np.asarray(s) / 3600.0 - lon) % 24)


def refs_by_zt(a, sb, lon):
    """P1 reference per light phase (0 light, 1 dark), phase from ZT (the one disclosed change)."""
    r = {0: [], 1: []}
    for (s1, n1, _), (s2, _, _) in zip(sb[:-1], sb[1:]):
        w0, w1 = s1 + n1, s2
        if w1 - w0 < 600 or (a[w0:w1] == -1).any():
            continue
        mid = (w0 + w1) // 2; ref = float((a[mid - 150:mid + 150] == 1).mean())
        pw = int(zt_of(w0, lon) >= 12); po = int(zt_of(w1, lon) >= 12)
        for ph in (0, 1):
            if pw == ph or po == ph:
                r[ph].append(ref)
    return {ph: (float(np.median(v)) if len(v) >= 3 else np.nan) for ph, v in r.items()}


def wake_table(F, deaths, L):
    """Rows: fly, start_zt, T, k (number of Q windows among the first L/60), for fully observed wake bouts T >= L."""
    acc = {}
    for c, w, a in o1.fly_weeks(F, deaths, [2, 3]):
        acc.setdefault(c, []).append((w, a))
    rows = []; excl = 0
    for c, lst in acc.items():
        if sum((a >= 0).sum() for _, a in lst) < 3 * 86400:
            continue
        for w, a in lst:
            lon = o1.STAMP[w]; sb = pl.sleep_bouts(a); ref = refs_by_zt(a, sb, lon)
            for (s1, n1, c1), (s2, n2, c2) in zip(sb[:-1], sb[1:]):
                w0, w1 = s1 + n1, s2
                if c1 or c2 or (a[w0:w1] == -1).any() or w1 - w0 < L:
                    continue
                k = 0; bad = False
                for i in range(L // 60):
                    st = w0 + 60 * i; r = ref[int(zt_of(st, lon) >= 12)]
                    if np.isnan(r): bad = True; break
                    k += (a[st:st + 60] == 1).mean() < THETA * r
                if bad: excl += 1; continue
                rows.append((c, float(zt_of(w0, lon)), w1 - w0, int(k)))
    return rows, excl


def logrank_perm(R, g, strata, nperm, rng):
    """Stratified log-rank Z for group g==1 (vs 0) and its within-stratum permutation distribution (all events)."""
    Zs = []; obs = None
    num_o = 0.0; num_p = np.zeros(nperm); V_o = 0.0; V_p = np.zeros(nperm)
    for s in np.unique(strata):
        m = strata == s; r = R[m]; gg = g[m]
        if gg.min() == gg.max():
            continue
        o = np.argsort(r, kind='stable'); r = r[o]; gg = gg[o]; n = len(r)
        G = np.vstack([gg] + [rng.permutation(gg) for _ in range(nperm)]).astype(float)   # row 0 observed
        ut, first, cnt = np.unique(r, return_index=True, return_counts=True)
        nrisk = n - first                                                                  # at risk at each distinct time
        cs = np.concatenate([np.zeros((G.shape[0], 1)), np.cumsum(G, 1)], 1)
        n1 = G.sum(1)[:, None] - cs[:, first]                                              # group-1 at risk
        d1 = cs[:, first + cnt] - cs[:, first]                                             # group-1 events
        E1 = cnt * n1 / nrisk
        V = np.where(nrisk > 1, cnt * (n1 / nrisk) * (1 - n1 / nrisk) * (nrisk - cnt) / np.maximum(nrisk - 1, 1), 0.0)
        diff = (d1 - E1).sum(1); var = V.sum(1)
        num_o += diff[0]; V_o += var[0]; num_p += diff[1:]; V_p += var[1:]
    Z = num_o / np.sqrt(V_o); Zp = num_p / np.sqrt(V_p)
    return float(Z), Zp


def analyse(rows, L, rich, poor, nperm, seed):
    fly = np.array([r[0] for r in rows]); zt = np.array([r[1] for r in rows]); T = np.array([r[2] for r in rows])
    k = np.array([r[3] for r in rows]); R = (T - L).astype(float)
    use = (k >= rich) | (k <= poor); g = (k >= rich).astype(int)
    strata = fly * 10 + np.floor(zt / 4).astype(int)
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    Z, Zp = logrank_perm(R[use], g[use], strata[use], nperm, rng)
    p = (1 + int((np.abs(Zp) >= abs(Z)).sum())) / (nperm + 1)
    per = {}
    for c in np.unique(fly):
        a_ = R[(fly == c) & (k >= rich)]; b_ = R[(fly == c) & (k <= poor)]
        if len(a_) >= 5 and len(b_) >= 5:
            per[int(c)] = dict(n_q=len(a_), n_act=len(b_), med_q=float(np.median(a_)), med_act=float(np.median(b_)))
    return dict(L=L, n_bouts=int(len(rows)), k_counts=np.bincount(k, minlength=L // 60 + 1).tolist(),
                n_qstart=int((k >= rich).sum()), n_activestart=int((k <= poor).sum()),
                median_R_qstart=float(np.median(R[k >= rich])) if (k >= rich).any() else None,
                median_R_activestart=float(np.median(R[k <= poor])), Z=Z, p=p, nperm=nperm,
                null_absZ_95=float(np.percentile(np.abs(Zp), 95)), per_fly=per,
                flies_qstart_shorter=int(sum(v['med_q'] < v['med_act'] for v in per.values())), flies_compared=len(per))


def main(outdir, mode):
    os.makedirs(outdir, exist_ok=True); F, deaths = o1.load_all()
    rows, excl = wake_table(F, deaths, 300)
    if mode == 'check':
        rng = np.random.default_rng(np.random.SeedSequence([12002])); ps = []
        k = np.array([r[3] for r in rows]); frac = float((k >= 3).mean() / ((k >= 3) | (k <= 1)).mean())
        for i in range(200):
            fake = [(r[0], r[1], r[2], int(5 if rng.random() < frac else 0)) for r in rows]
            ps.append(analyse(fake, 300, 3, 1, 999, [12003, i])['p'])
        ks = kstest(ps, 'uniform'); share = float(np.mean(np.array(ps) < 0.05))
        out = dict(draws=200, nperm=999, ks_p=float(ks.pvalue), share_p_below_005=share, excluded_bouts=excl)
        print('CHECK (random labels):', out); json.dump(dict(out, p=ps), open(os.path.join(outdir, 'q1_check.json'), 'w'))
        return
    res = dict(excluded_bouts_no_reference_L300=excl, primary=analyse(rows, 300, 3, 1, 9999, [12001]))
    rows6, excl6 = wake_table(F, deaths, 600)
    res['sensitivity_L600'] = analyse(rows6, 600, 6, 2, 9999, [12001]); res['excluded_bouts_no_reference_L600'] = excl6
    for kk in ('primary', 'sensitivity_L600'):
        r = res[kk]
        print(f"{kk}: L {r['L']} bouts {r['n_bouts']} k counts {r['k_counts']} | Q-START {r['n_qstart']} median R "
              f"{r['median_R_qstart']} s, ACTIVE-START {r['n_activestart']} median R {r['median_R_activestart']} s | "
              f"Z {r['Z']:.2f} p {r['p']:.4f} (null |Z| 95% {r['null_absZ_95']:.2f}) | flies Q-start shorter "
              f"{r['flies_qstart_shorter']}/{r['flies_compared']}")
    json.dump(res, open(os.path.join(outdir, 'q1_results.json'), 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
