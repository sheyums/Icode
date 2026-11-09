# PRESPEC 29 Amendment A (a133e8d): hidden test-retest with SEPARATE HMM fits per half (days {0,2} vs {1,3}),
# the other half's minutes masked as missing.
import json, numpy as np, sqw, sqw2, sqw_ref as R
from multiprocessing import Pool
from sqw15 import seqs
from sqw2_common import em, label3b, degenerate
from sqw29 import shm_seq, day_J, retest, holm, SHF, nph, PH
G = ['CS', 'per0LD', 'per0DD', 'fumin', 'sss', 'shmns']
def one(a):
    g, i = a
    y, ph = shm_seq(i) if g == 'shmns' else (lambda s: (np.asarray(s[0]), np.asarray(s[1])))(sqw2.seqs(g, i, 0.5))
    y = np.asarray(y).astype(int); ph = np.asarray(ph); nd = len(y) // 1440; day = np.arange(len(y)) // 1440; out = dict(g=g, fly=i, halves={})
    for h, keep in (('A', (0, 2)), ('B', (1, 3))):
        mask = np.isin(day, keep)
        if not mask.any(): out['halves'][h] = None; continue
        yh = np.where(mask, y, -1)
        try:
            ll, pi, P, B, it = em(yh, ph, nph(g), 3, sqw2.ST, 9000 + 100 * G.index(g) + i + (0 if h == 'A' else 50000))
            P, B, pi = label3b(np.array(P), np.array(B), np.array(pi))
            if degenerate(B): out['halves'][h] = 'deg'; continue
            Z = R.ffbs(pi, P, B, yh, ph, n_paths=200, seed=7000 + i); Jd, _ = day_J(Z, ph, nph(g))
            out['halves'][h] = Jd.tolist()
        except Exception as e:
            out['halves'][h] = 'fail:' + str(e)[:60]
    return out
if __name__ == '__main__':
    jobs = [(g, i) for g in G for i in range(len(SHF) if g == 'shmns' else len(sqw.files(g)))]
    with Pool(4) as pool: HB = pool.map(one, jobs)
    res = {}
    for g in G:
        rr = [r for r in HB if r['g'] == g]; ok = [r for r in rr if isinstance(r['halves'].get('A'), list) and isinstance(r['halves'].get('B'), list)]
        bad = {k: sum(1 for r in rr for h in ('A', 'B') if isinstance(r['halves'].get(h), str) and r['halves'][h].startswith(k)) for k in ('deg', 'fail')}
        for p in range(nph(g)):
            # combine: per fly, a pseudo day-array with A days from fit A and B days from fit B
            recs = []
            for r in ok:
                A = np.array(r['halves']['A']); B = np.array(r['halves']['B']); Jd = np.where(np.isin(np.arange(len(A)), (0, 2))[:, None], A, B)
                recs.append(dict(Jd=Jd.tolist()))
            k = f"{g}|{PH[p] if nph(g) == 2 else 'whole'}"; res[k] = dict(T1=retest(recs, p, scored=False), n_ok=len(ok), **bad)
            print(k, res[k], flush=True)
    for k, a in zip(['sss|light', 'shmns|light'], holm([res[k]['T1']['p'] for k in ['sss|light', 'shmns|light']])): res[k]['T1']['holm'] = float(a)
    json.dump(dict(res=res, flies=HB), open(sqw.OUTD + 's29b_hidden_split.json', 'w'))
    for k in ('sss|light', 'shmns|light'): print('PRIMARY', k, res[k])
