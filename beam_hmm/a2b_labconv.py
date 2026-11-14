"""POST HOC (not pre-specified): A2b on the pilot with the LAB activity conversion (labconv.py) instead of
pipeline.activity_1s. Same flies, valid seconds, ZT and statistics as a2b.py; only the 1-s activity differs.
Usage: python -I a2b_labconv.py <outdir> <analysis weeks, comma-separated> <all supplied weeks...>"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.stats import wilcoxon
import pipeline as pl, p1e, a2b, labconv

outdir = sys.argv[1]; aw = [int(w) for w in sys.argv[2].split(',')]; weeks = [int(w) for w in sys.argv[3:]]
os.makedirs(outdir, exist_ok=True)
F = {}
for w in weeks:
    d = np.load(f'data/week{w}/week{w}.npz'); Xall = d['X'][:, 16:32]
    X, t = pl.load_week(f'data/week{w}/week{w}.npz'); P, info = p1e.preprocess(X, t)
    F[w] = dict(P=P, t=t, info=info, stamp=p1e.phase_stamp(P, t, info) if w in aw else None,
                lab=labconv.lab_activity_per_second(Xall) if w in aw else None)
res = {}; cmp = {}
for st, (L, b) in a2b.SETTINGS.items():
    rows = []
    for j, c in enumerate(c + 1 for c in pl.CH):
        chg = [p1e.change_times(F[w]['P'][:, j], F[w]['t']) for w in weeks]
        death = p1e.death_times(chg, [F[w]['t'][-1] for w in weeks])
        D, blk = [], []
        for fi, w in enumerate(weeks):
            if w not in aw or not F[w]['info'][j]['included']:
                continue
            x = F[w]['P'][:, j]; t = F[w]['t']; v = p1e.valid_mask(t, w, death, fi)
            if v.sum() == 0 or (x[v].max() - x[v].min()) < 8:
                continue
            a_ours = pl.activity_1s(x, t, v); cnt = F[w]['lab'][:, j]; n = min(len(a_ours), len(cnt))
            a = np.where(a_ours[:n] == -1, -1, (cnt[:n] > 0).astype(np.int8)).astype(np.int8)
            if st == 'L300':
                ok = a_ours[:n] >= 0
                cmp.setdefault(c, []).append(dict(week=w, seconds_differ=float((a[ok] != a_ours[:n][ok]).mean()),
                    sleep_min_day_ours=sum(nn for _, nn, _ in pl.sleep_bouts(a_ours[:n])) / 60 / (ok.sum() / 86400),
                    sleep_min_day_lab=sum(nn for _, nn, _ in pl.sleep_bouts(a)) / 60 / (ok.sum() / 86400)))
            act = (a == 1).astype(float); sb = pl.sleep_bouts(a); lo = F[w]['stamp']['lights_on_h']
            zt = lambda s: ((s / 3600.0) - lo) % 24
            for k in range(len(sb) - 1):
                w0 = sb[k][0] + sb[k][1]; w1 = sb[k + 1][0]
                if w1 - w0 < 2 * L or (a[w0:w1] == -1).any():
                    continue
                z0, z1 = zt(w0), zt(w1 - 1); dh = (w1 - w0) / 3600
                if not ((1 <= z0 < 11 and 1 <= z1 < 11 and dh < 10) or (13 <= z0 < 23 and 13 <= z1 < 23 and dh < 10)):
                    continue
                D.append(act[w0:w0 + L].reshape(-1, b).mean(1) - act[w1 - L:w1][::-1].reshape(-1, b).mean(1))
                blk.append(w * 100 + ((w0 + w1) // 2) // 86400)
        if len(D) >= 5:
            r = a2b.fly_stats(np.array(D), np.array(blk), np.random.default_rng(27_000 + j)); r['fly'] = c; r['mean_d'] = np.mean(D, 0).tolist(); rows.append(r)
    gp = a2b.group(rows, 'z', 'p'); gb = a2b.group(rows, 'zb', 'pb'); md = np.median([r['mean_d'] for r in rows], 0)
    res[st] = dict(per_bout=gp, per_fly_block=gb, median_profile=md.tolist(), bouts=[r['n'] for r in rows])
    print(f"{st}: flies {len(rows)} bouts {[r['n'] for r in rows]} | z med {gp['zmed']:+.2f} Wilcoxon {gp['wil']:.3g} flies p<.05 {gp['frac']:.2f} (binom {gp['binom']:.3g}) | "
          f"day-block z med {gb['zmed']:+.2f} p {gb['wil']:.3g}\n   d(tau): " + ' '.join(f'{v:+.3f}' for v in md[:12]))
for c, L_ in cmp.items():
    print(f"ch{c}: " + '; '.join(f"wk{r['week']} seconds differ {r['seconds_differ']:.3f}, sleep ours {r['sleep_min_day_ours']:.0f} vs lab {r['sleep_min_day_lab']:.0f} min/day" for r in L_))
json.dump(dict(analysis_weeks=aw, results=res, conversion_comparison=cmp), open(os.path.join(outdir, f"a2b_labconv_weeks{'-'.join(map(str, aw))}.json"), 'w'), indent=1, default=float)
