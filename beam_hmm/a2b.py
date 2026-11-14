"""A2b — PRESPEC 21 (+ amendment B) ported to the pilot, LD-adapted (PRESPEC A2b).
Statistic and nulls are line-by-line ports of hmm_audit/planB_21.py (origin/claude/trusting-allen-7nxrsf).
Usage: python -I a2b.py <outdir> <analysis weeks, comma-separated> <all supplied weeks...>"""
import sys, os, json, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.stats import wilcoxon, binom
import pipeline as pl, p1e

SETTINGS = {'L300': (300, 10), 'L120': (120, 5)}; NFLIP = 2000


def week_bouts(x, t, v, lights_on_h, week, L, b, switch_exclusion=True):
    """Per-bout binned start-minus-end differences (activity and count intensity), block keys, levels."""
    a = pl.activity_1s(x, t, v)                                     # 1 active, 0 inactive, -1 excluded
    n = len(a); sec = np.floor(t).astype(int)
    chg = np.r_[False, (x[1:] != x[:-1]) & v[1:] & v[:-1]]
    cnt = np.bincount(sec[chg], minlength=n)[:n].astype(float)
    act = (a == 1).astype(float); sb = pl.sleep_bouts(a)
    zt = lambda s: ((s / 3600.0) - lights_on_h) % 24
    D, Dc, blk, lev = [], [], [], []
    for k in range(len(sb) - 1):
        w0 = sb[k][0] + sb[k][1]; w1 = sb[k + 1][0]
        if w1 - w0 < 2 * L or (a[w0:w1] == -1).any():
            continue
        z0, z1 = zt(w0), zt(w1 - 1); dur_h = (w1 - w0) / 3600
        if switch_exclusion and not ((1 <= z0 < 11 and 1 <= z1 < 11 and dur_h < 10) or (13 <= z0 < 23 and 13 <= z1 < 23 and dur_h < 10)):
            continue
        for arr, lst in ((act, D), (cnt, Dc)):
            st = arr[w0:w0 + L].reshape(-1, b).mean(1); en = arr[w1 - L:w1][::-1].reshape(-1, b).mean(1); lst.append(st - en)
        lev.append(act[w0:w0 + 60].mean() - act[w1 - 60:w1].mean())
        blk.append(week * 100 + ((w0 + w1) // 2) // 86400)
    return D, Dc, blk, lev


def S_of(C, signs):
    return np.sum(((signs @ C) / C.shape[0]) ** 2, axis=1)


def fly_stats(D, blk, rng):
    C = D - D.mean(1, keepdims=True); n = len(C); S = float(np.sum(C.mean(0) ** 2))
    sg = rng.choice([-1.0, 1.0], size=(NFLIP, n)); Sn = S_of(C, sg)
    p = (1 + np.sum(Sn >= S)) / (NFLIP + 1); z = (S - Sn.mean()) / Sn.std()
    ub = np.unique(blk); cfg = np.array(list(itertools.product([-1.0, 1.0], repeat=len(ub))))
    sgb = cfg[:, np.searchsorted(ub, blk)]; Sb = S_of(C, sgb)
    pb = np.sum(Sb >= S) / len(Sb); zb = (S - Sb.mean()) / Sb.std() if Sb.std() > 0 else np.nan
    return dict(S=S, p=float(p), z=float(z), pb=float(pb), zb=float(zb), nblocks=int(len(ub)), nconf=int(len(cfg)), n=int(n))


def group(rows, zk, pk):
    z = np.array([r[zk] for r in rows]); p = np.array([r[pk] for r in rows]); z = z[np.isfinite(z)]
    w = wilcoxon(z).pvalue if len(z) and np.any(z != 0) else 1.0; k = int(np.sum(p < 0.05))
    return dict(n=len(rows), zmed=float(np.median(z)), wil=float(w), frac=k / len(p), binom=float(binom.sf(k - 1, len(p), 0.05)))


if __name__ == '__main__':
    SWX = '--no-switch-exclusion' not in sys.argv          # POST HOC option (ledger 2026-11-14): all bouts
    args = [a for a in sys.argv[1:] if a != '--no-switch-exclusion']
    outdir = args[0]; aw = [int(w) for w in args[1].split(',')]; weeks = [int(w) for w in args[2:]]
    os.makedirs(outdir, exist_ok=True)
    F = {}
    for w in weeks:
        X, t = pl.load_week(f'data/week{w}/week{w}.npz'); P, info = p1e.preprocess(X, t)
        F[w] = dict(P=P, t=t, info=info, stamp=p1e.phase_stamp(P, t, info) if w in aw else None)
    chans = [c + 1 for c in pl.CH]; res = {}
    for st, (L, b) in SETTINGS.items():
        rows = []; excl = []; per_fly = {}
        for j, c in enumerate(chans):
            chg = [p1e.change_times(F[w]['P'][:, j], F[w]['t']) for w in weeks]
            death = p1e.death_times(chg, [F[w]['t'][-1] for w in weeks])
            D, Dc, blk, lev = [], [], [], []
            for fi, w in enumerate(weeks):
                if w not in aw or not F[w]['info'][j]['included']:
                    continue
                x = F[w]['P'][:, j]; t = F[w]['t']; v = p1e.valid_mask(t, w, death, fi)
                if v.sum() == 0 or (x[v].max() - x[v].min()) < 8:
                    continue
                d_, dc_, b_, l_ = week_bouts(x, t, v, F[w]['stamp']['lights_on_h'], w, L, b, SWX)
                D += d_; Dc += dc_; blk += b_; lev += l_
            if not D:
                continue
            if len(D) < 5:
                excl.append(c); continue
            D, Dc, blk, lev = np.array(D), np.array(Dc), np.array(blk), np.array(lev)
            r = fly_stats(D, blk, np.random.default_rng(27_000 + j)); r['fly'] = c
            rc = fly_stats(Dc, blk, np.random.default_rng(28_000 + j)); r['z_counts'] = rc['z']
            sg = np.random.default_rng(29_000 + j).choice([-1.0, 1.0], size=(NFLIP, len(lev)))
            ln = (sg @ lev) / len(lev); r['level'] = float(lev.mean()); r['p_level'] = float((1 + np.sum(np.abs(ln) >= abs(lev.mean()))) / (NFLIP + 1))
            r['mean_d'] = D.mean(0).tolist(); r['_C'] = D - D.mean(1, keepdims=True); r['_blk'] = blk; rows.append(r)
        gr = np.random.default_rng(30_000); Sobs = sum(r['S'] for r in rows); Sn = np.zeros(NFLIP)
        for r in rows:
            ub = np.unique(r['_blk']); sg = gr.choice([-1.0, 1.0], size=(NFLIP, len(ub)))[:, np.searchsorted(ub, r['_blk'])]
            Sn += S_of(r['_C'], sg)
        for r in rows:
            del r['_C'], r['_blk']
        gp = group(rows, 'z', 'p') if rows else None; gb = group(rows, 'zb', 'pb') if rows else None
        gbp = float((1 + np.sum(Sn >= Sobs)) / (NFLIP + 1)) if rows else None
        va = bool(gp and gp['wil'] < 0.05 and gp['zmed'] > 0 and gp['binom'] < 0.05)
        vb = bool(gb and gb['wil'] < 0.05 and gb['zmed'] > 0); vg = bool(gbp is not None and gbp < 0.05)
        verdict = ('TIME-ASYMMETRIC SHAPE' if va and vb and vg else
                   'not robust to dependence between bouts within a day (no claim)' if va else 'no shape asymmetry detected')
        lv = np.array([r['level'] for r in rows]); zc = np.array([r['z_counts'] for r in rows])
        md = np.median([r['mean_d'] for r in rows], 0).tolist() if rows else []
        res[st] = dict(rows=rows, excluded_lt5=excl, per_bout=gp, per_fly_block=gb, group_block_p=gbp, verdict=verdict,
                       level_median=float(np.median(lv)) if len(lv) else None, level_wilcoxon=float(wilcoxon(lv).pvalue) if len(lv) > 4 else None,
                       counts_shape_zmed=float(np.median(zc)) if len(zc) else None, counts_shape_wilcoxon=float(wilcoxon(zc).pvalue) if len(zc) > 4 else None,
                       median_profile_start_minus_end=md)
        print(f"\n== {st}: flies {len(rows)} (excluded <5 bouts: {excl}); bouts/fly {[r['n'] for r in rows]}")
        if rows:
            print(f"   per-bout: z med {gp['zmed']:+.2f} Wilcoxon {gp['wil']:.3g} flies p<.05 {gp['frac']:.2f} (binom {gp['binom']:.3g}) | "
                  f"per-fly day-block: z med {gb['zmed']:+.2f} Wilcoxon {gb['wil']:.3g} | group block p {gbp:.3g}")
            print(f"   VERDICT: {verdict}{' [secondary setting]' if st == 'L120' else ''}")
            print(f"   level (first 60 s, start - end) median {res[st]['level_median']:+.3f} Wilcoxon {res[st]['level_wilcoxon']} | counts-shape z med {res[st]['counts_shape_zmed']:+.2f} Wilcoxon {res[st]['counts_shape_wilcoxon']}")
            print('   median d(tau) start - end, bins from the sleep boundary:', ' '.join(f'{v:+.3f}' for v in md[:12]))
    json.dump(dict(analysis_weeks=aw, switch_exclusion=SWX, results=res), open(os.path.join(outdir, f"a2b_weeks{'-'.join(map(str, aw))}{'' if SWX else '_allbouts'}.json"), 'w'),
              indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
