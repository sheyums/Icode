"""PRESPEC 21 (+ amendment B, superseding A): 1-s wake-bout shape asymmetry beyond level; per-bout null, per-fly whole-day
block null (Wilcoxon of z only) and group-level fly-day block null.
Usage: python3 -I planB_21.py <data root> <CS DD folder>"""
import sys, os, glob, json, itertools, numpy as np
from scipy.stats import wilcoxon, binom
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R
SETS = {'per0DD': 'per0DD_activity_files_1sec_res', 'CSDD': None, 'CS': 'CS_activity_files_1sec_res',
        'per0LD': 'per0LD_activity_1sec_res', 'fumin': 'Fmn_activity_files_1sec', 'sss': 'sleepless_activity_1sec_res'}
GIDX = {g: k for k, g in enumerate(SETS)}
PRIMARY = ('per0DD', 'CSDD'); SETTINGS = {'L300': (300, 10), 'L120': (120, 5)}; NFLIP = 2000

def files(d): return sorted(f for f in glob.glob(os.path.join(d, '**', '*'), recursive=True) if os.path.isfile(f))

def traces(root, csdd):
    out = {}
    import planB_20u as U
    _, _, names, tl = U.load_ref(); up = files(os.path.join(root, SETS['per0DD'])); X = []
    for i, nm in enumerate(names):
        m = [f for f in up if os.path.basename(f).startswith(nm)][0]; x = R.read_counts(m)[:tl[i]]; X.append(x[:(len(x) // 86400) * 86400])
    out['per0DD'] = X
    for g, sub in SETS.items():
        if g == 'per0DD': continue
        d = csdd if g == 'CSDD' else os.path.join(root, sub)
        out[g] = [x[:(len(x) // 86400) * 86400] for x in (R.read_counts(f) for f in files(d))]
    return out

def windows(x, L, b):
    a = (np.asarray(x) > 0).astype(float); c = np.asarray(x, float); sb = R.sleep_bouts(x, True)
    D = []; Dc = []; blk = []; lev = []
    for k in range(len(sb) - 1):
        w0 = sb[k][0] + sb[k][1]; w1 = sb[k + 1][0]
        if w1 - w0 < 2 * L: continue
        for arr, lst in ((a, D), (c, Dc)):
            st = arr[w0:w0 + L].reshape(-1, b).mean(1); en = arr[w1 - L:w1][::-1].reshape(-1, b).mean(1); lst.append(st - en)
        lev.append(a[w0:w0 + 60].mean() - a[w1 - 60:w1].mean())
        blk.append(((w0 + w1) // 2) // 86400)
    return np.array(D), np.array(Dc), np.array(blk), np.array(lev)

def S_of(Dc_centred, signs):
    return np.sum(((signs @ Dc_centred) / Dc_centred.shape[0])**2, axis=1)

def fly_stats(D, blk, rng):
    C = D - D.mean(1, keepdims=True); n = len(C); S = float(np.sum(C.mean(0)**2))
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

def holm(d):
    ps = sorted((v, k) for k, v in d.items()); out = {}; run = 0
    for i, (v, k) in enumerate(ps): run = max(run, min(1, (len(ps) - i) * v)); out[k] = run
    return out

def main(root, csdd):
    T = traces(root, csdd); res = {}
    for st, (L, b) in SETTINGS.items():
        res[st] = {}
        for g, X in T.items():
            rows = []; excl = 0; prof_s = []; prof_e = []
            for i, x in enumerate(X):
                D, Dc, blk, lev = windows(x, L, b)
                if len(D) < 5: excl += 1; continue
                rng = np.random.default_rng(27_000 + 1000 * GIDX[g] + i)
                r = fly_stats(D, blk, rng); r['fly'] = i
                rc = fly_stats(Dc, blk, np.random.default_rng(28_000 + 1000 * GIDX[g] + i)); r['z_counts'] = rc['z']; r['p_counts'] = rc['p']
                sg = np.random.default_rng(29_000 + 1000 * GIDX[g] + i).choice([-1.0, 1.0], size=(NFLIP, len(lev)))
                ln = (sg @ lev) / len(lev); r['level'] = float(lev.mean()); r['p_level'] = float((1 + np.sum(np.abs(ln) >= abs(lev.mean()))) / (NFLIP + 1))
                r['mean_d'] = D.mean(0).tolist(); r['_C'] = D - D.mean(1, keepdims=True); r['_blk'] = blk; rows.append(r)
            # group-level block null (amendment B (b)): each fly-day flipped independently, 2000 configs, sum of S
            gr = np.random.default_rng(30_000 + GIDX[g]); Sobs = sum(r['S'] for r in rows); Sn = np.zeros(NFLIP)
            for r in rows:
                ub = np.unique(r['_blk']); sg = gr.choice([-1.0, 1.0], size=(NFLIP, len(ub)))[:, np.searchsorted(ub, r['_blk'])]
                Sn += S_of(r['_C'], sg)
            for r in rows: del r['_C'], r['_blk']
            res[st][g] = dict(rows=rows, excluded=excl, group_block_p=float((1 + np.sum(Sn >= Sobs)) / (NFLIP + 1)))
    for st in SETTINGS:
        print(f"\n===== setting {st} (L {SETTINGS[st][0]} s, bins {SETTINGS[st][1]} s) =====")
        gp = {g: group(res[st][g]['rows'], 'z', 'p') for g in PRIMARY}; gb = {g: group(res[st][g]['rows'], 'zb', 'pb') for g in PRIMARY}
        hp = holm({g: gp[g]['wil'] for g in PRIMARY}); hb = holm({g: gb[g]['wil'] for g in PRIMARY})
        for g in SETS:
            rows = res[st][g]['rows']; nb = [r['n'] for r in rows]
            line = (f"{g:7s} flies {len(rows)} (excluded <5 bouts: {res[st][g]['excluded']}) | bouts/fly median {np.median(nb):.0f} "
                    f"range {min(nb)}-{max(nb)} | ")
            if g in PRIMARY:
                a, c = gp[g], gb[g]
                va = hp[g] < 0.05 and a['zmed'] > 0 and a['binom'] < 0.05
                vb = hb[g] < 0.05 and c['zmed'] > 0; vg = res[st][g]['group_block_p'] < 0.05
                verdict = ('TIME-ASYMMETRIC SHAPE' if va and vb and vg else
                           'not robust to dependence between bouts within a day (no claim)' if va else 'no shape asymmetry detected')
                line += (f"per-bout: z med {a['zmed']:+.2f} Wilcoxon {a['wil']:.3g} Holm {hp[g]:.3g} flies p<.05 {a['frac']:.2f} (binom {a['binom']:.3g}) | "
                         f"per-fly day-block: z med {c['zmed']:+.2f} Holm {hb[g]:.3g} | group fly-day block p {res[st][g]['group_block_p']:.3g} -> "
                         f"{verdict if st == 'L300' else verdict + ' [secondary setting]'}")
                res[st][g]['verdict'] = verdict
            else:
                a = group(rows, 'z', 'p'); line += f"[descriptive] per-bout z med {a['zmed']:+.2f} p {a['wil']:.3g} flies p<.05 {a['frac']:.2f}"
            lv = np.array([r['level'] for r in rows]); zc = np.array([r['z_counts'] for r in rows])
            line += (f"\n          secondary: level (first 60 s, start - end) median {np.median(lv):+.3f} Wilcoxon {wilcoxon(lv).pvalue:.3g} | "
                     f"counts-intensity shape z med {np.median(zc):+.2f} Wilcoxon {wilcoxon(zc).pvalue:.3g}")
            md = np.median([r['mean_d'] for r in rows], 0)
            line += "\n          median per-fly d(tau) start-minus-end, bins from boundary: " + ' '.join(f"{v:+.3f}" for v in md[:12])
            print(line)
            print(f"          bouts per fly: {nb}")
    json.dump(res, open(os.path.join(HERE, 'planB_21.json'), 'w'), default=float)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
