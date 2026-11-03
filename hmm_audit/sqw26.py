# PRESPEC 26 (9a8ff19): persistence of loop direction. Per fly S = fraction of consecutive loop pairs with the same
# direction; excess E = S_obs - mean S_null under a clock-matched (4-h windows, primary) and a 12-h (secondary) chain,
# 500 simulated records per fly. Wilcoxon on E, Holm within primary (per0 DD, per0 LD) and reference (CS LD, CS DD).
import json, numpy as np, sqw
from multiprocessing import Pool
from sqw15 import seqs
from scipy.stats import wilcoxon, mannwhitneyu
NR = 500; MINPAIRS = 10; WIN = 240
G = ['per0DD', 'per0LD', 'CS', 'CSDD', 'fumin', 'sss']
FWD = {(1, 0, 2), (0, 2, 1), (2, 1, 0)}
def loop_dirs(y):
    """Directions (+1/-1) of loops completed, in order, with a segment id (missing minutes start a new segment)."""
    y = np.asarray(y); ch = np.flatnonzero(np.diff(y)) + 1; seq = np.r_[y[0], y[ch]]
    out = []; seg = 0; stack = []
    for x in seq.tolist():
        if x < 0: stack = []; seg += 1; continue
        if x in stack:
            k = stack.index(x); loop = stack[k:]
            if len(loop) == 3: out.append((seg, 1 if tuple(loop) in FWD else -1))
            del stack[k + 1:]
        else: stack.append(x)
    return out
def S_stat(L):
    same = n = 0
    for (s1, d1), (s2, d2) in zip(L[:-1], L[1:]):
        if s1 == s2: n += 1; same += d1 == d2
    return (same / n if n else np.nan), n
def runs(L):
    r = []; cur = None; ln = 0
    for s, d in L:
        if cur is not None and (s, d) == cur: ln += 1
        else:
            if cur is not None: r.append(ln)
            cur = (s, d); ln = 1
    if cur is not None: r.append(ln)
    return r
def chain(y, idx, nm, fallback=None):
    C = np.zeros((nm, 3, 3)); ok = (y[:-1] >= 0) & (y[1:] >= 0); np.add.at(C, (idx[1:][ok], y[:-1][ok], y[1:][ok]), 1); P = np.empty_like(C)
    for m in range(nm):
        for s in range(3):
            r = C[m, s].sum(); P[m, s] = C[m, s] / r if r > 0 else (fallback[m][s] if fallback is not None else np.full(3, 1 / 3))
    return P
def simulate(P, idx, y0, T, R, rng):
    Cc = np.cumsum(P, 2); Cc[..., -1] = 1.0; S = np.empty((R, T), np.int8); cur = np.full(R, y0); S[:, 0] = cur
    for t in range(1, T):
        cur = (rng.random(R)[:, None] > Cc[idx[t], cur]).sum(1); S[:, t] = cur
    return S
def fly(a):
    g, i = a; rng = np.random.default_rng(2600 + 1000 * G.index(g) + i)
    y, ph, _ = seqs(g, i); y = np.asarray(y).astype(int); ph = np.asarray(ph); T = len(y); miss = y < 0
    L = loop_dirs(y); S_obs, npairs = S_stat(L); dirs = np.array([d for _, d in L])
    rec = dict(g=g, fly=i, n_loops=len(L), n_pairs=npairs, frac_fwd=float(np.mean(dirs == 1)) if len(L) else None, S_obs=S_obs, runs_obs=runs(L))
    nP = int(ph.max()) + 1; Pph = chain(y, ph, nP); wpos = (np.arange(T) % 1440) // WIN; nw = 1440 // WIN
    fb = np.array([Pph[ph[np.argmax(wpos == m)]] for m in range(nw)]); Pclk = chain(y, wpos, nw, fb)
    y0 = int(y[np.argmax(y >= 0)]) if np.any(y >= 0) else 0
    for nm, P, idx in (('clock', Pclk, wpos), ('phase', Pph, ph)):
        Ssim = []; rl = []
        X = simulate(P, idx, y0, T, NR, rng); X[:, miss] = -1
        for r in range(NR):
            Lr = loop_dirs(X[r]); s, _ = S_stat(Lr); Ssim.append(s)
            if r < 100: rl += runs(Lr)
        Ssim = np.array(Ssim, float); rec[f'S_null_{nm}'] = float(np.nanmean(Ssim)); rec[f'E_{nm}'] = (S_obs - float(np.nanmean(Ssim))) if np.isfinite(S_obs) else None
        rec[f'S_null_sd_{nm}'] = float(np.nanstd(Ssim)); rec[f'runs_null_{nm}'] = rl
    # descriptive: within-fly shuffle of the observed loop order (segments kept)
    if npairs:
        sh = []; segs = [s for s, _ in L]
        for _ in range(2000):
            dd = rng.permutation(dirs); sh.append(S_stat(list(zip(segs, dd.tolist())))[0])
        rec['S_shuffle'] = float(np.nanmean(sh))
    return rec
if __name__ == '__main__':
    jobs = [(g, i) for g in G for i in range(len(sqw.files(g)))]
    with Pool() as pool: recs = pool.map(fly, jobs)
    res = {}
    for g in G:
        R = [r for r in recs if r['g'] == g and r['n_pairs'] >= MINPAIRS and r['E_clock'] is not None]
        e = {nm: np.array([r[f'E_{nm}'] for r in R]) for nm in ('clock', 'phase')}
        out = dict(n_flies=len(R), n_excluded=sum(1 for r in recs if r['g'] == g) - len(R), loops=int(sum(r['n_loops'] for r in R)),
                   S_obs_med=float(np.median([r['S_obs'] for r in R])), S_shuffle_med=float(np.median([r['S_shuffle'] for r in R])),
                   frac_fwd_med=float(np.median([r['frac_fwd'] for r in R])))
        for nm in ('clock', 'phase'):
            out[nm] = dict(S_null_med=float(np.median([r[f'S_null_{nm}'] for r in R])), E_med=float(np.median(e[nm])), E_mean=float(e[nm].mean()),
                           frac_E_pos=float(np.mean(e[nm] > 0)), p=float(wilcoxon(e[nm]).pvalue))
        res[g] = out
    for fam in (['per0DD', 'per0LD'], ['CS', 'CSDD']):
        for nm in ('clock', 'phase'):
            ps = [res[g][nm]['p'] for g in fam]; o = np.argsort(ps); adj = [0, 0]; run = 0
            for j, k in enumerate(o): run = max(run, min(1, (2 - j) * ps[k])); adj[k] = run
            for g, a in zip(fam, adj): res[g][nm]['holm'] = float(a)
    sec = {}
    for a, b in (('per0DD', 'CSDD'), ('per0LD', 'CS')):
        ea = [r['E_clock'] for r in recs if r['g'] == a and r['n_pairs'] >= MINPAIRS]; eb = [r['E_clock'] for r in recs if r['g'] == b and r['n_pairs'] >= MINPAIRS]
        sec[f'{a}_vs_{b}'] = dict(diff_med=float(np.median(ea) - np.median(eb)), p=float(mannwhitneyu(ea, eb).pvalue))
    json.dump(dict(genotypes=res, secondary=sec, flies=recs), open(sqw.OUTD + 's26_loopseq.json', 'w'))
    print('genotype flies(excl) loops | frac fwd | S obs | S shuffle | CLOCK null S, E med (mean), E>0, p, Holm | PHASE null S, E med, p, Holm')
    for g, v in res.items():
        c, p = v['clock'], v['phase']
        print(f"{g:7s} {v['n_flies']:3d}({v['n_excluded']}) {v['loops']:5d} | {v['frac_fwd_med']:.2f} | {v['S_obs_med']:.3f} | {v['S_shuffle_med']:.3f} | "
              f"{c['S_null_med']:.3f}, {c['E_med']:+.3f} ({c['E_mean']:+.3f}), {c['frac_E_pos']:.2f}, p {c['p']:.4f}, Holm {c.get('holm', float('nan')):.3f} | "
              f"{p['S_null_med']:.3f}, {p['E_med']:+.3f}, p {p['p']:.4f}, Holm {p.get('holm', float('nan')):.3f}")
    print('secondary', json.dumps(sec))
