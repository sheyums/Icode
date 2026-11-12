"""D1 per-fly procedure (PRESPEC D1 + D1-A). Input: 0.4-s readings (filtered), times, validity.
Returns a JSON-able dict. Truth-blind: nothing here sees the hidden state."""
import numpy as np
import pipeline as pl
import hmm_cat as hc
import sim_positions as sp

M = 12; NP = 2; KS = list(range(2, 9)); WPD = 8640          # windows per day
NB = 19


def to_days(y):
    """Split window symbols into fly-day sequences with phases."""
    nd = len(y) // WPD; seqs = []
    w = np.arange(WPD); ph = ((w * 10) % 86400 >= 43200).astype(np.int64)
    for d in range(nd):
        seqs.append((y[d*WPD:(d+1)*WPD].copy(), ph.copy()))
    return seqs


def sqw_per_window(x, t, valid):
    a = pl.activity_1s(x, t, valid)
    ymin, refs, bouts = pl.sqw_minutes(a)
    return np.repeat(ymin, 6)          # 6 windows per minute


def type_of(B, G, lab):
    """D1-A item 1. B: emission row (12), G: posterior mass per window (T), lab: S/Q/W per window."""
    Bm = B.reshape(4, 3)
    still = Bm[0].sum(); mv = Bm[1:].sum(1)
    band = 'QUIESCENT' if still >= 0.9 else ('MIXED' if still >= 0.5 else 'ACTIVE')
    move = '-' if band == 'QUIESCENT' else ['FLICK', 'LOCAL', 'TRAVEL'][int(np.argmax(mv))]
    zone = ['BEYOND1', 'END', 'MIDDLE'][int(np.argmax(Bm.sum(0)))]
    share = np.array([G[lab == L].sum() for L in (pl.S, pl.Q, pl.W)]); share = share / max(share.sum(), 1e-300)
    label = ['S', 'Q', 'W'][int(np.argmax(share))]
    return (band, move, zone, label), share


def describe(model, seqs, lab_days):
    G = hc.posterior(model, seqs)
    Gc = np.concatenate(G); lab = np.concatenate(lab_days)
    pi, P, B = model; K = B.shape[0]
    types, shares = [], []
    for k in range(K):
        tp, sh = type_of(B[k], Gc[:, k], lab); types.append(tp); shares.append(sh)
    dwell = [float(np.mean([10.0 / max(1 - P[p, k, k], 1e-12) for p in range(NP)])) for k in range(K)]
    argmax = Gc.argmax(1)
    return dict(types=types, shares=shares, dwell=dwell, argmax=argmax)


def candidates(desc):
    out = []
    labs = [tp[3] for tp in desc['types']]
    for k, (tp, sh) in enumerate(zip(desc['types'], desc['shares'])):
        splits = labs.count(tp[3]) > 1
        straddles = sh.max() < 0.6
        if splits or straddles:
            out.append(k)
    return out


def degenerate(B):
    K = B.shape[0]
    return any(hc.js_bits(B[i], B[j]) < 0.05 for i in range(K) for j in range(i + 1, K))


def heldout_gain_curve(seqs, fly_idx, Kset, base, restarts, odd, even):
    """Two-fold by day. Returns {K: dict(heldout per window, fits)}."""
    res = {}
    for K in Kset:
        hl, fits = [], []
        for fold, (tr, te) in enumerate(((odd, even), (even, odd))):
            tr_s = [seqs[i] for i in tr]; te_s = [seqs[i] for i in te]
            if K == 1:
                ys = np.concatenate([y for y, _ in tr_s]); ys = ys[ys >= 0]
                f = np.bincount(ys, minlength=M) + 1e-9; f = f / f.sum()
                model = (np.ones(1), np.ones((NP, 1, 1)), f[None, :]); lls = []
            else:
                ss = np.random.SeedSequence([base, fly_idx, K, fold])
                model, ll, lls, it = hc.fit(tr_s, K, M, NP, ss, n_restarts=restarts)
            nobs = sum((y >= 0).sum() for y, _ in te_s)
            hl.append(hc.loglik(model, te_s) / nobs); fits.append((model, lls))
        res[K] = dict(heldout=float(np.mean(hl)), fits=fits)
    return res


def run_fly(x, t, valid, fly_idx, seed_base=1001, boot_base=2001, extra=None, log=print):
    days_y = sp.window_symbols(x, valid)
    lab = sqw_per_window(x, t, valid)
    n = min(len(days_y), len(lab)) // WPD * WPD
    y = days_y[:n]; lab = lab[:n]
    seqs = to_days(y); lab_days = [lab[d*WPD:(d+1)*WPD] for d in range(len(seqs))]
    nd = len(seqs); odd = list(range(0, nd, 2)); even = list(range(1, nd, 2))
    out = dict(n_days=nd)
    # K sweep
    curve = heldout_gain_curve(seqs, fly_idx, KS, seed_base, 10, odd, even)
    out['heldout'] = {K: curve[K]['heldout'] for K in KS}
    out['restart_lls'] = {K: [list(map(float, f[1])) for f in curve[K]['fits']] for K in KS}
    deg = {K: any(degenerate(f[0][2]) for f in curve[K]['fits']) for K in KS}
    out['degenerate'] = deg
    elig = [K for K in KS if not deg[K]]
    Kstar = max(elig, key=lambda K: curve[K]['heldout']) if elig else None
    out['Kstar'] = Kstar
    if Kstar is None:
        return out
    log(f'fly {fly_idx}: K* = {Kstar}')
    # full fit at K*, BIC beside
    ss = np.random.SeedSequence([seed_base, fly_idx, Kstar, 9])
    full, ll_full, lls_full, _ = hc.fit(seqs, Kstar, M, NP, ss, n_restarts=10)
    nobs = int((y >= 0).sum())
    out['bic'] = {}
    desc = describe(full, seqs, lab_days)
    out['types'] = [list(tp) for tp in desc['types']]; out['shares'] = [list(map(float, s)) for s in desc['shares']]
    out['dwell_s'] = desc['dwell']
    cand = candidates(desc); out['candidates'] = cand
    # C1: observed gain vs bootstrap from K*-1
    if Kstar - 1 >= 2:
        prev = heldout_gain_curve(seqs, fly_idx, [Kstar - 1], seed_base, 10, odd, even)[Kstar - 1]
    else:
        prev = heldout_gain_curve(seqs, fly_idx, [1], seed_base, 10, odd, even)[1]
    gain_obs = curve[Kstar]['heldout'] - prev['heldout']
    # null model = K*-1 fitted on all days
    if Kstar - 1 >= 2:
        null, _, _, _ = hc.fit(seqs, Kstar - 1, M, NP, np.random.SeedSequence([seed_base, fly_idx, Kstar - 1, 9]), n_restarts=10)
    else:
        null = (np.ones(1), np.ones((NP, 1, 1)), prev['fits'][0][0][2])
    rng = np.random.default_rng(np.random.SeedSequence([boot_base, fly_idx]))
    gains = []
    for b in range(NB):
        sim = hc.simulate(null, [(len(yy), ph) for yy, ph in seqs], [yy < 0 for yy, _ in seqs], rng)
        Kset = [Kstar - 1, Kstar] if Kstar - 1 >= 2 else [1, Kstar]
        cb = heldout_gain_curve(sim, fly_idx * 1000 + b, Kset, boot_base, 3, odd, even)
        gains.append(cb[Kstar]['heldout'] - cb[Kset[0]]['heldout'])
    thr = float(np.percentile(gains, 95))
    out['C1'] = dict(gain_obs=float(gain_obs), null_gains=list(map(float, gains)), thr95=thr, pass_=bool(gain_obs > thr))
    # C2, C3, C4 per candidate
    halves = [curve[Kstar]['fits'][0][0], curve[Kstar]['fits'][1][0]]     # odd-days fit, even-days fit
    hdesc = [describe(halves[0], [seqs[i] for i in odd], [lab_days[i] for i in odd]),
             describe(halves[1], [seqs[i] for i in even], [lab_days[i] for i in even])]
    refits = {}
    for nm, kw in (('a_flick_merged', dict(merge_flick=True)), ('b_beyond_merged', dict(merge_beyond=True))):
        ym = sp.window_symbols(x, valid, **kw)[:n]; sm = to_days(ym)
        mod, _, _, _ = hc.fit(sm, Kstar, M, NP, np.random.SeedSequence([seed_base, fly_idx, Kstar, 7 if nm[0] == 'a' else 8]), n_restarts=10)
        refits[nm] = np.concatenate(hc.posterior(mod, sm)).argmax(1)
    per = {}
    for k in range(Kstar):
        tp = desc['types'][k]; Bk = full[2][k]
        c2 = desc['dwell'][k] >= 30.0
        match = []
        for h in (0, 1):
            ks = [j for j, tj in enumerate(hdesc[h]['types']) if tuple(tj) == tuple(tp)]
            match.append(min(ks, key=lambda j: hc.js_bits(halves[h][2][j], Bk)) if ks else None)
        if None in match:
            c3 = False
        else:
            a, b = match; Ba, Bb = halves[0][2][a], halves[1][2][b]
            dab = hc.js_bits(Ba, Bb)
            oth_a = min([hc.js_bits(Ba, halves[0][2][j]) for j in range(Kstar) if j != a] or [np.inf])
            oth_b = min([hc.js_bits(Bb, halves[1][2][j]) for j in range(Kstar) if j != b] or [np.inf])
            c3 = bool(dab < oth_a and dab < oth_b)
        win = desc['argmax'] == k
        c4 = {}
        for nm, am in refits.items():
            best = 0.0
            for j in range(Kstar):
                wj = am == j; inter = (win & wj).sum(); uni = (win | wj).sum()
                best = max(best, inter / uni if uni else 0.0)
            c4[nm] = float(best)
        per[k] = dict(type=list(tp), candidate=k in cand, dwell_s=desc['dwell'][k], C2=bool(c2), C3=bool(c3),
                      C4_jaccard=c4, C4_pass=bool(all(v >= 0.5 for v in c4.values())))
    out['states'] = per
    out['argmax_full'] = desc['argmax']        # removed before saving; used by the gate
    return out
