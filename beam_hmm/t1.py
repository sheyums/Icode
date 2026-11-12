"""T1 per-fly procedures: clock-model selection at K*, and the D1-B procedure with M* transitions (C6/T-G3).
Truth-blind (the gate adds truth-side evaluation)."""
import numpy as np
import hmm_cat as hc, iohmm as io, d1, sim_positions as sp

M = 4; KS = d1.KS; NB = d1.NB; MODELS = (0, 1, 2, 3)


def k1_model(seqs):
    ys = np.concatenate([y for y, _ in seqs]); ys = ys[ys >= 0]
    f = np.bincount(ys, minlength=M) + 1e-9; f = f / f.sum()
    return (np.ones(1), np.ones((io.WPD, 1, 1)), f[None, :])


def folds(nd):
    odd = list(range(0, nd, 2)); even = list(range(1, nd, 2))
    return ((odd, even), (even, odd))


def clock_selection(seqs_w, K, fly_idx, restarts=10):
    """Held-out logL per window for M0..M3 at K (two-fold by day). M1-M3 warm-start from that fold's M0 fit."""
    nd = len(seqs_w); res = {m: [] for m in MODELS}; fits = {m: [] for m in MODELS}
    for fold, (tr, te) in enumerate(folds(nd)):
        tr_s = [seqs_w[i] for i in tr]; te_s = [seqs_w[i] for i in te]
        nobs = sum((y >= 0).sum() for y, _ in te_s); m0 = None
        for m in MODELS:
            ss = np.random.SeedSequence([5001, fly_idx, K, m, fold])
            init = None if m == 0 else (m0[0][0], m0[1], m0[0][2])
            model, beta, ll, lls = io.fit(tr_s, K, M, m, ss, n_restarts=restarts, init=init)
            if m == 0:
                m0 = (model, beta)
            res[m].append(hc.loglik(model, te_s) / nobs); fits[m].append((model, beta, lls))
    held = {m: float(np.mean(res[m])) for m in MODELS}
    return held, fits


def heldout_curve_clock(seqs_w, fly_idx, Kset, model_id, base, restarts):
    res = {}
    for K in Kset:
        hl, fits = [], []
        for fold, (tr, te) in enumerate(folds(len(seqs_w))):
            tr_s = [seqs_w[i] for i in tr]; te_s = [seqs_w[i] for i in te]
            if K == 1:
                model = k1_model(tr_s); lls = []
            else:
                model, beta, ll, lls = io.fit(tr_s, K, M, model_id, np.random.SeedSequence([base, fly_idx, K, model_id, fold]),
                                              n_restarts=restarts)
            nobs = sum((y >= 0).sum() for y, _ in te_s)
            hl.append(hc.loglik(model, te_s) / nobs); fits.append((model, lls))
        res[K] = dict(heldout=float(np.mean(hl)), fits=fits)
    return res


def dwell_clock(P, k):
    return float(np.mean(10.0 / np.maximum(1 - P[:, k, k], 1e-12)))


def run_d1_clock(x, t, valid, fly_idx, model_id, seed_base=5001, boot_base=7001):
    """The D1-B procedure (K sweep, candidates, C1 bootstrap, C2, C3) with transitions of clock model model_id."""
    y = sp.window_symbols(x, valid); lab = d1.sqw_per_window(x, t, valid)
    n = min(len(y), len(lab)) // io.WPD * io.WPD; y = y[:n]; lab = lab[:n]
    seqs = io.to_days_w(y); lab_days = [lab[d*io.WPD:(d+1)*io.WPD] for d in range(len(seqs))]
    out = dict(model_id=model_id)
    curve = heldout_curve_clock(seqs, fly_idx, KS, model_id, seed_base, 10)
    out['heldout'] = {K: curve[K]['heldout'] for K in KS}
    deg = {K: any(d1.degenerate(f[0][2]) for f in curve[K]['fits']) for K in KS}; out['degenerate'] = deg
    elig = [K for K in KS if not deg[K]]
    Kstar = max(elig, key=lambda K: curve[K]['heldout']) if elig else None; out['Kstar'] = Kstar
    if Kstar is None:
        return out
    full, _, _, _ = io.fit(seqs, Kstar, M, model_id, np.random.SeedSequence([seed_base, fly_idx, Kstar, model_id, 9]), n_restarts=10)
    desc = d1.describe(full, seqs, lab_days)
    desc['dwell'] = [dwell_clock(full[1], k) for k in range(Kstar)]
    out['types'] = [list(tp) for tp in desc['types']]; out['dwell_s'] = desc['dwell']
    cand = d1.candidates(desc); out['candidates'] = cand
    Kn = Kstar - 1
    prev = heldout_curve_clock(seqs, fly_idx, [Kn], model_id, seed_base, 10)[Kn]
    gain_obs = curve[Kstar]['heldout'] - prev['heldout']
    null = (io.fit(seqs, Kn, M, model_id, np.random.SeedSequence([seed_base, fly_idx, Kn, model_id, 9]), n_restarts=10)[0]
            if Kn >= 2 else k1_model(seqs))
    rng = np.random.default_rng(np.random.SeedSequence([boot_base, fly_idx, model_id])); gains = []
    for b in range(NB):
        sim = hc.simulate(null, [(len(yy), w) for yy, w in seqs], [yy < 0 for yy, _ in seqs], rng)
        cb = heldout_curve_clock(sim, fly_idx * 1000 + b, [Kn, Kstar], model_id, boot_base, 3)
        gains.append(cb[Kstar]['heldout'] - cb[Kn]['heldout'])
    thr = float(np.percentile(gains, 95))
    out['C1'] = dict(gain_obs=float(gain_obs), null_gains=list(map(float, gains)), thr95=thr, pass_=bool(gain_obs > thr))
    halves = [curve[Kstar]['fits'][0][0], curve[Kstar]['fits'][1][0]]
    nd = len(seqs); odd = list(range(0, nd, 2)); even = list(range(1, nd, 2))
    hdesc = [d1.describe(halves[0], [seqs[i] for i in odd], [lab_days[i] for i in odd]),
             d1.describe(halves[1], [seqs[i] for i in even], [lab_days[i] for i in even])]
    per = {}
    for k in range(Kstar):
        tp = desc['types'][k]; Bk = full[2][k]
        match = []
        for h in (0, 1):
            ks = [j for j, tj in enumerate(hdesc[h]['types']) if tuple(tj) == tuple(tp)]
            match.append(min(ks, key=lambda j: hc.js_bits(halves[h][2][j], Bk)) if ks else None)
        if None in match:
            c3 = False
        else:
            a, b = match; Ba, Bb = halves[0][2][a], halves[1][2][b]; dab = hc.js_bits(Ba, Bb)
            oa = min([hc.js_bits(Ba, halves[0][2][j]) for j in range(Kstar) if j != a] or [np.inf])
            ob = min([hc.js_bits(Bb, halves[1][2][j]) for j in range(Kstar) if j != b] or [np.inf])
            c3 = bool(dab < oa and dab < ob)
        per[k] = dict(type=list(tp), candidate=k in cand, dwell_s=desc['dwell'][k], C2=bool(desc['dwell'][k] >= 30), C3=c3)
    out['states'] = per
    return out
