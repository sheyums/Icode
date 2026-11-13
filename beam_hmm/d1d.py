"""D1-D per-fly procedure: sample-level (0.4-s) behavioural HMM (PRESPEC D1-D). Truth-blind."""
import numpy as np
import pipeline as pl, hmm_cat as hc, d1

M = 4; NP = 2; KS = list(range(2, 9)); SPD = 216000          # samples per day
NB = 19; TOL = 0.01; MAXIT = 300; R_OBS = 5


def step_symbols(x, valid=None):
    x = np.asarray(x).astype(np.int64)
    y = np.r_[0, np.minimum(np.abs(np.diff(x)), 3)].astype(np.int64)
    if valid is not None:
        v = np.asarray(valid); bad = ~v | ~np.r_[True, v[:-1]]; y[bad] = -1
    return y


def remove_flicker(x):
    x = np.asarray(x).copy()
    i = np.nonzero((x[:-2] == x[2:]) & (np.abs(x[1:-1] - x[:-2]) == 1))[0] + 1
    x[i] = x[i - 1]
    return x


def to_days(y):
    nd = len(y) // SPD
    ph = ((np.arange(SPD) * 0.4) % 86400 >= 43200).astype(np.int64)
    return [(y[d*SPD:(d+1)*SPD].copy(), ph.copy()) for d in range(nd)]


def labels_per_sample(x, t, valid, n):
    a = pl.activity_1s(x, t, valid); ymin, _, _ = pl.sqw_minutes(a)
    lab = np.repeat(ymin, 150)
    if len(lab) < n:
        lab = np.r_[lab, np.full(n - len(lab), -1)]
    return lab[:n]


def _fit(seqs, K, fly_idx, tag, fold, restarts=R_OBS, inits=None):
    return hc.fit(seqs, K, M, NP, np.random.SeedSequence([1001, fly_idx, K, tag, fold]), n_restarts=restarts,
                  max_iter=MAXIT, tol=TOL, inits=inits)


def k1_model(seqs):
    ys = np.concatenate([y for y, _ in seqs]); ys = ys[ys >= 0]
    f = np.bincount(ys, minlength=M) + 1e-9; f = f / f.sum()
    return (np.ones(1), np.ones((NP, 1, 1)), f[None, :])


def curve(seqs, fly_idx, Kset, restarts, tag=0, inits_by_K=None):
    nd = len(seqs); odd = list(range(0, nd, 2)); even = list(range(1, nd, 2)); res = {}
    for K in Kset:
        hl, fits = [], []
        for fold, (tr, te) in enumerate(((odd, even), (even, odd))):
            tr_s = [seqs[i] for i in tr]; te_s = [seqs[i] for i in te]
            if K == 1:
                model = k1_model(tr_s); lls = []
            else:
                ini = inits_by_K.get(K) if inits_by_K else None
                model, ll, lls, it = _fit(tr_s, K, fly_idx, tag, fold, restarts, ini)
            nobs = sum((y >= 0).sum() for y, _ in te_s)
            hl.append(hc.loglik(model, te_s) / nobs); fits.append((model, lls))
        res[K] = dict(heldout=float(np.mean(hl)), fits=fits)
    return res


def type_of(B, G, lab):
    p0 = B[0]
    band = 'QUIESCENT' if p0 >= 0.999 else ('LOW' if p0 >= 0.95 else 'ACTIVE')
    move = '-' if band == 'QUIESCENT' else ['1-beam', '2-beam', '>=3-beam'][int(np.argmax(B[1:]))]
    share = np.array([G[lab == L].sum() for L in (pl.S, pl.Q, pl.W)]); share = share / max(share.sum(), 1e-300)
    return (band, move, ['S', 'Q', 'W'][int(np.argmax(share))]), share


def describe(model, seqs, lab, zone=None):
    pi, P, B = model; K = B.shape[0]
    Gs = hc.posterior(model, seqs); Gc = np.concatenate(Gs); del Gs
    types, shares = [], []
    for k in range(K):
        tp, sh = type_of(B[k], Gc[:, k], lab[:len(Gc)]); types.append(tp); shares.append(sh)
    dwell = [float(np.mean([0.4 / max(1 - P[p, k, k], 1e-12) for p in range(NP)])) for k in range(K)]
    place = None
    if zone is not None:
        z = zone[:len(Gc)]; place = []
        for k in range(K):
            tot = max(Gc[:, k].sum(), 1e-300)
            place.append([float(Gc[z == zz, k].sum() / tot) for zz in (0, 1, 2)])
    return dict(types=types, shares=shares, dwell=dwell, argmax=Gc.argmax(1).astype(np.int8), place=place)


def zones(x):
    x = np.asarray(x)
    return np.where(x == 0, 0, np.where((x <= 2) | (x >= 14), 1, 2))


def run_fly(x, t, valid, fly_idx, n_spike_fixed=None):
    y = step_symbols(x, valid); n = len(y) // SPD * SPD; y = y[:n]
    seqs = to_days(y); lab = labels_per_sample(x, t, valid, n); zone = zones(x)[:n]
    nd = len(seqs); odd = list(range(0, nd, 2)); even = list(range(1, nd, 2))
    out = dict(n_days=nd, n_spike_fixed=n_spike_fixed)
    cv = curve(seqs, fly_idx, KS, R_OBS)
    out['heldout'] = {K: cv[K]['heldout'] for K in KS}
    out['restart_lls'] = {K: [list(map(float, f[1])) for f in cv[K]['fits']] for K in KS}
    deg = {K: any(d1.degenerate(f[0][2]) for f in cv[K]['fits']) for K in KS}; out['degenerate'] = deg
    elig = [K for K in KS if not deg[K]]
    Kstar = max(elig, key=lambda K: cv[K]['heldout']) if elig else None; out['Kstar'] = Kstar
    if Kstar is None:
        return out
    full, llf, _, _ = _fit(seqs, Kstar, fly_idx, 0, 9)
    nobs = int((y >= 0).sum()); kpar = (Kstar - 1) + NP * Kstar * (Kstar - 1) + Kstar * (M - 1)
    out['bic_Kstar'] = float(-2 * llf + kpar * np.log(nobs))
    desc = describe(full, seqs, lab, zone)
    out['types'] = [list(tp) for tp in desc['types']]; out['shares'] = [list(map(float, s)) for s in desc['shares']]
    out['dwell_s'] = desc['dwell']; out['place'] = desc['place']
    cand = d1.candidates(desc); out['candidates'] = cand
    # C1
    Kn = Kstar - 1
    prev = cv[Kn] if Kn in cv else curve(seqs, fly_idx, [Kn], R_OBS)[Kn]
    gain_obs = cv[Kstar]['heldout'] - prev['heldout']
    null = _fit(seqs, Kn, fly_idx, 0, 9)[0] if Kn >= 2 else k1_model(seqs)
    rng = np.random.default_rng(np.random.SeedSequence([2001, fly_idx])); gains = []
    for b in range(NB):
        sim = hc.simulate_fast(null, [(len(yy), ph) for yy, ph in seqs], [yy < 0 for yy, _ in seqs], rng)
        ini = {Kstar: [full]}
        if Kn >= 2:
            ini[Kn] = [null]
        cb = curve(sim, fly_idx * 1000 + b, [Kn, Kstar], 2, tag=1, inits_by_K=ini)
        gains.append(cb[Kstar]['heldout'] - cb[Kn]['heldout'])
    thr = float(np.percentile(gains, 95))
    out['C1'] = dict(gain_obs=float(gain_obs), null_gains=list(map(float, gains)), thr95=thr, pass_=bool(gain_obs > thr))
    # C3 halves = fold fits
    halves = [cv[Kstar]['fits'][0][0], cv[Kstar]['fits'][1][0]]
    lab_odd = np.concatenate([lab[i*SPD:(i+1)*SPD] for i in odd]); lab_even = np.concatenate([lab[i*SPD:(i+1)*SPD] for i in even])
    hdesc = [describe(halves[0], [seqs[i] for i in odd], lab_odd), describe(halves[1], [seqs[i] for i in even], lab_even)]
    # C4(a) flicker removed
    yf = step_symbols(remove_flicker(x), valid)[:n]; sf = to_days(yf)
    mf, _, _, _ = _fit(sf, Kstar, fly_idx, 2, 9)
    amf = np.concatenate(hc.posterior(mf, sf)).argmax(1)
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
        win = desc['argmax'] == k
        best = 0.0
        for j in range(Kstar):
            wj = amf == j; uni = (win | wj).sum(); best = max(best, (win & wj).sum() / uni if uni else 0.0)
        b1 = desc['place'][k][0]
        per[k] = dict(type=list(tp), candidate=k in cand, dwell_s=desc['dwell'][k], C2=bool(desc['dwell'][k] >= 30),
                      C3=c3, C4a_jaccard=float(best), beyond1_share=b1, C4_pass=bool(best >= 0.5 and b1 < 0.5))
    out['states'] = per
    out['argmax_full'] = desc['argmax']
    return out
