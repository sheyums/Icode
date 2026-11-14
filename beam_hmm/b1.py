"""B1 (= D1-E): direction-aware autoregressive HMM at 0.4 s (PRESPEC B1). Truth-blind per-fly procedure.

Symbols: 0 STAY, 1 SAME, 2 REV, 3 BIG (relative to the direction of the last non-zero step).
Emission P(y_t | k, y_{t-1}) = B[k, y_{t-1}, y_t]; implemented through composite symbols c_t = 4 y_{t-1} + y_t
(y_{-1} := STAY at the start of a fly-day and after a missing sample), so hmm_cat's compiled forward-backward
is reused; only the M-step differs (normalisation within each (k, prev) row).
"""
import numpy as np
from numba import njit
import hmm_cat as hc, d1d, d1

STAY, SAME, REV, BIG = 0, 1, 2, 3
M = 4; MC = 16; NP = 2; KS = list(range(2, 9)); SPD = d1d.SPD
NB = 19; TOL = 0.01; MAXIT = 300; R_OBS = 5


@njit(cache=True)
def _symbols(x, valid):
    n = x.shape[0]; y = np.empty(n, np.int64); last = 0
    y[0] = STAY if valid[0] else -1
    for t in range(1, n):
        if not (valid[t] and valid[t-1]):
            y[t] = -1; last = 0; continue
        d = x[t] - x[t-1]
        if d == 0:
            y[t] = STAY
        elif abs(d) >= 2:
            y[t] = BIG; last = 1 if d > 0 else -1
        else:
            s = 1 if d > 0 else -1
            y[t] = SAME if (last == 0 or s == last) else REV
            last = s
    return y


def symbols(x, valid=None):
    x = np.asarray(x).astype(np.int64)
    v = np.ones(len(x), np.bool_) if valid is None else np.asarray(valid, np.bool_)
    return _symbols(x, v)


def composite(y):
    """c_t = 4 y_{t-1} + y_t per fly-day sequence (prev := STAY at index 0 and after a missing sample)."""
    prev = np.r_[STAY, y[:-1]]; prev[prev < 0] = STAY
    c = 4 * prev + y; c[y < 0] = -1
    return c


def to_days(y):
    ph = ((np.arange(SPD) * 0.4) % 86400 >= 43200).astype(np.int64)
    return [(composite(y[d*SPD:(d+1)*SPD]), ph.copy()) for d in range(len(y) // SPD)]


def normalise(Eb, floor=1e-12):
    K = Eb.shape[0]; B3 = Eb.reshape(K, 4, 4) + 1e-300
    B3 = np.maximum(B3 / B3.sum(2, keepdims=True), floor); B3 /= B3.sum(2, keepdims=True)
    return B3.reshape(K, MC)


def fit(seqs, K, seed_seq, n_restarts=R_OBS, max_iter=MAXIT, tol=TOL, inits=None):
    cs = np.concatenate([c for c, _ in seqs]); cs = cs[cs >= 0]
    cf = np.bincount(cs, minlength=MC).astype(float).reshape(4, 4) + 1e-6
    cf = (cf / cf.sum(1, keepdims=True)).ravel()
    children = seed_seq.spawn(n_restarts); best, best_ll, lls = None, -np.inf, []
    for r in range(n_restarts):
        rng = np.random.default_rng(children[r])
        if inits is not None and r < len(inits):
            pi, P, B = (a.copy() for a in inits[r])
        else:
            pi, P, _ = hc.init_model(K, 4, NP, np.full(4, .25), rng)
            B = np.stack([normalise((0.5 * rng.dirichlet(np.ones(4), size=4).ravel() + 0.5 * cf)[None, :])[0] for _ in range(K)])
        prev = -np.inf
        for it in range(max_iter):
            Xi = np.zeros((NP, K, K)); Eb = np.zeros((K, MC)); G0 = np.zeros(K); ll = 0.0
            for c, ph in seqs:
                ll += hc._estep(pi, P, B, c, ph, Xi, Eb, G0)
            pi = G0 / G0.sum()
            for p in range(NP):
                rs = Xi[p].sum(1, keepdims=True)
                P[p] = np.where(rs > 0, Xi[p] / np.where(rs > 0, rs, 1), P[p])
                P[p] = np.maximum(P[p], 1e-12); P[p] /= P[p].sum(1, keepdims=True)
            B = normalise(Eb)
            if ll - prev < tol:
                break
            prev = ll
        ll = hc.loglik((pi, P, B), seqs); lls.append(ll)
        if ll > best_ll:
            best, best_ll = (pi, P.copy(), B.copy()), ll
    return best, best_ll, lls


def k1_model(seqs):
    cs = np.concatenate([c for c, _ in seqs]); cs = cs[cs >= 0]
    cf = np.bincount(cs, minlength=MC).astype(float)
    return (np.ones(1), np.ones((NP, 1, 1)), normalise(cf[None, :]))


@njit(cache=True)
def _sim(pi, P, B, ph, u):
    T = ph.shape[0]; K = pi.shape[0]
    c = np.empty(T, np.int64); s = 0; acc = 0.0
    for j in range(K):
        acc += pi[j]
        if u[0] < acc:
            s = j; break
    prev = 0
    for t in range(T):
        if t > 0:
            acc = 0.0; ns = K - 1
            for j in range(K):
                acc += P[ph[t], s, j]
                if u[2*t] < acc:
                    ns = j; break
            s = ns
        acc = 0.0; y = 3
        for j in range(4):
            acc += B[s, 4 * prev + j]
            if u[2*t+1] < acc:
                y = j; break
        c[t] = 4 * prev + y; prev = y
    return c


def simulate(model, lengths_phases, missing_masks, rng):
    pi, P, B = model; out = []
    for (T, ph), miss in zip(lengths_phases, missing_masks):
        c = _sim(pi, P, B, ph, rng.random(2 * T)); c[miss] = -1; out.append((c, ph))
    return out


# --------------------------------------------------------------------------------- state description
def stationary_emission(Bk):
    """Stationary distribution of the symbol chain B[k] (4 x 4)."""
    T = Bk.reshape(4, 4); w, v = np.linalg.eig(T.T)
    s = np.real(v[:, np.argmin(np.abs(w - 1))]); s = np.abs(s) / np.abs(s).sum()
    return s


def dwell(P, k):
    return np.array([0.4 / max(1 - P[p, k, k], 1e-12) for p in range(NP)])


def degenerate(model):
    pi, P, B = model; K = B.shape[0]
    for i in range(K):
        for j in range(i + 1, K):
            close_em = hc.js_bits(stationary_emission(B[i]), stationary_emission(B[j])) < 0.05
            r = dwell(P, i) / dwell(P, j)
            close_dyn = bool(np.all((r <= 1.5) & (r >= 1 / 1.5)))
            if close_em and close_dyn:
                return True
    return False


def type_of(Bk, G, lab):
    s = stationary_emission(Bk); p0 = s[STAY]
    band = 'QUIESCENT' if p0 >= 0.999 else ('LOW' if p0 >= 0.95 else 'ACTIVE')
    T = Bk.reshape(4, 4)
    move = '-' if band == 'QUIESCENT' else ('STEP' if s[SAME] + s[REV] >= s[BIG] else 'BIG')
    rho = np.mean([T[pv, REV] / max(T[pv, SAME] + T[pv, REV], 1e-300) for pv in (SAME, REV)])
    pattern = '-' if band == 'QUIESCENT' else ('FLICKER-LIKE' if rho >= 0.8 else ('PERSISTENT' if rho <= 0.3 else 'RANDOM'))
    share = np.array([G[lab == L].sum() for L in (2, 0, 1)]); share = share / max(share.sum(), 1e-300)
    return (band, move, pattern, ['S', 'Q', 'W'][int(np.argmax(share))]), share, float(rho)


def describe(model, seqs, lab, zone=None):
    pi, P, B = model; K = B.shape[0]
    Gc = np.concatenate(hc.posterior(model, seqs)); lab = lab[:len(Gc)]
    types, shares, rhos = [], [], []
    for k in range(K):
        tp, sh, rho = type_of(B[k], Gc[:, k], lab); types.append(tp); shares.append(sh); rhos.append(rho)
    place = None
    if zone is not None:
        z = zone[:len(Gc)]
        place = [[float(Gc[z == zz, k].sum() / max(Gc[:, k].sum(), 1e-300)) for zz in (0, 1, 2)] for k in range(K)]
    return dict(types=types, shares=shares, rho=rhos, dwell=[float(dwell(P, k).mean()) for k in range(K)],
                argmax=Gc.argmax(1).astype(np.int8), place=place)


# --------------------------------------------------------------------------------- per-fly procedure
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
                model, ll, lls = fit(tr_s, K, np.random.SeedSequence([1001, fly_idx, K, tag, fold]), restarts, inits=ini)
            nobs = sum((c >= 0).sum() for c, _ in te_s)
            hl.append(hc.loglik(model, te_s) / nobs); fits.append((model, lls))
        res[K] = dict(heldout=float(np.mean(hl)), fits=fits)
    return res


def run_fly(x, t, valid, fly_idx, n_spike_fixed=None):
    y = symbols(x, valid); n = len(y) // SPD * SPD; y = y[:n]
    seqs = to_days(y); lab = d1d.labels_per_sample(x, t, valid, n); zone = d1d.zones(x)[:n]
    nd = len(seqs); odd = list(range(0, nd, 2)); even = list(range(1, nd, 2))
    out = dict(n_days=nd, n_spike_fixed=n_spike_fixed)
    cv = curve(seqs, fly_idx, KS, R_OBS)
    out['heldout'] = {K: cv[K]['heldout'] for K in KS}
    deg = {K: any(degenerate(f[0]) for f in cv[K]['fits']) for K in KS}; out['degenerate'] = deg
    elig = [K for K in KS if not deg[K]]
    Kstar = max(elig, key=lambda K: cv[K]['heldout']) if elig else None; out['Kstar'] = Kstar
    if Kstar is None:
        return out
    full, llf, _ = fit(seqs, Kstar, np.random.SeedSequence([1001, fly_idx, Kstar, 0, 9]))
    desc = describe(full, seqs, lab, zone)
    out['types'] = [list(tp) for tp in desc['types']]; out['shares'] = [list(map(float, s)) for s in desc['shares']]
    out['rho'] = desc['rho']; out['dwell_s'] = desc['dwell']; out['place'] = desc['place']
    cand = d1.candidates(desc); out['candidates'] = cand
    Kn = Kstar - 1
    prev = cv[Kn] if Kn in cv else curve(seqs, fly_idx, [Kn], R_OBS)[Kn]
    gain_obs = cv[Kstar]['heldout'] - prev['heldout']
    null = fit(seqs, Kn, np.random.SeedSequence([1001, fly_idx, Kn, 0, 9]))[0] if Kn >= 2 else k1_model(seqs)
    rng = np.random.default_rng(np.random.SeedSequence([2001, fly_idx])); gains = []
    for b in range(NB):
        sim = simulate(null, [(len(c), ph) for c, ph in seqs], [c < 0 for c, _ in seqs], rng)
        ini = {Kstar: [full]}
        if Kn >= 2:
            ini[Kn] = [null]
        cb = curve(sim, fly_idx * 1000 + b, [Kn, Kstar], 2, tag=1, inits_by_K=ini)
        gains.append(cb[Kstar]['heldout'] - cb[Kn]['heldout'])
    thr = float(np.percentile(gains, 95))
    out['C1'] = dict(gain_obs=float(gain_obs), null_gains=list(map(float, gains)), thr95=thr, pass_=bool(gain_obs > thr))
    halves = [cv[Kstar]['fits'][0][0], cv[Kstar]['fits'][1][0]]
    lab_o = np.concatenate([lab[i*SPD:(i+1)*SPD] for i in odd]); lab_e = np.concatenate([lab[i*SPD:(i+1)*SPD] for i in even])
    hdesc = [describe(halves[0], [seqs[i] for i in odd], lab_o), describe(halves[1], [seqs[i] for i in even], lab_e)]
    yf = symbols(d1d.remove_flicker(x), valid)[:n]; sf = to_days(yf)
    mf, _, _ = fit(sf, Kstar, np.random.SeedSequence([1001, fly_idx, Kstar, 2, 9]))
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
        win = desc['argmax'] == k; best = 0.0
        for j in range(Kstar):
            wj = amf == j; uni = (win | wj).sum(); best = max(best, (win & wj).sum() / uni if uni else 0.0)
        b1share = desc['place'][k][0]
        per[k] = dict(type=list(tp), candidate=k in cand, dwell_s=desc['dwell'][k], C2=bool(desc['dwell'][k] >= 30),
                      C3=c3, C4a_jaccard=float(best), beyond1_share=b1share, C4_pass=bool(best >= 0.5 and b1share < 0.5))
    out['states'] = per
    out['argmax_full'] = desc['argmax']
    return out
