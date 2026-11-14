"""T1-C: cohort clock test on OBSERVED sleep onsets and wake-ups (PRESPEC T1-C). Poisson log-linear model with
per-fly x phase intercepts and shared harmonics; null = rotation of each fly's bins within each 12-h phase."""
import numpy as np
from scipy.optimize import minimize
import pipeline as pl

NB = 144; HALF = 72


def harmonics(model_id):
    zt = (np.arange(NB) + 0.5) / 6.0; cols = []
    for k in range(1, {0: 0, 1: 1, 2: 2}[model_id] + 1):
        cols += [np.cos(2 * np.pi * k * zt / 24), np.sin(2 * np.pi * k * zt / 24)]
    return np.column_stack(cols) if cols else np.zeros((NB, 0))


PH = (np.arange(NB) >= HALF).astype(int)


def fly_bins(a, lights_on_h):
    """a: 1-s labels from the scoring (1 wake, 0 sleep, -1 excluded). Returns per-bin sums over all days:
    wake-ups, onsets, seconds asleep, seconds awake (arrays of 144)."""
    a = np.asarray(a); s = np.arange(len(a))
    b = (((s / 3600.0 - lights_on_h) % 24) * 6).astype(int) % NB
    ok = a >= 0; okp = np.r_[False, ok[1:] & ok[:-1]]
    prev = np.r_[-1, a[:-1]]
    wu = okp & (prev == 0) & (a == 1); on = okp & (prev == 1) & (a == 0)
    f = lambda m: np.bincount(b[m], minlength=NB).astype(float)
    return dict(wu=f(wu), on=f(on), es=f(ok & (a == 0)), ew=f(ok & (a == 1)))


def labels_from_activity(act):
    """P1 scoring: 1-s activity (1/0/-1) -> labels 1 wake / 0 sleep / -1 excluded."""
    lab = np.where(np.asarray(act) == -1, -1, 1).astype(np.int8)
    for s0, n0, _ in pl.sleep_bouts(act):
        lab[s0:s0 + n0] = 0
    return lab


def keep_bins(excl_h=0.0):
    """T1-C1: indices of the 10-min bins kept when the first excl_h hours after each light switch are dropped."""
    k = int(round(excl_h * 6)); b = np.arange(NB)
    return np.nonzero(~((b < k) | ((b >= HALF) & (b < HALF + k))))[0]


def _nll(theta, Y, E, Z, nf, ph=None):
    """Poisson negative log-likelihood (without the y! term) for one event type; Y, E: nf x bins."""
    ph = PH if ph is None else ph
    c = theta[:2 * nf].reshape(nf, 2); h = theta[2 * nf:]
    hz = Z @ h if len(h) else np.zeros(Y.shape[1])
    eta = c[:, ph] + hz[None, :]
    m = E > 0
    mu = np.exp(eta) * E
    nll = -(Y[m] * eta[m]).sum() + mu[m].sum()
    r = np.where(m, mu - Y, 0.0)
    gc = np.stack([r[:, ph == 0].sum(1), r[:, ph == 1].sum(1)], 1).ravel()
    gh = Z.T @ r.sum(0) if len(h) else np.zeros(0)
    return nll, np.r_[gc, gh]


def fit(Y, E, model_id, keep=None):
    keep = np.arange(NB) if keep is None else keep
    Y = Y[:, keep]; E = E[:, keep]; ph = PH[keep]
    nf = Y.shape[0]; Z = harmonics(model_id)[keep]
    c0 = np.zeros((nf, 2))
    for p_ in (0, 1):
        ys = Y[:, ph == p_].sum(1); es = E[:, ph == p_].sum(1)
        c0[:, p_] = np.log(np.maximum(ys, 0.5) / np.maximum(es, 1.0))
    th0 = np.r_[c0.ravel(), np.zeros(Z.shape[1])]
    r = minimize(_nll, th0, args=(Y, E, Z, nf, ph), jac=True, method='L-BFGS-B', options=dict(maxiter=1000))
    return float(-r.fun), r.x[2 * nf:]


def _profile(h, Ysum, ZY, E, Zk, ph):
    """Profile negative log-likelihood over the per-fly x phase intercepts (closed form), up to a constant:
    f(h) = sum_{f,p} Y_fp log S_fp(h) - h . sum_{f,b} y_fb z_b,  S_fp(h) = sum_{b in p} E_fb exp(z_b . h).
    Exactly equivalent to the full fit; LR = 2 (f(0) - min_h f(h))."""
    w = np.exp(Zk @ h)                                   # bins
    EW = E * w[None, :]
    S = np.stack([EW[:, ph == 0].sum(1), EW[:, ph == 1].sum(1)], 1)        # flies x 2
    m = Ysum > 0
    f = (Ysum[m] * np.log(S[m])).sum() - h @ ZY
    G = np.zeros_like(h)
    for p_ in (0, 1):
        a = Ysum[:, p_] / np.where(S[:, p_] > 0, S[:, p_], 1.0)          # Y_fp / S_fp
        G += ((a[:, None] * EW[:, ph == p_]).sum(0)) @ Zk[ph == p_]
    return f, G - ZY


def lr_profile(Y, E, model_id=2, keep=None):
    keep = np.arange(NB) if keep is None else keep
    Y = Y[:, keep]; E = E[:, keep]; ph = PH[keep]; Zk = harmonics(model_id)[keep]
    Ysum = np.stack([Y[:, ph == 0].sum(1), Y[:, ph == 1].sum(1)], 1); ZY = Zk.T @ Y.sum(0)
    f0, _ = _profile(np.zeros(Zk.shape[1]), Ysum, ZY, E, Zk, ph)
    r = minimize(_profile, np.zeros(Zk.shape[1]), args=(Ysum, ZY, E, Zk, ph), jac=True, method='BFGS')
    return 2 * (f0 - r.fun), r.x


def lr_stat_fast(B, model_id=2, keep=None):
    Ywu = np.stack([b['wu'] for b in B]); Ees = np.stack([b['es'] for b in B])
    Yon = np.stack([b['on'] for b in B]); Eew = np.stack([b['ew'] for b in B])
    la, h = lr_profile(Ywu, Ees, model_id, keep); lb, g = lr_profile(Yon, Eew, model_id, keep)
    return la + lb, h, g


def lr_stat(B, model_id=2, l0=None, keep=None):
    """B: list of fly_bins dicts. Returns LR (both event types), h (wake-ups), g (onsets), and the C0 logLs.
    The C0 fit has only per-fly x phase intercepts, so it is unchanged by a within-phase rotation: pass l0 to reuse it."""
    Ywu = np.stack([b['wu'] for b in B]); Ees = np.stack([b['es'] for b in B])
    Yon = np.stack([b['on'] for b in B]); Eew = np.stack([b['ew'] for b in B])
    if l0 is None:
        l0 = (fit(Ywu, Ees, 0, keep)[0], fit(Yon, Eew, 0, keep)[0])
    l2a, h = fit(Ywu, Ees, model_id, keep); l2b, g = fit(Yon, Eew, model_id, keep)
    return 2 * ((l2a - l0[0]) + (l2b - l0[1])), h, g, l0


def rotate(b, oL, oD, keep=None):
    """Rotate within each phase; with keep (T1-C1), only among the kept bins of each phase."""
    if keep is None:
        idx = np.r_[(np.arange(HALF) + oL) % HALF, HALF + (np.arange(HALF) + oD) % HALF]
        return {k: v[idx] for k, v in b.items()}
    L = keep[keep < HALF]; D = keep[keep >= HALF]; idx = np.arange(NB)
    idx[L] = L[(np.arange(len(L)) + oL) % len(L)]; idx[D] = D[(np.arange(len(D)) + oD) % len(D)]
    return {k: v[idx] for k, v in b.items()}


def rotation_test(B, seed_seq, n=999, model_id=2, excl_h=0.0):
    keep = keep_bins(excl_h) if excl_h > 0 else None
    nL = HALF if keep is None else int((keep < HALF).sum()); nD = HALF if keep is None else int((keep >= HALF).sum())
    LR, h, g = lr_stat_fast(B, model_id, keep=keep)       # profile likelihood: exact, much faster (T1-C1 note)
    rng = np.random.default_rng(seed_seq); null = []
    for _ in range(n):
        Br = [rotate(b, int(rng.integers(0, nL)), int(rng.integers(0, nD)), keep) for b in B]
        null.append(lr_stat_fast(Br, model_id, keep)[0])
    p = (1 + sum(x >= LR for x in null)) / (n + 1)
    return dict(LR=float(LR), p=float(p), h=h.tolist(), g=g.tolist(), null_median=float(np.median(null)))


def within_phase_peak(vec, model_id=2):
    prof = harmonics(model_id) @ np.asarray(vec)
    zt = (np.arange(NB) + 0.5) / 6.0
    return dict(light_peak_zt=float(zt[:HALF][np.argmax(prof[:HALF])]), dark_peak_zt=float(zt[HALF:][np.argmax(prof[HALF:])]),
                dark_range=float(prof[HALF:].max() - prof[HALF:].min()))
