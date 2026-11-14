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


def _nll(theta, Y, E, Z, nf):
    """Poisson negative log-likelihood (without the y! term) for one event type; Y, E: nf x 144."""
    c = theta[:2 * nf].reshape(nf, 2); h = theta[2 * nf:]
    hz = Z @ h if len(h) else np.zeros(NB)
    eta = c[:, PH] + hz[None, :]
    m = E > 0
    mu = np.exp(eta) * E
    nll = -(Y[m] * eta[m]).sum() + mu[m].sum()
    r = np.where(m, mu - Y, 0.0)
    gc = np.stack([r[:, PH == 0].sum(1), r[:, PH == 1].sum(1)], 1).ravel()
    gh = Z.T @ r.sum(0) if len(h) else np.zeros(0)
    return nll, np.r_[gc, gh]


def fit(Y, E, model_id):
    nf = Y.shape[0]; Z = harmonics(model_id)
    c0 = np.zeros((nf, 2))
    for ph in (0, 1):
        ys = Y[:, PH == ph].sum(1); es = E[:, PH == ph].sum(1)
        c0[:, ph] = np.log(np.maximum(ys, 0.5) / np.maximum(es, 1.0))
    th0 = np.r_[c0.ravel(), np.zeros(Z.shape[1])]
    if Z.shape[1] == 0:
        r = minimize(_nll, th0, args=(Y, E, Z, nf), jac=True, method='L-BFGS-B')
        return float(-r.fun), r.x[2 * nf:]
    r = minimize(_nll, th0, args=(Y, E, Z, nf), jac=True, method='L-BFGS-B', options=dict(maxiter=1000))
    return float(-r.fun), r.x[2 * nf:]


def lr_stat(B, model_id=2, l0=None):
    """B: list of fly_bins dicts. Returns LR (both event types), h (wake-ups), g (onsets), and the C0 logLs.
    The C0 fit has only per-fly x phase intercepts, so it is unchanged by a within-phase rotation: pass l0 to reuse it."""
    Ywu = np.stack([b['wu'] for b in B]); Ees = np.stack([b['es'] for b in B])
    Yon = np.stack([b['on'] for b in B]); Eew = np.stack([b['ew'] for b in B])
    if l0 is None:
        l0 = (fit(Ywu, Ees, 0)[0], fit(Yon, Eew, 0)[0])
    l2a, h = fit(Ywu, Ees, model_id); l2b, g = fit(Yon, Eew, model_id)
    return 2 * ((l2a - l0[0]) + (l2b - l0[1])), h, g, l0


def rotate(b, oL, oD):
    idx = np.r_[(np.arange(HALF) + oL) % HALF, HALF + (np.arange(HALF) + oD) % HALF]
    return {k: v[idx] for k, v in b.items()}


def rotation_test(B, seed_seq, n=999, model_id=2):
    LR, h, g, l0 = lr_stat(B, model_id)
    rng = np.random.default_rng(seed_seq); null = []
    for _ in range(n):
        Br = [rotate(b, *rng.integers(0, HALF, 2)) for b in B]
        null.append(lr_stat(Br, model_id, l0)[0])
    p = (1 + sum(x >= LR for x in null)) / (n + 1)
    return dict(LR=float(LR), p=float(p), h=h.tolist(), g=g.tolist(), null_median=float(np.median(null)))


def within_phase_peak(vec, model_id=2):
    prof = harmonics(model_id) @ np.asarray(vec)
    zt = (np.arange(NB) + 0.5) / 6.0
    return dict(light_peak_zt=float(zt[:HALF][np.argmax(prof[:HALF])]), dark_peak_zt=float(zt[HALF:][np.argmax(prof[HALF:])]),
                dark_range=float(prof[HALF:].max() - prof[HALF:].min()))
