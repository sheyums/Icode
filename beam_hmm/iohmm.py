"""Input-output HMM with clock-dependent transitions (PRESPEC T1).

Transitions vary with window-of-day w (0..8639, ZT = 10 w s) through a multinomial logit
P(j | i, z_w) = softmax_j(b_ij . z_w), b_ii = 0. Emissions are categorical and time-invariant.
Reuses hmm_cat's numba forward-backward by passing P with one matrix per window-of-day and ph = w.
"""
import numpy as np
from scipy.optimize import minimize
import hmm_cat as hc

WPD = 8640


def covariates(model_id):
    zt = np.arange(WPD) * 10.0 / 3600.0
    cols = [np.ones(WPD), (zt >= 12).astype(float)]
    for h in range({0: 0, 1: 1, 2: 2, 3: 4}[model_id]):
        k = h + 1
        cols += [np.cos(2 * np.pi * k * zt / 24), np.sin(2 * np.pi * k * zt / 24)]
    Z = np.column_stack(cols)
    pen = np.r_[0.0, 0.0, np.ones(Z.shape[1] - 2)]          # no ridge on intercept and dark step
    return Z, pen


def trans_matrices(beta, Z):
    """beta: K x K x p (beta[i,i] = 0). Returns P: W x K x K."""
    E = np.einsum('ijp,wp->wij', beta, Z)
    E -= E.max(2, keepdims=True)
    P = np.exp(E); P /= P.sum(2, keepdims=True)
    return P


def _mstep_beta(beta, Xi, Z, pen, lam):
    """Weighted multinomial logit per origin state, L-BFGS, warm-started."""
    K = beta.shape[0]; p = Z.shape[1]; out = beta.copy()
    for i in range(K):
        others = [j for j in range(K) if j != i]
        Ci = Xi[:, i, :]                                         # W x K expected counts from i
        if Ci.sum() <= 0 or K == 1:
            continue
        x0 = beta[i, others].ravel()

        def f(v):
            b = np.zeros((K, p)); b[others] = v.reshape(K - 1, p)
            E = Z @ b.T; m = E.max(1, keepdims=True); lse = m[:, 0] + np.log(np.exp(E - m).sum(1))
            ll = (Ci * E).sum() - (Ci.sum(1) * lse).sum()
            Pw = np.exp(E - lse[:, None])
            G = Ci - Ci.sum(1, keepdims=True) * Pw               # W x K
            g = (G.T @ Z)[others]                                # (K-1) x p
            bo = b[others]
            obj = -ll + 0.5 * lam * (pen * bo ** 2).sum()
            grad = -g + lam * pen * bo
            return obj, grad.ravel()
        r = minimize(f, x0, jac=True, method='L-BFGS-B', options=dict(maxiter=200))
        out[i, others] = r.x.reshape(K - 1, p)
    return out


def penalized(beta, pen, lam):
    return 0.5 * lam * float((pen * beta ** 2).sum())


def fit(seqs, K, M, model_id, seed_seq, n_restarts=10, max_iter=500, tol=1e-4, lam=0.01, init=None, floor=1e-12):
    """seqs: list of (y, w) with w = window-of-day index. init: optional (pi, beta, B) to warm-start every
    restart (perturbed). Returns (model_tuple(pi, P, B), beta, best penalized ll, all lls)."""
    Z, pen = covariates(model_id); p = Z.shape[1]
    ys = np.concatenate([y for y, _ in seqs]); ys = ys[ys >= 0]
    freq = np.bincount(ys, minlength=M).astype(float); freq = (freq + 1e-6) / (freq + 1e-6).sum()
    children = seed_seq.spawn(n_restarts)
    best, best_beta, best_ll, lls = None, None, -np.inf, []
    for r in range(n_restarts):
        rng = np.random.default_rng(children[r])
        if init is None:
            pi0, P0, B = hc.init_model(K, M, 2, freq, rng)
            beta = np.zeros((K, K, p))
            for i in range(K):
                for j in range(K):
                    if i != j:
                        beta[i, j, 0] = np.log(P0[0, i, j] / P0[0, i, i])
                        beta[i, j, 1] = np.log(P0[1, i, j] / P0[1, i, i]) - beta[i, j, 0]
            pi = pi0
        else:
            pi, beta, B = init[0].copy(), init[1].copy(), init[2].copy()
            if beta.shape[2] < p:
                beta = np.concatenate([beta, np.zeros((K, K, p - beta.shape[2]))], 2)
            beta[:, :, 2:] += rng.normal(0, 0.1, size=beta[:, :, 2:].shape)
            for i in range(K):
                beta[i, i] = 0
        prev = -np.inf
        for it in range(max_iter):
            P = trans_matrices(beta, Z)
            Xi = np.zeros((WPD, K, K)); Eb = np.zeros((K, M)); G0 = np.zeros(K); ll = 0.0
            for y, w in seqs:
                ll += hc._estep(pi, P, B, y, w, Xi, Eb, G0)
            pll = ll - penalized(beta, pen, lam)
            pi = G0 / G0.sum()
            B = np.maximum(Eb / np.maximum(Eb.sum(1, keepdims=True), 1e-300), floor); B /= B.sum(1, keepdims=True)
            beta = _mstep_beta(beta, Xi, Z, pen, lam)
            if pll - prev < tol:
                break
            prev = pll
        P = trans_matrices(beta, Z)
        ll = hc.loglik((pi, P, B), seqs) - penalized(beta, pen, lam)
        lls.append(ll)
        if ll > best_ll:
            best, best_beta, best_ll = (pi, P, B), beta.copy(), ll
    return best, best_beta, best_ll, lls


def to_days_w(y):
    nd = len(y) // WPD; w = np.arange(WPD, dtype=np.int64)
    return [(y[d*WPD:(d+1)*WPD].copy(), w.copy()) for d in range(nd)]
