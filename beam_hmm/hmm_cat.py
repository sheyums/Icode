"""Categorical HMM with phase-specific transitions and missing observations (D1, D1-A).

y: int array, symbol 0..M-1, or -1 = missing (emission likelihood 1).
ph: int array, phase of each step; the step INTO t uses P[ph[t]] (sqw_ref convention).
Multiple sequences (fly-days) each start from pi.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def _fb(pi, P, B, y, ph):
    T = y.shape[0]; K = pi.shape[0]
    A = np.empty((T, K)); c = np.empty(T); Bk = np.empty((T, K))
    for k in range(K):
        A[0, k] = pi[k] * (B[k, y[0]] if y[0] >= 0 else 1.0)
    c[0] = A[0].sum(); A[0] /= c[0]
    for t in range(1, T):
        Pt = P[ph[t]]
        for j in range(K):
            s = 0.0
            for i in range(K):
                s += A[t-1, i] * Pt[i, j]
            A[t, j] = s * (B[j, y[t]] if y[t] >= 0 else 1.0)
        c[t] = A[t].sum()
        for j in range(K):
            A[t, j] /= c[t]
    for k in range(K):
        Bk[T-1, k] = 1.0
    for t in range(T-2, -1, -1):
        Pn = P[ph[t+1]]
        for i in range(K):
            s = 0.0
            for j in range(K):
                s += Pn[i, j] * (B[j, y[t+1]] if y[t+1] >= 0 else 1.0) * Bk[t+1, j]
            Bk[t, i] = s / c[t+1]
    ll = 0.0
    for t in range(T):
        ll += np.log(c[t])
    return A, Bk, c, ll


@njit(cache=True)
def _estep(pi, P, B, y, ph, Xi, Eb, G0):
    A, Bk, c, ll = _fb(pi, P, B, y, ph)
    T = y.shape[0]; K = pi.shape[0]
    for k in range(K):
        G0[k] += A[0, k] * Bk[0, k]
    for t in range(T):
        if y[t] >= 0:
            for k in range(K):
                Eb[k, y[t]] += A[t, k] * Bk[t, k]
    for t in range(1, T):
        Pt = P[ph[t]]; p = ph[t]
        for i in range(K):
            for j in range(K):
                Xi[p, i, j] += A[t-1, i] * Pt[i, j] * (B[j, y[t]] if y[t] >= 0 else 1.0) * Bk[t, j] / c[t]
    return ll


def loglik(model, seqs):
    pi, P, B = model
    return float(sum(_fb(pi, P, B, y, ph)[3] for y, ph in seqs))


def posterior(model, seqs):
    pi, P, B = model; out = []
    for y, ph in seqs:
        A, Bk, c, ll = _fb(pi, P, B, y, ph)
        G = A * Bk; G /= G.sum(1, keepdims=True); out.append(G)
    return out


def init_model(K, M, nP, freq, rng):
    B = 0.5 * rng.dirichlet(np.ones(M), size=K) + 0.5 * freq[None, :]
    P = np.empty((nP, K, K))
    for p in range(nP):
        for i in range(K):
            if K == 1:
                P[p, i, 0] = 1.0; continue
            off = 0.1 * rng.dirichlet(np.ones(K - 1))
            P[p, i] = np.insert(off, i, 0.9)
    return np.full(K, 1.0 / K), P, B


def fit(seqs, K, M, nP=2, seed_seq=None, n_restarts=10, max_iter=500, tol=1e-4, floor=1e-12):
    """Baum-Welch from n_restarts random starts. Returns (best_model, best_ll, all_lls, n_iters_best)."""
    ys = np.concatenate([y for y, _ in seqs]); ys = ys[ys >= 0]
    freq = np.bincount(ys, minlength=M).astype(float); freq = (freq + 1e-6) / (freq + 1e-6).sum()
    children = seed_seq.spawn(n_restarts)
    best, best_ll, lls, best_it = None, -np.inf, [], 0
    for r in range(n_restarts):
        rng = np.random.default_rng(children[r])
        pi, P, B = init_model(K, M, nP, freq, rng)
        prev = -np.inf
        for it in range(max_iter):
            Xi = np.zeros((nP, K, K)); Eb = np.zeros((K, M)); G0 = np.zeros(K); ll = 0.0
            for y, ph in seqs:
                ll += _estep(pi, P, B, y, ph, Xi, Eb, G0)
            pi = G0 / G0.sum()
            for p in range(nP):
                rs = Xi[p].sum(1, keepdims=True)
                P[p] = np.where(rs > 0, Xi[p] / np.where(rs > 0, rs, 1), P[p])
                P[p] = np.maximum(P[p], floor); P[p] /= P[p].sum(1, keepdims=True)
            B = np.maximum(Eb / np.maximum(Eb.sum(1, keepdims=True), 1e-300), floor); B /= B.sum(1, keepdims=True)
            if ll - prev < tol:
                break
            prev = ll
        ll = loglik((pi, P, B), seqs)
        lls.append(ll)
        if ll > best_ll:
            best, best_ll, best_it = (pi, P.copy(), B.copy()), ll, it + 1
    return best, best_ll, lls, best_it


def simulate(model, lengths_phases, missing_masks, rng):
    """Window-level simulation from a fitted model, same lengths, phases and missing pattern."""
    pi, P, B = model; K, M = B.shape; out = []
    for (T, ph), miss in zip(lengths_phases, missing_masks):
        s = np.empty(T, int); y = np.empty(T, np.int64)
        s[0] = rng.choice(K, p=pi)
        for t in range(1, T):
            s[t] = rng.choice(K, p=P[ph[t], s[t-1]])
        u = rng.random(T); cum = np.cumsum(B, 1)
        for t in range(T):
            y[t] = min(np.searchsorted(cum[s[t]], u[t]), M - 1)
        y[miss] = -1
        out.append((y, ph))
    return out


def js_bits(p, q):
    m = 0.5 * (p + q)
    def kl(a, b):
        a = np.maximum(a, 1e-300); return float((a * np.log2(a / np.maximum(b, 1e-300))).sum())
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)
