"""T1-B: cohort-level clock test on sleep entry and exit (PRESPEC T1-B).

Per fly f: sample-level (0.4-s) HMM on step symbols (D1-D), K_f states, sleep state s_f.
Transitions vary with the 10-s bin of the day w (0..8639; ZT = 10 w s):
  out of sleep:   P(leave | w) = sigmoid(a_{f,ph(w)} + Z_w . h);  destinations split by d_{f,ph}(j)
  wake state i:   P(i -> sleep | w) = sigmoid(b_{f,i,ph(w)} + Z_w . g);  rest split by r_{f,i,ph}(j)
h, g shared by the cohort. C0: h = g = 0 (identical to the per-fly D1-D light/dark HMM). C1: 1st harmonic.
C2: 1st + 2nd harmonics.
"""
import numpy as np
from scipy.optimize import minimize
import hmm_cat as hc, d1d

W = 8640; SPB = 25          # bins per day; samples per bin
NHARM = {0: 0, 1: 1, 2: 2}


def harmonics(model_id):
    zt = np.arange(W) * 10.0 / 3600.0; cols = []
    for k in range(1, NHARM[model_id] + 1):
        cols += [np.cos(2 * np.pi * k * zt / 24), np.sin(2 * np.pi * k * zt / 24)]
    return np.column_stack(cols) if cols else np.zeros((W, 0))


PHASE = (np.arange(W) >= W // 2).astype(int)
_bin_of_sample = (np.arange(d1d.SPD) // SPB).astype(np.int64)


def bin_seqs(seqs):
    """Replace the per-sample phase index by the 10-s bin index (sequences are whole fly-days)."""
    return [(y, _bin_of_sample[:len(y)]) for y, _ in seqs]


def sig(x):
    return 1.0 / (1.0 + np.exp(-x))


def logit(p):
    p = np.clip(p, 1e-12, 1 - 1e-12); return np.log(p / (1 - p))


# ------------------------------------------------------------------ per-fly parameter blocks
def from_c0(model, s):
    """Convert a per-phase HMM (pi, P[2,K,K], B) into T1-B parameters for sleep state s."""
    pi, P, B = model; K = B.shape[0]; wake = [i for i in range(K) if i != s]
    a = np.array([logit(1 - P[p, s, s]) for p in (0, 1)])
    d = np.array([[P[p, s, j] / max(1 - P[p, s, s], 1e-300) if j != s else 0.0 for j in range(K)] for p in (0, 1)])
    b = np.array([[logit(P[p, i, s]) for p in (0, 1)] for i in wake])                    # nwake x 2
    r = np.array([[[P[p, i, j] / max(1 - P[p, i, s], 1e-300) if j != s else 0.0 for j in range(K)]
                   for p in (0, 1)] for i in wake])                                     # nwake x 2 x K
    return dict(pi=pi.copy(), B=B.copy(), s=s, wake=wake, a=a, d=d, b=b, r=r)


def build_P(F, h, g, Z):
    """Transition matrices per bin (W x K x K) for one fly."""
    K = F['B'].shape[0]; s = F['s']; P = np.zeros((W, K, K))
    hz = Z @ h if len(h) else np.zeros(W); gz = Z @ g if len(g) else np.zeros(W)
    e = sig(F['a'][PHASE] + hz)
    P[:, s, s] = 1 - e
    P[:, s, :] += e[:, None] * F['d'][PHASE]
    for n, i in enumerate(F['wake']):
        q = sig(F['b'][n][PHASE] + gz)
        P[:, i, :] = (1 - q)[:, None] * F['r'][n][PHASE]
        P[:, i, s] = q
    return P


# ------------------------------------------------------------------ joint logistic M-step
def _fit_logistic(groups, Z, theta0, nint):
    """groups: list of (succ[W], fail[W], intercept index per bin [W]); shared slope vector of length Z.shape[1].
    Parameters = nint intercepts + shared slopes. Minimises the binomial negative log-likelihood."""
    p = Z.shape[1]

    def f(theta):
        c, beta = theta[:nint], theta[nint:]
        zb = Z @ beta if p else 0.0
        obj = 0.0; gc = np.zeros(nint); gb = np.zeros(p)
        for S, Fl, idx in groups:
            eta = c[idx] + zb
            sp = np.logaddexp(0, eta)
            obj += ((S + Fl) * sp - S * eta).sum()
            resid = (S + Fl) * sig(eta) - S
            gc += np.bincount(idx, weights=resid, minlength=nint)
            if p:
                gb += Z.T @ resid
        return obj, np.r_[gc, gb]
    res = minimize(f, theta0, jac=True, method='L-BFGS-B', options=dict(maxiter=500))
    return res.x


def fit_cohort(flies, model_id, starts, max_iter=300, tol=0.01):
    """flies: list of dict(F=initial params (from_c0), seqs=bin-indexed sequences).
    starts: list of (h0, g0) initial shared vectors. Returns best dict(ll, h, g, F list)."""
    Z = harmonics(model_id); p = Z.shape[1]
    best = None
    for h0, g0 in starts:
        Fs = [{k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in fl['F'].items()} for fl in flies]
        h, g = np.array(h0, float), np.array(g0, float)
        prev = -np.inf
        for it in range(max_iter):
            ll = 0.0; ex_groups = []; en_groups = []; ex_c0 = []; en_c0 = []
            for f_i, (Fd, fl) in enumerate(zip(Fs, flies)):
                K = Fd['B'].shape[0]; s = Fd['s']
                P = build_P(Fd, h, g, Z)
                Xi = np.zeros((W, K, K)); Eb = np.zeros((K, Fd['B'].shape[1])); G0 = np.zeros(K)
                for y, w in fl['seqs']:
                    ll += hc._estep(Fd['pi'], P, Fd['B'], y, w, Xi, Eb, G0)
                Fd['pi'] = G0 / G0.sum()
                Bn = np.maximum(Eb / np.maximum(Eb.sum(1, keepdims=True), 1e-300), 1e-12); Fd['B'] = Bn / Bn.sum(1, keepdims=True)
                for ph in (0, 1):
                    m = PHASE == ph
                    row = Xi[m, s, :].sum(0); row[s] = 0; Fd['d'][ph] = row / max(row.sum(), 1e-300)
                    for n, i in enumerate(Fd['wake']):
                        row = Xi[m, i, :].sum(0); row[s] = 0; Fd['r'][n][ph] = row / max(row.sum(), 1e-300)
                off = len(ex_c0)
                ex_groups.append((Xi[:, s, :].sum(1) - Xi[:, s, s], Xi[:, s, s], off + PHASE)); ex_c0 += list(Fd['a'])
                for n, i in enumerate(Fd['wake']):
                    off = len(en_c0)
                    en_groups.append((Xi[:, i, s], Xi[:, i, :].sum(1) - Xi[:, i, s], off + PHASE)); en_c0 += list(Fd['b'][n])
            if ll - prev < tol:
                break
            prev = ll
            th = _fit_logistic(ex_groups, Z, np.r_[ex_c0, h], len(ex_c0)); a_all, h = th[:len(ex_c0)], th[len(ex_c0):]
            th = _fit_logistic(en_groups, Z, np.r_[en_c0, g], len(en_c0)); b_all, g = th[:len(en_c0)], th[len(en_c0):]
            ka = kb = 0
            for Fd in Fs:
                Fd['a'] = a_all[ka:ka + 2]; ka += 2
                for n in range(len(Fd['wake'])):
                    Fd['b'][n] = b_all[kb:kb + 2]; kb += 2
        # final log-likelihood at the returned parameters
        llf = sum(hc.loglik((Fd['pi'], build_P(Fd, h, g, Z), Fd['B']), fl['seqs']) for Fd, fl in zip(Fs, flies))
        if best is None or llf > best['ll']:
            best = dict(ll=float(llf), h=h.copy(), g=g.copy(), F=Fs, iters=it + 1)
    return best


def sleep_state(model, seqs_ph, lab):
    """QUIESCENT/S state with the longest dwell; fallback: the state with the highest P(no change)."""
    pi, P, B = model; K = B.shape[0]
    G = np.concatenate(hc.posterior(model, seqs_ph)); lab = lab[:len(G)]
    cands = []
    for k in range(K):
        sh = np.array([G[lab == L, k].sum() for L in (2, 0, 1)])           # S, Q, W (pipeline codes)
        if B[k, 0] >= 0.999 and sh.argmax() == 0:
            cands.append((np.mean([1 / max(1 - P[p, k, k], 1e-12) for p in (0, 1)]), k))
    if cands:
        return max(cands)[1], True
    return int(np.argmax(B[:, 0])), False


def harmonic_profile(vec, model_id):
    return harmonics(model_id) @ vec if len(vec) else np.zeros(W)
