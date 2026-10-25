"""
sqw_ref.py -- reference implementation of the S/Q/W hidden-Markov analysis (PRESPEC 19, 20, 20b-e).

One function per step, no file paths, no global state. This module is the SPECIFICATION for the MATLAB
port: every function here reproduces the numbers in RESULTS_LEDGER.md (checked by verify_ref.py against the
analysis scripts that produced them), and golden/*.json holds fixed inputs with expected outputs.

Conventions used everywhere (see PORTING_GUIDE.md for the full list):
  * time: count traces are 1 sample per second starting at ZT0; minute m covers seconds [60m, 60m+60).
  * symbols: Q = 0, W = 1, S = 2, missing = -1   (MATLAB export adds 1: 1 = Q, 2 = W, 3 = S).
  * hidden states after labelling: index 0 = H1 (quiet/regulator), 1 = H2 (wake), 2 = H3 (sleep).
  * phase arrays ph[t] give the transition matrix used for the step INTO minute t (ph[0] unused for transitions).
  * every count of a transition or route is attributed to the phase / window of its ARRIVAL minute.
  * sign: positive cycle current = forward cycle A -> Q -> S -> A, i.e. H2 -> H1 -> H3 -> H2.
"""
import numpy as np

DAY = 86400          # seconds
HALF = 43200         # lights off at ZT12 (LD 12:12)
MIN_SLEEP = 300      # 5-minute sleep criterion, seconds
Q, W, S = 0, 1, 2

# ---------------------------------------------------------------------------------------------------------
# 1. Data and sleep scoring
# ---------------------------------------------------------------------------------------------------------
def read_counts(path):
    """One integer per line = movement counts in that second. Returns int array (length = seconds)."""
    return np.loadtxt(path).astype(int)

def sleep_bouts(x, bridge=True, min_len=MIN_SLEEP):
    """Sleep bouts = runs of zero counts lasting >= min_len seconds.
    Bridging (lab convention, option b): a single active second is reset to inactive ONLY if the inactive runs on
    BOTH sides are each >= min_len seconds in the UNBRIDGED trace (301 + 1 + 300 -> one 602 s bout;
    301 + 1 + 100 -> not bridged). Returns list of (start_second, length_seconds), 0-based starts."""
    a = (np.asarray(x) > 0).astype(int)
    if bridge:
        z0 = 1 - a
        ch = np.r_[0, np.flatnonzero(np.diff(z0)) + 1, len(z0)]
        runs = [(s, e - s, z0[s]) for s, e in zip(ch[:-1], ch[1:])]
        for i in range(1, len(runs) - 1):
            s_, n_, zz = runs[i]
            if (zz == 0 and n_ == 1 and runs[i-1][2] == 1 and runs[i+1][2] == 1
                    and runs[i-1][1] >= min_len and runs[i+1][1] >= min_len):
                a[s_] = 0
    z = 1 - a
    ch = np.r_[0, np.flatnonzero(np.diff(z)) + 1, len(z)]
    return [(int(s), int(b - s)) for s, b in zip(ch[:-1], ch[1:]) if z[s] == 1 and b - s >= min_len]

def wake_bout_events(x, bouts, min_wake=600):
    """For every wake bout between two consecutive sleep bouts that lasts >= min_wake s:
       ref   = active fraction in the bout's middle 300 s ([mid-150, mid+150)),
       wake  = (phase of the waking second, active fraction in seconds 61-300 after waking)   [w0+60, w0+300)
       onset = (phase of the onset second,  active fraction in seconds 61-300 before onset)   [w1-300, w1-60)
    phase = 0 light (ZT 0-12), 1 dark, from the second's clock time. (PRESPEC 19)"""
    mv = (np.asarray(x) > 0).astype(float); ev = []
    for (a1, n1), (a2, _) in zip(bouts[:-1], bouts[1:]):
        w0 = a1 + n1; w1 = a2
        if w1 - w0 < min_wake:
            continue
        mid = (w0 + w1) // 2
        ev.append(dict(ref=mv[mid-150:mid+150].mean(),
                       wake=(int((w0 % DAY) >= HALF), mv[w0+60:w0+300].mean()),
                       onset=(int((w1 % DAY) >= HALF), mv[w1-300:w1-60].mean())))
    return ev

def wake_reference(events, phase=None, min_events=3):
    """Fly reference activity r_f for one block: median of `ref` over events whose waking OR onset falls in
    `phase` (None = whole record). NaN if fewer than min_events."""
    bb = [e['ref'] for e in events if phase is None or e['wake'][0] == phase or e['onset'][0] == phase]
    return float(np.median(bb)) if len(bb) >= min_events else np.nan

# ---------------------------------------------------------------------------------------------------------
# 2. Observable S/Q/W sequence (PRESPEC 20 section 1)
# ---------------------------------------------------------------------------------------------------------
def sqw_symbols(x, theta=0.5, light_dark=True, full_days=True):
    """1-minute symbols. S if >= 30 of the minute's 60 s lie in a (bridged) sleep bout; otherwise a wake minute,
    Q if its active fraction < theta * r_f(block of the minute), else W. r_f is computed per light/dark block
    (light_dark=True) or for the whole record (False, e.g. per0 in DD). A wake minute whose block has no
    reference gets -1 (missing). If full_days, the sequence is truncated to whole days (1440 min).
    Returns (y, refs) with refs = {'light': r, 'dark': r} or {'whole': r}."""
    x = np.asarray(x); nm = len(x) // 60; x = x[:nm * 60]
    sl = np.zeros(len(x), bool)
    bouts = sleep_bouts(x, True)
    for s, n in bouts:
        sl[s:s+n] = True
    ssec = sl.reshape(nm, 60).sum(1); act = (x > 0).reshape(nm, 60).mean(1)
    ev = wake_bout_events(x, bouts)
    refs = ({'light': wake_reference(ev, 0), 'dark': wake_reference(ev, 1)} if light_dark
            else {'whole': wake_reference(ev, None)})
    y = np.full(nm, W); y[ssec >= 30] = S
    for m in np.flatnonzero(ssec < 30):
        r = refs['whole'] if not light_dark else (refs['light'] if (m * 60) % DAY < HALF else refs['dark'])
        if np.isnan(r):
            y[m] = -1
        elif act[m] < theta * r:
            y[m] = Q
    if full_days:
        y = y[:(len(y) // 1440) * 1440]
    return y, refs

def phases(T, scheme='12h'):
    """Phase index of each minute (transition INTO minute t uses P[ph[t]]).
    '12h': 0 = ZT 0-12, 1 = ZT 12-24.  '6h': 0..3 = ZT 0-6, 6-12, 12-18, 18-24.  'none': all 0."""
    t = np.arange(T) % 1440
    if scheme == '12h': return (t >= 720).astype(int)
    if scheme == '6h': return t // 360
    if scheme == 'none': return np.zeros(T, int)
    raise ValueError(scheme)

# ---------------------------------------------------------------------------------------------------------
# 3. HMM: shared emissions B (K x 3), phase-specific transitions P (nP x K x K), initial pi (K)
# ---------------------------------------------------------------------------------------------------------
def emission_lik(B, y):
    """E[t, k] = B[k, y_t]; missing symbols (y < 0) contribute 1."""
    E = B[:, np.clip(y, 0, 2)].T.copy(); E[np.asarray(y) < 0] = 1.0; return E

def forward(pi, P, B, y, ph):
    """Scaled forward pass. Returns (alpha_hat (T x K), c (T)), with log-likelihood = sum(log c)."""
    E = emission_lik(B, y); T, K = E.shape; A = np.empty((T, K)); c = np.empty(T)
    a = pi * E[0]; c[0] = a.sum(); A[0] = a / c[0]
    for t in range(1, T):
        a = (A[t-1] @ P[ph[t]]) * E[t]; c[t] = a.sum(); A[t] = a / c[t]
    return A, c

def loglik(pi, P, B, y, ph):
    return float(np.log(forward(pi, P, B, y, ph)[1]).sum())

def em_step(pi, P, B, y, ph):
    """ONE Baum-Welch update (deterministic). Returns (pi', P', B', loglik of the INPUT parameters).
    Transition counts are split by the phase of the arrival minute; emissions pooled over phases."""
    y = np.asarray(y); ph = np.asarray(ph); nP = P.shape[0]
    E = emission_lik(B, y); A, c = forward(pi, P, B, y, ph); T, K = E.shape
    Bk = np.empty_like(A); Bk[-1] = 1.0
    for t in range(T - 2, -1, -1):
        Bk[t] = P[ph[t+1]] @ (E[t+1] * Bk[t+1]) / c[t+1]
    G = A * Bk; G /= G.sum(1, keepdims=True); R = E * Bk / c[:, None]
    P2 = np.empty_like(P)
    for p in range(nP):
        idx = np.flatnonzero(ph[1:] == p) + 1
        Xi = (A[idx-1].T @ R[idx]) * P[p]
        P2[p] = Xi / np.maximum(Xi.sum(1, keepdims=True), 1e-300)
    onehot = np.stack([(y == v) for v in range(3)], 1).astype(float)
    Bn = G.T @ onehot; B2 = Bn / np.maximum(Bn.sum(1, keepdims=True), 1e-300)
    return G[0].copy(), P2, B2, float(np.log(c).sum())

def em_fit(y, ph, nP, K, starts=20, seed=20, maxit=1000, tol=1e-7):
    """Maximum likelihood by Baum-Welch from `starts` random initialisations run in parallel (batched);
    stops when every start's relative log-likelihood gain < tol, or at maxit. Returns the best start:
    (loglik, pi, P (nP x K x K), B (K x 3), iterations).
    Initialisation: pi, rows of P, rows of B ~ Dirichlet(1); the first half of the starts have
    P <- 0.5 P + 0.5 I (diagonal-heavy). The random draws are numpy-specific; a port cannot reproduce them and
    should be verified with em_step (exact) and by reaching the same optimum (likelihood) on real data."""
    rng = np.random.default_rng(seed); Sn = starts; y = np.asarray(y); ph = np.asarray(ph)
    pi = rng.dirichlet(np.ones(K), Sn); P = rng.dirichlet(np.ones(K), (Sn, nP, K))
    half = Sn // 2; P[:half] = 0.5 * P[:half] + 0.5 * np.eye(K)
    B = rng.dirichlet(np.ones(3), (Sn, K)); ll_old = np.full(Sn, -np.inf); T = len(y)
    idx = [np.flatnonzero(ph[1:] == p) + 1 for p in range(nP)]
    onehot = np.stack([(y == v) for v in range(3)], 1).astype(float)
    def Eb(B):
        E = B[:, :, np.clip(y, 0, 2)].transpose(2, 0, 1).copy(); E[y < 0] = 1.0; return E
    def fwd(pi, P, E):
        A = np.empty((T, Sn, K)); c = np.empty((T, Sn)); a = pi * E[0]; c[0] = a.sum(1); A[0] = a / c[0][:, None]
        for t in range(1, T):
            a = np.einsum('sk,skl->sl', A[t-1], P[:, ph[t]]) * E[t]; c[t] = a.sum(1); A[t] = a / c[t][:, None]
        return A, c
    for it in range(maxit):
        E = Eb(B); A, c = fwd(pi, P, E); ll = np.log(c).sum(0)
        Bk = np.empty_like(A); Bk[-1] = 1.0
        for t in range(T - 2, -1, -1):
            Bk[t] = np.einsum('skl,sl->sk', P[:, ph[t+1]], E[t+1] * Bk[t+1]) / c[t+1][:, None]
        G = A * Bk; G /= G.sum(2, keepdims=True); R = E * Bk / c[:, :, None]
        for p in range(nP):
            Xi = np.einsum('tsk,tsl->skl', A[idx[p]-1], R[idx[p]]) * P[:, p]
            P[:, p] = Xi / np.maximum(Xi.sum(2, keepdims=True), 1e-300)
        Bn = np.einsum('tsk,tv->skv', G, onehot); B = Bn / np.maximum(Bn.sum(2, keepdims=True), 1e-300)
        pi = G[0]
        if np.all(ll - ll_old < tol * np.abs(ll)):
            break
        ll_old = ll
    E = Eb(B); _, c = fwd(pi, P, E); ll = np.log(c).sum(0); ll[~np.isfinite(ll)] = -np.inf
    s = int(np.argmax(ll)); return float(ll[s]), pi[s], P[s], B[s], it + 1

def n_params(K, nP):
    """Free parameters: (K-1) initial + nP*K*(K-1) transitions + 2K emissions (3 symbols)."""
    return (K - 1) + nP * K * (K - 1) + 2 * K

def bic(ll, K, nP, n):
    """BIC = -2 logL + k ln n. Primary n = number of minutes; sensitivity n_eff = number of symbol changes."""
    return -2 * ll + n_params(K, nP) * np.log(n)

def n_eff(y):
    return max(2, int(np.sum(np.diff(y) != 0)))

def label_states(P, B, pi):
    """Reorder a K = 3 fit to (H1, H2, H3): H3 = largest P(S); of the other two, H2 = larger P(W); H1 = rest.
    The same permutation is applied to every phase's matrix. P: (nP x 3 x 3)."""
    h3 = int(np.argmax(B[:, S])); rest = [i for i in range(3) if i != h3]
    h2 = rest[int(np.argmax(B[rest, W]))]; h1 = [i for i in rest if i != h2][0]; o = [h1, h2, h3]
    return P[:, o][:, :, o], B[o], pi[o]

def is_degenerate(B, thr=0.05):
    """H1 is not identified if its emission row lies within total variation thr of H2's or H3's."""
    tv = lambda a, b: 0.5 * np.abs(a - b).sum()
    return min(tv(B[0], B[1]), tv(B[0], B[2])) < thr

# ---------------------------------------------------------------------------------------------------------
# 4. Routing predicted from a transition matrix (one phase, labelled H1, H2, H3 = indices 0, 1, 2)
# ---------------------------------------------------------------------------------------------------------
def xn_pred(P):
    """Sleep -> wake: Xn = (r - 1)/(r + 1), r = P(H3->H1) * P(H1->H2)/(P(H1->H2)+P(H1->H3)) / P(H3->H2)."""
    with np.errstate(all='ignore'):
        kX = P[0, 1] / (P[0, 1] + P[0, 2]); a = P[2, 0] * kX; return (a - P[2, 1]) / (a + P[2, 1])

def yn_pred(P):
    """Wake -> sleep: Yn with r = P(H2->H1) * P(H1->H3)/(P(H1->H2)+P(H1->H3)) / P(H2->H3)."""
    with np.errstate(all='ignore'):
        kY = P[0, 2] / (P[0, 1] + P[0, 2]); a = P[1, 0] * kY; return (a - P[1, 2]) / (a + P[1, 2])

def c_pred(P):
    """Signed cycle current (Yn - Xn)/2. sign(c_pred) = sign(forward affinity) = sign(J_model)."""
    return (yn_pred(P) - xn_pred(P)) / 2

def stationary(P):
    w, v = np.linalg.eig(P.T); p = np.real(v[:, np.argmin(np.abs(w - 1))]); return p / p.sum()

def model_energetics(P, zero=1e-10, per_hour=60.0):
    """For one phase's 3x3 matrix at its stationary distribution pi, per hour (steps are minutes):
       J   = pi_H2 P(H2->H1) - pi_H1 P(H1->H2)    (forward-positive; equal on every edge of the cycle)
       EPR = sum_{i<j} (pi_i P_ij - pi_j P_ji) ln(pi_i P_ij / (pi_j P_ji))
       A   = ln[P(H2->H1) P(H1->H3) P(H3->H2) / (P(H1->H2) P(H3->H1) P(H2->H3))]   (forward affinity)
    An edge with P < zero in exactly one direction makes EPR and A infinite (never floored)."""
    P = np.asarray(P); Z = P < zero; pi = stationary(P)
    J = per_hour * (pi[1] * P[1, 0] - pi[0] * P[0, 1])
    with np.errstate(all='ignore'):
        A = np.log(P[1, 0] * P[0, 2] * P[2, 1] / (P[0, 1] * P[2, 0] * P[1, 2]))
    if np.any(Z != Z.T):
        return dict(J=J, EPR=np.inf, A=A)
    epr = sum((pi[i] * P[i, j] - pi[j] * P[j, i]) * np.log(pi[i] * P[i, j] / (pi[j] * P[j, i]))
              for i in range(3) for j in range(i + 1, 3) if not Z[i, j])
    return dict(J=J, EPR=per_hour * epr, A=A)

# ---------------------------------------------------------------------------------------------------------
# 5. Decoding: FFBS posterior paths and Viterbi
# ---------------------------------------------------------------------------------------------------------
def ffbs(pi, P, B, y, ph, n_paths=200, seed=0):
    """Forward filtering, backward sampling: n_paths hidden paths (n_paths x T) from the posterior."""
    A, _ = forward(pi, P, B, y, ph); T = len(y); rng = np.random.default_rng(seed)
    Z = np.empty((n_paths, T), int); c = np.cumsum(A[-1]); Z[:, -1] = np.searchsorted(c / c[-1], rng.random(n_paths))
    for t in range(T - 2, -1, -1):
        pr = np.cumsum(A[t][:, None] * P[ph[t+1]][:, Z[:, t+1]], 0); pr /= pr[-1]
        Z[:, t] = (rng.random(n_paths)[None, :] > pr).sum(0)
    return Z

def viterbi(pi, P, B, y, ph):
    """Most probable hidden path (log domain). Ties go to the lowest state index."""
    E = emission_lik(B, y); T = len(y)
    with np.errstate(divide='ignore'):
        lE = np.log(E); lP = np.log(P); d = np.log(pi) + lE[0]; bp = np.empty((T, P.shape[1]), int)
        for t in range(1, T):
            m = d[:, None] + lP[ph[t]]; bp[t] = m.argmax(0); d = m.max(0) + lE[t]
    z = np.empty(T, int); z[-1] = d.argmax()
    for t in range(T - 1, 0, -1):
        z[t-1] = bp[t, z[t]]
    return z

# ---------------------------------------------------------------------------------------------------------
# 6. Counting on a hidden path (run-compressed)
# ---------------------------------------------------------------------------------------------------------
ROUTES = ['N312', 'N32', 'N213', 'N23']

def route_events(z):
    """Routes between sleep and wake on the run-compressed path, with the minute of ARRIVAL in the destination:
       N312 = H3 -> H1 -> H2, N32 = H3 -> H2, N213 = H2 -> H1 -> H3, N23 = H2 -> H3.
    Excursions H3 -> H1 -> H3 and H2 -> H1 -> H2 are not routes. Returns list of (name, t)."""
    z = np.asarray(z); ch = np.r_[0, np.flatnonzero(np.diff(z)) + 1]; st = z[ch]; out = []
    for k in range(1, len(st)):
        a, b = st[k-1], st[k]
        if (a, b) == (2, 1): out.append(('N32', int(ch[k])))
        elif (a, b) == (1, 2): out.append(('N23', int(ch[k])))
        elif k >= 2 and st[k-1] == 0:
            p = st[k-2]
            if (p, b) == (2, 1): out.append(('N312', int(ch[k])))
            elif (p, b) == (1, 2): out.append(('N213', int(ch[k])))
    return out

def route_counts(z, ph, nph):
    """(nph x 4) counts in ROUTES order, by phase of arrival."""
    c = np.zeros((nph, 4))
    for r, t in route_events(z):
        c[ph[t], ROUTES.index(r)] += 1
    return c

def pair_counts(z, ph, nph):
    """(nph x 3 x 3) single-step counts N[i, j] = number of i -> j changes on the run-compressed path, by phase of
    the arrival minute. Same object as compute_HMM_energetics.m's N12..N32 (GLOBAL, no window seams)."""
    z = np.asarray(z); ch = np.flatnonzero(np.diff(z)) + 1; c = np.zeros((nph, 3, 3))
    np.add.at(c, (np.asarray(ph)[ch], z[ch-1], z[ch]), 1); return c

def xy_counts(c, min_routes=3):
    """Xn, Yn from route counts (last axis in ROUTES order); NaN if a denominator < min_routes."""
    c = np.asarray(c, float)
    with np.errstate(all='ignore'):
        X = (c[..., 0] - c[..., 1]) / (c[..., 0] + c[..., 1]); Y = (c[..., 2] - c[..., 3]) / (c[..., 2] + c[..., 3])
    X = np.where(c[..., 0] + c[..., 1] >= min_routes, X, np.nan); Y = np.where(c[..., 2] + c[..., 3] >= min_routes, Y, np.nan)
    return X, Y

def j_cycle_count(N, hours):
    """J_cycle = (N32 - N23)/hours, forward-positive (H3 -> H2 is the S -> A leg). Linear in counts, so the mean
    over posterior paths equals J of the mean counts (no Jensen issue)."""
    return (N[2][1] - N[1][2]) / hours

def epr_count_bound(N, hours):
    """Schnakenberg EPR from directed counts, nats per hour, with the rule-of-three 95% LOWER BOUND for a pair
    observed in one direction only (identical to epr_count_based.m, rule 'bound'):
        both > 0      : (a - b) ln(a / b)
        one-way, k    : (k - 3) ln(k / 3) if k > 3, else 0     [b* = min(3, k) enters prefactor AND log]
        both zero     : 0
    Returns (EPR, substituted?). Asserts EPR >= 0 (every term is a product of two same-sign factors)."""
    e = 0.0; sub = False
    for i in range(3):
        for j in range(i + 1, 3):
            a, b = float(N[i][j]), float(N[j][i])
            if a == 0 and b == 0:
                continue
            if b == 0: b = min(3.0, a); sub = True
            elif a == 0: a = min(3.0, b); sub = True
            e += (a - b) * np.log(a / b)
    assert e >= 0, 'EPR bound negative: impossible'
    return e / hours, sub

def posterior_summaries(Z, ph, nph, hours):
    """Per-path statistics averaged over posterior paths (E[f(N)], respecting Jensen), per phase.
    Returns dict of arrays (nph): Xn, Yn, c (means over paths where defined, NaN if defined on < 50% of paths),
    J (mean), EPR_bound (mean), frac_substituted, and the mean route / pair counts."""
    rp = np.stack([route_counts(z, ph, nph) for z in Z]); pp = np.stack([pair_counts(z, ph, nph) for z in Z])
    X, Y = xy_counts(rp); C = (Y - X) / 2
    def m(v):
        ok = np.isfinite(v).mean(0) >= 0.5
        with np.errstate(all='ignore'):
            return np.where(ok, np.nanmean(np.where(np.isfinite(v), v, np.nan), 0), np.nan)
    eb = np.array([[epr_count_bound(pp[k, p], hours[p]) for p in range(nph)] for k in range(len(Z))], dtype=object)
    return dict(Xn=m(X), Yn=m(Y), c=m(C),
                J=np.array([np.mean([j_cycle_count(pp[k, p], hours[p]) for k in range(len(Z))]) for p in range(nph)]),
                EPR_bound=np.array([np.mean([eb[k, p][0] for k in range(len(Z))]) for p in range(nph)]),
                frac_substituted=np.array([np.mean([eb[k, p][1] for k in range(len(Z))]) for p in range(nph)]),
                route_counts=rp.mean(0), pair_counts=pp.mean(0))

# ---------------------------------------------------------------------------------------------------------
# 7. Observable-level currents (no HMM)
# ---------------------------------------------------------------------------------------------------------
def c_obs_minutes(y, ph, nph):
    """Observable current on the minute symbols (PRESPEC 20d): P(last wake minute before a sleep onset is Q)
    - P(first wake minute after a waking is Q), per phase of the transition's second minute."""
    y = np.asarray(y); t = np.arange(1, len(y))
    on = t[(y[:-1] >= 0) & (y[:-1] < 2) & (y[1:] == 2)]; wk = t[(y[:-1] == 2) & (y[1:] >= 0) & (y[1:] < 2)]; out = []
    for p in range(nph):
        o = on[ph[on] == p]; w = wk[ph[wk] == p]
        out.append(np.mean(y[o-1] == Q) - np.mean(y[w] == Q) if len(o) and len(w) else np.nan)
    return np.array(out)

def cycle_current_model_free(events, phase, theta, hours):
    """PRESPEC 19 signed cycle current from the 1-s counts (+-5 min window, no HMM). Transit activity of each
    onset / waking (wake_bout_events) is Q if < theta * r_f, else A. c = P(Q | onset) - P(Q | waking);
    J = [(nQS - nSQ) + (nSA - nAS)] / 2 / hours. None if < 3 events or no onset / waking in the phase."""
    bb = [e for e in events if phase is None or e['wake'][0] == phase or e['onset'][0] == phase]
    if len(bb) < 3:
        return None
    thr = theta * np.median([e['ref'] for e in bb])
    on = [e['onset'][1] < thr for e in events if phase is None or e['onset'][0] == phase]
    wk = [e['wake'][1] < thr for e in events if phase is None or e['wake'][0] == phase]
    if not on or not wk:
        return None
    nQS = sum(on); nAS = len(on) - nQS; nSQ = sum(wk); nSA = len(wk) - nSQ
    return dict(c=nQS / len(on) - nSQ / len(wk), J=((nQS - nSQ) + (nSA - nAS)) / 2 / hours,
                nQS=nQS, nAS=nAS, nSQ=nSQ, nSA=nSA)

# ---------------------------------------------------------------------------------------------------------
# 8. Simulation (for validation, bootstrap and nulls)
# ---------------------------------------------------------------------------------------------------------
def simulate(pi, P, B, ph, rng):
    """Hidden path z and symbols y of length len(ph) from (pi, P[ph[t]], B)."""
    T = len(ph); K = len(pi); z = np.empty(T, int); z[0] = rng.choice(K, p=pi); u = rng.random(T)
    for t in range(1, T):
        z[t] = min(np.searchsorted(np.cumsum(P[ph[t]][z[t-1]]), u[t]), K - 1)
    cb = np.cumsum(B, 1); v = rng.random(T)
    y = np.array([min(np.searchsorted(cb[z[t]], v[t]), 2) for t in range(T)])
    return z, y

# ---------------------------------------------------------------------------------------------------------
# 9. Observable-level counts and round trips around sleep bouts (PRESPEC 20f, 20g)
# ---------------------------------------------------------------------------------------------------------
ROUND_TRIPS = ['forward', 'reverse', 'direct', 'quiet_both']   # (entry, exit) = (Q, A), (A, Q), (A, A), (Q, Q)

def _rt_index(q_in, q_out):
    return {(1, 0): 0, (0, 1): 1, (0, 0): 2, (1, 1): 3}[(int(q_in), int(q_out))]

def observable_pair_counts(y, ph, nph):
    """(nph x 3 x 3) single-step counts on the run-compressed S/Q/W symbols (Q = 0, W = 1, S = 2 share the H1, H2,
    H3 indices), by phase of the arrival minute. A change that touches a missing minute (-1) is not counted.
    Feed to j_cycle_count / epr_count_bound for the observable-level J and EPR (PRESPEC 20f)."""
    y = np.asarray(y); t = np.flatnonzero(y >= 0); ys = y[t]; ps = np.asarray(ph)[t]
    ch = np.flatnonzero(np.diff(ys)) + 1; ok = (t[ch] - t[ch - 1]) == 1
    N = np.zeros((nph, 3, 3)); np.add.at(N, (ps[ch[ok]], ys[ch - 1][ok], ys[ch][ok]), 1); return N

def round_trips_observed(x, theta=0.5, light_dark=True, min_wake=600):
    """(nph x 4) counts of round-trip types from the 1-s counts. A sleep bout (bridged scoring) is used if the wake
    bouts before AND after it each last >= min_wake s. Entry is Q if the active fraction in seconds 61-300 before the
    onset is < theta * r_f; exit is Q if the active fraction in seconds 61-300 after waking is < theta * r_f; r_f is
    the reference of the onset's phase (wake_reference). Phase of a round trip = phase of its sleep onset."""
    nph = 2 if light_dark else 1; bouts = sleep_bouts(x, True); mv = (np.asarray(x) > 0).astype(float)
    ev = wake_bout_events(x, bouts)
    ref = [wake_reference(ev, 0), wake_reference(ev, 1)] if light_dark else [wake_reference(ev, None)]
    C = np.zeros((nph, 4))
    for k in range(1, len(bouts) - 1):
        s, n = bouts[k]; e = s + n
        if s - (bouts[k-1][0] + bouts[k-1][1]) < min_wake or bouts[k+1][0] - e < min_wake:
            continue
        p = int((s % DAY) >= HALF) if light_dark else 0; r = ref[p]
        if not np.isfinite(r):
            continue
        C[p, _rt_index(mv[s-300:s-60].mean() < theta * r, mv[e+60:e+300].mean() < theta * r)] += 1
    return C

def round_trips_hidden(z, ph, nph, min_sleep=0, min_gap=0):
    """(nph x 4) counts of round-trip types on one hidden path (run-compressed); H1 plays Q, H2 plays A.
    Defaults (PRESPEC 20g primary): every H3 run with a neighbour on both sides. Matched definition (20g amendment B):
    min_sleep=5, min_gap=10 -- H3 runs >= 5 min with >= 10 min free of H3 before and after (the previous and the next
    H3 run must exist). Phase = phase of the H3 run's first minute."""
    z = np.asarray(z); ch = np.r_[0, np.flatnonzero(np.diff(z)) + 1, len(z)]; st = z[ch[:-1]]; ln = np.diff(ch)
    s_idx = np.flatnonzero(st == 2); C = np.zeros((nph, 4))
    for a, k in enumerate(s_idx):
        if k == 0 or k == len(st) - 1 or ln[k] < min_sleep:
            continue
        if min_gap > 0:
            if a == 0 or a + 1 == len(s_idx):
                continue
            if ch[k] - ch[s_idx[a-1] + 1] < min_gap or ch[s_idx[a+1]] - ch[k+1] < min_gap:
                continue
        C[ph[ch[k]], _rt_index(st[k-1] == 0, st[k+1] == 0)] += 1
    return C

def round_trip_stats(counts, hours):
    """Summary of round trips for one fly and phase. `counts` is either (4,) for one record (observed level) or
    (paths x 4) for posterior paths (hidden level: every statistic is computed per path, then averaged; paths with
    no round trip are skipped). Returns n (mean count), fractions (4), rates per hour (4), P(Q entry), P(Q exit),
    and dependence = fraction - product of marginals (4), the independence check of PRESPEC 20g amendment A."""
    C = np.atleast_2d(np.asarray(counts, float)); n = C.sum(1); ok = n > 0
    if not ok.any():
        return None
    Fr = C[ok] / n[ok, None]; qi = Fr[:, 0] + Fr[:, 3]; qo = Fr[:, 1] + Fr[:, 3]
    prod = np.stack([qi * (1 - qo), (1 - qi) * qo, (1 - qi) * (1 - qo), qi * qo], 1)
    return dict(n=float(n.mean()), frac=Fr.mean(0), rate=C.mean(0) / hours, p_q_entry=float(qi.mean()),
                p_q_exit=float(qo.mean()), dependence=(Fr - prod).mean(0), frac_paths_used=float(ok.mean()))
