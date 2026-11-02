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
WAKE_UNKNOWN = -2   # a wake minute whose Q/W type is unknown (no reference r_f): partial observation, see PORTING_GUIDE

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
def sqw_symbols(x, theta=0.5, light_dark=True, full_days=True, unknown_wake=-1):
    """1-minute symbols. S if >= 30 of the minute's 60 s lie in a (bridged) sleep bout; otherwise a wake minute,
    Q if its active fraction < theta * r_f(block of the minute), else W. r_f is computed per light/dark block
    (light_dark=True) or for the whole record (False, e.g. per0 in DD). A wake minute whose block has no
    reference gets `unknown_wake`: -1 (missing, as in all analyses run so far) or WAKE_UNKNOWN = -2 (recommended:
    'wake of unknown type', likelihood P(Q|H) + P(W|H); see PORTING_GUIDE section 4). If full_days, the sequence is truncated to whole days (1440 min).
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
            y[m] = unknown_wake
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
    """E[t, k] = B[k, y_t]; missing symbols (-1) contribute 1; WAKE_UNKNOWN (-2) contributes B[k, Q] + B[k, W]."""
    y = np.asarray(y); E = B[:, np.clip(y, 0, 2)].T.copy(); E[y == -1] = 1.0
    E[y == WAKE_UNKNOWN] = (B[:, Q] + B[:, W])[None, :]; return E

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
    Bn = G.T @ onehot
    u = y == WAKE_UNKNOWN
    if u.any():                    # split each unknown-wake minute between Q and W in proportion to B
        g = G[u].sum(0); qw = B[:, Q] + B[:, W]; qw = np.where(qw > 0, qw, 1.0)
        Bn[:, Q] += g * B[:, Q] / qw; Bn[:, W] += g * B[:, W] / qw
    B2 = Bn / np.maximum(Bn.sum(1, keepdims=True), 1e-300)
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
    u = y == WAKE_UNKNOWN
    def Eb(B):
        E = B[:, :, np.clip(y, 0, 2)].transpose(2, 0, 1).copy(); E[y == -1] = 1.0
        if u.any(): E[u] = (B[:, :, Q] + B[:, :, W])[None]
        return E
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
        Bn = np.einsum('tsk,tv->skv', G, onehot)
        if u.any():
            gq = G[u].sum(0); qw = B[:, :, Q] + B[:, :, W]; qw = np.where(qw > 0, qw, 1.0)
            Bn[:, :, Q] += gq * B[:, :, Q] / qw; Bn[:, :, W] += gq * B[:, :, W] / qw
        B = Bn / np.maximum(Bn.sum(2, keepdims=True), 1e-300)
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

# ---------------------------------------------------------------------------------------------------------
# 10. Dwell times (PRESPEC 20h): runs with censoring, Kaplan-Meier, life-table hazard
# ---------------------------------------------------------------------------------------------------------
DWELL_EDGES = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 721]   # minutes, fixed in PRESPEC 20h

def runs_censored(seq, ph, values, nph, censor_at_phase=True):
    """Runs of any value in `values` (e.g. {2} for S, {0, 1} for wake) on a 1-min sequence. Returns, per phase of the
    run's first minute, arrays (durations in minutes, event flags). A run is censored (event 0) where it is cut by a
    phase boundary (if censor_at_phase), a missing minute (-1) or the end of the record; event 1 = it ended in a
    different valid value."""
    s = np.asarray(seq); ph = np.asarray(ph); inr = np.isin(s, list(values)); out = [([], []) for _ in range(nph)]
    T = len(s); t = 0
    while t < T:
        if not inr[t]: t += 1; continue
        p = ph[t]; u = t + 1
        while u < T and inr[u] and not (censor_at_phase and ph[u] != p): u += 1
        event = int(u < T and s[u] >= 0 and not inr[u])
        out[p][0].append(u - t); out[p][1].append(event); t = u
        if u < T and inr[u]:            # cut at the phase boundary: the remainder starts a new (left-truncated) run
            while u < T and inr[u]: u += 1
            t = u                       # the remainder is not counted as a new run (its start is not observed)
    return [(np.array(d, int), np.array(e, int)) for d, e in out]

def km_survival(d, e, w=None):
    """Kaplan-Meier survival of durations d (integer minutes) with event flags e and optional weights w.
    Returns (t, S(t)) at the distinct durations; S(t) = P(duration > t)."""
    d = np.asarray(d); e = np.asarray(e); w = np.ones(len(d)) if w is None else np.asarray(w, float)
    ts = np.unique(d); S = []; s = 1.0
    for t in ts:
        risk = w[d >= t].sum(); ev = w[(d == t) & (e == 1)].sum()
        if risk > 0: s *= 1 - ev / risk
        S.append(s)
    return ts, np.array(S)

def life_table_hazard(d, e, edges=DWELL_EDGES, w=None):
    """Per-minute hazard in bins [lo, hi): events / minutes at risk. A run of length d is at risk at minutes 1..d.
    Returns (events, at_risk_minutes, n_runs_at_risk) per bin."""
    d = np.asarray(d); e = np.asarray(e); w = np.ones(len(d)) if w is None else np.asarray(w, float)
    ev, ex, nr = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        ex.append((w * np.clip(np.minimum(d, hi - 1) - lo + 1, 0, None)).sum())
        ev.append((w * ((e == 1) & (d >= lo) & (d < hi))).sum()); nr.append((w * (d >= lo)).sum())
    return np.array(ev), np.array(ex), np.array(nr)

def geometric_survival(p_stay, t):
    """P(dwell > t) for a hidden state with self-transition probability p_stay (per minute): p_stay ** t."""
    return np.asarray(p_stay, float) ** np.asarray(t, float)

# ---------------------------------------------------------------------------------------------------------
# 11. EPR vs J_cycle across flies: linear vs quadratic, Onsager conductance (PRESPEC 20j; Fig. 11, 12)
# ---------------------------------------------------------------------------------------------------------
def onsager_fit(J, EPR):
    """EPR = J^2 / L by ordinary least squares THROUGH THE ORIGIN (user's Methods Eq. 25). Returns dict with
    slope s = 1/L, its SE (df = n - 1), L, SE_L = L^2 SE_s, the 95% interval of L (t interval of s inverted; upper
    bound inf if the lower slope bound <= 0) and R^2 through the origin (1 - RSS / sum EPR^2)."""
    from scipy.stats import t as _t
    x = np.asarray(J, float) ** 2; y = np.asarray(EPR, float); n = len(x)
    s = (x @ y) / (x @ x); res = y - s * x; se = np.sqrt((res @ res) / (n - 1) / (x @ x))
    q = _t.ppf(0.975, n - 1); lo, hi = s - q * se, s + q * se
    return dict(s=s, se_s=se, L=1 / s, se_L=se / s ** 2, ci_L=(1 / hi, 1 / lo if lo > 0 else np.inf), n=n,
                r2=1 - (res @ res) / (y @ y))

def onsager_bootstrap(points, n_boot=2000, seed=20):
    """Fly bootstrap of L. points: list of (fly_id, J, EPR), several per fly (e.g. light and dark); flies are
    resampled with replacement and all their points kept together. Returns the array of bootstrap L values."""
    rng = np.random.default_rng(seed); flies = sorted({p[0] for p in points}); by = {f: [p for p in points if p[0] == f] for f in flies}
    out = []
    for _ in range(n_boot):
        rr = [p for f in rng.choice(flies, len(flies)) for p in by[f]]
        out.append(onsager_fit([p[1] for p in rr], [p[2] for p in rr])['L'])
    return np.array(out)

def bootstrap_difference_p(boot_a, boot_b):
    """Difference of two independent bootstrap distributions: (95% interval of a - b, two-sided p =
    2 min(P(diff <= 0), P(diff >= 0)))."""
    d = np.asarray(boot_a) - np.asarray(boot_b)
    return np.percentile(d, [2.5, 97.5]), min(1.0, 2 * min(np.mean(d <= 0), np.mean(d >= 0)))

def epr_j_models(J, EPR, offset=0.02, j_min=0.05, epr_min=0.05):
    """PRESPEC 20j descriptions with a free intercept b (absorbs the estimator offset; not physics):
    M1 constant affinity EPR = b + a|J|; M2 linear response EPR = b + c J^2 (c = 1/L); compared by residual sum of
    squares (both 2 parameters). Exponent n: slope of ln(EPR - offset) on ln|J| over points with |J| >= j_min and
    EPR >= epr_min (NaN if fewer than 5)."""
    J = np.asarray(J, float); E = np.asarray(EPR, float)
    def ls(x):
        A = np.c_[np.ones_like(x), x]; coef = np.linalg.lstsq(A, E, rcond=None)[0]; return coef[1], coef[0], float(((A @ coef - E) ** 2).sum())
    a, b1, r1 = ls(np.abs(J)); c, b2, r2 = ls(J ** 2)
    m = (np.abs(J) >= j_min) & (E >= epr_min)
    n = np.polyfit(np.log(np.abs(J[m])), np.log(E[m] - offset), 1)[0] if m.sum() >= 5 else np.nan
    return dict(a=a, b_linear=b1, rss_linear=r1, c=c, b_quadratic=b2, rss_quadratic=r2, exponent=n, n_points=len(J), n_exponent=int(m.sum()))

# ---------------------------------------------------------------------------------------------------------
# 12. Route vs bout duration (PRESPEC 20k)
# ---------------------------------------------------------------------------------------------------------
def route_duration_observed(x, theta=0.5, light_dark=True, min_wake=600):
    """One record per bridged sleep bout with a neighbour on both sides: dict(phase_on, phase_wake_end,
    sleep_min, wake_before_min, wake_after_min, q_in, q_out). q_in (q_out) is None unless the preceding (following)
    wake bout lasts >= min_wake s; otherwise True/False for a quiet (Q) route, with r_f of the onset's phase."""
    bouts = sleep_bouts(x, True); mv = (np.asarray(x) > 0).astype(float); ev = wake_bout_events(x, bouts)
    ref = [wake_reference(ev, 0), wake_reference(ev, 1)] if light_dark else [wake_reference(ev, None)]
    out = []
    for k in range(1, len(bouts) - 1):
        s, n = bouts[k]; e = s + n; wb = s - (bouts[k-1][0] + bouts[k-1][1]); wa = bouts[k+1][0] - e
        p = int((s % DAY) >= HALF) if light_dark else 0; r = ref[p]
        if not np.isfinite(r): continue
        out.append(dict(phase_on=p, phase_wake_end=p, sleep_min=n / 60, wake_before_min=wb / 60, wake_after_min=wa / 60,
                        q_in=(mv[s-300:s-60].mean() < theta * r) if wb >= min_wake else None,
                        q_out=(mv[e+60:e+300].mean() < theta * r) if wa >= min_wake else None))
    return out

def route_duration_hidden(z, ph, min_sleep=5, min_gap=10):
    """Same records on one hidden path: H3 runs >= min_sleep min with a previous and a next H3 run; entry (exit)
    classified only if the gap since the previous (until the next) H3 run is >= min_gap min; quiet = H1."""
    z = np.asarray(z); ch = np.r_[0, np.flatnonzero(np.diff(z)) + 1, len(z)]; st = z[ch[:-1]]; ln = np.diff(ch)
    s_idx = np.flatnonzero(st == 2); out = []
    for a in range(1, len(s_idx) - 1):
        k = s_idx[a]
        if ln[k] < min_sleep: continue
        gb = ch[k] - ch[s_idx[a-1] + 1]; ga = ch[s_idx[a+1]] - ch[k+1]
        out.append(dict(phase_on=int(ph[ch[k]]), phase_wake_end=int(ph[ch[k]]), sleep_min=float(ln[k]),
                        wake_before_min=float(gb), wake_after_min=float(ga),
                        q_in=(st[k-1] == 0) if gb >= min_gap else None, q_out=(st[k+1] == 0) if ga >= min_gap else None))
    return out

RELATIONS = {'R1_exit_vs_sleep': ('q_out', 'sleep_min'), 'R2_entry_vs_wake_before': ('q_in', 'wake_before_min'),
             'R3_entry_vs_sleep': ('q_in', 'sleep_min')}

def route_duration_D(records, relation, phase, min_each=2):
    """D = mean ln(duration | quiet) - mean ln(duration | active) for one relation and phase; NaN if fewer than
    min_each of either kind."""
    lab, dur = RELATIONS[relation]
    rr = [r for r in records if r['phase_on'] == phase and r[lab] is not None]
    q = [np.log(r[dur]) for r in rr if r[lab]]; a = [np.log(r[dur]) for r in rr if not r[lab]]
    return np.mean(q) - np.mean(a) if len(q) >= min_each and len(a) >= min_each else np.nan

# ---------------------------------------------------------------------------------------------------------
# 13. Energy landscape (user's Methods Eqs. 5-16), per phase matrix (PRESPEC 20l)
# ---------------------------------------------------------------------------------------------------------
def landscape(P):
    """For one labelled 3x3 matrix (H1, H2, H3 = indices 0, 1, 2): state energies relative to H1,
    E2 - E1 = ln(pi_1/pi_2), E3 - E1 = ln(pi_1/pi_3); branching b12 = P12/(P12+P13), b13 = 1 - b12;
    branching-corrected barriers E31c = -ln(P31 b12), E21c = -ln(P21 b13); direct barriers E32 = -ln P32,
    E23 = -ln P23; predicted routing F_SW = tanh((E32 - E31c)/2), F_WS = tanh((E23 - E21c)/2), identical to
    xn_pred / yn_pred."""
    P = np.asarray(P, float); pi = stationary(P)
    with np.errstate(all='ignore'):
        b12 = P[0, 1] / (P[0, 1] + P[0, 2]); b13 = 1 - b12
        E31c = -np.log(P[2, 0] * b12); E21c = -np.log(P[1, 0] * b13); E32 = -np.log(P[2, 1]); E23 = -np.log(P[1, 2])
        return dict(E2_minus_E1=np.log(pi[0] / pi[1]), E3_minus_E1=np.log(pi[0] / pi[2]), b12=b12, b13=b13,
                    E31c=E31c, E21c=E21c, E32=E32, E23=E23, F_SW=np.tanh((E32 - E31c) / 2), F_WS=np.tanh((E23 - E21c) / 2))

# ---------------------------------------------------------------------------------------------------------
# 14. Symbols on an arbitrary time grid (PRESPEC 20m)
# ---------------------------------------------------------------------------------------------------------
def sqw_symbols_bin(x, bin_s=60, theta=0.5, light_dark=True, full_days=True, unknown_wake=-1):
    """As sqw_symbols but with bins of bin_s seconds (bin_s must divide 86400): S if >= bin_s/2 seconds of the bin
    lie in a sleep bout; else Q if the bin's active fraction < theta * r_f (r_f from the 1-s data, per phase),
    else W. bin_s = 60 reproduces sqw_symbols exactly. Returns (y, refs)."""
    x = np.asarray(x); nb = len(x) // bin_s; x = x[:nb * bin_s]
    sl = np.zeros(len(x), bool); bouts = sleep_bouts(x, True)
    for s, n in bouts: sl[s:s+n] = True
    ssec = sl.reshape(nb, bin_s).sum(1); act = (x > 0).reshape(nb, bin_s).mean(1)
    ev = wake_bout_events(x, bouts)
    refs = ({'light': wake_reference(ev, 0), 'dark': wake_reference(ev, 1)} if light_dark else {'whole': wake_reference(ev, None)})
    y = np.full(nb, W); y[ssec * 2 >= bin_s] = S
    for m in np.flatnonzero(ssec * 2 < bin_s):
        r = refs['whole'] if not light_dark else (refs['light'] if (m * bin_s) % DAY < HALF else refs['dark'])
        if np.isnan(r): y[m] = unknown_wake
        elif act[m] < theta * r: y[m] = Q
    if full_days: y = y[:(len(y) // (DAY // bin_s)) * (DAY // bin_s)]
    return y, refs

def phases_bin(T, bin_s=60, scheme='12h'):
    """Phase per step on a bin_s grid (transition INTO step t uses P[ph[t]])."""
    spd = DAY // bin_s; t = np.arange(T) % spd
    if scheme == '12h': return (t >= spd // 2).astype(int)
    if scheme == 'none': return np.zeros(T, int)
    raise ValueError(scheme)

# ---------------------------------------------------------------------------------------------------------
# 15. PLAN A: the observed S/Q/W chain and model-free probes (PRESPEC 20n, 20o, 20p, 20q, 20s-A)
#     Copied from the scripts that produced the ledger numbers (sqw13-sqw17); verify_ref.py checks equality.
# ---------------------------------------------------------------------------------------------------------
def mc_fit(y, ph, nP, tol=1e-10, maxit=500):
    """First-order Markov chain on the symbols = the HMM with B fixed to the identity (Q, W, S = states 0, 1, 2),
    fitted by em_step so that missing minutes (-1) are marginalised exactly. Start: transition counts + 1 per phase of
    the arrival minute; stop when the relative log-likelihood gain < tol. On complete data equals the count MLE.
    Returns (loglik, pi, P (nP x 3 x 3)). PRESPEC 20n."""
    y = np.asarray(y); I3 = np.eye(3); P = np.full((nP, 3, 3), 1.0)
    for t in range(1, len(y)):
        if y[t-1] >= 0 and y[t] >= 0: P[ph[t], y[t-1], y[t]] += 1
    P /= P.sum(2, keepdims=True); pi = np.full(3, 1/3); old = -np.inf
    for it in range(maxit):
        pi, P, B, ll = em_step(pi, P, I3, y, ph)
        occ = [v for v in range(3) if np.any(y == v)]; assert np.allclose(B[occ], I3[occ], atol=1e-12)
        if ll - old < tol * abs(ll): break
        old = ll
    return loglik(pi, P, I3, y, ph), pi, P

B_PAIR = np.zeros((9, 3)); B_PAIR[np.arange(9), np.arange(9) % 3] = 1          # pair state k = 3*y_{t-1} + y_t emits y_t
ALLOW_PAIR = np.array([[k % 3 == l // 3 for l in range(9)] for k in range(9)])    # (a,b) -> (b,c) only

def mc2_fit(y, ph, nP, tol=1e-10, maxit=1000):
    """Second-order chain P(y_t | y_t-2, y_t-1) as a 9-state chain on pairs with fixed deterministic emissions and
    structural zeros; missing minutes marginalised. Returns (loglik, pi (9), P (nP x 9 x 9)). PRESPEC 20o."""
    y = np.asarray(y); P = np.where(ALLOW_PAIR, 1.0, 0.0)[None].repeat(nP, 0)
    for t in range(2, len(y)):
        if min(y[t-2], y[t-1], y[t]) >= 0: P[ph[t], 3*y[t-2] + y[t-1], 3*y[t-1] + y[t]] += 1
    P /= P.sum(2, keepdims=True); pi = np.full(9, 1/9); old = -np.inf
    for it in range(maxit):
        pi, P, _, ll = em_step(pi, P, B_PAIR, y, ph); assert np.all(P[:, ~ALLOW_PAIR] == 0)
        if ll - old < tol * abs(ll): break
        old = ll
    return loglik(pi, P, B_PAIR, y, ph), pi, P

def sim_rows(P, allow=None):
    """FOR SIMULATION ONLY (20p amendment B): a row with zero expected count is replaced by its count-start row,
    uniform over the allowed transitions, so a simulated record cannot get stuck in a state the data never leave."""
    P = np.array(P, float); A = np.ones(P.shape[1:], bool) if allow is None else allow
    for p in range(P.shape[0]):
        for k in range(P.shape[1]):
            if P[p, k].sum() < 0.5: P[p, k] = A[k] / A[k].sum()
    return P

def block_kl(y, ph, nph, qw=False, nmax=4):
    """D_n (n = 1..nmax) per phase: KL divergence between n-minute blocks and their time reverses,
    D_n = sum_b p(b) ln[p(b)/p(reverse b)], pseudocount 0.5 on every block of the alphabet ({Q,W,S}, or {Q,W} if qw).
    Blocks lie within one phase and contain no missing minute. d_n = D_n - D_{n-1} (x 60 for nats/h); for a
    stationary first-order chain d_3 = d_2 (PRESPEC 20p). Returns (nph x nmax), column n-1 = D_n, D_1 = 0."""
    y = np.asarray(y); ph = np.asarray(ph); T = len(y); out = np.zeros((nph, nmax))
    ok = (y >= 0) & ((y < 2) if qw else True)
    for n in range(2, nmax + 1):
        m = T - n + 1; good = np.ones(m, bool); code = np.zeros(m, int); rc = np.zeros(m, int)
        for k in range(n):
            good &= ok[k:k+m] & (ph[k:k+m] == ph[:m]); yk = np.clip(y[k:k+m], 0, 2)
            code = code * 3 + yk; rc += yk * 3**k
        allc = np.arange(3**n); digits = np.array([[(c // 3**k) % 3 for k in range(n)] for c in allc])
        alph = (digits < 2).all(1) if qw else np.ones(3**n, bool)
        rev = np.array([sum(digits[c, n-1-k] * 3**k for k in range(n)) for c in allc])
        for p in range(nph):
            c = np.bincount(code[good & (ph[:m] == p)], minlength=3**n).astype(float)
            pr = np.where(alph, c + 0.5, 0.0); pr /= pr.sum(); s = alph
            out[p, n-1] = float(np.sum(pr[s] * np.log(pr[s] / pr[rev[s]])))
    return out

def block_kl_stats(y, ph, nph):
    """Per phase: d2, d3, d4 (nats/h), Delta3 = d3 - d2, Delta4 = d4 - d3, and Delta3 on wake-only blocks (D3QW).
    The excess E of PRESPEC 20p is Delta3 minus the median Delta3 of records simulated from the fly's own mc_fit."""
    D = block_kl(y, ph, nph); Dq = block_kl(y, ph, nph, True); res = []
    for p in range(nph):
        d = np.diff(D[p]) * 60; dq = np.diff(Dq[p]) * 60
        res.append(dict(d2=d[0], d3=d[1], d4=d[2], D3=d[1] - d[0], D4=d[2] - d[1], D3QW=dq[1] - dq[0]))
    return res

D_CYCLE = np.zeros((3, 3))
for _a, _b in ((W, Q), (Q, S), (S, W)): D_CYCLE[_a, _b] = 1; D_CYCLE[_b, _a] = -1   # forward cycle W->Q->S->W = +1

def tur_estimate(y, ph, nph, taus=(180, 120, 360)):
    """Per phase: J (counted forward-cycle increments per hour) and, per block length tau (minutes), D = Var(N_tau)/(2 tau)
    (per hour) from non-overlapping blocks lying in one phase with no missing minute (>= 4 blocks), and the TUR bound
    B = J^2 / D (k_B per hour). PRESPEC 20q (amendment C: B = J^2/D, not 2 J^2/D)."""
    y = np.asarray(y); ph = np.asarray(ph); T = len(y)
    val = np.zeros(T, bool); val[1:] = (y[1:] >= 0) & (y[:-1] >= 0)
    inc = np.zeros(T); inc[1:] = np.where(val[1:], D_CYCLE[np.clip(y[:-1], 0, 2), np.clip(y[1:], 0, 2)], 0)
    seg = np.flatnonzero(np.diff(ph) != 0) + 1; starts = np.r_[0, seg]; ends = np.r_[seg, T]; out = []
    for p in range(nph):
        m = (ph == p) & val; hours = m.sum() / 60; J = inc[m].sum() / hours if hours > 0 else np.nan; r = dict(J=J)
        for tau in taus:
            Nb = []
            for s, e in zip(starts, ends):
                if ph[s] != p: continue
                for b0 in range(s, e - tau + 1, tau):
                    if b0 == 0 or not val[b0:b0+tau].all(): continue
                    Nb.append(inc[b0:b0+tau].sum())
            D = np.var(Nb, ddof=1) / (2 * tau / 60) if len(Nb) >= 4 else np.nan
            r[f'D{tau}'] = D; r[f'B{tau}'] = J**2 / D if D and np.isfinite(D) and D > 0 else np.nan; r[f'n{tau}'] = len(Nb)
        out.append(r)
    return out

def _log_perron(M):
    return float(np.log(np.max(np.linalg.eigvals(M).real)))

def tur_exact(T, d=D_CYCLE, h=1e-3):
    """Exact J, D, B = J^2/D (per hour) of the current with increments d on the per-minute chain T, from the tilted matrix
    T_ij exp(s d_ij): lambda(s) = log Perron root, J = lambda'(0), D = lambda''(0)/2 (central differences)."""
    lp = _log_perron(T * np.exp(h * d)); lm = _log_perron(T * np.exp(-h * d)); l0 = _log_perron(T)
    J = (lp - lm) / (2 * h) * 60; D = (lp - 2 * l0 + lm) / h**2 / 2 * 60; return J, D, (J**2 / D if D > 0 else np.nan)

def tur_exact_ct(T, d=D_CYCLE, h=1e-3):
    """The same for the continuous-time generator T - I (off-diagonal tilt). The steady-state TUR B <= sigma holds here;
    the 1-min discrete chain can exceed it slightly (kappa = B_discrete / B_continuous, PRESPEC 20q amendment C)."""
    def lam(s):
        M = T * np.exp(s * d); np.fill_diagonal(M, np.diag(T)); return float(np.max(np.linalg.eigvals(M).real)) - 1.0
    lp, lm = lam(h), lam(-h); J = (lp - lm) / (2 * h) * 60; D = (lp + lm) / h**2 / 2 * 60; return J, D, (J**2 / D if D > 0 else np.nan)

def hmm_observed_chain(P, B):
    """9-state chain on (hidden z, symbol y): T[(z,y),(z',y')] = P[z,z'] B[z',y'], with the cycle increments on the symbols,
    so that tur_exact(T, d) gives the HMM's exact J, D, B for the OBSERVED current."""
    T = np.einsum('ab,bc->abc', P, B)[:, None, :, :].repeat(3, 1).reshape(9, 9); d = np.tile(D_CYCLE, (3, 3)); return T, d

def routes_observed(y, ph, nph):
    """Route counts (ROUTES order) on the observed symbols, each run of non-missing minutes separately, by phase of the
    arrival minute (PRESPEC 20s-A). Equals route_counts(y, ph, nph) on a complete record."""
    y = np.asarray(y); c = np.zeros((nph, 4)); ok = np.r_[False, y >= 0, False]; d = np.diff(ok.astype(int))
    for s, e in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)):
        for r, t in route_events(y[s:e]): c[ph[s + t], ROUTES.index(r)] += 1
    return c

def observable_quantities(y, ph, nph, P_chain):
    """PLAN A per phase (PRESPEC 20s-A): c_obs, count J (cycles/h), count EPR (rule-of-three bound) and whether the rule
    was used, chain J on mc_fit's matrix; chain affinity A / EPR and landscape on the rule-of-three matrix (section 16;
    2026-11-02 amendment, finite where an edge is one-way), predicted routing, counted routing Xn, Yn on the observed symbols. hours = minutes in the phase / 60."""
    N = observable_pair_counts(y, ph, nph); c = c_obs_minutes(y, ph, nph); X, Y = xy_counts(routes_observed(y, ph, nph)); out = []
    for p in range(nph):
        h = np.sum(np.asarray(ph) == p) / 60; e = model_energetics_r3(P_chain[p], N[p]); ls = landscape(rule_of_three_matrix(P_chain[p], N[p])[0]); epr, sub = epr_count_bound(N[p], h)
        J_fit = model_energetics(P_chain[p])['J']
        out.append(dict(c_obs=float(c[p]), J=float(j_cycle_count(N[p], h)), EPR=float(epr), rule_of_three=bool(sub),
                        J_chain=float(J_fit), A_chain=float(e['A']), EPR_chain=float(e['EPR']),
                        E21=float(ls['E2_minus_E1']), E31=float(ls['E3_minus_E1']), F_SW=float(ls['F_SW']), F_WS=float(ls['F_WS']),
                        Xn_cnt=float(X[p]), Yn_cnt=float(Y[p])))
    return out


# ---------------------------------------------------------------------------------------------------------
# 16. Rule of three for model-based quantities (2026-11-02, Sheyum's decision; replaces 'infinite')
# ---------------------------------------------------------------------------------------------------------
def rule_of_three_matrix(P, N, zero=1e-10):
    """Transition matrix with every ONE-WAY edge filled by the rule of three. P: fitted 3x3 matrix (one phase);
    N: off-diagonal transition counts of the same phase (N[i,j], i != j). For an edge with P[a,b] < zero and
    P[b,a] >= zero, the unseen count is set to min(3, N[b,a]). Row a's total number of steps is recovered once,
    from the ORIGINAL row, as n_a = sum_{j != a} N[a,j] / (1 - P[a,a]); all of the row's fills are made against
    that original row and the row is renormalized ONCE, so the result does not depend on the order in which the
    one-way edges are visited (2026-11-02 revision; the first version renormalized after each fill).
    A row with no observed exit (sum_{j != a} N[a,j] = 0) is left unchanged: the state is never left, there is no
    cycle through it, and such a fly-phase is excluded as no-cycle. Pairs unseen in BOTH directions are left at
    zero. Returns (P_r3, changed)."""
    P0 = np.array(P, float); P = P0.copy(); N = np.asarray(N, float); Z = P0 < zero; changed = False
    for a in range(3):
        off = sum(N[a, j] for j in range(3) if j != a)
        if off <= 0 or P0[a, a] >= 1.0: continue
        fills = {b: min(3.0, N[b, a]) for b in range(3)
                 if b != a and Z[a, b] and not Z[b, a] and N[b, a] > 0}
        if not fills: continue
        row = P0[a] * (off / (1.0 - P0[a, a]))
        for b, c in fills.items(): row[b] = c
        P[a] = row / row.sum(); changed = True
    return P, changed

def model_energetics_r3(P, N, zero=1e-10, per_hour=60.0):
    """model_energetics on the rule-of-three matrix: finite 95% lower bounds on |A| and EPR where an edge is
    one-way; identical to model_energetics otherwise. Adds 'r3' (bool)."""
    Pr, ch = rule_of_three_matrix(P, N, zero)
    e = model_energetics(Pr, zero, per_hour); e['r3'] = ch
    # No-cycle convention: a chain whose stationary flux vanishes (e.g. a state that is never left, so the stationary
    # mass sits in it) carries no cycle; EPR is 0 and the affinity is undefined, and the fly-phase is excluded.
    # Applied only where the value would otherwise be non-finite; finite round-off values are left as computed.
    if abs(e['J']) < 1e-12 and not (np.isfinite(e['EPR']) and np.isfinite(e['A'])): e['EPR'] = 0.0; e['A'] = np.nan
    return e
