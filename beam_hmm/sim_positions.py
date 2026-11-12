"""Synthetic 0.4-s positions with a known hidden state (D1-A item 7) and the 10-s symbol coder (D1)."""
import numpy as np
from numba import njit

DT = 0.4
SLEEP, QUIET, TRAVEL, BOUNDARY = 0, 1, 2, 3
DWELL_MIN = {0: np.array([20, 3, 0.67, 3.0]), 1: np.array([40, 2, 0.5, 3.0])}   # light, dark
JUMP = np.array([[0, .6, .3, .1],
                 [.4, 0, .4, .2],
                 [.1, .6, 0, .3],
                 [.5, .3, .2, 0]])


def truth(scenario):
    """Per-phase per-sample exit probabilities (2 x S) and jump matrix (S x S)."""
    S = 4 if scenario == 'G1' else 3
    J = JUMP[:S, :S].copy(); J /= J.sum(1, keepdims=True)
    exitp = np.array([DT / (60 * DWELL_MIN[p][:S]) for p in (0, 1)])
    return exitp, J


@njit(cache=True)
def _sim(n, exitp, J, mod_amp, u, x0):
    S = J.shape[0]
    state = np.empty(n, np.int8); pos = np.empty(n, np.int16); read = np.empty(n, np.int16)
    s = 0; p = x0; d = 1; flick = 0; last = x0
    cumJ = np.empty((S, S))
    for i in range(S):
        acc = 0.0
        for j in range(S):
            acc += J[i, j]; cumJ[i, j] = acc
    k = 0
    for t in range(n):
        tt = t * 0.4
        ph = 0 if (tt % 86400.0) < 43200.0 else 1
        e = exitp[ph, s]
        if s == 0 and mod_amp > 0:
            e = e * (1.0 + mod_amp * np.cos(2 * np.pi * tt / 21600.0))
        if u[k] < e:
            r = u[k+1]; ns = S - 1
            for j in range(S):
                if r < cumJ[s, j]:
                    ns = j; break
            s = ns
        k += 2
        # movement
        if s == 1:
            if u[k] < 0.03:
                step = 1 if u[k+1] < 0.5 else -1
                p = p + step
        elif s == 2:
            if u[k] < 0.5:
                if u[k+1] > 0.85:
                    d = -d
                step = 2 if u[k+2] < 0.1 else 1
                p = p + d * step
        if p < 0:
            p = -p; d = 1
        if p > 16:
            p = 32 - p; d = -1
        k += 3
        rd = p
        if s == 3 and u[k] < 0.4 / 60.0:
            rd = p + 1 if (u[k+1] < 0.5 and p < 16) or p == 0 else p - 1
        k += 2
        if rd == 15:
            rd = last           # dead beam 15 reads as the previous reading
        state[t] = s; pos[t] = p; read[t] = rd; last = rd
    return state, pos, read


def simulate_fly(scenario, days, seed_seq):
    rng = np.random.default_rng(seed_seq)
    exitp, J = truth(scenario)
    if scenario == 'G2':
        exitp = exitp * np.exp(rng.normal(0, 0.3, size=exitp.shape))
        amp = 0.5
    else:
        amp = 0.0
    n = int(days * 86400 / DT)
    u = rng.random(7 * n + 10)
    state, pos, read = _sim(n, exitp, J, amp, u, int(rng.integers(0, 17)))
    return state, read, np.arange(n) * DT


# ------------------------------------------------------------------------------------------- symbols
STILL, FLICK, LOCAL, TRAVELC = 0, 1, 2, 3
BEYOND1, END, MIDDLE = 0, 1, 2
WIN = 25    # samples per 10-s window


def window_zones(x):
    """Zone of each 10-s window (BEYOND1 / END / MIDDLE), D1 definitions. Place layer of D1-B."""
    x = np.asarray(x).astype(np.int64); n = (len(x) // WIN) * WIN; X = x[:n].reshape(n // WIN, WIN)
    n0 = (X == 0).sum(1)
    mode = np.array([np.bincount(r, minlength=17).argmax() for r in X])
    return np.where(n0 >= 13, BEYOND1, np.where((mode <= 2) | (mode >= 14), END, MIDDLE))


def window_symbols(x, valid=None, merge_flick=False, merge_beyond=False, activity_only=True):
    """10-s symbols; -1 for windows with any invalid sample. x: readings at 0.4 s.
    D1-B (default): activity class only (M = 4). activity_only=False gives D1's activity*3 + zone (M = 12)."""
    x = np.asarray(x).astype(np.int64); n = (len(x) // WIN) * WIN; x = x[:n]
    dx = np.r_[0, np.diff(x)]; chg = dx != 0
    fl = np.zeros(n, bool)
    i = np.nonzero((x[:-2] == x[2:]) & (np.abs(x[1:-1] - x[:-2]) == 1))[0] + 1
    fl[i] = True; fl[np.minimum(i + 1, n - 1)] = True
    W = n // WIN
    C = chg.reshape(W, WIN); F = (chg & fl).reshape(W, WIN); D = np.abs(dx).reshape(W, WIN)
    nc = C.sum(1); nf = F.sum(1); path = D.sum(1)
    act = np.where(nc == 0, STILL, np.where(nc == nf, FLICK, np.where(path <= 3, LOCAL, TRAVELC)))
    if merge_flick:
        act = np.where(act == FLICK, STILL, act)
    if activity_only:
        y = act.copy()
        if valid is not None:
            bad = (~np.asarray(valid)[:n]).reshape(W, WIN).any(1); y[bad] = -1
        return y.astype(np.int64)
    X = x.reshape(W, WIN)
    n0 = (X == 0).sum(1)
    mode = np.array([np.bincount(r, minlength=17).argmax() for r in X])
    zone = np.where(n0 >= 13, BEYOND1, np.where((mode <= 2) | (mode >= 14), END, MIDDLE))
    if merge_beyond:
        zone = np.where(zone == BEYOND1, END, zone)
    y = act * 3 + zone
    if valid is not None:
        bad = (~np.asarray(valid)[:n]).reshape(W, WIN).any(1); y[bad] = -1
    return y.astype(np.int64)
