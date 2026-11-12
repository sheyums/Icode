"""beam_hmm pilot pipeline (PRESPEC P1). Python prototype; MATLAB port later with golden tests.

Steps: load week file -> drop -1 lines -> spike filter -> inclusion (empty / artifact / handling / death)
-> 1-s activity -> sleep bouts and S/Q/W minutes (scoring copied from hmm_audit/sqw_ref.py @ 5ea02e7)
-> per-minute 0.4-s descriptors.
Positions: 0 = beyond beam 1, 1..16 = beams. Time step 0.4 s. Row 0 of a file = lights-on (ZT0).
"""
import numpy as np

DT = 0.4
DAY, HALF, MIN_SLEEP = 86400, 43200, 300
Q, W, S = 0, 1, 2
CH = list(range(16, 32))           # 0-based columns of channels 17-32


# ----------------------------------------------------------------------------------- loading
def load_week(npz):
    d = np.load(npz); X, T = d['X'], d['T']
    ok = ~np.isnan(T)
    assert ((X[~ok] == -1).all()), 'non -1 values on an untimed line'
    return X[ok][:, CH].astype(np.int16), T[ok]


def spike_filter(x, thr=4):
    """One pass: B with equal neighbours A and |B-A| >= thr is set to A. Returns (filtered, n_fixed)."""
    x = x.copy(); a, b, c = x[:-2], x[1:-1], x[2:]
    m = (a == c) & (np.abs(b - a) >= thr)
    idx = np.nonzero(m)[0] + 1
    x[idx] = x[idx - 1]
    return x, len(idx)


# ----------------------------------------------------------------------------------- inclusion
def change_counts_10min(P, t):
    w = (t // 600).astype(int); n = w.max() + 1
    M = np.vstack([np.zeros(P.shape[1], bool), P[1:] != P[:-1]])
    return np.array([np.bincount(w, weights=M[:, j], minlength=n) for j in range(P.shape[1])])


def handling_start(C, include):
    """First 10-min window where >= 50% of included channels exceed their own 99th percentile. Returns
    window index or None."""
    Ci = C[include]
    thr = np.percentile(Ci, 99, axis=1, keepdims=True)
    frac = (Ci > thr).mean(0)
    hit = np.nonzero(frac >= 0.5)[0]
    return int(hit[0]) if len(hit) else None


def longest_silences(x, t, valid):
    """Silences (no change) measured within valid samples. Returns list of (start_t, length_s)."""
    tv, xv = t[valid], x[valid]
    ch = np.nonzero(np.diff(xv))[0] + 1
    edges = np.r_[tv[0], tv[ch], tv[-1]]
    return [(edges[i], edges[i + 1] - edges[i]) for i in range(len(edges) - 1)]


# ----------------------------------------------------------------------------------- scoring (sqw_ref)
def sleep_bouts(a, min_len=MIN_SLEEP):
    """a: 1-s activity 0/1 with excluded seconds = -1. Bridging per lab rule (sqw_ref). Runs of inactivity are
    broken by excluded seconds. Returns list of (start, length, censored) where censored = touches an
    excluded second or the record edge."""
    a = np.asarray(a).copy(); n = len(a)
    z0 = (a == 0).astype(int)
    ch = np.r_[0, np.flatnonzero(np.diff(z0)) + 1, n]
    runs = [(s, e - s, z0[s]) for s, e in zip(ch[:-1], ch[1:])]
    for i in range(1, len(runs) - 1):
        s_, n_, zz = runs[i]
        if (zz == 0 and n_ == 1 and a[s_] == 1 and runs[i-1][2] == 1 and runs[i+1][2] == 1
                and runs[i-1][1] >= min_len and runs[i+1][1] >= min_len):
            a[s_] = 0
    z = (a == 0).astype(int)
    ch = np.r_[0, np.flatnonzero(np.diff(z)) + 1, n]
    out = []
    for s, e in zip(ch[:-1], ch[1:]):
        if z[s] == 1 and e - s >= min_len:
            cens = (s == 0) or (e == n) or a[s - 1] == -1 or (e < n and a[e] == -1)
            out.append((int(s), int(e - s), bool(cens)))
    return out


def wake_reference_events(a, bouts, min_wake=600):
    """sqw_ref.wake_bout_events, restricted to wake bouts with no excluded second."""
    ev = []
    for (a1, n1, _), (a2, _, _) in zip(bouts[:-1], bouts[1:]):
        w0, w1 = a1 + n1, a2
        if w1 - w0 < min_wake or (a[w0:w1] == -1).any():
            continue
        mid = (w0 + w1) // 2
        ev.append(dict(ref=(a[mid-150:mid+150] == 1).mean(),
                       wake=int((w0 % DAY) >= HALF), onset=int((w1 % DAY) >= HALF)))
    return ev


def sqw_minutes(a, theta=0.5):
    """S/Q/W per minute; -1 = minute containing an excluded second."""
    a = np.asarray(a); nm = len(a) // 60; a = a[:nm * 60]
    bouts = sleep_bouts(a)
    sl = np.zeros(len(a), bool)
    for s, n, _ in bouts:
        sl[s:s+n] = True
    ev = wake_reference_events(a, bouts)
    refs = {}
    for ph in (0, 1):
        r = [e['ref'] for e in ev if e['wake'] == ph or e['onset'] == ph]
        refs[ph] = float(np.median(r)) if len(r) >= 3 else np.nan
    A = a.reshape(nm, 60); ssec = sl.reshape(nm, 60).sum(1); act = (A == 1).mean(1)
    y = np.full(nm, W); y[ssec >= 30] = S
    for m in np.flatnonzero(ssec < 30):
        r = refs[int((m * 60) % DAY >= HALF)]
        if np.isnan(r): y[m] = -2
        elif act[m] < theta * r: y[m] = Q
    y[(A == -1).any(1)] = -1
    return y, refs, bouts


# ----------------------------------------------------------------------------------- descriptors
def activity_1s(x, t, valid):
    """Second s active if a change occurs at a sample in [s, s+1). Excluded seconds = -1."""
    n = int(np.floor(t[-1])) + 1
    sec = np.floor(t).astype(int)
    chg = np.r_[False, x[1:] != x[:-1]]
    a = np.zeros(n, np.int8)
    np.maximum.at(a, sec[chg & valid], 1)
    bad = np.zeros(n, bool); bad[sec[~valid]] = True
    a[bad] = -1
    return a


def minute_descriptors(x, t, valid):
    """Per-minute 0.4-s descriptors (NaN for minutes with any excluded sample)."""
    m = (t // 60).astype(int); nm = m.max() + 1
    dx = np.r_[0, np.diff(x)]; chg = dx != 0
    flick = np.zeros(len(x), bool)
    i = np.nonzero((x[:-2] == x[2:]) & (np.abs(x[1:-1] - x[:-2]) == 1))[0] + 1
    flick[i] = True; flick[i + 1] = True        # the excursion and the return
    def per_min(v): return np.bincount(m, weights=v.astype(float), minlength=nm)
    nsamp = per_min(np.ones(len(x)))
    out = dict(changes=per_min(chg), dist=per_min(np.abs(dx)),
               at0=per_min(x == 0) / nsamp, ends=per_min((x <= 1) | (x >= 15)) / nsamp,
               flick=per_min(flick & chg))
    # seconds whose only changes are edge flicker
    sec = np.floor(t).astype(int); ns = sec.max() + 1
    c_all = np.bincount(sec, weights=chg.astype(float), minlength=ns)
    c_fl = np.bincount(sec, weights=(chg & flick).astype(float), minlength=ns)
    act = c_all > 0; onlyfl = act & (c_all == c_fl)
    sm = np.arange(ns) // 60
    out['act_s'] = np.bincount(sm, weights=act, minlength=nm)[:nm]
    out['flickonly_s'] = np.bincount(sm, weights=onlyfl, minlength=nm)[:nm]
    # longest still run inside the minute (s)
    ls = np.zeros(nm)
    run = 0.0
    for k in range(len(x)):
        run = run + DT if not chg[k] else 0.0
        if run > ls[m[k]]: ls[m[k]] = run
    out['longest_still'] = np.minimum(ls, 60)
    bad = np.bincount(m, weights=(~valid).astype(float), minlength=nm) > 0
    for k in out: out[k] = np.where(bad[:len(out[k])], np.nan, out[k][:nm])
    return out
