"""Port of the lab's raw-beam -> 1-s activity conversion: PreProcess.m + multibeam_locomotion.m (binner = 5).
Operates on ALL data lines of a raw week file (the -1 lines included, as dlmread reads them).
PreProcess per column: leading 0/1 values set to the first other value; then, row by row, a value equal to 0, -1,
or 1 with the previous (processed) value != 2 is replaced by the previous value. multibeam_locomotion: displacement
|A[k+1] - A[k]| per row (last row copies the one before), summed in consecutive blocks of 5 rows = 1 s.
(PreProcess's while-loop has an operator-precedence quirk with the row bound; it only matters for a column made
entirely of 0/1, which raises an index error in MATLAB; here such a column is left unchanged.)"""
import numpy as np
from numba import njit


@njit(cache=True)
def _preprocess_col(a):
    a = a.copy(); n = a.shape[0]
    i = 0
    while i < n and (a[i] == 0 or a[i] == 1):
        i += 1
    if i < n:
        for k in range(i + 1):
            a[k] = a[i]
    for k in range(1, n):
        if a[k] == 0 or a[k] == -1 or (a[k] == 1 and a[k - 1] != 2):
            a[k] = a[k - 1]
    return a


def lab_activity_per_second(Xall):
    """Xall: all data lines (rows) x channels, -1 lines included. Returns counts per second (beams travelled)."""
    out = []
    for j in range(Xall.shape[1]):
        a = _preprocess_col(Xall[:, j].astype(np.int64))
        d = np.empty_like(a); d[:-1] = np.abs(a[1:] - a[:-1]); d[-1] = d[-2]
        n = (len(d) // 5) * 5
        out.append(d[:n].reshape(-1, 5).sum(1))
    return np.stack(out, 1)


def sleeptrace(activity, dt=1.0):
    """Port of sleeptrace.m: 1 = wake, 0 = sleep; inactivity runs (activity <= 0) of N samples with N*dt >= 300 s."""
    act = np.asarray(activity); out = np.ones(len(act))
    inact = (act <= 0).astype(np.int8); tr = np.diff(np.r_[0, inact, 0])
    st = np.nonzero(tr == 1)[0]; en = np.nonzero(tr == -1)[0] - 1
    for s, e in zip(st, en):
        if (e - s + 1) * dt >= 300:
            out[s:e + 1] = 0
    return out


def ironout(wave, max_wake_s=1.0, dt=1.0):
    """Port of ironout.m: interior wake runs (1s) of <= floor(max_wake_s/dt) samples become sleep (0); one pass on
    the original wave. Returns (wave_out, n_bridged, n_samples)."""
    w = np.asarray(wave, float).copy(); mb = int(np.floor(max_wake_s / dt))
    if mb < 1:
        return w, 0, 0
    iw = (w == 1).astype(np.int8); tr = np.diff(np.r_[0, iw, 0])
    st = np.nonzero(tr == 1)[0]; en = np.nonzero(tr == -1)[0] - 1; ln = en - st + 1
    q = (st > 0) & (en < len(w) - 1) & (ln <= mb)
    for s, e in zip(st[q], en[q]):
        w[s:e + 1] = 0
    return w, int(q.sum()), int(ln[q].sum())


def beam_counts_1s(x, t, keep_beyond_beam1=True):
    """beam_hmm mode, the reference for RawToActivity1s.m: x = positions on TIMED lines only (already spike-filtered
    and hub-repaired if required). Counts per second = beams travelled (sum of |dx|) over position changes whose
    sample time falls in [s, s+1); seconds 0 .. floor(t[-1]). If keep_beyond_beam1 is False, a reading of 0 is
    replaced by the previous reading first (the PreProcess convention for 0 only)."""
    x = np.asarray(x).astype(np.int64).copy()
    if not keep_beyond_beam1:
        for k in range(1, len(x)):
            if x[k] == 0:
                x[k] = x[k - 1]
    d = np.r_[0, np.abs(np.diff(x))]; sec = np.floor(np.asarray(t)).astype(int)
    n = int(np.floor(t[-1])) + 1
    return np.bincount(sec, weights=d, minlength=n)[:n].astype(np.int64)


def loading_check(x, t, hours=2.0, start_hour=0.0):
    """Python twin of LoadingCheck.m for one channel: x positions on timed lines, t times (s). The window is
    [start_hour, start_hour + hours) after the first reading (LoadingCheck's 'StartHour' and 'Hours')."""
    t = np.asarray(t); rel = t - t[0]
    m = (rel >= start_hour * 3600) & (rel < (start_hour + hours) * 3600); x = np.asarray(x)[m].astype(np.int64)
    d = np.diff(x); nchg = int((d != 0).sum()); big = np.nonzero(np.abs(d) >= 4)[0] + 1
    span = int(x.max() - x.min()) if len(x) else 0
    hub = None; share = 0.0; ret = 0.0
    if len(big):
        v, c = np.unique(x[big], return_counts=True); hub = int(v[np.argmax(c)]); share = float(c.max() / len(big))
        rs = []
        for i in big[x[big] == hub]:
            k = i
            while k + 1 < len(x) and x[k + 1] == hub: k += 1
            rs.append(k + 1 < len(x) and x[k + 1] == x[i - 1])
        ret = float(np.mean(rs)) if rs else 0.0
    flags = []
    if len(big) >= 10 and share >= 0.30 and ret >= 0.20:
        flags.append(f"BLOCKED BEAM {hub} ({'end' if hub in (1, 2, 15, 16) else 'mid-tube'})")
    if nchg >= 50 and span < 8: flags.append('NARROW RANGE')
    if len(x):
        vals, cnts = np.unique(x, return_counts=True); fv = int(vals[np.argmax(cnts)]); fs = cnts.max() / len(x)
        if fs >= 0.99:
            flags.append(f"FROZEN AT {'END BEAM' if fv in (1, 2, 15, 16) else 'BEAM'} {fv}" if fv != 0 else 'FROZEN AT 0 (no beam broken)')
    if nchg < 20: flags.append('LOW ACTIVITY')
    to0 = int(((x[1:] == 0) & (x[:-1] != 0)).sum()) if len(x) else 0      # entries into 0 (information only)
    return dict(changes=nchg, big_jumps=int(len(big)), hub=hub, hub_share=share, hub_returns=ret, span=span, flags=flags,
                entries_to_zero=to0)


def loading_change(prev_flags, now_flags):
    """LoadingCheck's 'Previous' label for one channel; flags are lists (empty = OK)."""
    b, a = bool(prev_flags), bool(now_flags)
    return 'RESOLVED' if b and not a else 'STILL FLAGGED' if b and a else 'NEW FLAG' if a else 'OK'


def at_zero_1s(x, t):
    """Reference for RawToActivity1s out.AtZero: per second s = 0 .. floor(t[-1]), the share of the timed readings
    with time in [s, s+1) that read 0 (NaN for a second without a reading). x: positions after the same repairs as
    the activity (spike filter + hub repair when Repair is on)."""
    x = np.asarray(x); sec = np.floor(np.asarray(t)).astype(int); n = int(np.floor(t[-1])) + 1
    tot = np.bincount(sec, minlength=n)[:n].astype(float); z = np.bincount(sec, weights=(x == 0), minlength=n)[:n]
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(tot > 0, z / tot, np.nan)


def entries_to_zero_per_day(x, t):
    """Reference for RawToActivity1s QC EntriesToZeroPerDay."""
    x = np.asarray(x); return float(((x[1:] == 0) & (x[:-1] != 0)).sum() / ((t[-1] - t[0]) / 86400))
