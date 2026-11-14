"""P1-E: behavioural ZT per file, first-cycle / transfer exclusions, death rule across files (PRESPEC P1-E)."""
import numpy as np
import pipeline as pl

BIN = 600  # s


def preprocess(X, t):
    """Spike filter + P1-C hub repair per channel; returns positions, flags per channel (P1-D)."""
    P = np.empty_like(X); info = []
    for j in range(X.shape[1]):
        x, _ = pl.spike_filter(X[:, j]); xr, hub = pl.hub_repair(x)
        d0 = np.abs(np.diff(x.astype(int))); c0 = d0 != 0
        d1 = np.abs(np.diff(xr.astype(int))); c1 = d1 != 0
        live = (np.diff(xr[t < 86400]) != 0).sum() >= 100
        f0 = float((d0[c0] >= 4).mean()) if c0.any() else 1.0
        f1 = float((d1[c1] >= 4).mean()) if c1.any() else 1.0
        P[:, j] = xr
        info.append(dict(live=bool(live), frac_big_raw=f0, frac_big_repaired=f1, hub=hub,
                         included=bool(live and f1 <= 0.05), never_artifact=bool(live and f0 <= 0.05)))
    return P, info


def phase_stamp(P, t, info):
    """Lights-on time in hours into the file (P1-E rule 1), with the guard."""
    inc = [j for j, i in enumerate(info) if i['included']]
    last_whole = int(t[-1] // 86400) * 86400
    keep = (t >= 86400) & (t < last_whole)
    b = ((t[keep] % 86400) // BIN).astype(int); nb = 86400 // BIN
    prof = []
    for j in inc:
        c = np.r_[False, P[1:, j] != P[:-1, j]][keep].astype(float)
        p = np.bincount(b, weights=c, minlength=nb) / np.maximum(np.bincount(b, minlength=nb), 1)
        prof.append(p / max(p.mean(), 1e-12))
    pr = np.mean(prof, 0); rise = pr - np.roll(pr, 1)
    off_h = int(np.argmax(rise)) * BIN / 3600
    on_h = (off_h - 12) % 24
    dist = min(on_h, 24 - on_h)
    return dict(lights_on_h=float(on_h), lights_off_h=float(off_h), n_flies=len(inc), bin_min=BIN / 60,
                flagged=bool(dist > 1.5))


def change_times(x, t):
    return t[1:][x[1:] != x[:-1]]


def death_times(files_changes, files_tend):
    """files_changes: list over files of change-time arrays (file-local seconds) for one fly; files_tend: file
    lengths (s). Returns (file index, death time in that file) or None; 'possible' if in the last file."""
    nf = len(files_changes)
    def busy(tc):
        # some 24-h window (starting at a change) holds >= 100 changes
        return len(tc) >= 100 and (np.searchsorted(tc, tc + 86400) - np.arange(len(tc))).max() >= 100
    def active_after(fi, s):
        tc = files_changes[fi]
        if busy(tc[tc > s]):
            return True
        return any(busy(files_changes[fj]) for fj in range(fi + 1, nf))
    for fi in range(nf):
        tc = files_changes[fi]
        edges = np.r_[0.0, tc, files_tend[fi]]
        gaps = np.diff(edges)
        for k in np.nonzero(gaps >= 43200)[0]:
            s = edges[k]
            if not active_after(fi, s + gaps[k] - 1):
                return dict(file=fi, t=float(s), possible=(fi == nf - 1))
    return None


def valid_mask(t, week, death=None, file_index=None):
    """P1-E rules 2-4 for one fly in one file."""
    v = t >= 86400                                  # first cycle
    if week == 1:
        v &= t < 150 * 3600                         # week 1: tubes out of the monitor from ~153 h
    if death is not None:
        if file_index > death['file']:
            v[:] = False
        elif file_index == death['file']:
            v &= t < death['t'] - 86400
    return v
