"""PRESPEC 20u (+ amendment A): reversible-null control of the per0 DD wake-hazard rise.
Stage 0 here: G1 (extractor port vs per0DD_bouts_20r_replication.mat, per fly, exact) and G2 (BoutHazard RiseRatio
port vs 1.416/1.402/1.385/1.435). Usage: python3 -I planB_20u.py stage0 <per0DD data dir>"""
import sys, os, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import sqw_ref as R
REF = os.path.join(ROOT, 'analysis_20r', 'per0DD_bouts_20r_replication.mat')
NB = (12, 16, 20, 24)

def mround(x):
    """MATLAB round: half away from zero."""
    x = np.asarray(x, float); return np.sign(x) * np.floor(np.abs(x) + 0.5)

# ------------------------------------------------------------------ extractor port (amendment A item 3)
def extract(x, trace_len):
    x = np.asarray(x)[:trace_len]; T = len(x)
    sb = R.sleep_bouts(x, bridge=True)                       # (start, length), all >= 300 s, bridged
    wake = [sb[k + 1][0] - (sb[k][0] + sb[k][1]) for k in range(len(sb) - 1)]
    sleep = [n for s, n in sb if s > 0 and s + n < T]
    return np.array(wake, float), np.array(sleep, float)

# ------------------------------------------------------------------ BoutHazard RiseRatio port (verbatim logic)
def build_edges(t, xmin, dt, nbins):
    lo = xmin - dt; hi = np.max(t); first = max(lo, dt / 2)
    e = np.exp(np.linspace(np.log(first), np.log(hi), nbins + 1))
    edges = np.r_[lo, e[1:]]
    return np.unique(mround(edges / dt) * dt)

def rise_ratio(t, xmin=2.0, dt=1.0, nbins=12, min_at_risk=10):
    t = np.asarray(t, float); t = t[np.isfinite(t) & (t >= xmin)]
    edges = build_edges(t, xmin, dt, nbins); nb = len(edges) - 1
    ts = np.sort(t)
    atR = len(ts) - np.searchsorted(ts, edges[:-1], side='right')          # sum(t > edges(j))
    ev = atR - (len(ts) - np.searchsorted(ts, edges[1:], side='right'))    # t > e_j & t <= e_j+1
    sat = ev >= atR; keep = (atR >= min_at_risk) & ~sat
    w = np.diff(edges)[keep]; p = ev[keep] / atR[keep]
    h = -np.log1p(-np.minimum(p, 1)) / w
    hh = np.where(np.isfinite(h), h, np.inf)
    iT = int(np.argmin(hh[:-1])) if len(hh) >= 2 else 0
    if iT < len(hh) - 1:
        tail = hh[iT + 1:].copy(); tail[~np.isfinite(tail)] = -np.inf; iP = iT + 1 + int(np.argmax(tail))
    else:
        iP = iT
    if not np.isfinite(h[iP]): iP = iT
    return h[iP] / h[iT] if h[iT] > 0 else np.inf

# ------------------------------------------------------------------ stage 0
def load_ref():
    import scipy.io as sio
    d = sio.loadmat(REF)
    files = [str(f[0]) if isinstance(f, np.ndarray) else str(f) for f in np.ravel(d['files'])]
    return d['Wake'], d['Sleep'], files, np.ravel(d['traceLen']).astype(int)

def stage0(datadir):
    W, Sl, files, tl = load_ref()
    up = sorted(f for f in glob.glob(os.path.join(datadir, '**', '*'), recursive=True) if os.path.isfile(f))
    ok_all = True; nW = nS = 0; maxsil = 0; rows = []
    for i, name in enumerate(files):
        m = [f for f in up if os.path.basename(f).startswith(name)]
        assert len(m) == 1, (name, m)
        x = R.read_counts(m[0]); assert len(x) >= tl[i]
        w, s = extract(x, tl[i])
        rw = np.sort(W[i][np.isfinite(W[i])]); rs = np.sort(Sl[i][np.isfinite(Sl[i])])
        okw = len(rw) == len(w) and np.array_equal(rw, np.sort(w))
        oks = len(rs) == len(s) and np.array_equal(rs, np.sort(s))
        a = (x[:tl[i]] > 0).astype(int); ch = np.flatnonzero(np.diff(np.r_[1, a, 1]))
        runs0 = [ch[k + 1] - ch[k] for k in range(0, len(ch) - 1, 2)]      # inactive runs
        maxsil = max(maxsil, max(runs0) if runs0 else 0)
        nW += len(w); nS += len(s); ok_all &= okw and oks
        rows.append(dict(fly=i, file=name, ref_wake=len(rw), port_wake=len(w), ref_sleep=len(rs), port_sleep=len(s),
                         wake_identical=bool(okw), sleep_identical=bool(oks)))
        if not (okw and oks):
            print(f"  G1 MISMATCH fly {i} {name}: wake ref {len(rw)} port {len(w)}, sleep ref {len(rs)} port {len(s)}")
    print(f"G1: port wake {nW} (ref {int(np.isfinite(W).sum())}), sleep {nS} (ref {int(np.isfinite(Sl).sum())}); "
          f"flies identical {sum(r['wake_identical'] and r['sleep_identical'] for r in rows)}/{len(rows)} -> {'PASS' if ok_all else 'FAIL'}")
    print(f"     longest inactive run on real traces {maxsil} s = {maxsil/3600:.2f} h (amendment A: must be < 12 h)")
    x = W[np.isfinite(W)]
    rr = [rise_ratio(x, 2.0, 1.0, k) for k in NB]; target = [1.416, 1.402, 1.385, 1.435]
    g2 = all(abs(round(a, 3) - b) < 1e-9 for a, b in zip(rr, target))
    print("G2: RiseRatio " + ' '.join(f"{k} bins {a:.4f}" for k, a in zip(NB, rr)) + f" (target {target}) -> {'PASS' if g2 else 'FAIL'}")
    json.dump(dict(G1=bool(ok_all), G2=bool(g2), rr=rr, maxsil=int(maxsil), flies=rows),
              open(os.path.join(HERE, 'planB_20u_stage0.json'), 'w'))
    return ok_all and g2

if __name__ == '__main__':
    if sys.argv[1] == 'stage0':
        stage0(sys.argv[2])
