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

# ================================================================== stage 1: reversible aggregated HMM (numba)
import math
from numba import njit

def runs_of(a):
    a = np.asarray(a).astype(np.int64); ch = np.flatnonzero(np.diff(a)) + 1
    st = np.r_[0, ch]; en = np.r_[ch, len(a)]
    return a[st].astype(np.int64), (en - st).astype(np.int64)

@njit(cache=True)
def _kernel(sym, L, P, Ka, K, lam1, V1, Vi1, lam0, V0, Vi0, pi):
    nR = len(sym); F = np.zeros((nR, 3)); U = np.zeros((nR, 3)); Bk = np.zeros((nR, 3)); G = np.zeros((nR, 3))
    mm = (K - Ka, Ka); oo = (Ka, 0)
    b = sym[0]; m = mm[b]; off = oo[b]
    f = np.zeros(3); s = 0.0
    for i in range(m):
        f[i] = pi[off + i]; s += f[i]
    for i in range(m): f[i] /= s
    ll = math.log(s)
    for r in range(nR):
        b = sym[r]; m = mm[b]; off = oo[b]
        lam = lam1 if b == 1 else lam0; V = V1 if b == 1 else V0; Vi = Vi1 if b == 1 else Vi0
        for i in range(m): F[r, i] = f[i]
        at = np.zeros(3)
        for k in range(m):
            for i in range(m): at[k] += f[i] * V[i, k]
        lmx = 0.0
        for k in range(m): lmx = max(lmx, lam[k])
        if lmx <= 0.0: lmx = 1e-300
        ll += (L[r] - 1.0) * math.log(lmx)
        u = np.zeros(3)
        for j in range(m):
            for k in range(m): u[j] += at[k] * math.pow(lam[k] / lmx, L[r] - 1.0) * Vi[k, j]
        if r == nR - 1:
            s = 0.0
            for j in range(m): s += u[j]
            if not (s > 0.0): return -np.inf, np.zeros((K, K))
            ll += math.log(s)
            for j in range(m): U[r, j] = u[j] / s
            break
        b2 = sym[r + 1]; m2 = mm[b2]; off2 = oo[b2]
        fn = np.zeros(3); s = 0.0
        for j in range(m2):
            for i in range(m): fn[j] += u[i] * P[off + i, off2 + j]
            s += fn[j]
        if not (s > 0.0): return -np.inf, np.zeros((K, K))
        ll += math.log(s)
        for j in range(m2): f[j] = fn[j] / s
        for i in range(m): U[r, i] = u[i] / s
    bl = sym[nR - 1]
    for i in range(mm[bl]): Bk[nR - 1, i] = 1.0
    for r in range(nR - 2, -1, -1):
        b = sym[r]; m = mm[b]; off = oo[b]; b2 = sym[r + 1]; m2 = mm[b2]; off2 = oo[b2]
        lam = lam1 if b2 == 1 else lam0; V = V1 if b2 == 1 else V0; Vi = Vi1 if b2 == 1 else Vi0
        bt = np.zeros(3)
        for k in range(m2):
            for j in range(m2): bt[k] += Vi[k, j] * Bk[r + 1, j]
        lmx2 = 0.0
        for k in range(m2): lmx2 = max(lmx2, lam[k])
        if lmx2 <= 0.0: lmx2 = 1e-300
        g = np.zeros(3)
        for j in range(m2):
            for k in range(m2): g[j] += V[j, k] * math.pow(lam[k] / lmx2, L[r + 1] - 1.0) * bt[k]
        for j in range(m2): G[r + 1, j] = g[j]
        s = 0.0; bb = np.zeros(3)
        for i in range(m):
            for j in range(m2): bb[i] += P[off + i, off2 + j] * g[j]
            s += bb[i]
        if not (s > 0.0): return -np.inf, np.zeros((K, K))
        for i in range(m): Bk[r, i] = bb[i] / s
    C = np.zeros((K, K))
    for r in range(nR):
        b = sym[r]; m = mm[b]; off = oo[b]
        lam = lam1 if b == 1 else lam0; V = V1 if b == 1 else V0; Vi = Vi1 if b == 1 else Vi0
        if L[r] >= 2:
            at = np.zeros(3); bt = np.zeros(3)
            for k in range(m):
                for i in range(m):
                    at[k] += F[r, i] * V[i, k]; bt[k] += Vi[k, i] * Bk[r, i]
            lmx = 0.0
            for k in range(m): lmx = max(lmx, lam[k])
            if lmx <= 0.0: lmx = 1e-300
            Z = 0.0
            for k in range(m): Z += at[k] * math.pow(lam[k] / lmx, L[r] - 1.0) * bt[k]
            Gm = np.zeros((3, 3))
            for k in range(m):
                for l in range(m):
                    if abs(lam[k] - lam[l]) < 1e-12:
                        Gm[k, l] = (L[r] - 1.0) * math.pow(lam[k] / lmx, L[r] - 2.0) / lmx
                    else:
                        Gm[k, l] = (math.pow(lam[k] / lmx, L[r] - 1.0) - math.pow(lam[l] / lmx, L[r] - 1.0)) / (lam[k] - lam[l])
            for i in range(m):
                for j in range(m):
                    Mij = 0.0
                    for k in range(m):
                        for l in range(m): Mij += Vi[k, i] * at[k] * Gm[k, l] * bt[l] * V[j, l]
                    C[off + i, off + j] += P[off + i, off + j] * Mij / Z
        if r < nR - 1:
            b2 = sym[r + 1]; m2 = mm[b2]; off2 = oo[b2]; s = 0.0
            for i in range(m):
                for j in range(m2): s += U[r, i] * P[off + i, off2 + j] * G[r + 1, j]
            for i in range(m):
                for j in range(m2): C[off + i, off2 + j] += U[r, i] * P[off + i, off2 + j] * G[r + 1, j] / s
    return ll, C

def _block_eig(X, idx):
    x = X.sum(1); sub = X[np.ix_(idx, idx)]; d = x[idx]
    S = sub / np.sqrt(np.outer(d, d)); lam, Uq = np.linalg.eigh(S)
    V = Uq / np.sqrt(d)[:, None]; Vi = Uq.T * np.sqrt(d)[None, :]
    lp = np.zeros(3); Vp = np.zeros((3, 3)); Vip = np.zeros((3, 3)); m = len(idx)
    lp[:m] = lam; Vp[:m, :m] = V; Vip[:m, :m] = Vi
    return lp, Vp, Vip

def loglik_counts(X, Ka, sym, L):
    K = X.shape[0]; x = X.sum(1); P = X / x[:, None]
    l1, V1, Vi1 = _block_eig(X, np.arange(Ka)); l0, V0, Vi0 = _block_eig(X, np.arange(Ka, K))
    return _kernel(sym, L, P, Ka, K, l1, V1, Vi1, l0, V0, Vi0, x / x.sum())

def reversible_mle(C, X, iters=2000, tol=1e-13):
    c = C.sum(1)
    for _ in range(iters):
        x = X.sum(1); Xn = (C + C.T) / (c[:, None] / x[:, None] + c[None, :] / x[None, :]); Xn /= Xn.sum()
        if np.max(np.abs(Xn - X)) < tol: return Xn
        X = Xn
    return X

def random_start(K, rng):
    pi = rng.dirichlet(np.ones(K)); F = np.zeros((K, K))
    for i in range(K):
        for j in range(i + 1, K):
            F[i, j] = F[j, i] = pi[i] * pi[j] * math.exp(rng.uniform(math.log(1e-3), math.log(0.3)))
    for i in range(K): F[i, i] = pi[i] - F[i].sum()
    return F / F.sum()

def fit_order(sym, L, Ka, Ki, rng, starts=10, maxit=2000, tol=1e-9):
    K = Ka + Ki; best = (-np.inf, None, 0)
    for s in range(starts):
        X = random_start(K, rng); old = -np.inf; it = 0
        for it in range(maxit):
            ll, C = loglik_counts(X, Ka, sym, L)
            if not np.isfinite(ll): break
            if it > 0 and ll - old < tol * abs(ll): break
            old = ll; X = reversible_mle(C, X)
        if ll > best[0]: best = (ll, X, it)
    return best

def fit_fly(args):
    i, a = args
    sym, L = runs_of(a); rng = np.random.default_rng(3000 + i); res = []
    for Ka in (1, 2, 3):
        for Ki in (1, 2, 3):
            ll, X, it = fit_order(sym, L, Ka, Ki, rng); K = Ka + Ki; k = K * (K + 1) // 2 - 1
            res.append(dict(Ka=Ka, Ki=Ki, ll=float(ll), k=k, bic=float(-2 * ll + k * math.log(len(L))), X=X.tolist(), iters=it))
    best = min(res, key=lambda r: r['bic'])
    X = np.array(best['X']); x = X.sum(1); S = X / np.sqrt(np.outer(x, x))
    best['min_eig'] = float(np.linalg.eigvalsh(S).min()); best['n_runs'] = int(len(L))
    return dict(fly=i, best=best, all=[{k: v for k, v in r.items() if k != 'X'} for r in res])

# ================================================================== simulation and scoring
@njit(cache=True)
def _seed(s): np.random.seed(s)

@njit(cache=True)
def _sim(P, pi, Ka, T):
    K = P.shape[0]; out = np.zeros(T, np.int8); u = np.random.random(); c = 0.0; z = 0
    for k in range(K):
        c += pi[k]
        if u <= c: z = k; break
    t = 0
    while t < T:
        stay = P[z, z]
        n = np.random.geometric(1.0 - stay) if stay < 1.0 else T
        e = min(T, t + n)
        if z < Ka:
            for q in range(t, e): out[q] = 1
        t = e
        u = np.random.random() * (1.0 - stay); c = 0.0
        for k in range(K):
            if k == z: continue
            c += P[z, k]
            if u <= c: z = k; break
    return out

def longest_silence(a):
    a = np.asarray(a); ch = np.flatnonzero(np.diff(np.r_[1, a, 1]))
    return int(np.max(ch[1::2] - ch[::2])) if len(ch) > 1 else 0

DEATH = 12 * 3600

# ================================================================== driver: stage 1, 1b, 2, verdict
def load_flies(datadir):
    W, Sl, files, tl = load_ref()
    up = sorted(f for f in glob.glob(os.path.join(datadir, '**', '*'), recursive=True) if os.path.isfile(f))
    A = []
    for i, name in enumerate(files):
        m = [f for f in up if os.path.basename(f).startswith(name)]; x = R.read_counts(m[0])[:tl[i]]
        A.append((x > 0).astype(np.int8))
    return A, W, Sl, files, tl

def score(a):
    w, s = extract(a, len(a)); return w, s

def run(datadir):
    from multiprocessing import Pool
    import time
    A, W, Sl, files, tl = load_flies(datadir); n = len(A); t0 = time.time()
    xObs = W[np.isfinite(W)]; nObs = len(xObs); rrObs = [rise_ratio(xObs, 2.0, 1.0, k) for k in NB]
    obs_w = np.array([np.isfinite(W[i]).sum() for i in range(n)]); obs_s = np.array([np.nansum(Sl[i]) for i in range(n)])
    # ---- stage 1
    ck = os.path.join(HERE, 'planB_20u_fits'); os.makedirs(ck, exist_ok=True)
    todo = [(i, A[i]) for i in range(n) if not os.path.exists(os.path.join(ck, f'fly{i:02d}.json'))]
    print(f"stage 1: {n - len(todo)} flies already checkpointed, fitting {len(todo)}", flush=True)
    with Pool(4) as pool:
        for f in pool.imap_unordered(fit_fly, todo, chunksize=1):
            json.dump(f, open(os.path.join(ck, f"fly{f['fly']:02d}.json"), 'w'))
            print(f"  fitted fly {f['fly']} ({time.time()-t0:.0f} s)", flush=True)
    fits = [json.load(open(os.path.join(ck, f'fly{i:02d}.json'))) for i in range(n)]
    out = dict(stage1=fits, rrObs=rrObs, nObs=int(nObs), t_stage1=time.time() - t0)
    print(f"stage 1 done {time.time()-t0:.0f} s")
    for f in fits:
        b = f['best']; print(f"  fly {f['fly']:2d} runs {b['n_runs']:6d} | selected Ka {b['Ka']} Ki {b['Ki']} | logL {b['ll']:.1f} | "
                            f"min eig {b['min_eig']:.4f} | iters {b['iters']}")
    Ps = []; pis = []; Kas = []
    for f in fits:
        X = np.array(f['best']['X']); x = X.sum(1); Ps.append(X / x[:, None]); pis.append(x / x.sum()); Kas.append(f['best']['Ka'])
    nonlazy = [f['fly'] for f in fits if f['best']['min_eig'] <= 0]
    print(f"  laziness: flies with min eigenvalue <= 0: {nonlazy}")
    # ---- stage 1b adequacy
    rw = np.zeros(n); rs = np.zeros(n); dead_ad = 0
    for i in range(n):
        _seed(4000 + i); ws = []; ss = []
        for k in range(20):
            a = _sim(Ps[i], pis[i], Kas[i], int(tl[i]))
            if longest_silence(a) >= DEATH: dead_ad += 1; continue
            w, s = score(a); ws.append(len(w)); ss.append(s.sum())
        rw[i] = np.mean(ws) / obs_w[i]; rs[i] = np.mean(ss) / obs_s[i]
    gate = 0.85 <= np.median(rw) <= 1.15 and 0.85 <= np.median(rs) <= 1.15
    print(f"stage 1b adequacy: wake-bout count ratio median {np.median(rw):.3f} [IQR {np.percentile(rw,25):.3f}, {np.percentile(rw,75):.3f}] | "
          f"total sleep ratio median {np.median(rs):.3f} [IQR {np.percentile(rs,25):.3f}, {np.percentile(rs,75):.3f}] | "
          f"death-excluded sims {dead_ad}/{20*n} -> {'PASS' if gate else 'FAIL (NULL INADEQUATE)'}")
    out.update(adequacy=dict(rw=rw.tolist(), rs=rs.tolist(), dead=dead_ad, PASS=bool(gate)))
    json.dump(out, open(os.path.join(HERE, 'planB_20u_run.json'), 'w'))
    if not gate:
        print('VERDICT: NULL INADEQUATE (stage 2 not run, per PRESPEC 20u)'); return
    # ---- stage 2
    rr = np.zeros((199, 4)); nsyn = np.zeros(199, int); dead = 0
    for b in range(1, 200):
        _seed(5000 + b); g = np.random.default_rng(5000 + b); xs = []
        for i in range(n):
            a = _sim(Ps[i], pis[i], Kas[i], int(tl[i]))
            if longest_silence(a) >= DEATH: dead += 1; continue
            w, _ = score(a); xs.append(w)
        xs = np.concatenate(xs); nsyn[b - 1] = len(xs)
        if len(xs) > nObs: xs = xs[g.permutation(len(xs))[:nObs]]
        rr[b - 1] = [rise_ratio(xs, 2.0, 1.0, k) for k in NB]
        if b % 20 == 0: print(f"  replicate {b} ({time.time()-t0:.0f} s)", flush=True)
    p = [(1 + np.sum(rr[:, j] >= rrObs[j])) / 200 for j in range(4)]
    print('stage 2: observed RiseRatio ' + ' '.join(f'{v:.3f}' for v in rrObs))
    print('         null median ' + ' '.join(f'{v:.3f}' for v in np.median(rr, 0)) + ' | 95th ' + ' '.join(f'{v:.3f}' for v in np.percentile(rr, 95, axis=0)))
    print('         p_rev ' + ' '.join(f'{v:.3f}' for v in p) + f' | pooled synthetic wake bouts median {np.median(nsyn):.0f} (obs {nObs}), min {nsyn.min()}')
    frac_dead = dead / (199 * n)
    nlo = sum(v < 0.05 for v in p); nhi = sum(v >= 0.05 for v in p)
    verdict = ('RISE BEYOND REVERSIBLE DYNAMICS' if nlo >= 3 else 'EXPLAINED BY REVERSIBLE DYNAMICS PLUS SCORING' if nhi >= 3
               else 'INCONCLUSIVE (bin-sensitive)')
    lim = ' [limitation: death-criterion exclusions > 5%]' if frac_dead > 0.05 else ''
    print(f"         death-excluded fly-replicates {dead}/{199*n} ({frac_dead:.3%})")
    print(f"VERDICT: {verdict}{lim}")
    out.update(stage2=dict(rr=rr.tolist(), p=p, nsyn=nsyn.tolist(), dead=dead, verdict=verdict + lim), seconds=time.time() - t0)
    json.dump(out, open(os.path.join(HERE, 'planB_20u_run.json'), 'w'))

if __name__ == '__main__':
    if sys.argv[1] == 'stage0':
        stage0(sys.argv[2])
    elif sys.argv[1] == 'run':
        run(sys.argv[2])
