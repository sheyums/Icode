"""
planB_f1f5.py -- PRESPEC 20t, F1 (synthetic recovery of clock a_ij and homeostat b_ij) and F5 (recovery of the
housekeeping / driving split of entropy production). Synthetic only; no data.

Model (PRESPEC 20t): per-minute competing-risks chain on S/Q/W with
  log k_ij(t) = c_ij + alpha_ij cos(w t) + beta_ij sin(w t) + b_ij H(t)
H = Process S from the chain's own sleep. Rows fitted independently given H; (tau_r, tau_d) profiled on a fixed grid.
Run:  python3 planB_f1f5.py  [n_flies]   -> writes planB_f1f5.json and prints the report lines.
"""
import json, sys, os, time
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm
from multiprocessing import Pool

Q, W, S = 0, 1, 2
HERE = os.path.dirname(os.path.abspath(__file__))
DAYMIN = 1440
WCLK = 2 * np.pi / DAYMIN                      # tau = 24 h, known
GRID = [1, 2, 4, 8, 16, 32]                   # hours, tau_r and tau_d
TRANS = [(Q, W), (Q, S), (W, Q), (W, S), (S, Q), (S, W)]   # 6 off-diagonal transitions; row i -> its 2 destinations
DEST = {i: [j for j in range(3) if j != i] for i in range(3)}

# ------------------------------------------------------------------------------------------------ baseline rates
def baseline_rates():
    g = json.load(open(os.path.join(HERE, 'golden', 'G12_planA.json')))
    P = np.array(g['mc_fit']['P'])[0]                       # light phase
    k = np.zeros((3, 3))
    for i in range(3):
        K = -np.log(P[i, i])
        for j in DEST[i]:
            k[i, j] = K * P[i, j] / (1 - P[i, i])
    return k, P

# ------------------------------------------------------------------------------------------------ Process S
def h_update(h, sleep, tr, td):
    return h * np.exp(-1 / td) if sleep else 1 - (1 - h) * np.exp(-1 / tr)

def h_path(y, tr_h, td_h):
    """H_t for t = 0..T-1, H_t from y_{t-1} (PRESPEC rule: initialise by running over the first 1440 min from 0.5)."""
    tr, td = tr_h * 60.0, td_h * 60.0
    er, ed = np.exp(-1 / tr), np.exp(-1 / td)
    sl = (np.asarray(y) == S)
    h = 0.5
    for t in range(min(DAYMIN, len(y))):
        h = h * ed if sl[t] else 1 - (1 - h) * er
    H = np.empty(len(y)); H[0] = h
    for t in range(1, len(y)):
        h = h * ed if sl[t - 1] else 1 - (1 - h) * er
        H[t] = h
    return H

# ------------------------------------------------------------------------------------------------ simulation
def rates_at(coef, cs, sn, h, step, dark):
    """coef: dict (i,j) -> [c, alpha, beta, b]; step: dict (i,j) -> dark offset or None. Returns 3x3 rate matrix."""
    k = np.zeros((3, 3))
    for (i, j), v in coef.items():
        e = v[0] + v[1] * cs + v[2] * sn + v[3] * h
        if step is not None and dark:
            e += step.get((i, j), 0.0)
        k[i, j] = np.exp(e)
    return k

def comp_risk_P(k):
    P = np.zeros((3, 3))
    for i in range(3):
        K = k[i].sum(); P[i, i] = np.exp(-K)
        for j in DEST[i]:
            P[i, j] = -np.expm1(-K) * k[i, j] / K
    return P

def simulate(truth, days, rng, burn_days=2):
    """Returns y (T,), H_true (T,), with H centred inside the coefficient c (truth['Hc'])."""
    coef, tr, td, step, Hc = truth['coef'], truth['tr'] * 60, truth['td'] * 60, truth.get('step'), truth['Hc']
    T = days * DAYMIN; Tb = burn_days * DAYMIN
    y = np.empty(T + Tb, int); H = np.empty(T + Tb)
    y[0] = S; h = 0.5
    t0 = -Tb                                         # clock time of the first burn-in minute (whole days, so ZT aligned)
    for n in range(T + Tb):
        H[n] = h
        if n > 0:
            tt = t0 + n - 1                          # rates at minute t-1
            k = rates_at(coef, np.cos(WCLK * tt), np.sin(WCLK * tt), h - Hc, step, (tt % DAYMIN) >= 720)
            P = comp_risk_P(k)
            y[n] = rng.choice(3, p=P[y[n - 1]])
            h = h * np.exp(-1 / td) if y[n - 1] == S else 1 - (1 - h) * np.exp(-1 / tr)
            H[n] = h
    # H[n] here is the value after the update from y[n-1], i.e. H_n; the rate for the step into n used H_{n-1}. Consistent.
    return y[Tb:], H[Tb:]

# ------------------------------------------------------------------------------------------------ fitting
def row_nll(theta, X, o, p):
    th = theta.reshape(2, p)
    eta = np.clip(X @ th.T, -30, 8); k = np.exp(eta); K = k.sum(1)
    stay = (o == 0)
    ll = np.where(stay, -K, np.log(-np.expm1(-K)) - np.log(K))
    mv = ~stay
    ll[mv] += eta[mv, o[mv] - 1]
    g = 1.0 / np.expm1(K) - 1.0 / K
    d = np.where(stay[:, None], -k, g[:, None] * k)
    d[mv, o[mv] - 1] += 1.0
    grad = -(d.T @ X).ravel()
    return -ll.sum(), grad

def fit_row(X, o, cols, x0=None):
    Xs = X[:, cols]; p = len(cols)
    if x0 is None:
        n = len(o); c = [max(np.sum(o == m), 0.5) / n for m in (1, 2)]
        x0 = np.zeros((2, p)); x0[:, 0] = np.log(c)
        x0 = x0.ravel()
    r = minimize(row_nll, x0, args=(Xs, o, p), jac=True, method='L-BFGS-B',
                 options=dict(maxiter=500, gtol=1e-7))
    return -r.fun, r.x.reshape(2, p)

def design(T, H, t0=0):
    tt = np.arange(T) - 1 + t0                      # transition into t uses minute t-1
    X = np.column_stack([np.ones(T), np.cos(WCLK * tt), np.sin(WCLK * tt), np.r_[H[0], H[:-1]]])
    return X                                         # row t = covariates for the step INTO minute t

def row_data(y):
    out = {}
    y = np.asarray(y)
    for i in range(3):
        idx = np.flatnonzero((y[:-1] == i) & (y[1:] >= 0)) + 1
        dest = y[idx]; o = np.zeros(len(idx), int)
        for m, j in enumerate(DEST[i]):
            o[dest == j] = m + 1
        out[i] = (idx, o)
    return out

def fit_fly(y, clock=True):
    """Full grid-profiled fit. Returns dict with LLs, per-row tests, best grid estimates."""
    T = len(y); rd = row_data(y)
    cols_full = [0, 1, 2, 3] if clock else [0, 3]
    cols_noclk = [0, 3]; cols_noh = [0, 1, 2] if clock else [0]
    X0 = design(T, np.zeros(T))
    const = {}; noh = {}
    for i in range(3):
        idx, o = rd[i]
        const[i] = fit_row(X0[idx], o, [0])
        noh[i] = fit_row(X0[idx], o, cols_noh, None)
    G = [(a, b) for a in GRID for b in GRID]
    full = np.zeros((len(G), 3)); nock = np.zeros((len(G), 3)); thetas = [None] * len(G)
    for gi, (tr, td) in enumerate(G):
        H = h_path(y, tr, td); X = design(T, H); th = {}
        for i in range(3):
            idx, o = rd[i]
            full[gi, i], th[i] = fit_row(X[idx], o, cols_full)
            if clock:
                nock[gi, i], _ = fit_row(X[idx], o, cols_noclk)
        thetas[gi] = th
    tot = full.sum(1); gbest = int(np.argmax(tot)); LLfull = tot[gbest]
    res = dict(LL_full=float(LLfull), LL_const=float(sum(const[i][0] for i in range(3))),
               tau_r=G[gbest][0], tau_d=G[gbest][1], LR_clock=[], LR_homeo=[])
    for i in range(3):
        oth = tot - full[:, i]
        if clock:
            res['LR_clock'].append(float(2 * (LLfull - np.max(oth + nock[:, i]))))
        res['LR_homeo'].append(float(2 * (LLfull - (np.max(oth) + noh[i][0]))))
    th = thetas[gbest]; Hhat = h_path(y, *G[gbest])
    coef = {}
    for i in range(3):
        for m, j in enumerate(DEST[i]):
            v = np.zeros(4); v[cols_full] = th[i][m]; coef[(i, j)] = v
    res['coef'] = {f'{i}{j}': coef[(i, j)].tolist() for (i, j) in coef}
    res['sd_Hhat'] = float(np.std(Hhat))
    return res, coef, Hhat

# ------------------------------------------------------------------------------------------------ thermodynamics (F5)
def thermo(coef, H, T, step=None, Hc=0.0):
    """Housekeeping and non-adiabatic EPR (nats/h), averaged over the record, with H along the given path."""
    tt = np.arange(T)
    Kt = np.zeros((T, 3, 3))
    for (i, j), v in coef.items():
        e = v[0] + v[1] * np.cos(WCLK * tt) + v[2] * np.sin(WCLK * tt) + v[3] * (H - Hc)
        if step is not None:
            e = e + step.get((i, j), 0.0) * ((tt % DAYMIN) >= 720)
        Kt[:, i, j] = np.exp(e)
    for i in range(3):
        Kt[:, i, i] = -Kt[:, i, :].sum(1)
    A = np.transpose(Kt, (0, 2, 1)).copy(); A[:, 2, :] = 1.0
    pi = np.linalg.solve(A, np.broadcast_to(np.array([0, 0, 1.0]), (T, 3))[..., None])[..., 0]
    p = pi[0].copy(); hk = 0.0; na = 0.0
    for t in range(T):
        k = Kt[t]; s = pi[t]
        for i in range(3):
            for j in range(3):
                if i != j:
                    f = p[i] * k[i, j]
                    hk += f * np.log(s[i] * k[i, j] / (s[j] * k[j, i]))
                    na += f * np.log(p[j] * s[i] / (p[i] * s[j]))
        p = p @ expm(k); p = np.clip(p, 1e-300, None); p /= p.sum()
    return 60 * hk / T, 60 * na / T

def multiple_r(H):
    T = len(H); tt = np.arange(T)
    X = np.column_stack([np.ones(T), np.cos(WCLK * tt), np.sin(WCLK * tt)])
    b, *_ = np.linalg.lstsq(X, H, rcond=None); r2 = 1 - np.var(H - X @ b) / np.var(H)
    return float(np.sqrt(max(r2, 0)))

# ------------------------------------------------------------------------------------------------ scenarios
SCEN = {  # name: (s, days, a, Eb, clock_on, homeo_on, step)
    'N0': (0, 4, 0.0, 0.0, False, False, False),
    'DD-CH': (1, 4, 0.40, 0.20, True, True, False),
    'DD-CH-half': (2, 4, 0.20, 0.10, True, True, False),
    'DD-C': (3, 4, 0.40, 0.0, True, False, False),
    'DD-H': (4, 4, 0.0, 0.20, False, True, False),
    'LD2-CH': (5, 2, 0.40, 0.20, True, True, False),
    'LD3-step-C': (6, 3, 0.40, 0.0, True, False, True),
    'N0-2d': (7, 2, 0.0, 0.0, False, False, False),
}
PEAK = {(Q, S): 18 * 60, (S, W): 6 * 60}
BSIGN = {(Q, S): +1, (S, W): -1}
STEP = {(Q, S): +0.4, (S, W): -0.4}

def make_truth(k0, rng, a, b, Hc, clock_on, homeo_on, step, jitter=True):
    coef = {}
    for (i, j) in TRANS:
        c = np.log(k0[i, j]) + (rng.normal(0, 0.25) if jitter else 0.0)
        al = be = bb = 0.0
        if (i, j) in PEAK:
            amp = rng.uniform(0.75, 1.25) if jitter else 1.0
            if clock_on:
                tp = PEAK[(i, j)] + (rng.normal(0, 60) if jitter else 0.0)
                al, be = a * amp * np.cos(WCLK * tp), a * amp * np.sin(WCLK * tp)
            if homeo_on:
                bb = BSIGN[(i, j)] * b * amp
        coef[(i, j)] = np.array([c, al, be, bb])
    tr = 6 * (np.exp(rng.normal(0, 0.2)) if jitter else 1.0)
    td = 3 * (np.exp(rng.normal(0, 0.2)) if jitter else 1.0)
    return dict(coef=coef, tr=tr, td=td, Hc=Hc, step=STEP if step else None)

def calibrate(k0, a, clock_on, days):
    rng = np.random.default_rng(19_999); sds = []; ms = []
    for f in range(10):
        tru = make_truth(k0, rng, a, 0.0, 0.0, clock_on, False, False)
        y, H = simulate(tru, days, rng); sds.append(np.std(H)); ms.append(np.mean(H))
    return float(np.mean(sds)), float(np.mean(ms))

def run_fly(args):
    name, f, truth, days, seed = args
    rng = np.random.default_rng(seed)
    y, Htrue = simulate(truth, days, rng)
    T = len(y)
    res, coef_hat, Hhat = fit_fly(y)
    out = dict(scen=name, fly=f, tr_true=truth['tr'], td_true=truth['td'], **res)
    out['coef_true'] = {f'{i}{j}': truth['coef'][(i, j)].tolist() for (i, j) in truth['coef']}
    out['sd_Htrue'] = float(np.std(Htrue))
    out['r_true'] = multiple_r(Htrue); out['r_hat'] = multiple_r(Hhat)
    out['frac_S'] = float(np.mean(y == S))
    if name in ('N0', 'DD-CH', 'LD2-CH', 'DD-CH-half'):
        hk, na = thermo(truth['coef'], Htrue, T, truth['step'], truth['Hc'])
        hk2, na2 = thermo(coef_hat, Hhat, T)
        out.update(hk_true=hk, na_true=na, hk_hat=hk2, na_hat=na2)
    return out

def main(nfly=25):
    k0, _ = baseline_rates()
    jobs = []; meta = {}
    for name, (s, days, a, Eb, con, hon, step) in SCEN.items():
        if hon:
            sd, mH = calibrate(k0, a, con, days); b = Eb / sd
        else:
            sd, mH, b = np.nan, 0.5, 0.0
        meta[name] = dict(days=days, a=a, Eb=Eb, b=b, sdH_calib=sd, Hc=mH)
        rng = np.random.default_rng(20_000 + s)
        for f in range(nfly):
            tru = make_truth(k0, rng, a, b, mH, con, hon, step)
            jobs.append((name, f, tru, days, int(rng.integers(2**31))))
    t0 = time.time()
    with Pool(4) as pool:
        out = pool.map(run_fly, jobs, chunksize=1)
    json.dump(dict(meta=meta, flies=out, seconds=time.time() - t0), open(os.path.join(HERE, 'planB_f1f5.json'), 'w'))
    print('done', time.time() - t0)

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 25)
