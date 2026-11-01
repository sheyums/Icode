"""PRESPEC 20w-B: clock-face cosine test with a per-fly circular-rotation null and a synthetic calibration gate.
Usage: python3 -I planB_20wB.py <data root> <CS DD folder>
Per fly, transitions are stored by arrival minute-of-day (1440 x 3 x 3) and onset/waking events (1440 x 4); rotating a
whole-day record by r minutes moves arrival minute-of-day m to (m - r) mod 1440, i.e. np.roll(.., -r). The original
record has no transition into its first minute, so nothing crosses the wrap junction."""
import sys, os, glob, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R
SETS = {'CS': ('CS_activity_files_1sec_res', True), 'CSDD': (None, False), 'per0LD': ('per0LD_activity_1sec_res', True),
        'per0DD': ('per0DD_activity_files_1sec_res', False), 'fumin': ('Fmn_activity_files_1sec', True),
        'sss': ('sleepless_activity_1sec_res', True)}
FAMILY = ('CS', 'CSDD', 'per0LD'); B = 999; NB = 12; I3 = np.eye(3)
bb = np.arange(NB); XC = np.column_stack([np.ones(NB), np.cos(2 * np.pi * bb / NB), np.sin(2 * np.pi * bb / NB)])

def files(d): return sorted(f for f in glob.glob(os.path.join(d, '**', '*'), recursive=True) if os.path.isfile(f))

def minute_arrays(y):
    y = np.asarray(y); t = np.arange(1, len(y)); ok = (y[:-1] >= 0) & (y[1:] >= 0); m = t % 1440
    N = np.zeros((1440, 3, 3)); np.add.at(N, (m[ok], y[:-1][ok], y[1:][ok]), 1)
    on = ok & (y[:-1] < 2) & (y[1:] == 2); wk = ok & (y[:-1] == 2) & (y[1:] < 2); E = np.zeros((1440, 4))
    np.add.at(E[:, 0], m[on], 1); np.add.at(E[:, 1], m[on & (y[:-1] == 0)], 1)
    np.add.at(E[:, 2], m[wk], 1); np.add.at(E[:, 3], m[wk & (y[1:] == 0)], 1)
    return N, E

def J_bins(Nb):
    """Vectorised model_energetics J for (12,3,3) pooled counts; NaN where a row is empty."""
    rs = Nb.sum(2); good = np.all(rs > 0, 1); J = np.full(NB, np.nan)
    if not good.any(): return J
    P = Nb[good] / rs[good][:, :, None]
    A = np.transpose(P, (0, 2, 1)) - I3; A[:, 2, :] = 1.0
    pi = np.linalg.solve(A, np.broadcast_to([0, 0, 1.0], (len(P), 3))[..., None])[..., 0]
    J[good] = 60 * (pi[:, 1] * P[:, 1, 0] - pi[:, 0] * P[:, 0, 1]); return J

def c_bins(Eb):
    with np.errstate(all='ignore'):
        return np.where((Eb[:, 0] > 0) & (Eb[:, 2] > 0), Eb[:, 1] / Eb[:, 0] - Eb[:, 3] / Eb[:, 2], np.nan)

def cosfit(v):
    g = np.isfinite(v)
    if g.sum() < 4: return np.nan, np.nan
    (m, a, c), *_ = np.linalg.lstsq(XC[g], v[g], rcond=None)
    return float(np.hypot(a, c)), float((np.arctan2(c, a) % (2 * np.pi)) / (2 * np.pi) * 24 + 1)   # peak clock time, bin centre

def stat(Ns, Es, rolls=None):
    if rolls is None: N = Ns.sum(0); E = Es.sum(0)
    else:
        N = sum(np.roll(Ns[i], -r, 0) for i, r in enumerate(rolls)); E = sum(np.roll(Es[i], -r, 0) for i, r in enumerate(rolls))
    Nb = N.reshape(NB, 120, 3, 3).sum(1); Eb = E.reshape(NB, 120, 4).sum(1)
    J = J_bins(Nb); c = c_bins(Eb); AJ, phJ = cosfit(J); Ac, phc = cosfit(c)
    return AJ, phJ, Ac, phc, J

def test(Ns, Es, rng):
    AJ, phJ, Ac, phc, J = stat(Ns, Es); nJ = np.zeros(B); nc = np.zeros(B)
    for k in range(B):
        r = rng.integers(0, 1440, size=len(Ns)); a, _, ac, _, _ = stat(Ns, Es, r); nJ[k] = a; nc[k] = ac
    return dict(AJ=AJ, phJ=phJ, pJ=float((1 + np.sum(nJ >= AJ)) / (B + 1)), nullJ95=float(np.percentile(nJ, 95)),
                Ac=Ac, phc=phc, pc=float((1 + np.sum(nc >= Ac)) / (B + 1)), nullc95=float(np.percentile(nc, 95)), J=J.tolist())

def test_seeded(Ns, Es, seed_of_k):
    """Real-data version: replicate k uses default_rng(20_900 + k) (PRESPEC); offsets drawn per fly from it."""
    AJ, phJ, Ac, phc, J = stat(Ns, Es); nJ = np.zeros(B); nc = np.zeros(B)
    for k in range(B):
        r = np.random.default_rng(seed_of_k(k)).integers(0, 1440, size=len(Ns)); a, _, ac, _, _ = stat(Ns, Es, r); nJ[k] = a; nc[k] = ac
    return dict(AJ=AJ, phJ=phJ, pJ=float((1 + np.sum(nJ >= AJ)) / (B + 1)), nullJ95=float(np.percentile(nJ, 95)),
                Ac=Ac, phc=phc, pc=float((1 + np.sum(nc >= Ac)) / (B + 1)), nullc95=float(np.percentile(nc, 95)), J=J.tolist())

def symbols(path, ld):
    y, _ = R.sqw_symbols(R.read_counts(path), theta=0.5, light_dark=ld); return np.asarray(y)

def sim_rows(P):
    P = np.array(P, float)
    for p in range(P.shape[0]):
        for k in range(3):
            if P[p, k].sum() < 0.5: P[p, k] = 1 / 3
    return P

def main(root, csdd):
    t0 = time.time(); Y = {}
    for g, (sub, ld) in SETS.items():
        Y[g] = [symbols(f, ld) for f in files(csdd if g == 'CSDD' else os.path.join(root, sub))]
    # ---- check: vectorised J equals sqw_ref.model_energetics on CS LD's pooled bins
    Ns = np.array([minute_arrays(y)[0] for y in Y['CS']]); Nb = Ns.sum(0).reshape(NB, 120, 3, 3).sum(1)
    ref = np.array([R.model_energetics(Nb[k] / Nb[k].sum(1)[:, None])['J'] for k in range(NB)])
    print(f"check: max |J_vec - model_energetics| = {np.max(np.abs(J_bins(Nb) - ref)):.2e}")
    # ---- calibration gate (per0 DD matched, constant chains)
    fits = []
    for y in Y['per0DD']:
        ph = np.zeros(len(y), int); _, pi, P = R.mc_fit(y, ph, 1); fits.append((pi, sim_rows(P), ph, y == -1))
    pc = []
    for d in range(200):
        g = np.random.default_rng(21_000 + d); arr = []
        for pi, P, ph, miss in fits:
            _, ys = R.simulate(pi, P, I3, ph, g); ys[miss] = -1; arr.append(minute_arrays(ys))
        r = test(np.array([a[0] for a in arr]), np.array([a[1] for a in arr]), np.random.default_rng(22_000 + d))
        pc.append(r['pJ'])
        if (d + 1) % 50 == 0: print(f"  calibration {d+1}/200 ({time.time()-t0:.0f} s)", flush=True)
    rate = float(np.mean(np.array(pc) < 0.05)); gate = 0.02 <= rate <= 0.09
    print(f"CALIBRATION: fraction p < 0.05 = {rate:.3f} over 200 per0-DD-matched constant-chain datasets (gate [0.02, 0.09]) -> "
          f"{'PASS' if gate else 'FAIL'}")
    out = dict(calibration=dict(p=pc, rate=rate, PASS=bool(gate)))
    json.dump(out, open(os.path.join(HERE, 'planB_20wB.json'), 'w'))
    if not gate:
        print('STOP: calibration gate failed; no real-data verdict (PRESPEC 20w-B).'); return
    # ---- optional power: per0 DD size, CS LD's pooled 12-bin chains
    PB = sim_rows(Nb / Nb.sum(2, keepdims=True)); pw = []
    for d in range(100):
        g = np.random.default_rng(23_000 + d); arr = []
        for y in Y['per0DD']:
            T = len(y); ph = (np.arange(T) % 1440) // 120; pi0 = R.stationary(PB[0])
            _, ys = R.simulate(pi0, PB, I3, ph, g); ys[y == -1] = -1; arr.append(minute_arrays(ys))
        r = test(np.array([a[0] for a in arr]), np.array([a[1] for a in arr]), np.random.default_rng(24_000 + d)); pw.append(r['pJ'])
    power = float(np.mean(np.array(pw) < 0.05))
    print(f"POWER at per0 DD size vs CS LD's pooled modulation (no heterogeneity, optimistic): {power:.2f} (100 datasets)")
    out['power'] = dict(p=pw, power=power)
    # ---- real data
    res = {}
    for g in SETS:
        A = [minute_arrays(y) for y in Y[g]]
        res[g] = test_seeded(np.array([a[0] for a in A]), np.array([a[1] for a in A]), lambda k: 20_900 + k)
        res[g]['n'] = len(A)
    ps = sorted((res[g]['pJ'], g) for g in FAMILY); holm = {}; run = 0
    for k, (pv, g) in enumerate(ps):
        run = max(run, min(1.0, (len(ps) - k) * pv)); holm[g] = run
    flag = res['per0DD']['pJ'] < 0.05
    for g in SETS:
        r = res[g]; tag = (f"Holm {holm[g]:.3f} -> {'CLOCK-LOCKED MODULATION' if holm[g] < 0.05 else 'not detected'}" if g in FAMILY
                           else ('NEGATIVE CONTROL' if g == 'per0DD' else 'descriptive'))
        print(f"{g:7s} n {r['n']:2d} | J amplitude {r['AJ']:.4f} (null 95th {r['nullJ95']:.4f}) peak {r['phJ']:.1f} h | p {r['pJ']:.3f} | {tag} | "
              f"c_obs amplitude {r['Ac']:.4f} (null 95th {r['nullc95']:.4f}) peak {r['phc']:.1f} h p {r['pc']:.3f}")
    print('FLAG: per0 DD p < 0.05 -> every positive carries the miscalibration/residual-rhythm caveat' if flag else
          'per0 DD null (p >= 0.05): no flag')
    out.update(real=res, holm=holm, flag=bool(flag), seconds=time.time() - t0)
    json.dump(out, open(os.path.join(HERE, 'planB_20wB.json'), 'w'))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
