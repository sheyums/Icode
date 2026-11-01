"""PRESPEC 20v amendment F: the scoring-matched null.
  python3 -I planB_20vF.py f1 <per0DD dir> <CS DD dir>     # F1: time symmetry of sqw_symbols in the code
  python3 -I planB_20vF.py f2 <per0DD dir> <CS DD dir>     # F2/F3 (needs 20u stage-1 fits; fits CS DD itself)"""
import sys, os, glob, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R

def files(d): return sorted(f for f in glob.glob(os.path.join(d, '**', '*'), recursive=True) if os.path.isfile(f))

def per0_traces(d):
    """per0 DD retained traces: first traceLen seconds (20r reference), floored to whole days; reference fly order."""
    import planB_20u as U
    W, Sl, names, tl = U.load_ref(); up = files(d); out = []
    for i, name in enumerate(names):
        m = [f for f in up if os.path.basename(f).startswith(name)]; assert len(m) == 1
        x = R.read_counts(m[0])[:tl[i]]; out.append(x[:(len(x) // 86400) * 86400])
    return out, names

def csdd_traces(d):
    fs = files(d); out = []
    for f in fs:
        x = R.read_counts(f); out.append(x[:(len(x) // 86400) * 86400])
    return out, [os.path.basename(f) for f in fs]

def f1(per0dir, csdddir):
    res = {}; allpass = True
    for g, (X, names) in (('per0DD', per0_traces(per0dir)), ('CSDD', csdd_traces(csdddir))):
        rows = []
        for x, nm in zip(X, names):
            y1 = R.sqw_symbols(x, theta=0.5, light_dark=False)[0][::-1]
            y2 = R.sqw_symbols(x[::-1], theta=0.5, light_dark=False)[0]
            assert len(y1) == len(y2)
            d = np.flatnonzero(y1 != y2); frac = 1 - len(d) / len(y1)
            kinds = {}
            for k in d: kinds[f'{y1[k]}->{y2[k]}'] = kinds.get(f'{y1[k]}->{y2[k]}', 0) + 1
            rows.append(dict(file=nm, minutes=len(y1), differ=len(d), agree=frac, kinds=kinds))
        ok = all(r['agree'] >= 0.999 for r in rows); allpass &= ok; res[g] = rows
        agg = {}
        for r in rows:
            for k, v in r['kinds'].items(): agg[k] = agg.get(k, 0) + v
        print(f"F1 {g:7s} flies {len(rows)} | minutes differing total {sum(r['differ'] for r in rows)} of "
              f"{sum(r['minutes'] for r in rows)} | worst fly agreement {min(r['agree'] for r in rows):.5f} | "
              f"kinds (forward symbol->reversed-scored symbol; Q0 W1 S2 missing-1) {agg} -> {'PASS' if ok else 'FAIL'}")
        for r in rows:
            if r['differ']: print(f"     {r['file'][:45]:45s} differ {r['differ']:4d} agree {r['agree']:.5f} {r['kinds']}")
    print(f"F1 overall: {'PASS' if allpass else 'FAIL'}")
    json.dump(dict(F1=res, PASS=bool(allpass)), open(os.path.join(HERE, 'planB_20vF_f1.json'), 'w'))

# ======================================================================== F2 / F3
import planB_20v as V
from planB_20u import _seed, _sim, longest_silence, DEATH, extract, fit_fly
from scipy.stats import wilcoxon, binom

CELLS = [(s,) + ex for s in V.PAIRS for ex in V.PAIRS[s]]

def cell_bins(y):
    """Per (s, e, x) cell: list of dwell-bin indices (whole record, one phase)."""
    y = np.asarray(y); T = len(y); ch = np.flatnonzero(np.diff(y)) + 1; st = np.r_[0, ch]; en = np.r_[ch, T]
    out = {c: [] for c in CELLS}; nrun = np.zeros(3)
    for a, b in zip(st, en):
        s = y[a]
        if s >= 0: nrun[s] += 1
        if s not in V.PAIRS or a == 0 or b == T: continue
        c = (s, y[a - 1], y[b])
        if c in out: out[c].append(np.searchsorted(V.EDGES, b - a, side='right') - 1)
    return out, nrun

def sigma_from_cells(cells, hours):
    H = {c: np.bincount(np.asarray(cells[c], int), minlength=6).astype(float) for c in CELLS}
    return V.sigma_wt(H, hours)

def rep_job(args):
    g, b, seed, models, tls, obs_n, hours, fly_seed_base = args
    _seed(seed); out = []
    for i, (P, pi, Ka) in enumerate(models):
        a = _sim(P, pi, Ka, int(tls[i]))
        if longest_silence(a) >= DEATH: out.append(None); continue
        nd = (int(tls[i]) // 86400) * 86400; y = R.sqw_symbols(a[:nd], theta=0.5, light_dark=False)[0]
        cells, nrun = cell_bins(y); rng = np.random.default_rng(fly_seed_base + 100 * b + i); flag = False
        for c in CELLS:
            k = obs_n[i][c]
            if len(cells[c]) > k: cells[c] = list(rng.choice(cells[c], size=k, replace=False))
            elif len(cells[c]) < k: flag = True
        out.append(dict(sigma=float(sigma_from_cells(cells, hours[i])), flag=flag, nrun=nrun.tolist()))
    return b, out

def models_from(fits):
    M = []
    for f in fits:
        X = np.array(f['best']['X']); x = X.sum(1); M.append((X / x[:, None], x / x.sum(), int(f['best']['Ka'])))
    return M

def block(g, X, tls, models, seed_of_b, fly_seed_base):
    from multiprocessing import Pool
    obs_n = []; hours = []; obs_sigma = []; obs_runs = []
    for x in X:
        y = R.sqw_symbols(x, theta=0.5, light_dark=False)[0]; cells, nrun = cell_bins(y)
        obs_n.append({c: len(cells[c]) for c in CELLS}); hours.append(len(y) / 60)
        obs_sigma.append(sigma_from_cells(cells, len(y) / 60)); obs_runs.append(nrun)
    jobs = [(g, b, seed_of_b(b), models, tls, obs_n, hours, fly_seed_base) for b in range(1, 200)]
    with Pool(4) as pool:
        reps = dict(pool.map(rep_job, jobs, chunksize=4))
    n = len(X); E = np.zeros(n); p = np.zeros(n); kept = np.zeros(n, int); dead = 0; flags = 0
    nullmed = np.zeros(n); runratio = np.zeros((n, 3))
    for i in range(n):
        v = [reps[b][i] for b in range(1, 200)]; dead += sum(r is None for r in v); v = [r for r in v if r is not None]
        sg = np.array([r['sigma'] for r in v]); kept[i] = len(sg); flags += sum(r['flag'] for r in v)
        nullmed[i] = np.median(sg); E[i] = obs_sigma[i] - nullmed[i]; p[i] = (1 + np.sum(sg >= obs_sigma[i])) / (len(sg) + 1)
        runratio[i] = np.mean([r['nrun'] for r in v], 0) / np.maximum(obs_runs[i], 1)
    return dict(E=E.tolist(), p=p.tolist(), kept=kept.tolist(), dead=int(dead), flags=int(flags), obs_sigma=list(map(float, obs_sigma)),
                nullmed=nullmed.tolist(), runratio=runratio.tolist())

def f2(per0dir, csdddir):
    import planB_20u as U
    t0 = time.time(); out = {}
    # ---- per0 DD: 20u fits and adequacy
    run = json.load(open(os.path.join(HERE, 'planB_20u_run.json')))
    fits = [json.load(open(os.path.join(HERE, 'planB_20u_fits', f'fly{i:02d}.json'))) for i in range(25)]
    Xp, _ = per0_traces(per0dir); _, _, _, tl = U.load_ref()
    cp = os.path.join(HERE, 'planB_20vF_per0_cache.json')
    if os.path.exists(cp): out['per0DD'] = json.load(open(cp))
    else:
        out['per0DD'] = block('per0DD', Xp, tl, models_from(fits), lambda b: 5000 + b, 25_000); json.dump(out['per0DD'], open(cp, 'w'))
    out['per0DD']['adequate'] = bool(run['adequacy']['PASS'])
    print(f"per0 DD null done ({time.time()-t0:.0f} s)", flush=True)
    # ---- CS DD: fits (20u stage 1 procedure, seeds 3000+100+fly), adequacy (4000+100+fly), nulls (5000+100+b)
    Xc, names = csdd_traces(csdddir); tlc = [len(x) for x in Xc]
    ck = os.path.join(HERE, 'planB_20vF_csdd_fits'); os.makedirs(ck, exist_ok=True)
    from multiprocessing import Pool
    todo = [(100 + i, (Xc[i] > 0).astype(np.int8)) for i in range(len(Xc)) if not os.path.exists(os.path.join(ck, f'fly{i:02d}.json'))]
    with Pool(4) as pool:
        for f in pool.imap_unordered(fit_fly, todo, chunksize=1):
            json.dump(f, open(os.path.join(ck, f"fly{f['fly'] - 100:02d}.json"), 'w')); print(f"  CS DD fitted fly {f['fly']-100} ({time.time()-t0:.0f} s)", flush=True)
    cfits = [json.load(open(os.path.join(ck, f'fly{i:02d}.json'))) for i in range(len(Xc))]; Mc = models_from(cfits)
    rw = []; rs = []
    for i, x in enumerate(Xc):
        ow, os_ = extract((x > 0).astype(np.int8), len(x)); _seed(4000 + 100 + i); ws = []; ss = []
        for k in range(20):
            a = _sim(*Mc[i][:2], Mc[i][2], tlc[i])
            if longest_silence(a) >= DEATH: continue
            w, s_ = extract(a, len(a)); ws.append(len(w)); ss.append(s_.sum())
        rw.append(np.mean(ws) / len(ow)); rs.append(np.mean(ss) / os_.sum())
    adq = 0.85 <= np.median(rw) <= 1.15 and 0.85 <= np.median(rs) <= 1.15
    print(f"CS DD adequacy: wake-bout ratio median {np.median(rw):.3f}, total sleep ratio median {np.median(rs):.3f} -> {'PASS' if adq else 'FAIL (NULL INADEQUATE)'}")
    out['CSDD_adequacy'] = dict(rw=list(map(float, rw)), rs=list(map(float, rs)), PASS=bool(adq),
                                selected=[(f['best']['Ka'], f['best']['Ki']) for f in cfits], min_eig=[f['best']['min_eig'] for f in cfits])
    if adq:
        out['CSDD'] = block('CSDD', Xc, tlc, Mc, lambda b: 5000 + 100 + b, 26_000); out['CSDD']['adequate'] = True
    # ---- F3 decision
    v20 = json.load(open(os.path.join(HERE, 'planB_20v.json')))
    blocks = [g for g in ('per0DD', 'CSDD') if g in out and out[g].get('adequate')]
    wil = {}
    for g in blocks:
        E = np.array(out[g]['E']); wil[g] = wilcoxon(E).pvalue
    ps = sorted((wil[g], g) for g in blocks); holm = {}; r_ = 0
    for k, (pv, g) in enumerate(ps): r_ = max(r_, min(1, (len(ps) - k) * pv)); holm[g] = r_
    for g in ('per0DD', 'CSDD'):
        if g not in blocks:
            print(f"{g}: NULL INADEQUATE (no verdict)" if g in out or g == 'CSDD' else f"{g}: not run"); continue
        o = out[g]; E = np.array(o['E']); p = np.array(o['p']); k = int(np.sum(p < 0.05)); bp = binom.sf(k - 1, len(p), 0.05)
        sup = holm[g] < 0.05 and np.median(E) > 0 and bp < 0.05
        chain_null = [r['obs'][0] - r['E_A'][0] for r in v20 if r['g'] == g and r['ok'][0]]
        orig = [r['obs'][0] for r in v20 if r['g'] == g and r['ok'][0]]
        rr = np.array(o['runratio'])
        print(f"{g:7s} n {len(E)} | sigma_WT obs (matched recomputation) median {np.median(o['obs_sigma']):.4f} (20v original {np.median(orig):.4f}) | "
              f"E_score median {np.median(E):+.4f} | Wilcoxon p {wil[g]:.3g} Holm {holm[g]:.3g} | flies p<.05 {k/len(p):.2f} (binom {bp:.3g}) -> "
              f"{'SUPPORTED' if sup else 'NOT SUPPORTED'}")
        print(f"        null median sigma_WT: scored reversible {np.median(o['nullmed']):.4f} vs 20v chain null {np.median(chain_null):.4f} | "
              f"death-excluded fly-replicates {o['dead']} | count-match flags {o['flags']} | run ratios sim/obs median Q {np.median(rr[:,0]):.2f} "
              f"W {np.median(rr[:,1]):.2f} S {np.median(rr[:,2]):.2f}")
        o['verdict'] = 'SUPPORTED' if sup else 'NOT SUPPORTED'; o['holm'] = holm[g]; o['wilcoxon'] = float(wil[g])
    out['seconds'] = time.time() - t0
    json.dump(out, open(os.path.join(HERE, 'planB_20vF_f2.json'), 'w'))

if __name__ == '__main__':
    if sys.argv[1] == 'f1': f1(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == 'f2': f2(sys.argv[2], sys.argv[3])
