"""T1-B synthetic gate (PRESPEC T1-B + T1-B1). 12-fly cohorts x 14 days; B-P (rhythm) and B-N (no rhythm).
Phase 1: one job per cohort. Phase 2: bootstrap replicates (B = 99) for cohorts 0-4 of each scenario.
Resumable. Usage: python -I gate_t1b.py <outdir> [nproc]  |  python -I gate_t1b.py <outdir> --report"""
import sys, os, json, pickle, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.stats import chi2
import pipeline as pl, sim_positions as sp, hmm_cat as hc, d1d, t1b

N_COH, N_FLY, DAYS, NBOOT, NBC = 20, 12, 14, 99, 5
SC = {'TBP': 1, 'TBN': 2}
DARK = t1b.PHASE == 1
ZT = np.arange(t1b.W) * 10 / 3600


def fly_data(scen, c, f):
    code = SC[scen]
    ss = np.random.SeedSequence([9001, code, c]).spawn(N_FLY)[f]
    st, rd, t = sp.simulate_fly(scen, DAYS, ss)
    x, _ = pl.spike_filter(rd); valid = np.ones(len(x), bool)
    y = d1d.step_symbols(x, valid); n = len(y) // d1d.SPD * d1d.SPD; y = y[:n]
    return d1d.to_days(y), d1d.labels_per_sample(x, t, valid, n)


def c2_starts(c, code):
    out = [(np.zeros(4), np.zeros(4))]
    for r in range(2):
        rng = np.random.default_rng(np.random.SeedSequence([8001, code, c, 2, r]))
        out.append((rng.normal(0, 0.3, 4), rng.normal(0, 0.3, 4)))
    return out


def cohort(args):
    scen, c, outdir = args; code = SC[scen]
    f = os.path.join(outdir, f'{scen}_coh{c:02d}.json')
    if os.path.exists(f):
        return f
    flies, c0_models, ll0, fallback = [], [], 0.0, 0
    for k in range(N_FLY):
        seqs, lab = fly_data(scen, c, k)
        m, ll, _, _ = hc.fit(seqs, 3, 4, 2, np.random.SeedSequence([9101, code, c, k]), n_restarts=3, max_iter=d1d.MAXIT, tol=d1d.TOL)
        s, ok = t1b.sleep_state(m, seqs, lab); fallback += (not ok)
        flies.append(dict(F=t1b.from_c0(m, s), seqs=t1b.bin_seqs(seqs))); c0_models.append((m, s)); ll0 += ll
    fit = t1b.fit_cohort(flies, 2, c2_starts(c, code))
    LR = 2 * (fit['ll'] - ll0); p = float(chi2.sf(max(LR, 0), 8))
    h = t1b.harmonic_profile(fit['h'], 2); hd = h[DARK]
    r = dict(scenario=scen, cohort=c, ll_C0=ll0, ll_C2=fit['ll'], LR=LR, p_chi2=p, reject=bool(p < 0.05),
             h=list(fit['h']), g=list(fit['g']), peak_dark_zt=float(ZT[DARK][np.argmax(hd)]),
             range_ratio=float(np.exp(hd.max() - hd.min()) / 3.0), sleep_state_fallbacks=fallback, iters=fit['iters'])
    if c < NBC:
        with open(os.path.join(outdir, f'{scen}_coh{c:02d}_c0.pkl'), 'wb') as fh:
            pickle.dump(c0_models, fh)
    json.dump(r, open(f + '.tmp', 'w')); os.replace(f + '.tmp', f)
    return f


def boot_rep(args):
    scen, c, b, outdir = args; code = SC[scen]
    f = os.path.join(outdir, f'boot_{scen}_coh{c:02d}_r{b:02d}.json')
    if os.path.exists(f):
        return f
    c0_models = pickle.load(open(os.path.join(outdir, f'{scen}_coh{c:02d}_c0.pkl'), 'rb'))
    rng = np.random.default_rng(np.random.SeedSequence([9201, code, c, b]))
    ph = ((np.arange(d1d.SPD) * 0.4) % 86400 >= 43200).astype(np.int64)
    flies, ll0 = [], 0.0
    for m, s in c0_models:
        seqs = hc.simulate_fast(m, [(d1d.SPD, ph)] * DAYS, [np.zeros(d1d.SPD, bool)] * DAYS, rng)
        m2, ll, _, _ = hc.fit(seqs, 3, 4, 2, np.random.SeedSequence([9201, code, c, b, 1]), n_restarts=1,
                              max_iter=d1d.MAXIT, tol=d1d.TOL, inits=[m])
        flies.append(dict(F=t1b.from_c0(m2, s), seqs=t1b.bin_seqs(seqs))); ll0 += ll
    fit = t1b.fit_cohort(flies, 2, [(np.zeros(4), np.zeros(4))])
    json.dump(dict(LR=2 * (fit['ll'] - ll0)), open(f + '.tmp', 'w')); os.replace(f + '.tmp', f)
    return f


def report(outdir):
    R = {s: [json.load(open(os.path.join(outdir, f'{s}_coh{c:02d}.json'))) for c in range(N_COH)
             if os.path.exists(os.path.join(outdir, f'{s}_coh{c:02d}.json'))] for s in SC}
    ok = {}
    P = R['TBP']
    if P:
        rej = sum(r['reject'] for r in P); pk = sum(abs(r['peak_dark_zt'] - 12.5) <= 2 for r in P)
        med = float(np.median([r['range_ratio'] for r in P]))
        print(f"B-P ({len(P)}/20): reject {rej} (need >= 16); dark-phase peak within 2 h of ZT 12.5: {pk} (need >= 16); "
              f"median range ratio {med:.2f} (need 0.67-1.5)")
        print('   p_chi2', [f"{r['p_chi2']:.2g}" for r in P]); print('   peak ZT', [round(r['peak_dark_zt'], 1) for r in P])
        ok['B-P'] = rej >= 16 and pk >= 16 and 0.67 <= med <= 1.5
    N = R['TBN']
    if N:
        rej = sum(r['reject'] for r in N)
        print(f"B-N ({len(N)}/20): reject {rej} (allowed <= 2)"); print('   p_chi2', [f"{r['p_chi2']:.2g}" for r in N])
        ok['B-N'] = rej <= 2
    print('sleep-state fallbacks per cohort:', [r['sleep_state_fallbacks'] for r in P + N])
    for s in SC:
        for r in R[s][:NBC]:
            reps = [json.load(open(os.path.join(outdir, f'boot_{s}_coh{r["cohort"]:02d}_r{b:02d}.json')))['LR']
                    for b in range(NBOOT) if os.path.exists(os.path.join(outdir, f'boot_{s}_coh{r["cohort"]:02d}_r{b:02d}.json'))]
            if reps:
                pb = (1 + sum(x >= r['LR'] for x in reps)) / (1 + len(reps))
                print(f"  bootstrap {s} coh{r['cohort']:02d}: {len(reps)}/{NBOOT} reps, p_boot {pb:.3f} vs p_chi2 {r['p_chi2']:.3g}")
    if all(len(R[s]) == N_COH for s in SC):
        print('T1-B GATE:', {k: 'PASS' if v else 'FAIL' for k, v in ok.items()}, '->', 'PASS' if all(ok.values()) else 'FAIL')


if __name__ == '__main__':
    outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
    if '--report' in sys.argv:
        report(outdir); sys.exit()
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    from multiprocessing import Pool
    with Pool(nproc) as p:
        for fn in p.imap_unordered(cohort, [(s, c, outdir) for c in range(N_COH) for s in SC]):
            print('done', os.path.basename(fn), flush=True)
        jobs = [(s, c, b, outdir) for c in range(NBC) for s in SC for b in range(NBOOT)]
        for fn in p.imap_unordered(boot_rep, jobs):
            print('done', os.path.basename(fn), flush=True)
    report(outdir)
