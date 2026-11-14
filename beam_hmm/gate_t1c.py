"""T1-C synthetic gate (PRESPEC T1-C): TBP/TBN cohorts (T1-B's truth and data seeds), P1 scoring, rotation null.
Usage: python -I gate_t1c.py <outdir> [nproc] | python -I gate_t1c.py <outdir> --report"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl, sim_positions as sp, t1c

N_COH, N_FLY, DAYS = 20, 12, 14; SC = {'TBP': 1, 'TBN': 2}


def cohort(args):
    scen, c, outdir = args; code = SC[scen]; f = os.path.join(outdir, f'{scen}_coh{c:02d}.json')
    if os.path.exists(f):
        return f
    B = []
    for k in range(N_FLY):
        st, rd, t = sp.simulate_fly(scen, DAYS, np.random.SeedSequence([9001, code, c]).spawn(N_FLY)[k])
        x, _ = pl.spike_filter(rd); v = t >= 86400
        B.append(t1c.fly_bins(t1c.labels_from_activity(pl.activity_1s(x, t, v)), 0.0))
    r = t1c.rotation_test(B, np.random.SeedSequence([11101, code, c]))
    r.update(scenario=scen, cohort=c, reject=bool(r['p'] < 0.05), peaks_h=t1c.within_phase_peak(r['h']), peaks_g=t1c.within_phase_peak(r['g']))
    json.dump(r, open(f + '.tmp', 'w')); os.replace(f + '.tmp', f)
    return f


def report(outdir):
    R = {s: [json.load(open(os.path.join(outdir, f'{s}_coh{c:02d}.json'))) for c in range(N_COH)
             if os.path.exists(os.path.join(outdir, f'{s}_coh{c:02d}.json'))] for s in SC}
    P, N = R['TBP'], R['TBN']; ok = {}
    if P:
        rej = sum(r['reject'] for r in P); pk = sum(abs(r['peaks_h']['dark_peak_zt'] - 12.5) <= 2 for r in P)
        print(f"B-P ({len(P)}/20): reject {rej} (need >= 16); wake-up dark peak within 2 h of ZT 12.5: {pk} (need >= 16)")
        print('   p', [r['p'] for r in P]); print('   dark peak ZT (h)', [round(r['peaks_h']['dark_peak_zt'], 1) for r in P])
        ok['B-P'] = rej >= 16 and pk >= 16
    if N:
        rej = sum(r['reject'] for r in N); print(f"B-N ({len(N)}/20): reject {rej} (allowed <= 2)"); print('   p', [r['p'] for r in N])
        ok['B-N'] = rej <= 2
    if len(P) == N_COH and len(N) == N_COH:
        print('T1-C GATE:', {k: 'PASS' if v else 'FAIL' for k, v in ok.items()}, '->', 'PASS' if all(ok.values()) else 'FAIL')


if __name__ == '__main__':
    outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
    if '--report' in sys.argv:
        report(outdir); sys.exit()
    from multiprocessing import Pool
    with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as p:
        for fn in p.imap_unordered(cohort, [(s, c, outdir) for c in range(N_COH) for s in SC]):
            print('done', os.path.basename(fn), flush=True)
    report(outdir)
