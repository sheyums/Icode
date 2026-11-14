"""T1-C1 synthetic gate (PRESPEC T1-C1): fresh TBP (20) / TBN (60) cohorts, 12 flies x 12 days, P1 scoring, rotation
null with the first E h after each light switch excluded; E = 2 (decision), 1 and 3 (sensitivity), same data.
Usage: python -I gate_t1c1.py <outdir> [nproc] | python -I gate_t1c1.py <outdir> --report"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl, sim_positions as sp, t1c

N_FLY, DAYS = 12, 12; COH = {'TBP': 20, 'TBN': 60}; SC = {'TBP': 1, 'TBN': 2}; ES = (2.0, 1.0, 3.0)


def cohort(args):
    scen, c, outdir = args; code = SC[scen]; f = os.path.join(outdir, f'{scen}_coh{c:02d}.json')
    if os.path.exists(f):
        return f
    B = []
    for k in range(N_FLY):
        st, rd, t = sp.simulate_fly(scen, DAYS, np.random.SeedSequence([9501, code, c]).spawn(N_FLY)[k])
        x, _ = pl.spike_filter(rd)
        B.append(t1c.fly_bins(t1c.labels_from_activity(pl.activity_1s(x, t, t >= 86400)), 0.0))
    out = dict(scenario=scen, cohort=c)
    for E in ES:
        r = t1c.rotation_test(B, np.random.SeedSequence([11201, code, c]), excl_h=E)
        r.update(reject=bool(r['p'] < 0.05), peaks_h=t1c.within_phase_peak(r['h']))
        out[f'E{E:g}'] = r
    json.dump(out, open(f + '.tmp', 'w')); os.replace(f + '.tmp', f)
    return f


def report(outdir):
    R = {s: [json.load(open(os.path.join(outdir, f'{s}_coh{c:02d}.json'))) for c in range(COH[s])
             if os.path.exists(os.path.join(outdir, f'{s}_coh{c:02d}.json'))] for s in SC}
    P, N = R['TBP'], R['TBN']; ok = {}
    for E in ES:
        k = f'E{E:g}'; tag = 'DECISION' if E == 2.0 else 'sensitivity'
        rp = sum(r[k]['reject'] for r in P); pk = sum(abs(r[k]['peaks_h']['dark_peak_zt'] - 12.5) <= 2 for r in P)
        rn = sum(r[k]['reject'] for r in N)
        print(f"E = {E:g} h [{tag}]: B-P ({len(P)}/20) reject {rp}, peak within 2 h {pk} | B-N ({len(N)}/60) reject {rn} "
              f"({rn / max(len(N), 1):.3f}); null p mean {np.mean([r[k]['p'] for r in N]) if N else float('nan'):.2f}")
        if E == 2.0:
            ok = {'B-P': rp >= 16 and pk >= 16, 'B-N': rn <= 5}
    if len(P) == 20 and len(N) == 60:
        print('T1-C1 GATE (E = 2 h):', {k: 'PASS' if v else 'FAIL' for k, v in ok.items()}, '->', 'PASS' if all(ok.values()) else 'FAIL')


if __name__ == '__main__':
    outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
    if '--report' in sys.argv:
        report(outdir); sys.exit()
    from multiprocessing import Pool
    jobs = [(s, c, outdir) for s in ('TBP', 'TBN') for c in range(COH[s])]
    with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as p:
        for fn in p.imap_unordered(cohort, jobs):
            print('done', os.path.basename(fn), flush=True)
    report(outdir)
