"""T1 synthetic gate (T-G1 recovery, T-G2 false positive, T-G3 spurious states). Resumable, one JSON per fly.
Usage: python3 -I gate_t1.py <outdir> [nproc] [--scenarios TG1,TG2]  |  python3 -I gate_t1.py <outdir> --report
T1-A: TG3 is deferred; pass --scenarios TG1,TG2 to run the scenarios that remain specified."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl, sim_positions as sp, d1, t1, iohmm as io, hmm_cat as hc

N_FLIES, DAYS = 20, 14
SCEN = {'TG1': 1, 'TG2': 2, 'TG3': 3}


def kstar_step(seqs, fly_idx):
    nd = len(seqs); odd = list(range(0, nd, 2)); even = list(range(1, nd, 2))
    curve = d1.heldout_gain_curve(seqs, fly_idx, d1.KS, 1001, 10, odd, even)
    elig = [K for K in d1.KS if not any(d1.degenerate(f[0][2]) for f in curve[K]['fits'])]
    return (max(elig, key=lambda K: curve[K]['heldout']) if elig else None), {K: curve[K]['heldout'] for K in d1.KS}


def circ_h(a, b):
    d = abs(a - b) % 24; return min(d, 24 - d)


def one(args):
    scen, i, outdir = args
    f = os.path.join(outdir, f'{scen}_fly{i:02d}.json')
    if os.path.exists(f):
        return f
    code = SCEN[scen]; fly_idx = 100 * code + i
    st, rd, t = sp.simulate_fly(scen, DAYS, np.random.SeedSequence([6001, code, i]))
    x, _ = pl.spike_filter(rd); valid = np.ones(len(x), bool)
    y = sp.window_symbols(x, valid); n = len(y) // io.WPD * io.WPD; y = y[:n]
    r = dict(scenario=scen, fly=i)
    if scen == 'TG3':
        r0 = d1.run_fly(x, t, valid, fly_idx=fly_idx, log=lambda s: None); r0.pop('argmax_full', None)
        Kstar = r0['Kstar']
        r['M0_d1'] = dict(Kstar=Kstar, C1=r0.get('C1', {}).get('pass_'),
                          spurious=bool(r0.get('C1', {}).get('pass_') and any(s['candidate'] and s['C2'] and s['C3'] for s in r0.get('states', {}).values())))
    else:
        Kstar, held_k = kstar_step(d1.to_days(y), fly_idx); r['heldout_K_step'] = held_k
    r['Kstar_M0'] = Kstar
    seqs_w = io.to_days_w(y)
    held, fits = t1.clock_selection(seqs_w, Kstar, fly_idx)
    Mstar = max(held, key=held.get); r['heldout_models'] = held; r['Mstar'] = int(Mstar)
    if scen == 'TG1':
        full, beta, ll, lls = io.fit(seqs_w, Kstar, 4, Mstar, np.random.SeedSequence([5001, fly_idx, Kstar, Mstar, 9]), n_restarts=10)
        am = np.concatenate(hc.posterior(full, seqs_w)).argmax(1)
        tw = (st[:n * 25].reshape(n, 25) == sp.SLEEP).mean(1) >= 0.5
        jac = [float((tw & (am == k)).sum() / max((tw | (am == k)).sum(), 1)) for k in range(Kstar)]
        ks = int(np.argmax(jac)); prof = 1 - full[1][:, ks, ks]; true = sp.true_sleep_exit_profile(scen)
        zt = np.arange(io.WPD) * 10 / 3600
        r['TG1'] = dict(sleep_state=ks, jaccard=jac[ks], peak_fit=float(zt[np.argmax(prof)]), peak_true=float(zt[np.argmax(true)]),
                        peak_dist_h=circ_h(zt[np.argmax(prof)], zt[np.argmax(true)]),
                        ptt_ratio=float((prof.max() / prof.min()) / (true.max() / true.min())))
    if scen == 'TG3':
        rc = t1.run_d1_clock(x, t, valid, fly_idx, Mstar)
        r['Mstar_d1'] = dict(Kstar=rc.get('Kstar'), C1=rc.get('C1', {}).get('pass_'),
                             spurious=bool(rc.get('C1', {}).get('pass_') and any(s['candidate'] and s['C2'] and s['C3'] for s in rc.get('states', {}).values())),
                             detail=rc)
    with open(f + '.tmp', 'w') as fh:
        json.dump(r, fh, default=lambda o: o.tolist() if hasattr(o, 'tolist') else (int(o) if isinstance(o, np.integer) else str(o)))
    os.replace(f + '.tmp', f)
    return f


def report(outdir):
    R = {s: [json.load(open(os.path.join(outdir, f'{s}_fly{i:02d}.json'))) for i in range(N_FLIES)
             if os.path.exists(os.path.join(outdir, f'{s}_fly{i:02d}.json'))] for s in SCEN}
    ok = {}
    g = R['TG1']
    if g:
        sel = sum(r['Mstar'] > 0 for r in g); pk = sum(r['TG1']['peak_dist_h'] <= 2 for r in g)
        med = float(np.median([r['TG1']['ptt_ratio'] for r in g]))
        print(f"TG1 ({len(g)}/20): M* {[r['Mstar'] for r in g]}; clock chosen {sel}; peak within 2 h {pk}; median ptt ratio {med:.2f}")
        print('   peak dist h', [round(r['TG1']['peak_dist_h'], 1) for r in g])
        ok['TG1'] = sel >= 16 and pk >= 16 and 0.67 <= med <= 1.5
    g = R['TG2']
    if g:
        sel = sum(r['Mstar'] > 0 for r in g)
        print(f"TG2 ({len(g)}/20): M* {[r['Mstar'] for r in g]}; clock chosen {sel} (PASS iff <= 3)"); ok['TG2'] = sel <= 3
    g = R['TG3']
    if g:
        s1 = sum(r['Mstar_d1']['spurious'] for r in g); s0 = sum(r['M0_d1']['spurious'] for r in g)
        print(f"TG3 ({len(g)}/20): M* {[r['Mstar'] for r in g]}; spurious promoted under M* {s1} (PASS iff <= 1); under M0 {s0}")
        ok['TG3'] = s1 <= 1
    if all(len(R[s]) == N_FLIES for s in SCEN if R[s] or s != 'TG3') and ok:
        if not R['TG3']:
            print('TG3 deferred (T1-A): verdict below covers TG1 and TG2 only')
        print('T1 GATE:', {k: ('PASS' if v else 'FAIL') for k, v in ok.items()}, '->', 'PASS' if all(ok.values()) else 'FAIL')


if __name__ == '__main__':
    outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
    if '--report' in sys.argv:
        report(outdir); sys.exit()
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 and not sys.argv[2].startswith('--') else 4
    scen = ('TG3', 'TG1', 'TG2')
    if '--scenarios' in sys.argv:
        scen = tuple(x for x in sys.argv[sys.argv.index('--scenarios') + 1].split(',') if x)
        assert all(x in SCEN for x in scen), scen
    jobs = [(s, i, outdir) for s in scen for i in range(N_FLIES)]    # longest first
    from multiprocessing import Pool
    with Pool(nproc) as p:
        for f in p.imap_unordered(one, jobs):
            print('done', os.path.basename(f), flush=True)
    report(outdir)
