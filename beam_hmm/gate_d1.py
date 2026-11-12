"""D1 synthetic gate (G1 planted state, G2 negative control). Resumable: one JSON per synthetic fly.
Usage: python3 -I gate_d1.py <outdir> [nproc]   then   python3 -I gate_d1.py <outdir> --report"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as pl, sim_positions as sp, d1

N_FLIES, DAYS = 20, 14


def one(args):
    scen, i, outdir = args
    f = os.path.join(outdir, f'{scen}_fly{i:02d}.json')
    if os.path.exists(f):
        return f
    code = 1 if scen == 'G1' else 2
    st, rd, t = sp.simulate_fly(scen, DAYS, np.random.SeedSequence([3001, code, i]))
    x, _ = pl.spike_filter(rd)
    valid = np.ones(len(x), bool)
    r = d1.run_fly(x, t, valid, fly_idx=code * 100 + i, log=lambda s: None)
    # truth-side evaluation (D1-A item 6)
    am = r.pop('argmax_full', None)
    if scen == 'G1' and am is not None:
        n = len(am)
        tw = (st[:n * 25].reshape(n, 25) == sp.BOUNDARY).mean(1) >= 0.5
        jac = []
        for k in range(r['Kstar']):
            w = am == k; inter = (w & tw).sum(); uni = (w | tw).sum(); jac.append(float(inter / uni) if uni else 0.0)
        kb = int(np.argmax(jac)); s = r['states'][kb]
        r['planted'] = dict(state=kb, jaccard=jac[kb], candidate=s['candidate'], C1=r['C1']['pass_'], C2=s['C2'], C3=s['C3'],
                            C4_pass=s['C4_pass'], promoted=bool(s['candidate'] and r['C1']['pass_'] and s['C2'] and s['C3']))
    if scen == 'G2' and 'states' in r:
        cands = [k for k, s in r['states'].items() if s['candidate']]
        r['spurious'] = dict(n_candidates=len(cands), C1=bool(r['C1']['pass_'] and len(cands) > 0),
                             promoted=bool(r['C1']['pass_'] and any(r['states'][k]['C2'] and r['states'][k]['C3'] for k in cands)))
    with open(f + '.tmp', 'w') as fh:
        json.dump(r, fh, default=lambda o: o.tolist() if hasattr(o, 'tolist') else (int(o) if isinstance(o, np.integer) else str(o)))
    os.replace(f + '.tmp', f)
    return f


def report(outdir):
    rows = {}
    for scen in ('G1', 'G2'):
        R = []
        for i in range(N_FLIES):
            f = os.path.join(outdir, f'{scen}_fly{i:02d}.json')
            if os.path.exists(f):
                R.append(json.load(open(f)))
        rows[scen] = R
    g1, g2 = rows['G1'], rows['G2']
    print(f'G1: {len(g1)}/{N_FLIES} flies done')
    if g1:
        ks = [r['Kstar'] for r in g1]; pl_ = [r.get('planted', {}) for r in g1]
        prom = sum(p.get('promoted', False) for p in pl_); k4 = sum(k == 4 for k in ks)
        print('  K*:', ks)
        for key in ('candidate', 'C1', 'C2', 'C3', 'C4_pass', 'promoted'):
            print(f'  planted {key}: {sum(p.get(key, False) for p in pl_)}/{len(g1)}')
        print('  planted jaccard:', [round(p.get('jaccard', 0), 2) for p in pl_])
        print(f'  G1 PASS iff promoted >= 16/20 and K*=4 >= 14/20: promoted {prom}, K*=4 {k4}')
        prom5 = sum(p.get('promoted', False) and p.get('jaccard', 0) >= 0.5 for p in pl_)
        print(f'  [D1-A2, secondary] promoted AND matched state Jaccard >= 0.5: {prom5}/{len(g1)}')
    print(f'G2: {len(g2)}/{N_FLIES} flies done')
    if g2:
        ks = [r['Kstar'] for r in g2]; sp_ = [r.get('spurious', {}) for r in g2]
        c1 = sum(s.get('C1', False) for s in sp_); pr = sum(s.get('promoted', False) for s in sp_)
        print('  K*:', ks)
        print(f'  G2 PASS iff spurious C1 <= 2/20 and promoted <= 1/20: C1 {c1}, promoted {pr}')
    if len(g1) == N_FLIES and len(g2) == N_FLIES:
        G1 = prom >= 16 and k4 >= 14; G2 = c1 <= 2 and pr <= 1
        print(f'GATE: G1 {"PASS" if G1 else "FAIL"}, G2 {"PASS" if G2 else "FAIL"} -> {"PASS" if G1 and G2 else "FAIL"}')


if __name__ == '__main__':
    outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
    if '--report' in sys.argv:
        report(outdir); sys.exit()
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    jobs = [(s, i, outdir) for i in range(N_FLIES) for s in ('G1', 'G2')]
    from multiprocessing import Pool
    with Pool(nproc) as p:
        for f in p.imap_unordered(one, jobs):
            print('done', os.path.basename(f), flush=True)
    report(outdir)
