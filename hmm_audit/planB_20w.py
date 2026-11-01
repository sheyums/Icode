"""PRESPEC 20w: circadian "clock face" (descriptive). 12 bins of 2 h of clock time (arrival minute), pooled transitions
over flies and days per genotype, pooled first-order chain per bin -> J, A, EPR (sqw_ref.model_energetics), E3-E1
(sqw_ref.landscape), pooled c_obs. 2000 fly-bootstrap replicates (seed 2026). Inf ranks highest (amendment A item 4).
Usage: python3 -I planB_20w.py <data root> <CS DD folder>"""
import sys, os, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sqw_ref as R
SETS = {'CS': ('CS_activity_files_1sec_res', True), 'per0LD': ('per0LD_activity_1sec_res', True),
        'per0DD': ('per0DD_activity_files_1sec_res', False), 'CSDD': (None, False)}
NB = 12

def files(d): return sorted(f for f in glob.glob(os.path.join(d, '**', '*'), recursive=True) if os.path.isfile(f))

def fly_counts(path, ld):
    y, _ = R.sqw_symbols(R.read_counts(path), theta=0.5, light_dark=ld); y = np.asarray(y)
    t = np.arange(1, len(y)); ok = (y[:-1] >= 0) & (y[1:] >= 0); b = (t % 1440) // 120
    N = np.zeros((NB, 3, 3)); np.add.at(N, (b[ok], y[:-1][ok], y[1:][ok]), 1)
    on = ok & (y[:-1] < 2) & (y[1:] == 2); wk = ok & (y[:-1] == 2) & (y[1:] < 2)
    E = np.zeros((NB, 4))                      # onsets, onsets from Q, wakings, wakings to Q
    np.add.at(E[:, 0], b[on], 1); np.add.at(E[:, 1], b[on & (y[:-1] == 0)], 1)
    np.add.at(E[:, 2], b[wk], 1); np.add.at(E[:, 3], b[wk & (y[1:] == 0)], 1)
    return N, E

def quantities(N, E):
    out = np.full((NB, 5), np.nan)            # J, A, EPR, E3-E1, c_obs
    for k in range(NB):
        rs = N[k].sum(1)
        if np.any(rs == 0): continue
        P = N[k] / rs[:, None]; e = R.model_energetics(P); ls = R.landscape(P)
        c = E[k, 1] / E[k, 0] - E[k, 3] / E[k, 2] if E[k, 0] > 0 and E[k, 2] > 0 else np.nan
        out[k] = [e['J'], e['A'], e['EPR'], ls['E3_minus_E1'], c]
    return out

def rng_(v):
    v = v[np.isfinite(v)] if not np.any(np.isinf(v)) else v
    return np.inf if np.any(np.isinf(v)) else float(np.nanmax(v) - np.nanmin(v))

def main(root, csdd):
    rng = np.random.default_rng(2026); res = {}
    for g, (sub, ld) in SETS.items():
        fl = files(csdd if g == 'CSDD' else os.path.join(root, sub))
        C = [fly_counts(f, ld) for f in fl]; Ns = np.array([c[0] for c in C]); Es = np.array([c[1] for c in C])
        pt = quantities(Ns.sum(0), Es.sum(0)); bs = []; rgJ = []; rgS = []
        for b in range(2000):
            idx = rng.integers(len(fl), size=len(fl)); q = quantities(Ns[idx].sum(0), Es[idx].sum(0)); bs.append(q)
            rgJ.append(rng_(q[:, 0])); rgS.append(rng_(q[:, 2]))
        bs = np.array(bs); lo = np.nanpercentile(bs, 2.5, axis=0); hi = np.nanpercentile(bs, 97.5, axis=0)
        res[g] = dict(n=len(fl), point=pt.tolist(), lo=lo.tolist(), hi=hi.tolist(),
                      rangeJ=[rng_(pt[:, 0])] + list(np.percentile(rgJ, [2.5, 97.5])),
                      rangeEPR=[rng_(pt[:, 2])] + list(np.percentile(rgS, [2.5, 97.5])),
                      maxbinJ=int(np.nanargmax(pt[:, 0])), maxbinEPR=int(np.nanargmax(pt[:, 2])))
        print(f"\n{g} (n {len(fl)} flies; narrowest block: one 2-h clock bin pooled over flies and days)")
        print("  bin(ZT/CT h)  J cyc/h [95%]            A nats         pooled-chain EPR [95%]    E3-E1     c_obs")
        for k in range(NB):
            print(f"  {2*k:2d}-{2*k+2:2d}   {pt[k,0]:+.3f} [{lo[k,0]:+.3f},{hi[k,0]:+.3f}]   {pt[k,1]:+.2f}   "
                  f"{pt[k,2]:.3f} [{lo[k,2]:.3f},{hi[k,2]:.3f}]   {pt[k,3]:+.2f}   {pt[k,4]:+.3f}")
        r = res[g]
        print(f"  range of J {r['rangeJ'][0]:.3f} [{r['rangeJ'][1]:.3f}, {r['rangeJ'][2]:.3f}], max bin {2*r['maxbinJ']}-{2*r['maxbinJ']+2} h | "
              f"range of EPR {r['rangeEPR'][0]:.3f} [{r['rangeEPR'][1]:.3f}, {r['rangeEPR'][2]:.3f}], max bin {2*r['maxbinEPR']}-{2*r['maxbinEPR']+2} h")
    ref = res['per0DD']
    for g in ('CS', 'per0LD', 'CSDD'):
        for key in ('rangeJ', 'rangeEPR'):
            above = res[g][key][1] > ref[key][2]
            print(f"{g:7s} {key}: band [{res[g][key][1]:.3f}, {res[g][key][2]:.3f}] vs per0DD [{ref[key][1]:.3f}, {ref[key][2]:.3f}] -> "
                  f"{'RHYTHMIC (band entirely above per0 DD)' if above else 'not distinguishable from per0 DD'}")
    json.dump(res, open(os.path.join(HERE, 'planB_20w.json'), 'w'))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
