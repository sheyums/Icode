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

if __name__ == '__main__':
    if sys.argv[1] == 'f1': f1(sys.argv[2], sys.argv[3])
