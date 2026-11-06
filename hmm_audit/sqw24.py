# PRESPEC 24 (c5d7349): 23E10>Kir2.1 (active) vs 23E10>Kir2.1-NC, same background. Readouts as 20y.
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from scipy.stats import mannwhitneyu
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R
D = '/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/data_23E10'
G = {'kir': D + '/kir', 'nc': D + '/nc'}; PH = ['light', 'dark']; TR = [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]; SYM = 'QWS'
def trim(x):
    """Sheyum's death rule: first run of >= 12 h without activity; cut 24 h before its start. Then first 8 whole days."""
    act = x > 0; dead_at = None; run = 0
    for i, a in enumerate(act):
        run = 0 if a else run + 1
        if run == 43200: dead_at = i - 43200 + 1; break
    if dead_at is not None: x = x[:max(0, dead_at - 86400)]
    nd = min(len(x) // 86400, 8); return x[:nd * 86400], nd, dead_at
def one(a):
    g, f = a; x0 = np.loadtxt(f).astype(int); x, nd, dead_at = trim(x0)
    out = dict(g=g, f=os.path.basename(f), days=nd, died=dead_at is not None, died_day=None if dead_at is None else dead_at / 86400)
    if nd < 2: out['excluded'] = True; return out
    out['excluded'] = False
    y, refs = R.sqw_symbols(x, 0.5, light_dark=True); ph = R.phases(len(y), '12h')
    for tag, s in (('main', 0), ('nod0', 1440)):
        yy, pp = y[s:], ph[s:]
        if len(yy) < 1440: out[tag] = None; continue
        _, _, P = R.mc_fit(yy, pp, 2); oq = R.observable_quantities(yy, pp, 2, P); res = []
        for p in range(2):
            ok = not np.isnan(refs[PH[p]]); m = pp == p; h = m.sum() / 60
            d = dict(oq[p], ok=ok, sleep=float((yy[m] == 2).sum() / h), P={f'{SYM[i]}{SYM[j]}': float(P[p][i, j]) for i, j in TR})
            d['cycle'] = abs(d['J_chain']) >= 1e-12; res.append(d)
        out[tag] = res
    return out
def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); adj = np.empty(len(ps)); run = 0
    for j, i in enumerate(o): run = max(run, min(1, (len(ps) - j) * ps[i])); adj[i] = run
    return adj
rng = np.random.default_rng(24)
def diffci(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    bs = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(2000)]
    return float(np.median(a) - np.median(b)), [float(x) for x in np.percentile(bs, [2.5, 97.5])]
def vals(F, g, tag, p, key, cyc=False):
    v = []
    for d in F:
        if d['g'] != g or d['excluded'] or d[tag] is None: continue
        r = d[tag][p]
        if not r['ok'] or (cyc and not r['cycle']): continue
        x = r['P'][key] if key in r['P'] else r[key]
        if key in ('E21', 'E31') and x is not None: x = -x   # E_W - E_Q, E_S - E_Q
        if x is not None and np.isfinite(x): v.append(x)
    return v
def test(F, tag, key, p, alt='two-sided', cyc=False):
    a, b = vals(F, 'kir', tag, p, key, cyc), vals(F, 'nc', tag, p, key, cyc); u = mannwhitneyu(a, b, alternative=alt)
    md, ci = diffci(a, b); return dict(key=key, phase=PH[p], n=(len(a), len(b)), med=(float(np.median(a)), float(np.median(b))), diff=md, ci=ci,
                                     p=float(u.pvalue), cliff=float(2 * u.statistic / (len(a) * len(b)) - 1))
if __name__ == '__main__':
    jobs = [(g, f) for g, d in G.items() for f in sorted(glob.glob(d + '/**/*.txt', recursive=True))]
    with Pool(4) as pool: F = pool.map(one, jobs)
    res = {}
    for g in G:
        Fg = [d for d in F if d['g'] == g]
        print(g, 'flies', len(Fg), '| died', sum(d['died'] for d in Fg), '| excluded (<2 days)', sum(d['excluded'] for d in Fg),
              '| days used', np.bincount([d['days'] for d in Fg]).tolist())
    for tag in ('main', 'nod0'):
        R_ = {}
        R_['primary'] = test(F, tag, 'J', 1, alt='greater')
        fa = [test(F, tag, k, p) for k in ('c_obs', 'J', 'EPR') for p in (0, 1)]
        for t, a in zip(fa, holm([t['p'] for t in fa])): t['holm'] = float(a)
        sl = [test(F, tag, 'sleep', p) for p in (0, 1)]
        for t, a in zip(sl, holm([t['p'] for t in sl])): t['holm'] = float(a)
        trs = [test(F, tag, f'{SYM[i]}{SYM[j]}', p) for i, j in TR for p in (0, 1)]
        for t, a in zip(trs, holm([t['p'] for t in trs])): t['holm'] = float(a)
        ls = [test(F, tag, k, p, cyc=True) for k in ('E21', 'E31') for p in (0, 1)]
        for t, a in zip(ls, holm([t['p'] for t in ls])): t['holm'] = float(a)
        R_.update(familyA=fa, sleep=sl, transitions=trs, landscape=ls); res[tag] = R_
        pr = R_['primary']
        print(f"\n[{tag}] PRIMARY J dark one-sided: Kir {pr['med'][0]:+.3f} vs NC {pr['med'][1]:+.3f} (n {pr['n']}) diff {pr['diff']:+.3f} [{pr['ci'][0]:+.3f},{pr['ci'][1]:+.3f}] p {pr['p']:.4g} Cliff {pr['cliff']:+.2f}")
        for name, L in (('familyA', fa), ('sleep', sl), ('transitions', trs), ('landscape', ls)):
            for t in L:
                print(f"  {name:11s} {t['key']:5s} {t['phase']:5s} Kir {t['med'][0]:+.4f} NC {t['med'][1]:+.4f} diff {t['diff']:+.4f} [{t['ci'][0]:+.4f},{t['ci'][1]:+.4f}] p {t['p']:.3g} Holm {t['holm']:.3g} d {t['cliff']:+.2f}")
    # across backgrounds: 20y effects (dfb, vnc vs empty) for the same readouts
    Y = [d for d in json.load(open('y_kir.json')) if not d['dead']]
    def yv(g, p, key):
        v = []
        for d in Y:
            if d['g'] != g: continue
            r = d['main'][p]
            if not r['ok']: continue
            x = r['P'][key] if key in r['P'] else r[key]
            if key in ('E21', 'E31'):
                if not r['cycle']: continue
                x = -x
            if x is not None and np.isfinite(x): v.append(x)
        return np.median(v)
    print('\nACROSS BACKGROUNDS (median difference vs own control; descriptive)')
    rows = []
    for t in res['main']['familyA'] + res['main']['sleep'] + res['main']['transitions'] + res['main']['landscape']:
        p = PH.index(t['phase']); e = yv('empty', p, t['key']); dd = yv('dfb', p, t['key']) - e; vv = yv('vnc', p, t['key']) - e
        rows.append(dict(key=t['key'], phase=t['phase'], d23E10=t['diff'], dfb=float(dd), vnc=float(vv)))
        print(f"  {t['key']:5s} {t['phase']:5s} 23E10 {t['diff']:+.4f} | dfb {dd:+.4f} | vnc {vv:+.4f} | sign 23E10=dfb {np.sign(t['diff']) == np.sign(dd)} =vnc {np.sign(t['diff']) == np.sign(vv)}")
    json.dump(dict(flies=F, results=res, across=rows), open(R.__file__.rsplit('/', 1)[0] + '/sqw_out/s24_23E10.json', 'w'))
