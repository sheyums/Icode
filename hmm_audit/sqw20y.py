# PRESPEC 20y: constitutive Kir2.1 silencing (empty / dfb / vnc split-GAL4), observable S/Q/W cycle. One run.
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from scipy.stats import mannwhitneyu
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R
from sqw16b import kind
D = '/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/data_kir'
G = ['empty', 'dfb', 'vnc']; MON = {715157: 'rc2', 715166: 'rc', 715185: 'wc'}
PH = ['light', 'dark']; TR = [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]; SYM = 'QWS'
def one(a):
    g, f = a; x = np.loadtxt(f).astype(int); mon = MON[len(x)]
    if not (x[-86400:] > 0).any(): return dict(g=g, f=os.path.basename(f), mon=mon, dead=True)
    y, refs = R.sqw_symbols(x, 0.5, light_dark=True); ph = R.phases(len(y), '12h')
    out = dict(g=g, f=os.path.basename(f), mon=mon, dead=False, refs=refs)
    for tag, s in (('main', 0), ('nod0', 1440)):
        yy, pp = y[s:], ph[s:]; _, _, P = R.mc_fit(yy, pp, 2); oq = R.observable_quantities(yy, pp, 2, P); res = []
        for p in range(2):
            ok = not np.isnan(refs[PH[p]]); m = pp == p; h = m.sum() / 60
            d = dict(oq[p], ok=ok, sleep=float((yy[m] == 2).sum() / h), P={f'{SYM[i]}{SYM[j]}': float(P[p][i, j]) for i, j in TR})
            d['cycle'] = abs(d['J_chain']) >= 1e-12
            d['kind'] = None if d['cycle'] else kind(__import__('sqw16b').counts(yy, pp, p), yy, pp, p)
            res.append(d)
        out[tag] = res
    return out
def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); adj = np.empty_like(ps); run = 0
    for r, i in enumerate(o): run = max(run, (len(ps) - r) * ps[i]); adj[i] = min(1, run)
    return adj
rng = np.random.default_rng(20)
def bmed(v):
    v = np.asarray(v, float); return np.median(v)
def diffci(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    bs = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(2000)]
    return np.median(a) - np.median(b), *np.percentile(bs, [2.5, 97.5])
def vals(F, g, tag, p, key, cyc=False, mon=None):
    out = []
    for d in F:
        if d['g'] != g or d['dead'] or (mon and d['mon'] != mon): continue
        r = d[tag][p]
        if key != 'sleep' and not r['ok']: continue
        if cyc and not r['cycle']: continue
        out.append(r['P'][key] if key in r['P'] else r[key])
    return out
def mw(a, b):
    u = mannwhitneyu(a, b, alternative='two-sided'); return u.pvalue, u.statistic / (len(a) * len(b))
def family(F, tag, keys, cyc_keys, label):
    rows = []
    for k in keys:
        for p in range(2):
            for g in ('dfb', 'vnc'):
                a = vals(F, g, tag, p, k, k in cyc_keys); b = vals(F, 'empty', tag, p, k, k in cyc_keys)
                pv, ps = mw(a, b); d, lo, hi = diffci(a, b)
                rows.append(dict(key=k, ph=PH[p], g=g, n=(len(a), len(b)), med=(bmed(a), bmed(b)), diff=d, ci=(lo, hi), p=pv, PS=ps))
    for r, a in zip(rows, holm([r['p'] for r in rows])): r['holm'] = a
    print(f'\n== {label} ({tag}; Holm over {len(rows)})')
    for r in rows:
        print(f"{r['key']:>6} {r['ph']:<5} {r['g']} vs empty  n {r['n'][0]}/{r['n'][1]}  med {r['med'][0]:+.4f} vs {r['med'][1]:+.4f}"
              f"  diff {r['diff']:+.4f} [{r['ci'][0]:+.4f}, {r['ci'][1]:+.4f}]  PS {r['PS']:.2f}  p {r['p']:.4g}  Holm {r['holm']:.4g}"
              + ('  -> CHANGES' if r['holm'] < 0.05 else ''))
    return rows
def onsager(F, g, B=2000):
    pts = [(d['f'], r['J'], r['EPR']) for d in F if d['g'] == g and not d['dead'] for r in d['main'] if r['ok'] and r['cycle'] and np.isfinite(r['EPR'])]
    L = lambda q: (np.array([t[1] for t in q]) ** 4).sum() / (np.array([t[1] ** 2 * t[2] for t in q])).sum()
    flies = sorted({t[0] for t in pts}); bs = []
    for _ in range(B):
        pick = rng.choice(flies, len(flies)); qq = [t for f in pick for t in pts if t[0] == f]; bs.append(L(qq))
    return L(pts), np.array(bs)
if __name__ == '__main__':
    jobs = [(g, f) for g in G for f in sorted(glob.glob(f'{D}/{g}/*'))]
    with Pool(8) as pool: F = pool.map(one, jobs)
    json.dump(F, open('y_kir.json', 'w'), default=float)
    print('flies:', {g: sum(d['g'] == g for d in F) for g in G}, 'dead:', [(d['g'], d['f']) for d in F if d['dead']])
    print('no Q/W reference:', {f'{g} {PH[p]}': sum(1 for d in F if d['g'] == g and not d['dead'] and not d['main'][p]['ok']) for g in G for p in range(2)})
    K = {}
    for d in F:
        if d['dead']: continue
        for p in range(2):
            r = d['main'][p]
            if r['ok'] and not r['cycle']: K[f"{d['g']} {PH[p]} {r['kind']}"] = K.get(f"{d['g']} {PH[p]} {r['kind']}", 0) + 1
    print('no-cycle fly-phases by kind (i one-way, ii symmetric, iii never observed / state absent):', K)
    family(F, 'main', ['sleep'], [], 'POSITIVE CONTROL: sleep (min/h)')
    family(F, 'main', ['c_obs', 'J', 'EPR'], ['J', 'EPR'], 'PRIMARY')
    family(F, 'main', [f'{SYM[i]}{SYM[j]}' for i, j in TR], [], 'S1 transitions')
    family(F, 'main', ['E21', 'E31'], ['E21', 'E31'], 'S2 landscape')
    print('\n== S3 Onsager L (light+dark pooled; fly bootstrap)')
    Ls = {g: onsager(F, g) for g in G}
    for g in G: print(f'{g}: L {Ls[g][0]:.3f} [{np.percentile(Ls[g][1], 2.5):.3f}, {np.percentile(Ls[g][1], 97.5):.3f}]')
    for g in ('dfb', 'vnc'):
        dd = Ls[g][1] - Ls['empty'][1]; p = 2 * min((dd <= 0).mean(), (dd >= 0).mean())
        print(f'{g} - empty: {Ls[g][0] - Ls["empty"][0]:+.3f}, bootstrap p {max(p, 1 / 2000):.4f}')
    print('\n== Monitor check (descriptive, uncorrected)')
    for g, m1, m2 in (('vnc', 'rc2', 'wc'), ('empty', 'rc2', 'rc')):
        for k in ('sleep', 'c_obs', 'J', 'EPR'):
            for p in range(2):
                a = vals(F, g, 'main', p, k, k in ('J', 'EPR'), m1); b = vals(F, g, 'main', p, k, k in ('J', 'EPR'), m2)
                pv = mw(a, b)[0] if len(a) and len(b) else np.nan
                print(f'{g} {k:>6} {PH[p]:<5} {m1} {bmed(a):+.4f} (n {len(a)}) vs {m2} {bmed(b):+.4f} (n {len(b)})  diff {bmed(a) - bmed(b):+.4f}  p {pv:.3g}')
    family(F, 'nod0', ['c_obs', 'J', 'EPR'], ['J', 'EPR'], 'SENSITIVITY: primary without day 0')
