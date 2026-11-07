# PRESPEC 29 (ee691ad): test-retest (days A = 0,2 vs B = 1,3) of per-fly J, observed and hidden; excess spread of
# per-fly J against a genotype-pooled chain (observed). Sh^mns HMMs fitted here with the 20b procedure.
import glob, json, numpy as np, sqw, sqw2, sqw_ref as R
from multiprocessing import Pool
from scipy.stats import rankdata
from sqw15 import seqs
from sqw2_common import em, label3b, degenerate
from sqw23 import trim, D as SHD
G = ['CS', 'CSDD', 'per0LD', 'per0DD', 'fumin', 'sss', 'shmns']; PH = ['light', 'dark']
def nph(g): return 1 if g in ('CSDD', 'per0DD') else 2
SHF = sorted(glob.glob(SHD + '/**/*.txt', recursive=True))
def shm_seq(i):
    x, nd, _ = trim(np.loadtxt(SHF[i]).astype(int)); y, _ = R.sqw_symbols(x, 0.5, light_dark=True); return np.asarray(y), R.phases(len(y), '12h')
def day_J(z, ph, npp):
    """(days x phases) net S->W minus W->S per 12 h; z may be (paths x T)."""
    z = np.atleast_2d(z); T = z.shape[1]; nd = T // 1440; t = np.arange(1, T); ok = (z[:, :-1] >= 0) & (z[:, 1:] >= 0)
    net = ((z[:, :-1] == 2) & (z[:, 1:] == 1) & ok).astype(int) - ((z[:, :-1] == 1) & (z[:, 1:] == 2) & ok).astype(int)
    d = t // 1440; out = np.zeros((z.shape[0], nd, npp)); scored = np.zeros((nd, npp))
    for k in range(nd):
        for p in range(npp):
            m = (d == k) & (ph[t] == p); out[:, k, p] = net[:, m].sum(1) / 12.0
            scored[k, p] = np.mean(ok[0, m]) if m.any() else 0
    return out.mean(0), scored
def obs_one(a):
    g, i = a
    y, ph = shm_seq(i) if g == 'shmns' else (lambda s: (np.asarray(s[0]), np.asarray(s[1])))(seqs(g, i))
    Jd, sc = day_J(y, ph, nph(g)); return dict(g=g, fly=i, Jd=Jd.tolist(), scored=sc.tolist(), T=len(y))
def hid_one(a):
    g, i = a; rng_seed = 2900 + 100 * G.index(g) + i
    if g == 'shmns':
        y, ph = shm_seq(i); ll, pi, P, B, it = em(y, ph, 2, 3, sqw2.ST, 20 + 1000 * G.index(g) + i)
        P, B, pi = label3b(np.array(P), np.array(B), np.array(pi)); deg = bool(degenerate(B))
    else:
        F = HF.get((g, i))
        if F is None or F['deg']: return dict(g=g, fly=i, skip=True)
        y, ph, _ = sqw2.seqs(g, i, 0.5); y = np.asarray(y); ph = np.asarray(ph); P, B, pi = np.array(F['P3']), np.array(F['B3']), np.array(F['pi3']); deg = False
    if deg: return dict(g=g, fly=i, skip=True, deg=True)
    Z = R.ffbs(pi, P, B, y, ph, n_paths=200, seed=rng_seed); Jd, _ = day_J(Z, ph, nph(g))
    return dict(g=g, fly=i, Jd=Jd.tolist(), skip=False)
HF = sqw2.load(0.5)
rng = np.random.default_rng(29)
def rk(v): return rankdata(np.round(np.asarray(v, float), 10))
def halves(Jd, sc, p):
    Jd = np.asarray(Jd); nd = Jd.shape[0]; A = [k for k in (0, 2) if k < nd and (sc is None or sc[k][p] > 0.5)]; B = [k for k in (1, 3) if k < nd and (sc is None or sc[k][p] > 0.5)]
    return (Jd[A, p].mean(), Jd[B, p].mean()) if A and B else None
def retest(recs, p, scored=True):
    pr = [h for r in recs for h in [halves(r['Jd'], r.get('scored') if scored else None, p)] if h is not None]
    a = np.array([x for x, _ in pr]); b = np.array([y for _, y in pr]); ra, rb = rk(a), rk(b); r0 = np.corrcoef(ra, rb)[0, 1]
    nl = np.array([np.corrcoef(ra, rng.permutation(rb))[0, 1] for _ in range(10000)])
    return dict(n=len(pr), rho=float(r0), p=float((1 + np.sum(nl >= r0 - 1e-12)) / 10001))
def spread(g, recs, p, NS=500):
    Ys = []
    for r in recs:
        y, ph = shm_seq(r['fly']) if g == 'shmns' else (lambda s: (np.asarray(s[0]), np.asarray(s[1])))(seqs(g, r['fly'])); Ys.append((y, ph))
    C = np.zeros((nph(g), 3, 3))
    for y, ph in Ys:
        ok = (y[:-1] >= 0) & (y[1:] >= 0); np.add.at(C, (ph[1:][ok], y[:-1][ok], y[1:][ok]), 1)
    P = C / C.sum(2, keepdims=True); Cc = np.cumsum(P, 2); Cc[..., -1] = 1
    def Jwhole(z, ph):
        ok = (z[:, :-1] >= 0) & (z[:, 1:] >= 0); m = (ph[1:] == p)
        net = (((z[:, :-1] == 2) & (z[:, 1:] == 1) & ok).astype(int) - ((z[:, :-1] == 1) & (z[:, 1:] == 2) & ok).astype(int))[:, m].sum(1)
        h = (ok[:, m].sum(1)) / 60.0; return net / np.maximum(h, 1e-9)
    obs = np.array([Jwhole(y[None], ph)[0] for y, ph in Ys]); sims = []
    for y, ph in Ys:
        T = len(y); Z = np.empty((NS, T), np.int8); cur = np.full(NS, int(y[np.argmax(y >= 0)]))
        Z[:, 0] = cur
        for t in range(1, T): cur = (rng.random(NS)[:, None] > Cc[ph[t], cur]).sum(1); Z[:, t] = cur
        Z[:, y < 0] = -1; sims.append(Jwhole(Z, ph))
    sd_null = np.std(np.array(sims), axis=0); sd_obs = float(np.std(obs))
    return dict(sd_obs=sd_obs, sd_null_med=float(np.median(sd_null)), p=float((1 + np.sum(sd_null >= sd_obs)) / (NS + 1)), J=obs.tolist())
def holm(ps):
    ps = np.asarray(ps); o = np.argsort(ps); a = np.empty(len(ps)); run = 0
    for j, i in enumerate(o): run = max(run, min(1, (len(ps) - j) * ps[i])); a[i] = run
    return a
if __name__ == '__main__':
    jobs = [(g, i) for g in G for i in range(len(SHF) if g == 'shmns' else len(sqw.files(g)))]
    with Pool(4) as pool: OB = pool.map(obs_one, jobs)
    out = {'obs': {}, 'hid': {}}
    for g in G:
        recs = [r for r in OB if r['g'] == g]
        for p in range(nph(g)):
            k = f'{g}|{PH[p] if nph(g) == 2 else "whole"}'; out['obs'][k] = dict(T1=retest(recs, p), T2=spread(g, recs, p)); print('obs', k, out['obs'][k]['T1'], {x: out['obs'][k]['T2'][x] for x in ('sd_obs', 'sd_null_med', 'p')}, flush=True)
    prim = ['sss|light', 'shmns|light']; ps = [out['obs'][k][t]['p'] for k in prim for t in ('T1', 'T2')]
    for (k, t), a in zip([(k, t) for k in prim for t in ('T1', 'T2')], holm(ps)): out['obs'][k][t]['holm'] = float(a)
    hjobs = [j for j in jobs if j[0] != 'CSDD']
    with Pool(4) as pool: HB = pool.map(hid_one, hjobs)
    for g in G:
        if g == 'CSDD': continue
        recs = [r for r in HB if r['g'] == g and not r['skip']]
        for p in range(nph(g)):
            k = f'{g}|{PH[p] if nph(g) == 2 else "whole"}'; out['hid'][k] = dict(T1=retest(recs, p, scored=False), n_skip=sum(1 for r in HB if r['g'] == g and r['skip'])); print('hid', k, out['hid'][k], flush=True)
    for k, a in zip(prim, holm([out['hid'][k]['T1']['p'] for k in prim])): out['hid'][k]['T1']['holm'] = float(a)
    json.dump(dict(out=out, obs_flies=OB, hid_flies=HB), open(sqw.OUTD + 's29_hetero.json', 'w'))
    for lev in ('obs', 'hid'):
        for k in prim: print('PRIMARY', lev, k, out[lev][k])
