# PRESPEC 23 (8571978): Sh^mns direction, C14 within genotype, Onsager L vs CS, landscape sign.
import glob, json, numpy as np, sqw, sqw_ref as R
from multiprocessing import Pool
from scipy.stats import wilcoxon, rankdata
D = '/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/data_shmns'; PH = ['light', 'dark']
def trim(x):
    run = 0; dead_at = None
    for i, a in enumerate(x > 0):
        run = 0 if a else run + 1
        if run == 43200: dead_at = i - 43200 + 1; break
    if dead_at is not None: x = x[:max(0, dead_at - 86400)]
    nd = min(len(x) // 86400, 4); return x[:nd * 86400], nd, dead_at is not None
def one(f):
    x, nd, died = trim(np.loadtxt(f).astype(int)); out = dict(f=f.rsplit('/', 1)[1], days=nd, died=died, excluded=nd < 2)
    if nd < 2: return out
    y, refs = R.sqw_symbols(x, 0.5, light_dark=True); ph = R.phases(len(y), '12h'); _, _, P = R.mc_fit(y, ph, 2)
    oq = R.observable_quantities(y, ph, 2, P); res = []
    for p in range(2):
        m = ph == p; d = dict(oq[p], ok=not np.isnan(refs[PH[p]]), sleep=float(60 * np.mean(y[m] == 2)))
        d['cycle'] = abs(d['J_chain']) >= 1e-12; res.append(d)
    out['main'] = res; return out
rng = np.random.default_rng(23)
def rk(v): return rankdata(np.round(np.asarray(v, float), 10))
def sp_perm(x, y, n=10000):
    rx, ry = rk(x), rk(y); r0 = np.corrcoef(rx, ry)[0, 1]; nl = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(n)])
    return float(r0), float((1 + np.sum(np.abs(nl) >= abs(r0) - 1e-12)) / (n + 1))
def holm(ps):
    ps = np.asarray(ps); o = np.argsort(ps); a = np.empty(len(ps)); run = 0
    for j, i in enumerate(o): run = max(run, min(1, (len(ps) - j) * ps[i])); a[i] = run
    return a
def Lfit(pts): J = np.array([q[1] for q in pts]); E = np.array([q[2] for q in pts]); return (J ** 4).sum() / (J ** 2 * E).sum()
def Lboot(pts, B=2000):
    flies = sorted({q[0] for q in pts}); by = {f: [q for q in pts if q[0] == f] for f in flies}; out = []
    for _ in range(B):
        qq = [q for f in rng.choice(flies, len(flies)) for q in by[f]]; out.append(Lfit(qq))
    return np.array(out)
if __name__ == '__main__':
    with Pool(4) as pool: F = pool.map(one, sorted(glob.glob(D + '/**/*.txt', recursive=True)))
    U = [d for d in F if not d['excluded']]
    print('flies', len(F), '| died', sum(d['died'] for d in F), '| excluded', len(F) - len(U), '| days', np.bincount([d['days'] for d in F]).tolist())
    out = dict(n=len(U)); ps = []
    for p in range(2):
        R_ = [d['main'][p] for d in U if d['main'][p]['ok']]; J = np.array([r['J_chain'] for r in R_ if r['cycle']])
        w = wilcoxon(np.round(J, 10)); out[PH[p]] = dict(n=len(R_), n_nocycle=sum(not r['cycle'] for r in R_), med_J=float(np.median(J)), frac_neg=float(np.mean(J < 0)),
                p=float(w.pvalue), mean_J_all=float(np.mean([r['J_chain'] if r['cycle'] else 0 for r in R_])), sleep_med=float(np.median([r['sleep'] for r in R_])),
                c_obs_med=float(np.nanmedian([r['c_obs'] for r in R_])), E31_med=float(np.median([r['E31'] for r in R_ if r['cycle']])))
        ps.append(w.pvalue)
    for p, a in zip(range(2), holm(ps)): out[PH[p]]['holm'] = float(a)
    sec = []
    for p in range(2):
        R_ = [d['main'][p] for d in U if d['main'][p]['ok']]
        a = [(r['sleep'], r['c_obs']) for r in R_ if r['c_obs'] is not None and np.isfinite(r['c_obs'])]; b = [(r['sleep'], r['J_chain']) for r in R_ if r['cycle']]
        for nm, v in (('c_obs', a), ('J', b)): r0, pp = sp_perm([t[0] for t in v], [t[1] for t in v]); sec.append(dict(phase=PH[p], y=nm, rho=r0, p=pp, n=len(v)))
    for s, a in zip(sec, holm([s['p'] for s in sec])): s['holm'] = float(a)
    out['C14'] = sec
    pts = [(d['f'], r['J'], r['EPR']) for d in U for r in d['main'] if r['ok'] and r['cycle'] and np.isfinite(r['EPR'])]
    O = json.load(open(sqw.OUTD + 's_observable.json')); cs = [(d['fly'], r['J'], r['EPR']) for d in O if d['g'] == 'CS' for r in d['main'] if r['ok'] and r['J_chain'] is not None and abs(r['J_chain']) >= 1e-12 and np.isfinite(r['EPR'])]
    Ls, Lc = Lfit(pts), Lfit(cs); bs, bc = Lboot(pts), Lboot(cs); dlt = bs - bc
    out['L'] = dict(shmns=float(Ls), shmns_ci=np.percentile(bs, [2.5, 97.5]).tolist(), CS=float(Lc), CS_ci=np.percentile(bc, [2.5, 97.5]).tolist(),
                    diff=float(Ls - Lc), p=float(min(1, 2 * min(np.mean(dlt <= 0), np.mean(dlt >= 0)))))
    for p in range(2):
        e = np.array([d['main'][p]['E31'] for d in U if d['main'][p]['ok'] and d['main'][p]['cycle']]); b = [np.median(rng.choice(e, len(e))) for _ in range(2000)]
        out[PH[p]]['E31_ci'] = np.percentile(b, [2.5, 97.5]).tolist()
    json.dump(dict(summary=out, flies=F), open(sqw.OUTD + 's23_shmns.json', 'w'), default=float)
    for p in PH:
        v = out[p]; print(f"{p}: n {v['n']} (no cycle {v['n_nocycle']}) | sleep {v['sleep_med']:.1f} min/h | J median {v['med_J']:+.3f}, mean(all) {v['mean_J_all']:+.3f}, J<0 in {v['frac_neg']:.0%} | Wilcoxon p {v['p']:.4g} Holm {v['holm']:.4g} | c_obs med {v['c_obs_med']:+.3f} | E_S-E_Q med {v['E31_med']:+.2f} [{v['E31_ci'][0]:+.2f},{v['E31_ci'][1]:+.2f}]")
    for s in sec: print(f"C14 {s['phase']} sleep vs {s['y']}: rho {s['rho']:+.2f} p {s['p']:.3g} Holm {s['holm']:.3g} (n {s['n']})")
    L = out['L']; print(f"Onsager L: Sh^mns {L['shmns']:.3f} [{L['shmns_ci'][0]:.3f},{L['shmns_ci'][1]:.3f}] vs CS {L['CS']:.3f} [{L['CS_ci'][0]:.3f},{L['CS_ci'][1]:.3f}] diff {L['diff']:+.3f} p {L['p']:.3f}")
