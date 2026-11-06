# PRESPEC 27 (58c6ec0): decomposition of the clock-bin cycle-flux rhythm (C10) by the state whose exits keep their
# daily modulation. Flatten / only-i counterfactuals on the pooled 2-h bin matrices; fly bootstrap (1000, seed 27).
import json, numpy as np, sqw, sqw_ref as R
from multiprocessing import Pool
NB = 12; I3 = np.eye(3); bb = np.arange(NB); XC = np.column_stack([np.ones(NB), np.cos(2 * np.pi * bb / NB), np.sin(2 * np.pi * bb / NB)])
G = {'CSDD': False, 'CS': True, 'per0LD': True, 'per0DD': False}; C10 = {'CS': 0.077, 'CSDD': 0.063}
NBOOT = 1000
def minute_N(y):
    y = np.asarray(y); t = np.arange(1, len(y)); ok = (y[:-1] >= 0) & (y[1:] >= 0); m = t % 1440
    N = np.zeros((1440, 3, 3)); np.add.at(N, (m[ok], y[:-1][ok], y[1:][ok]), 1); return N.reshape(NB, 120, 3, 3).sum(1)
def stat_pi(P):
    A = P.T - I3; A[2, :] = 1.0; return np.linalg.solve(A, np.array([0, 0, 1.0]))
def J_of(P):
    pi = stat_pi(P); return 60 * (pi[1] * P[1, 0] - pi[0] * P[0, 1]), pi[2]
def cosfit(v):
    (m, a, c), *_ = np.linalg.lstsq(XC, v, rcond=None)
    return float(np.hypot(a, c)), float((np.arctan2(c, a) % (2 * np.pi)) / (2 * np.pi) * 24 + 1), float(m)
def analyse(Nb):
    Pb = Nb / Nb.sum(2, keepdims=True); Pday = Nb.sum(0) / Nb.sum(0).sum(1, keepdims=True)
    def rhythm(keep):   # keep: rows that stay bin-specific
        J = []; S = []
        for b in range(NB):
            P = Pb[b].copy()
            for i in range(3):
                if i not in keep: P[i] = Pday[i]
            j, s = J_of(P); J.append(j); S.append(s)
        return np.array(J), np.array(S)
    out = {}
    Jf, Sf = rhythm({0, 1, 2}); A, ph, mu = cosfit(Jf); AS, phS, muS = cosfit(Sf)
    out['full'] = dict(A=A, peak=ph, mean=mu, J=Jf.tolist(), AS=AS, peakS=phS)
    for i, nm in enumerate('QWS'):
        Jx, Sx = rhythm({0, 1, 2} - {i}); Jo, So = rhythm({i})
        a1, p1, _ = cosfit(Jx); a2, p2, _ = cosfit(Jo); s1, _, _ = cosfit(Sx); s2, _, _ = cosfit(So)
        out[nm] = dict(R=1 - a1 / A, O=a2 / A, peak_flat=p1, peak_only=p2, RS=1 - s1 / AS, OS=s2 / AS, J_flat=Jx.tolist(), J_only=Jo.tolist())
    # secondary: rhythm of each off-diagonal exit probability, relative to its day value
    tp = {}
    for i in range(3):
        for j in range(3):
            if i != j:
                a, p, _ = cosfit(Pb[:, i, j]); tp['QWS'[i] + '->' + 'QWS'[j]] = dict(rel_amp=a / Pday[i, j], peak=p)
    out['tp'] = tp; return out
def boot(args):
    g, seed = args; Ns = np.load(f'/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit/sqw_out/s27_{g}_N.npy')
    rng = np.random.default_rng(seed); res = []
    for _ in range(NBOOT // 10):
        idx = rng.integers(0, len(Ns), len(Ns)); res.append(analyse(Ns[idx].sum(0)))
    return res
if __name__ == '__main__':
    out = {}
    for g, ld in G.items():
        Ns = np.array([minute_N(R.sqw_symbols(R.read_counts(f), theta=0.5, light_dark=ld)[0]) for f in sqw.files(g)])
        np.save(sqw.OUTD + f's27_{g}_N.npy', Ns); full = analyse(Ns.sum(0))
        if g in C10: print(f'GATE {g}: A {full["full"]["A"]:.4f} vs C10 {C10[g]}'); assert abs(round(full['full']['A'], 3) - C10[g]) <= 1e-3 + 1e-9, 'gate failed'
        with Pool(4) as pool: parts = pool.map(boot, [(g, 2700 + 10 * k + list(G).index(g)) for k in range(10)])
        B = [r for p in parts for r in p]
        def ci(f): v = np.array([f(r) for r in B]); return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
        summ = dict(n_flies=len(Ns), point=full, A_ci=ci(lambda r: r['full']['A']), n_boot_lowA=int(sum(r['full']['A'] < 0.01 for r in B)))
        for nm in 'QWS':
            summ[nm] = dict(R=full[nm]['R'], R_ci=ci(lambda r: r[nm]['R']), O=full[nm]['O'], O_ci=ci(lambda r: r[nm]['O']),
                            RS=full[nm]['RS'], RS_ci=ci(lambda r: r[nm]['RS']), OS=full[nm]['OS'], OS_ci=ci(lambda r: r[nm]['OS']))
        summ['tp'] = {k: dict(rel_amp=v['rel_amp'], peak=v['peak'], ci=ci(lambda r: r['tp'][k]['rel_amp'])) for k, v in full['tp'].items()}
        dom = [nm for nm in 'QWS' if summ[nm]['R_ci'][0] >= 0.5 and summ[nm]['O'] == max(summ[x]['O'] for x in 'QWS')]
        summ['verdict'] = f'exits from {dom[0]} carry the rhythm' if dom else 'shared across states'
        out[g] = summ
        print(f"{g}: A {full['full']['A']:.4f} [{summ['A_ci'][0]:.3f},{summ['A_ci'][1]:.3f}] peak {full['full']['peak']:.1f} h | low-A boots {summ['n_boot_lowA']} | " +
              ' | '.join(f"{nm}: R {summ[nm]['R']:+.2f} [{summ[nm]['R_ci'][0]:+.2f},{summ[nm]['R_ci'][1]:+.2f}] O {summ[nm]['O']:.2f} [{summ[nm]['O_ci'][0]:.2f},{summ[nm]['O_ci'][1]:.2f}]" for nm in 'QWS') + f" | {summ['verdict']}")
        print('   sleep occupancy: A_S %.4f peak %.1f | ' % (full['full']['AS'], full['full']['peakS']) + ' | '.join(f"{nm}: RS {summ[nm]['RS']:+.2f} [{summ[nm]['RS_ci'][0]:+.2f},{summ[nm]['RS_ci'][1]:+.2f}] OS {summ[nm]['OS']:.2f}" for nm in 'QWS'))
        print('   TP rel. amplitude: ' + ', '.join(f"{k} {v['rel_amp']:.2f} [{v['ci'][0]:.2f},{v['ci'][1]:.2f}] pk {v['peak']:.0f}" for k, v in summ['tp'].items()), flush=True)
    json.dump(out, open(sqw.OUTD + 's27_decomp.json', 'w'))
