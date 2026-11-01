# PRESPEC 20z: (1) export sleep bouts of the Kir2.1 flies for Part A (MATLAB, lClaude); (2) Part C, exit route vs
# time asleep (20k R1 statistic), per fly x phase; C1 within-fly permutation null; C2 Mann-Whitney vs empty, Holm over 4.
import os, sys, glob, csv, json, numpy as np
from multiprocessing import Pool
from scipy.stats import mannwhitneyu
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sqw_ref as R
D = '/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/data_kir'
G = ['empty', 'dfb', 'vnc']; MON = {715157: 'rc2', 715166: 'rc', 715185: 'wc'}; PH = ['light', 'dark']; DAY = 86400
def one(a):
    g, f = a; x = np.loadtxt(f).astype(int); mon = MON[len(x)]; xd = x[:(len(x) // DAY) * DAY]
    b = [(s, n) for s, n in R.sleep_bouts(xd, True) if s > 0 and s + n < len(xd)]
    rec = R.route_duration_observed(xd, 0.5, light_dark=True)
    return dict(g=g, fly=os.path.basename(f).replace('.txtcutdata_4.8.2025', '').replace('.txt', ''), mon=mon, bouts=b, rec=rec)
def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); adj = np.empty_like(ps); run = 0
    for r, i in enumerate(o): run = max(run, (len(ps) - r) * ps[i]); adj[i] = min(1, run)
    return adj
if __name__ == '__main__':
    F = Pool(4).map(one, [(g, f) for g in G for f in sorted(glob.glob(f'{D}/{g}/*'))])
    ids = [f"{d['g']}:{d['fly']}" for d in F]; assert len(set(ids)) == len(ids), 'fly ids not unique'
    with open('kir_sleep_bouts_20z.csv', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['group', 'fly', 'monitor', 'start_s', 'duration_s', 'start_ZT_h'])
        for d in F:
            for s, n in d['bouts']: w.writerow([d['g'], d['fly'], d['mon'], s, n, round((s % DAY) / 3600, 4)])
    print('sleep bouts exported:', {g: sum(len(d['bouts']) for d in F if d['g'] == g) for g in G},
          '| min duration', min(n for d in F for _, n in d['bouts']))
    rng = np.random.default_rng(20); Dv = {}
    print('\n== C1 group mean D (R1 exit vs sleep length) vs within-fly permutation null (2000)')
    for g in G:
        for p in range(2):
            fl = []
            for d in F:
                if d['g'] != g: continue
                rr = [r for r in d['rec'] if r['phase_on'] == p and r['q_out'] is not None]
                v = R.route_duration_D(rr, 'R1_exit_vs_sleep', p)
                if np.isfinite(v): fl.append((rr, v))
            Dv[(g, p)] = [v for _, v in fl]; obs = np.mean(Dv[(g, p)]); null = np.zeros(2000)
            for b in range(2000):
                vals = []
                for rr, _ in fl:
                    lab = rng.permutation([r['q_out'] for r in rr])
                    vals.append(R.route_duration_D([dict(r, q_out=bool(l)) for r, l in zip(rr, lab)], 'R1_exit_vs_sleep', p))
                null[b] = np.nanmean(vals)
            pv = (1 + np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean()))) / 2001
            print(f'{g:<5} {PH[p]:<5} flies {len(fl):>2} | mean D {obs:+.3f} | null mean {null.mean():+.3f} | p {pv:.4f}')
    print('\n== C2 per-fly D vs empty, two-sided Mann-Whitney, Holm over 4')
    rows = []
    for p in range(2):
        for g in ('dfb', 'vnc'):
            a, b = Dv[(g, p)], Dv[('empty', p)]; u = mannwhitneyu(a, b, alternative='two-sided')
            rows.append((g, PH[p], len(a), len(b), np.median(a), np.median(b), u.pvalue, u.statistic / (len(a) * len(b))))
    for r, h in zip(rows, holm([r[6] for r in rows])):
        print(f'{r[0]} vs empty {r[1]:<5} n {r[2]}/{r[3]} | median D {r[4]:+.3f} vs {r[5]:+.3f} | PS {r[7]:.2f} | p {r[6]:.4g} Holm {h:.4g}'
              + ('  -> CHANGES' if h < 0.05 else ''))
    json.dump({f'{g}|{PH[p]}': v for (g, p), v in Dv.items()}, open('z_partC.json', 'w'))
