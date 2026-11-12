"""P1-B item 4: R1 synthetic check of the hub repair. Usage: python3 -I r1_check.py"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pipeline as pl, sim_positions as sp

def add_hub(read, rng):
    n = len(read); hub = int(rng.integers(2, 15)); out = read.copy(); t = 0
    p_on = 0.4 / 60.0
    starts = np.nonzero(rng.random(n) < p_on)[0]
    for s in starts:
        L = max(1, int(round(rng.exponential(1.5) / 0.4)))
        out[s:s + L] = hub
    return out, hub

def score(x, t):
    v = np.ones(len(x), bool); a = pl.activity_1s(x, t, v); y, _, _ = pl.sqw_minutes(a)
    return y, sp.window_symbols(x)

ok = True; rows = []
for i in range(10):
    st, rd, t = sp.simulate_fly('G1', 7, np.random.SeedSequence([4001, i]))
    rng = np.random.default_rng(np.random.SeedSequence([4001, i, 1]))
    art, hub = add_hub(rd, rng)
    clean_x, _ = pl.spike_filter(rd); clean_rep, info_clean = pl.hub_repair(clean_x)
    y0, s0 = score(clean_x, t)
    ax, _ = pl.spike_filter(art); yb, sb = score(ax, t)
    rx, info = pl.hub_repair(ax); yr, sr = score(rx, t)
    m = (y0 >= 0) & (yr >= 0)
    a_min, a_win = (yr[m] == y0[m]).mean(), (sr == s0).mean()
    b_min, b_win = (yb[m] == y0[m]).mean(), (sb == s0).mean()
    chg_clean = (clean_rep != clean_x).mean()
    rows.append((i, hub, info['hub'], round(info['share'], 2), round(info['returns'], 2), b_min, a_min, b_win, a_win, chg_clean))
    ok &= a_min >= 0.95 and a_win >= 0.95 and chg_clean <= 0.001
    print(f'fly {i}: true hub {hub} found {info["hub"]} share {info["share"]:.2f} ret {info["returns"]:.2f} | '
          f'S/Q/W agree before {b_min:.3f} after {a_min:.3f} | symbols before {b_win:.3f} after {a_win:.3f} | '
          f'clean fly changed {chg_clean:.5f} (hub found {info_clean["hub"]}, share {info_clean["share"]:.2f}, ret {info_clean["returns"]:.2f})')
print('R1', 'PASS' if ok else 'FAIL')
