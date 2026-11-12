"""Smoke and recovery checks for hmm_cat / sim_positions (code correctness, not the gate)."""
import sys, os, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, hmm_cat as hc, sim_positions as sp, pipeline as pl
rng = np.random.default_rng(0)
# 1) HMM recovers a known 3-state categorical model
P = np.array([[[.98, .01, .01], [.02, .96, .02], [.03, .03, .94]]] * 2)
B = np.zeros((3, 12)); B[0, [0, 1]] = [.9, .1]; B[1, [3, 6, 4]] = [.5, .3, .2]; B[2, [9, 10, 11]] = [.4, .3, .3]
true = (np.ones(3) / 3, P, B)
ph = ((np.arange(8640) * 10) % 86400 >= 43200).astype(np.int64)
seqs = hc.simulate(true, [(8640, ph)] * 4, [np.zeros(8640, bool)] * 4, rng)
t0 = time.time(); m, ll, lls, it = hc.fit(seqs, 3, 12, 2, np.random.SeedSequence(5), n_restarts=3); dt = time.time() - t0
lt = hc.loglik(true, seqs)
print(f'fit K=3: {dt:.1f}s, iters {it}, ll {ll:.1f} vs truth {lt:.1f}')
assert ll >= lt - 1.0
perm = [int(np.argmax([-hc.js_bits(m[2][j], B[k]) for j in range(3)])) for k in range(3)]
assert sorted(perm) == [0, 1, 2]
assert max(abs(m[1][0][perm][:, perm].diagonal() - P[0].diagonal())) < 0.02
# monotone EM check (manual few steps)
# 2) simulator: 2 days of G1, symbols and scoring run, state fractions sane
st, rd, t = sp.simulate_fly('G1', 2, np.random.SeedSequence([1, 2]))
assert rd.min() >= 0 and rd.max() <= 16 and not (rd == 15).any()
y = sp.window_symbols(pl.spike_filter(rd)[0]); assert y.min() >= 0 and y.max() < 12
frac = np.bincount(st, minlength=4) / len(st); print('state fractions (S,Q,T,B):', np.round(frac, 3))
print('symbol freq:', np.round(np.bincount(y, minlength=12) / len(y), 3))
# 3) timing of one K=8 fit on 7 days
st, rd, t = sp.simulate_fly('G1', 7, np.random.SeedSequence([1, 3])); y = sp.window_symbols(rd)
seqs = [(y[d*8640:(d+1)*8640], ph) for d in range(7)]
t0 = time.time(); m, ll, lls, it = hc.fit(seqs, 8, 12, 2, np.random.SeedSequence(7), n_restarts=1); print(f'K=8, 7 days, 1 restart: {time.time()-t0:.1f}s, {it} iters')
print('all checks passed')
