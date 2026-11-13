"""Code checks for d1d.py (symbols, flicker removal, compiled simulator, warm start)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, d1d, hmm_cat as hc
x = np.array([5, 5, 6, 5, 9, 9, 16, 14])
assert list(d1d.step_symbols(x)) == [0, 0, 1, 1, 3, 0, 3, 2]
v = np.ones(8, bool); v[4] = False
assert list(d1d.step_symbols(x, v)) == [0, 0, 1, 1, -1, -1, 3, 2]
assert list(d1d.remove_flicker(x)) == [5, 5, 5, 5, 9, 9, 16, 14]
# compiled simulator reproduces emission and self-transition frequencies
P = np.array([[[.99, .01], [.02, .98]]] * 2); B = np.array([[.97, .02, .01, 0], [.5, .3, .1, .1]])
ph = np.zeros(200000, np.int64); y = hc.simulate_fast((np.array([1., 0]), P, B), [(200000, ph)], [np.zeros(200000, bool)], np.random.default_rng(0))[0][0]
f = np.bincount(y, minlength=4) / len(y); exp = (2/3) * B[0] + (1/3) * B[1]
assert np.allclose(f, exp, atol=0.01), (f, exp)
# warm start from the truth is never worse than the truth
seqs = [(y, ph)]; m, ll, lls, it = hc.fit(seqs, 2, 4, 2, np.random.SeedSequence(1), n_restarts=1, max_iter=50, tol=0.01, inits=[(np.array([1., 0]), P, B)])
assert ll >= hc.loglik((np.array([1., 0]), P, B), seqs) - 1e-6
print('d1d checks passed')
