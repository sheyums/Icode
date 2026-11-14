"""Code checks for b1.py: symbols, AR emission recovery, compiled simulator, degeneracy rule."""
import sys, os, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, b1, hmm_cat as hc
# symbols: 5 5 6 5 5 7 8 9 9 8  -> STAY STAY SAME REV STAY BIG SAME SAME STAY REV
x = np.array([5, 5, 6, 5, 5, 7, 8, 9, 9, 8])
assert list(b1.symbols(x)) == [0, 0, 1, 2, 0, 3, 1, 1, 0, 2], b1.symbols(x)
v = np.ones(10, bool); v[4] = False
assert list(b1.symbols(x, v)) == [0, 0, 1, 2, -1, -1, 1, 1, 0, 2], b1.symbols(x, v)   # direction resets after a gap
assert list(b1.composite(np.array([0, 1, 2, -1, 1]))) == [0, 1, 6, -1, 1]
# a 2-state AR truth: still (rare flicker pairs: step then REV) and persistent travel
B = np.zeros((2, 4, 4))
B[0] = [[.99, .005, .005, 0], [.01, 0, .99, 0], [.01, 0, .99, 0], [1, 0, 0, 0]]          # after a step: reverse at once
B[1] = [[.5, .3, .1, .1], [.4, .45, .05, .1], [.4, .45, .05, .1], [.5, .3, .1, .1]]       # persistent
B = B.reshape(2, 16); P = np.array([[[.999, .001], [.004, .996]]] * 2)
truth = (np.array([.5, .5]), P, B)
rng = np.random.default_rng(0); ph = np.zeros(100000, np.int64)
seqs = b1.simulate(truth, [(100000, ph)] * 4, [np.zeros(100000, bool)] * 4, rng)
t0 = time.time(); m, ll, lls = b1.fit(seqs, 2, np.random.SeedSequence(1), n_restarts=2)
print(f'fit {time.time()-t0:.0f}s, ll {ll:.1f} vs truth {hc.loglik(truth, seqs):.1f}')
assert ll >= hc.loglik(truth, seqs) - 1.0
k = int(np.argmax(m[2][:, 0]))                                   # the still state
assert abs(m[2][k].reshape(4, 4)[1, 2] - 0.99) < 0.03, m[2][k].reshape(4, 4)
tp, _, rho = b1.type_of(m[2][k], np.ones(1), np.array([0])); assert tp[2] == 'FLICKER-LIKE', (tp, rho)
tp, _, rho = b1.type_of(m[2][1 - k], np.ones(1), np.array([0])); assert tp[2] == 'PERSISTENT', (tp, rho)
# degeneracy: a state duplicated with the same dwell is degenerate; sleep vs quiet (close emissions,
# different dwells) is not
Bs = np.zeros((2, 4, 4)); Bs[0, :, 0] = 1; Bs[1] = [[.97, .015, .015, 0]] * 4
Pd = np.array([[[.9998, .0002], [.002, .998]]] * 2)
assert not b1.degenerate((np.ones(2) / 2, Pd, Bs.reshape(2, 16)))
assert b1.degenerate((np.ones(2) / 2, np.array([[[.998, .002], [.002, .998]]] * 2), np.stack([Bs[1], Bs[1]]).reshape(2, 16)))
print('b1 checks passed')
