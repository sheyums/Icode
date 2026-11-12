"""Code checks for iohmm.py (recovery of a known clock-modulated 2-state model; M0 matches hmm_cat)."""
import sys, os, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, iohmm as io, hmm_cat as hc
rng = np.random.default_rng(1)
Z, pen = io.covariates(1)
beta = np.zeros((2, 2, Z.shape[1])); beta[0, 1] = [-4.0, -0.7, 1.0, 0.0]; beta[1, 0] = [-3.0, 0.5, 0.0, 0.0]
B = np.array([[.9, .05, .03, .02], [.2, .2, .3, .3]])
P = io.trans_matrices(beta, Z); true = (np.array([.5, .5]), P, B)
w = np.arange(8640, dtype=np.int64)
seqs = hc.simulate(true, [(8640, w)] * 14, [np.zeros(8640, bool)] * 14, rng)
t0 = time.time()
m, b, ll, lls = io.fit(seqs, 2, 4, 1, np.random.SeedSequence(3), n_restarts=2, lam=0.0)
k = int(np.argmax(m[2][:, 0]))      # the 'still' state
print(f'{time.time()-t0:.1f}s; fitted beta from still state: {np.round(b[k, 1-k], 2)}  truth {beta[0, 1]}')
assert np.allclose(b[k, 1 - k], beta[0, 1], atol=0.35), b[k, 1 - k]
assert ll >= hc.loglik(true, seqs) - 1.0
# M0 IOHMM equals the light/dark step HMM: same maximum likelihood (within 1 nat)
ph = ((w * 10) % 86400 >= 43200).astype(np.int64)
s2 = [(y, ph) for y, _ in seqs]
_, l_step, _, _ = hc.fit(s2, 2, 4, 2, np.random.SeedSequence(4), n_restarts=2)
_, _, l_m0, _ = io.fit(seqs, 2, 4, 0, np.random.SeedSequence(4), n_restarts=2, lam=0.0)
print(f'step HMM ll {l_step:.1f}  IOHMM M0 ll {l_m0:.1f}'); assert abs(l_step - l_m0) < 1.0
print('iohmm checks passed')
