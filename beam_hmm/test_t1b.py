"""Code checks for t1b.py: C0 equivalence and recovery of a shared exit harmonic (direct HMM simulation)."""
import sys, os, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, hmm_cat as hc, t1b
rng = np.random.default_rng(0)
K, M = 3, 4
P0 = np.array([[[.9995, .0003, .0002], [.004, .99, .006], [.002, .01, .988]],
               [[.9998, .0001, .0001], [.006, .988, .006], [.004, .01, .986]]])
B = np.array([[1, 0, 0, 0], [.95, .04, .01, 0], [.5, .35, .1, .05]])
model = (np.ones(K) / K, P0, B)
F = t1b.from_c0(model, 0); Z = t1b.harmonics(2)
P = t1b.build_P(F, np.zeros(4), np.zeros(4), Z)
assert np.allclose(P[0], P0[0]) and np.allclose(P[-1], P0[1]), 'C0 equivalence'
h_true = np.array([0.6, 0.0, 0.4, 0.0])
flies = []
for f in range(12):
    Ff = t1b.from_c0(model, 0); Ff['a'] = Ff['a'] + rng.normal(0, 0.2, 2)
    Pt = t1b.build_P(Ff, h_true, np.zeros(4), Z)
    w = t1b._bin_of_sample.copy()
    seqs = hc.simulate_fast((np.ones(K) / K, Pt, B), [(len(w), w)] * 14, [np.zeros(len(w), bool)] * 14, rng)
    flies.append(dict(F=t1b.from_c0(model, 0), seqs=seqs))
t0 = time.time()
fit2 = t1b.fit_cohort(flies, 2, [(np.zeros(4), np.zeros(4))])
fit0 = t1b.fit_cohort(flies, 0, [(np.zeros(0), np.zeros(0))])
print(f'{time.time()-t0:.0f}s; iters C2 {fit2["iters"]}; h fitted {np.round(fit2["h"], 2)} truth {h_true}; g {np.round(fit2["g"], 2)}')
print(f'LR = {2*(fit2["ll"]-fit0["ll"]):.1f} (df 8)')
# only the within-phase shape is identified (sin of the 1st harmonic ~ the light/dark step): compare
# the fitted and true log-odds profiles after removing each phase's mean
def centred(v):
    p = t1b.harmonic_profile(v, 2); out = p.copy()
    for ph in (0, 1):
        m = t1b.PHASE == ph; out[m] -= p[m].mean()
    return out
dev = np.abs(centred(fit2['h']) - centred(h_true)).max()
print(f'max within-phase profile deviation {dev:.3f}')
assert dev < 0.15, dev
assert np.all(np.abs(fit2['g']) < 0.3)
assert fit2['ll'] >= fit0['ll'] - 1e-6
print('t1b checks passed')
