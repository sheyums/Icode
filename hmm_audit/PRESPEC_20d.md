# PRESPEC 20d — is the model's larger circulation real or manufactured? (before running)
Date 2026-10-23. Requested by the user after the 20b primary results: c_model is 1.5-6x the model-free c.
Coarse-graining can only LOSE irreversibility (Roldan & Parrondo 2010; Esposito 2012; Lynn, Holmes, Bialek &
Schwab, PRE 106, 034102 and PRL 129, 118101, 2022), so a hidden model SHOULD show more -- but c is not an
entropy production, the two c's have different definitions, and a fitted model can create circulation by
bias. Two checks separate these. Model, fitting and routing code as PRESPEC 20b (theta 0.5).

## (a) Reversible null (synthetic)
Truth: 3 states, the V1b emission matrix (distinct H1), transition matrices obeying detailed balance in each
phase, built from symmetric flows F_ij = pi_i P_ij:
- light: pi = (0.1, 0.4, 0.5); F_01 = F_02 = 0.0165, F_12 = 0.0115;
- dark:  pi = (0.1, 0.3, 0.6); F_01 = F_02 = 0.0165, F_12 = 0.0045.
True c = 0, J = 0, EPR = 0 in both phases. 20 flies x 5760 min, fit K = 3 (20 starts), FFBS routing.
As a positive reference, the same on the V1b truth (true c, EPR computed from its P).
Reported per phase: mean c_pred and c_ffbs with fly-bootstrap 95% intervals; distribution of fitted EPR
(the EPR "noise floor" at this record length).
**Pass:** |mean c| < 0.05 and its interval contains 0, for both c_pred and c_ffbs, in both phases.
A fail means the pipeline manufactures circulation, and the size of the bias is reported next to every fly result.
Fly EPRs are reported against the null's 95th percentile (fraction of flies above it).

## (b) Attenuation, per fly (posterior-predictive)
Observable-level current on the minute symbols, defined identically for real and simulated records:
c_obs_min = P(last wake minute before a sleep onset is Q) - P(first wake minute after a waking is Q),
onset = a wake symbol (Q or W) followed by S, waking = S followed by a wake symbol, phase = that of the second minute of the pair.
For every non-degenerate fly of the 20b fits: 20 records simulated from its fitted model (its length, phases).
On each simulated record: c_obs_min, and c_true from the TRUE hidden path's route counts.
Reported per genotype x phase:
- real mean c_obs_min vs the simulated distribution of the group mean (flies averaged within each replicate);
  **consistent** if the real value lies in the central 95% of that distribution;
- per fly: fraction of flies whose real c_obs_min lies outside its own simulated 95% band (expect ~5%);
- attenuation ratio: mean simulated c_true / mean simulated c_obs_min, beside the real ratio
  c_model(ffbs) / c_obs_min.
Reading, fixed now: if (a) passes and (b) is consistent, the gap between model and observable current is what the
fitted hidden dynamics predict from coarse-graining alone. If (b) is inconsistent for a genotype x phase (expected
candidate: per0LD), the model's circulation there is not supported by the data.

## Not included
EPR against the model-free irreversibility bound (E3) needs E3 per hour (TODO.md); a separate plan.

## Amendment (2026-10-23, before 20d runs; user remark that hidden and observed currents may legitimately differ)
Added to (b)'s report, no change to any test: per genotype x phase, the fraction of simulated records in which
c_true and c_obs_min have OPPOSITE signs, and the same per fly (fraction of its 20 replicates). Reading: the
model is judged only by whether its SIMULATED observable current matches the REAL observable current (same
level). A hidden current differing from the observable one in size or sign is a result about coarse-graining,
not a model failure.
