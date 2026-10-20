# PRESPEC 11 (DRAFT — not approved, not run) — flat count HMM, K = 2..6, as an assumption-free baseline
Date 2026-10-20. Written after PRESPEC 10 C3 results (K=3 spends states on second-scale movement;
K=4 beats K=3 on the held-out day in 98-100% of flies). Same data and symbols as PRESPEC 10 C3:
count 0,1,2,>=3 per second. REVISED BLOCK (user discussion, before any run): fly x 4 h ZT bin
(ZT 0-4, 4-8, 8-12, 12-16, 16-20, 20-24, cut at fixed clock times, NOT bout-adjusted), each day's
segment of that bin a separate sequence sharing one parameter set (CS 4 sequences, fumin 3; a
segment shorter than 90% of 14400 s is dropped). Rationale: 12 h halves contain within-phase drift,
which a time-homogeneous HMM represents with spurious slow states; hmmlearn CategoricalHMM, EM, tol 0.01, n_iter <= 500.

## Stage 1 — how many states? (cross-validation, no parameter counting)
- CV subset fixed now: 10 flies per genotype drawn with numpy default_rng(11).choice(n, 10, replace=False)
  over sorted filenames; all six ZT bins. (Full CV on all 78 flies would take ~18 h on 4 cores.)
- Leave-one-day-out: train on the other days' segments of the bin, score the held-out day's segment; K = 2..6; 3 restarts per fit
  (random_state 0..2), best training logL kept. Score = held-out log-likelihood per second, averaged
  over folds.
- Per fly x bin: K* = argmax score. Per genotype x bin: K_g = the mode of K* over the 10 flies
  (tie -> smaller K). Also reported: median score gain of each K over K-1.

## Stage 2 — full fits and what the states are (all flies)
- Every fly x bin fitted at its genotype x bin's K_g on all days, 10 restarts; best kept; number of
  restarts within 0.5 nats of the best reported (optimizer reliability, NOT identifiability).
- States reported sorted by P(count 0): emission rows, mean dwell 1/(1-P_ii), stationary occupancy.

## Stage 3 — do slow behavioural groupings emerge without being imposed? (metastability)
- Implied timescales of the fitted P: tau_i = -1/ln|lambda_i| (s), i = 2..K.
- "Slow mode" = tau >= 60 s (primary). Sensitivity thresholds fixed now: 30 s and 300 s, all reported.
- m = number of slow modes per fly x bin; distribution of m per genotype x bin, and pooled per genotype.
  m = 1 -> two metastable sets (e.g. sleep/wake); m = 2 -> three (e.g. sleep/transitional/wake).
- For m >= 1: PCCA+ into m+1 macrostates (inner-simplex solution, Deuflhard & Weber 2005, without
  the subsequent optimisation step), memberships chi; coarse transition matrix
  P_c = (chi' D chi)^-1 chi' D P chi, D = diag(stationary pi).
- Macrostate description: membership-weighted movement fraction, occupancy, mean dwell from P_c.
  Labels by movement fraction: lowest = S (sleep-like), highest = W (wake-like), others = T.
## Stage 4 (descriptive only) — routing through T, where m = 2
- Transition path theory on the microstate P (Metzner, Schuette & Vanden-Eijnden 2009), with S and W
  = microstates of membership > 0.5 in the S and W macrostates. Report the fraction of reactive
  S->W flux whose path visits T-microstates, and the same for W->S. The two fractions are the
  count-HMM analogues of Xn and Yn.
- NOT an identifiability-checked quantity: reported with the range over restarts within 0.5 nats
  (where >= 2 exist), and explicitly labelled descriptive. A profile-likelihood check is a separate plan.
## Read-out fixed now
- Primary: the distribution of m (60 s threshold) per genotype x bin, and whether m = 2 (three
  metastable sets) is the most common value in any genotype x bin.
- Day-to-day drift check: held-out-day score of the chosen K, by fold; a fold far below the others
  flags a non-stationary day (reported, not used to drop data). Everything else is description.
## Runtime estimate
Stage 1 ~ 3-4 h, Stage 2 ~ 5-6 h on 4 cores (6 fits per fly on 1/3 the data each: similar total) (K up to 6 is ~4x slower per iteration than K = 3).
