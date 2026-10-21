# PRESPEC 12 — hierarchical (two-timescale) count HMM: is there a transitional COARSE state? (before running)
Date 2026-10-21. Requested by the user, with the long-run aim of mapping onto the scored 3-state HMM.
Motivated by PRESPEC 11: K = 4 = ~2 slow groups x 2 fast states; no third slow group in flat models.
Same data, symbols (0,1,2,>=3), blocks (fly x 4 h ZT bin, days as sequences, fixed clock cuts) as PRESPEC 11.

## Model (structured flat HMM on coarse x fine states)
Coarse states c = 1..Nc, each with Nf = 2 fine states; emissions free per (c,f) (categorical, 4 symbols).
Per-second transitions from (c,f):
  to (c,f')  : (1 - eps_cf) * A_c[f,f']           (fine dynamics inside the coarse state)
  to (c',f') : eps_cf * G[c,c'] * Bin_c'[f'], c' != c
eps_cf = probability of leaving coarse c from fine f; G = coarse jump matrix (zero diagonal);
Bin_c' = entry distribution over fine states of c' (independent of origin). Fitted by EM with the exact
structured M-step (the expected-transition log-likelihood separates into eps, A, G, Bin).
Assumptions stated: Nf = 2 everywhere; entry fine state independent of where you came from; coarse
destination independent of the fine state you left from.

## Coarse labels and routing (the scored-pipeline formulas, applied to coarse states)
Coarse movement fraction mf_c = stationary-weighted mean of (1 - P(count 0)) over its fine states.
Nc = 3: S = lowest mf, W = highest, T = middle.
Xn_c = (G_ST*G_TW - G_SW)/(G_ST*G_TW + G_SW);  Yn_c = (G_WT*G_TS - G_WS)/(G_WT*G_TS + G_WS)
(identical in form to Xn_pred/Yn_pred of compute_HMM_energetics, with kX = G_TW, kY = G_TS).
Mean coarse dwell = Bin_c . (I - Q_c)^-1 . 1, Q_c = within-block transition matrix.

## A (primary): is a transitional coarse state supported?
- Every fly x bin (all 78 flies x 6 bins): leave-one-day-out, fit Nc = 2 and Nc = 3 (2 restarts each,
  best training logL kept), score the held-out day; gain = mean held-out logL difference (3 - 2), per hour.
- Supported in a genotype x bin iff Nc = 3 wins in more flies than expected by chance (one-sided exact
  binomial vs 0.5, p < 0.05) AND the median gain > 0. Twelve genotype x bin tests, all reported.
## B: full fits and coarse routing (all fly x bins)
- Nc = 3 on all days, 6 restarts, best kept; report mf, coarse dwell (S, T, W), Xn_c, Yn_c, fine emission
  rows; restarts within 0.5 nats reported.
- T is called "metastable" only if its mean coarse dwell >= 60 s (same threshold as PRESPEC 11).
## C: identifiability by parametric bootstrap (PRESPEC 11 CV subset: 10 flies/genotype, all 6 bins)
- 10 replicates: simulate from the best Nc = 3 fit with the observed sequence lengths, refit (2 restarts),
  recompute Xn_c, Yn_c (labels re-derived by mf in each replicate).
- Report replicate SD and min-max. Sign of Xn_c (Yn_c) is "determined" in a fly x bin if all 10 replicates
  share the fitted sign. Only sign-determined values are carried forward to any mapping.
## Validation before running on data
Synthetic recovery on a known 3-coarse model with a genuine T (code in hhmm.py, test printed to the ledger).

## Synthetic validation — pass criteria fixed before running (2026-10-21)
Truth P3: 3 coarse x 2 fine; S (silence + rare twitch), T (sparse movement), W (pause + move); 4 days x 4 h
of 1 s data from P3, plus truth P2 (same S and W, no T). Both fitted with hhmm.fit (6 restarts, best kept).
PASS requires all of:
 1. On P3 data: fitted Xn_c and Yn_c each within 0.15 of truth and of the same sign.
 2. On P3 data: fitted coarse dwells (S, T, W) each within a factor 1.5 of truth.
 3. On P3 data: held-out (train 3 days, test day 4) gain of Nc=3 over Nc=2 > 0.
 4. NEGATIVE CONTROL, P2 data: held-out gain of Nc=3 over Nc=2 <= +5 nats/hour (no false transitional state).
 5. Structured M-step never decreases training logL (monotone EM) on both datasets.
If any fails, the data analysis is NOT run; the failure is logged and the model is fixed first.

## Addendum 2 — replicate validation, written AFTER the single-replicate test FAILED criterion 1
The single test cannot separate bias from sampling noise (~100 coarse jumps per 4 x 4 h dataset).
12 independent datasets from P3 (seeds 101..112), each fitted exactly as before (6 restarts, best kept).
Reported: mean error (bias) and SD of Xn_c and Yn_c, and dwell ratios, across replicates.
Decision rule fixed now: the data analysis may proceed only if |mean error| <= 0.10 for both Xn_c and
Yn_c AND the 95% normal-approx. interval of the mean error contains 0. Otherwise the model or the fitting
(restarts) must be changed first. The per-fly bootstrap (part C) is what then quantifies noise in real data.
