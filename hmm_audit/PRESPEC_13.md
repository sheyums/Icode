# PRESPEC 13 — POOLED-FLY hierarchical count HMM (before running anything)
Date 2026-10-21. Follows PRESPEC 12: per fly-bin coarse routing is not estimable (replicate SD 0.6-0.7
with the true model). Here coarse AND fine parameters are SHARED across all flies of a genotype within
a 4 h ZT bin (each fly-day segment a separate sequence). Model, symbols, labels, Xn_c/Yn_c formulas and
dwell definition exactly as PRESPEC 12 (hhmm.py).

## Fitting (fit_pooled)
Warm start: 6 random-start fits (Nc given) on a random subset of 12 sequences (seeded), the one with the
best FULL-data logL is used to initialise a fit on all sequences (n_iter <= 300, tol 0.01).

## Part V — synthetic validation at pooled size (must PASS before Part D runs)
Size: 100 sequences x 14400 s (= CS: 25 flies x 4 days).
V1 homogeneous P3 (truth of PRESPEC 12), 6 replicates (seeds 201..206):
   PASS iff, for Xn_c and Yn_c: |mean error| <= 0.10 AND SD <= 0.15; and mean dwell ratio of S, T, W each
   in [0.8, 1.25].
V2 HETEROGENEOUS negative control (no T): P2 with per-fly variation (each fly: eps x exp(N(0,0.3)), each
   emission row mixed 10% with a Dirichlet(1) draw), 2 replicates (seeds 301, 302). Train on days 1-3 of
   all flies, test day 4. PASS iff held-out gain Nc=3 over Nc=2 <= +5 nats/hour per fly (median over flies)
   in both replicates — i.e. fly heterogeneity does not create a false transitional state.
V3 HETEROGENEOUS P3 (same perturbation), 2 replicates (seeds 401, 402): reported (errors), not a gate.
If V1 or V2 fails, Part D is not run; the failure is logged.

## Part D — data (only if V passes)
Per genotype x 4 h bin (12 blocks), all flies:
D1 Is a transitional coarse state supported? Fit Nc=2 and Nc=3 on all days but the last; score each fly's
   last day. Per-fly gain (nats/hour). Supported iff median gain > 0 AND one-sided sign test across flies
   p < 0.05. All 12 reported.
D2 Nc=3 on all days: mf, dwell (S,T,W), Xn_c, Yn_c, fine emission rows. T "metastable" iff dwell >= 60 s.
D3 Uncertainty: fly-level bootstrap, 20 resamples of flies (seed 13), each refitted warm-started from the
   D2 fit (n_iter <= 100); 95% percentile intervals for Xn_c, Yn_c, T dwell. Sign "determined" iff the
   interval excludes 0. Genotype comparison per bin: difference of bootstrap draws, interval excludes 0.
Only D-results with determined sign are carried to the planned mapping onto the scored 3-state HMM.

## Addendum (2026-10-21 15:27, before any validation value was viewed): optimizer-vs-data check
Run 1 of Part V was stopped after 1 of 10 jobs (values not viewed; kept as validate13_run1_aborted.*)
to add, for every V1 and V3 replicate: dLL = logL(fitted pooled model) - logL(true generating model)
on the replicate's own data (V3 truth = each fly's own perturbed model). Read-out: dLL < -1 nat
means the optimizer stopped BELOW the truth (optimizer failure); dLL >= -1 means the fit is at least
as good as the truth, so any error in Xn_c/Yn_c is statistical (data), not the optimizer.
Reported per replicate; not a pass/fail gate (the V1/V2 gates are unchanged).

## Amendment 2 (2026-10-21, user decision, BEFORE Part D has run): DAY COVARIATE (option 2)
Reason: the descriptive drift check (aging_check.log) showed systematic change across recording days.
Model for Part D becomes hhmm_day.HHMMDay: coarse exit probabilities eps_{c,f,d} = eps_{c,f} * kappa_{c,d},
with sum_d log kappa_{c,d} = 0 (one free day factor per coarse state per day). Shared across days:
fine dynamics A_c, emissions, coarse jump matrix G, entry distributions. So coarse DWELLS may change by
day; ROUTING (G, hence Xn_c, Yn_c) is assumed day-invariant and estimated from all days.
Changes to Part D:
- D1: a held-out DAY cannot be predicted (its kappa is unknown), so the test uses held-out FLIES: a fixed
  20% of flies per genotype (numpy default_rng(13), sorted filenames) are held out; Nc=2 and Nc=3 fitted
  on the rest (all days); per held-out fly gain (nats/hour, all its days). Read-out unchanged (median > 0
  and one-sided sign test p < 0.05).
- D2: as before on all flies, plus kappa per coarse state and day, and dwell per day.
- D3 unchanged (fly bootstrap), using the day model.
New validation gates for the day model (run after Part V; Part D needs V1, V2 AND V5, V6 to pass):
- V5: P3 with day effect on S exits (kappa_S = exp([-.4,-.13,.13,.4])), 25 flies x 4 days, 3 replicates
  (seeds 501-503): PASS iff |mean error| <= 0.10 and SD <= 0.15 for Xn_c and Yn_c, and fitted kappa_S
  within +-25% of truth on average for every day.
- V6: heterogeneous P2 (as V2) WITH the same day effect, 2 replicates (seeds 601, 602), held-out flies as
  in D1: PASS iff median held-out gain Nc=3 over Nc=2 <= +5 nats/hour per fly in both.
