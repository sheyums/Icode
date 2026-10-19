# PRESPEC 7 — two kinds of sleep bout within ONE fly at ONE time of day? (before running)
Date 2026-10-19. Written while PRESPEC 5b is running and before PRESPEC 6 has started;
no output of either has been viewed. User judged this the more relevant question.

## Data and block
- Sleep bouts as PRESPEC 5 (interior runs, xmin = 300 s, dt = 1 s). Bout time = start frame,
  empty columns counted as 14400 frames (as PRESPEC 5b); recording starts at ZT0.
- Block = fly x ZT bin, 6 bins of 4 h: [0,4) [4,8) [8,12) [12,16) [16,20) [20,24), pooling the
  fly's days. Blocks with >= 20 sleep bouts are analysed; the number excluded is reported.
- Day = floor(start_frame / 86400) + 1.

## Test per block: K=1 vs K=2 (discrete left-truncated hyperexponential, pooled.py fit())
- LR = 2[logL(K=2) - logL(K=1)]; K=2 fitted with 10 starts, K=1 with 3.
- Mixture order is non-regular, so p comes from a parametric bootstrap, 100 replicates, under
  two nulls, each replicate refitting both K exactly as the observed:
  N1 (single exponential): one rate for the block (its K=1 MLE).
  N2 (day heterogeneity): one rate PER DAY within the block (plug-in K=1 MLE per day, same n per
     day). N2 is the stricter null: it lets day-to-day drift produce the mixture.
- p = (1 + #null LR >= obs)/(1 + 100).
## Read-out per genotype (fixed now)
- Fraction of blocks with p < 0.05 under N1 and under N2; exact one-sided binomial test of that
  count against 5% (the rate expected if no block has two kinds of bout).
- PRIMARY: N2. "Two kinds of sleep bout within a fly at a time of day" is supported for a genotype
  only if its N2 binomial p < 0.05.
- For blocks significant under N2: median q of the shorter component and the two component means.
- Also reported, not tested: the same fractions per ZT bin (pooled over genotypes' flies).
- Seeds: SeedSequence([7, genotype, fly, bin]).

## Amendment (implementation only), 2026-10-19, before any within-fly result existed
The first run was stopped by a 2 h job limit having printed only block counts (CS 94, per0 113,
fumin 104, sss 4 analysed). Runtime with L-BFGS (13 fits per LR) was ~10 h. The K=2 MLE is now
computed by `em_fit.fit2g`: a deterministic profile grid over the two component means (40 x 40,
q profiled by EM, concave in q) with Nelder-Mead polish from the best 3 grid points; K=1 is the
closed form p = 1/mean(e). Same likelihood (a left-truncated discrete hyperexponential is a
mixture of geometrics on e = d - 299 with observed weights q). On 60 synthetic datasets it was
never below the better of L-BFGS and multi-start EM by > 1e-3 in LR and exceeded it in 3.
Everything else (blocks, nulls, replicates, read-out) unchanged.
