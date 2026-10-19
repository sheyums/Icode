# PRESPEC 5 — pooled, identifiable analysis of sleep/wake bout structure (before any result)
Date 2026-10-19. Structure inspected only: min interior inactive run = 300 frames, min active
run = 2 frames, every boundary between consecutive non-empty windows switches state (windows
were cut at bout boundaries). No duration distribution, fit or correlation has been computed.
Implemented in Python (no MATLAB here), following the Icode engine's discrete convention
P(N=n | N>=n_min) = [S((n-1)dt) - S(n dt)] / S((n_min-1)dt). Not a run of Icode's .m files.

## Data
- All 4 genotypes; ALL non-empty R windows (not only HMM-accepted ones: acceptance depended on
  the per-window model choice). dt = 1 s.
- Per fly, consecutive non-empty columns are treated as one contiguous segment (assumed
  physically contiguous, as compute_HMM_energetics_test does); an empty column ends a segment.
- Bouts = runs. Sleep = inactive (obs 1), xmin = 300 s. Wake = active (obs 2), xmin = 2 s.
- The first and last run of every segment are excluded (censored by the record edge).
  Interior window boundaries censor nothing (windows were cut at bout boundaries).
- Block: genotype x bout class, whole record. No ZT split (recording ZT start not supplied).

## P1. Number of exponential phases per bout class (identifiable, replaces HMM state count)
- Left-truncated discrete hyperexponential, K = 1..4, per genotype x class. Multi-start ML
  (30 starts per K). Report q (observed weights), component means, smallest-branch count q_min*n.
- Selection, primary: 5-fold cross-validation BY FLY (fold = (fly index - 1) mod 5),
  summed held-out log-likelihood; choose the K with the largest CV log-likelihood.
- Also reported (not used for selection): AICc and BIC with n = bouts, k = 2K - 1.
- Read-out: if sleep needs K >= 2, the 3-state HMM (1 sleep + 2 wake phases) is too small.

## P2. Coupling between consecutive bouts, and its direction (identifiable)
- Within each contiguous segment, pairs SW = (sleep_i, the wake that follows) and
  WS = (wake_i, the sleep that follows).
- Each duration is converted to a within-WINDOW centred rank, u = (rank - (n+1)/2)/n, among
  bouts of the same class in that window (removes fly and time-of-day level differences);
  bouts in windows with < 3 bouts of their class are dropped. Pairs may straddle windows.
- rho_SW = Pearson corr of u over SW pairs (pooled within genotype); rho_WS likewise;
  asymmetry A = rho_SW - rho_WS.
- Coupling test: within-window permutation of durations (each class separately), 2000
  permutations, two-sided p for rho_SW and rho_WS. Same n as observed by construction.
- Asymmetry: fly-level bootstrap (2000 resamples of flies), 95% percentile CI for A and for
  each rho; A is called nonzero only if its CI excludes 0. Genotype differences in A:
  bootstrap CI of the difference (independent resampling per genotype) excludes 0.
- Read-out: in a lumped Markov chain, consecutive-dwell correlation requires >= 2 phases in
  BOTH classes (Colquhoun & Hawkes; Fredkin & Rice); A != 0 indicates broken time-reversal
  symmetry (Song & Magleby 1994), the identifiable counterpart of sleep->wake vs wake->sleep
  asymmetry. A circadian trend is only partly removed by within-window ranking: say so.
- Seeds: numpy default_rng(5). All four genotypes reported, including nulls.
