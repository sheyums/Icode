# PRESPEC 8 — is sleep/wake time-irreversible, and by how much? (APPROVED by user 2026-10-19; unchanged from draft 1f4dd8e; not yet run)
Date 2026-10-19. Written after PRESPEC 5b results were seen (CS dark: A = +0.127, CI [0.024, 0.231]).
That result MOTIVATED this plan, so a CS-dark result here on the same flies is NOT an independent
confirmation; independent confirmation needs flies not used here.

## Principle
Under time-reversibility (detailed balance) a stationary sequence and its time reversal have the
same law. With alternating sleep (S) and wake (W) bouts, reversal maps each pair
(S_i, then W_{i+1}) onto a pair (W_j, then S_{j+1}) read backwards. So under H0 the joint law of
SW pairs (s, w) equals the joint law of WS pairs with coordinates swapped, (s = S_{j+1}, w = W_j).
This is model-free: it needs no HMM, no state count and no identifiable parametrization.

## Data and blocks
- Bouts as PRESPEC 5 (interior runs; sleep >= 300 s, wake >= 2 s). Segment = one 4 h WINDOW
  (narrowest block; limits circadian nonstationarity, which breaks exactness of the null).
- Blocks reported: genotype x {whole record, light ZT0-12, dark ZT12-24} (pair phase = phase of
  its first bout). All 12 reported. No other split.

## Statistic
- Each bout class is binned by its genotype-wide quantiles (reversal leaves both marginals
  unchanged, so the bins are not tuned to the asymmetry). Bin sets fixed now: 3x3, 4x4, 5x5;
  all three reported.
- D = KL( P_SW || P_WS,swapped ) in nats per pair, plug-in with +0.5 pseudocount per cell.
  By the data-processing inequality, D (bias-corrected, below) is a LOWER bound on the entropy
  production per S-W cycle (Roldan & Parrondo, PRL 105:150607, 2010); with the mean cycle
  duration it gives a lower bound in nats per hour.

## Null and bias correction
- Reversal null: each window's bout sequence is independently kept or time-reversed with
  probability 1/2, then pairs are re-formed and D recomputed. 2000 replicates, seed 8.
  Exact under stationary reversibility within a window; same n by construction.
- p = (1 + #null >= obs)/2001. Bias-corrected bound D* = D_obs - mean(D_null), floored at 0.
- Confidence interval for D*: fly-level bootstrap (2000), recomputing D* with the null mean
  re-estimated per resample (500 reversal replicates each).
## Read-out fixed now
- Irreversibility is demonstrated in a block if p < 0.05 at ALL three bin sets.
- 12 blocks x 3 bin sets are run; report the number of blocks expected to pass by chance
  under the global null (estimated from the reversal replicates themselves).
- Reported for every block: D*, its CI, D* in nats/hour, p for each bin set.
## Not in this plan
- The model-based route (aggregated Markov chain with the Hawkes-Jalali-Colquhoun exact
  missed-event correction for the 300 s sleep criterion, fitted with and without the
  Kolmogorov cycle constraint). Heavier; only if the model-free test finds irreversibility.

## Implementation notes fixed before running
- Reversing a window swaps its SW and WS pair-count matrices exactly (an SW pair read backwards
  is a WS pair with swapped coordinates), so each null replicate is a random swap per window.
- A pair's phase = phase of its earlier bout (original time), and is unchanged by reversal.
- The same random flips are used for the three bin sets, so the global-null pass rate (p < 0.05
  at all three) is estimated jointly from the replicates.
- Lower bound in nats/hour = D* / (mean sleep + mean wake duration of the block's SW pairs, in hours).
