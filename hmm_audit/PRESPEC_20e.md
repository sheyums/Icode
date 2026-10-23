# PRESPEC 20e — does a fly's cycle direction change with time? (before running; written while the 20b day split
# had not yet run, so its results were unseen)
Date 2026-10-23. User question. Statistic: c_pred = (Yn_pred - Xn_pred)/2 from a fitted K = 3 S/Q/W model
(PRESPEC 20b definitions, theta 0.5). Per-fly estimates are noisy (bootstrap SD ~0.1-0.15 for c at 2880 min, more
for shorter blocks), so a sign change in one fly is evidence of nothing unless it exceeds what noise produces
with NO change. Every comparison below is therefore against a no-change null built from the same fly.

## 1. Across days (from the 20b day split: first half of the days vs second half, per phase)
- Per genotype x phase, flies non-degenerate in both halves: Spearman(c_first, c_second); mean change
  c_second - c_first with 95% fly-bootstrap interval; fraction of flies whose sign flips.
- No-change null: for each fly, 8 records simulated from its FULL-record 20b K = 3 fit (constant parameters by
  construction), split at the same day, each half refitted (K = 3, 10 starts); per fly, the null flip
  probability p_f = fraction of the 8 replicates whose two halves differ in sign.
- Test: observed number of flipping flies vs the Poisson-binomial distribution with the p_f (one-sided, more
  flips than noise). Also the null distribution of the mean change (each replicate's flies averaged).
- A flip excess means some flies change direction across days beyond noise; no excess means the data cannot
  distinguish changing from stable direction at this resolution (NOT that direction is stable).
## 2. Light vs dark (20b, same fly; descriptive)
Per genotype (not per0DD): Spearman(c_light, c_dark), fraction of flies with opposite signs, with the expected
opposite-sign fraction from the 20b parametric bootstrap SDs (normal approximation around each fly's estimates).
## 3. Within the day (20c, CS / per0LD / fumin; descriptive)
Per fly, c in each 6 h block; fraction of flies whose sign differs between blocks; noise reference = 20c step 2
(realistic-truth) SDs. Reported with 20c's USABLE verdict per genotype.
## Order
After PRESPEC 20d. Section 1's null costs ~2-3 h (all flies, 8 replicates).
