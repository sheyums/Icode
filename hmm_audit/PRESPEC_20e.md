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

## Amendment A (2026-10-23, before any 20e or 20c result exists; user: ask the question at the HIDDEN level)
- 20e already works at the hidden level: c_pred is computed from the fitted hidden transition matrix. For a
  3-state stationary chain, sign(c_pred) = sign(forward cycle affinity ln[P21 P13 P32 / (P12 P31 P23)])
  = sign(hidden cycle flux J) (algebra: Yn > Xn iff P_{2->1} P_{1->3} P_{3->2} > P_{3->1} P_{1->2} P_{2->3}, H indices).
  Checked on the 20b fits: identical signs in 292 of 301 fly-phases; the other 9 have a numerically zero cycle
  edge (|J| ~ 1e-17, infinite affinity), i.e. no cycle. A sign flip of c_pred IS a flip of the hidden J.
  Fly-blocks with |J| < 1e-12 are classed "no cycle" and excluded from flip counts (counted and reported).
  J (cycles/h) is reported beside c.
- Section 3 becomes a TEST for each genotype that 20c step 2 judges USABLE (descriptive otherwise):
  no-change null within each phase: 8 records per fly simulated from its 20b fit (P_light for blocks 1-2,
  P_dark for 3-4), refitted with the 6 h model (K = 3, 10 starts); per fly null flip probability between blocks
  1-2 and between blocks 3-4; observed flip count vs Poisson-binomial, one-sided. Also reported, descriptively,
  flips across the light-dark boundaries (2-3, 4-1).
