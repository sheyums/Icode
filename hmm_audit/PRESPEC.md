# PRESPEC — HMM routing identifiability audit (written before any result was computed)

Date: 2026-10-19. Only data STRUCTURE has been inspected (shapes, model indices,
cell alignment, sequence lengths). No P/B summary, count or route has been computed.

## Data and window acceptance (mirrors compute_HMM_energetics_test)
- CS (wild type): C = CS_RealFinaldatacell (75x500 -> 25 flies), R = CSseq (25x25).
- fumin: C = RealFinaldatacell (159x500 -> 53 flies), R = fuminobservable (53x41).
- C column j pairs with R column j (offset 0; verified structurally). Only j <= size(R,2).
- Accept a window iff model index == 3, P is 3x3, B is 3x2 and R{i,j} is non-empty.
- Deviation from the pipeline: `smooth_observations` was not supplied, so analysis C
  below uses the RAW R sequences. Analyses A and B use only P and B and are unaffected.

## A. Emission degeneracy (descriptive)
- L1_w = sum |B(1,:) - B(2,:)| per accepted window, range [0, 2].
- Thresholds fixed now, ALL reported: L1 < 1e-6 (exactly degenerate), < 0.05, < 0.2.
- Report per genotype: fraction of windows below each threshold.

## B. Identified set of the routing prediction (primary)
- Exact only when L1_w < 1e-6: then every P' = S P S^-1 with S = blkdiag([a 1-a; b 1-b], 1),
  S invertible and P' >= 0 entrywise, gives the same observation likelihood (initial
  distribution transformed accordingly). Computed only for those windows.
- Search: grid a in [-3, 4], b in [-3, 4], 701 x 701 points, plus the identity.
- Quantities: Xn_pred = (P31 kX - P32)/(P31 kX + P32), kX = P12/(P12+P13);
  Yn_pred = (P21 kY - P23)/(P21 kY + P23), kY = P13/(P12+P13). Ungated (no P_min_bio).
- Report per window the interval [min, max] of Xn_pred and Yn_pred over the feasible set.
- Decision rule fixed now: a window's routing direction is UNDETERMINED if its interval
  contains 0. Report the fraction of undetermined windows per genotype, for X and Y.
- Spot check: on 3 windows per genotype (first accepted window of flies 1,2,3), verify
  the log-likelihood of the raw R sequence is equal under P and the extreme P'.
- The grid is a finite search: intervals are INNER bounds (true set can only be wider).

## C. Viterbi vs FFBS route counts (secondary)
- Per accepted window, independently ("within-window"), uniform initial distribution,
  fitted P/B as given, raw R. Viterbi with the pipeline's 1e-300 floor; FFBS exact
  (no floor), M = 200 draws, seed 1.
- Report per genotype, pooled over windows: N32, N312, N23, N213 by Viterbi and FFBS
  posterior mean; Xn, Yn from pooled counts for both; and the per-window paired
  difference (Viterbi - FFBS mean) for N312 and N213 (mean and median).
- Windows where FFBS errors (zero-probability data) are counted and reported, not dropped silently.
- Caveat fixed in advance: both condition on fitted parameters; in degenerate windows
  B shows FFBS counts are not identified either, so C measures decoder bias only.

## Block statement
Narrowest block for every number: a single fitted window (the C/R column); genotype
summaries pool windows across flies and days. No re-windowing or threshold change after results.
