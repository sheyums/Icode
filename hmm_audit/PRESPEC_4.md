# PRESPEC 4 — per-fly scatter of Xn vs Yn (before plotting)
Date 2026-10-19. Uses only results already computed (C_*.npy, B_*.npy; within-window FFBS M=200, raw R).
- Windows: exact H1==H2 windows only (L1_12 < 1e-6), the dominant class in every genotype.
  H1==H3 and non-degenerate windows are excluded so the two H1 classes are never pooled.
- Fly point (filled): Xn = (sum N312 - sum N32)/(sum N312 + sum N32), FFBS posterior mean counts
  summed over the fly's windows (count-based, as in the pipeline's global Xn). Yn likewise
  with N213, N23. Legacy Viterbi point (hollow) computed the same way from Viterbi counts.
- Fly range (bars): identified set of the fly's window-MEAN model Xn_pred / Yn_pred,
  [mean of per-window minima, mean of per-window maxima] from analysis B. Marginal, not joint;
  it is the model analogue of the count-based point, not the exact range of the pooled counts.
- Every fly with >= 1 qualifying window and a defined point is plotted; n reported per genotype.
- Axes fixed to [-1, 1] x [-1, 1]. Four panels (one per genotype), same scales.
