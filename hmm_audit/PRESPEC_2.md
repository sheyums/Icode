# PRESPEC 2 — can the Xn vs Yn ASYMMETRY survive non-identifiability? (before any result)

Date 2026-10-19. New data: sss (Sssseq 57x18, sss_RealFinaldatacell 171x500) and per0
(Per0seq 25x25, per0_RealFinaldatacell 75x500). Only shapes inspected.
The user keeps the Xn vs Yn analysis; this asks whether the COMPARISON is identified even
though Xn and Yn individually were not (PRESPEC 1, results in RESULTS_LEDGER.md).

Window acceptance, C/R pairing (offset 0) and raw-R deviation: as in PRESPEC 1, all 4 genotypes.

## A, B, C on sss and per0: exactly as PRESPEC 1.

## D (primary). Joint identified set of (Xn_pred, Yn_pred), all 4 genotypes
- Windows: exact H1==H2 (L1_12 < 1e-6), S on block {1,2}; and exact H1==H3 (L1_13 < 1e-6),
  S on block {1,3}. Same 701x701 grid on [-3,4]^2 plus identity. Feasible: P' >= 0, S invertible.
- Contrast: Delta = Xn_pred - Yn_pred, evaluated JOINTLY at each feasible S (not from the two
  marginal intervals). Also Sigma = Xn_pred + Yn_pred.
- Per window report the range of Delta and Sigma. A quantity is "invariant" in a window if
  its range width < 0.05. Sign of Delta is "determined" if its range excludes 0.
- Per genotype: fraction of windows with Delta invariant; with sign(Delta) determined;
  same for Sigma, Xn, Yn.
- Genotype level: windows are fitted independently, so the identified set of the
  window-MEAN of Delta is [mean of per-window minima, mean of per-window maxima]
  (windows where Delta is undefined everywhere excluded and counted). Report for each
  genotype and each H1-twin class separately (never pooled across the two classes).
- Decision rule fixed now: a genotype difference in mean Delta is identified only if the
  two genotype intervals do not overlap.
- Spot check: first accepted exact-H1==H2 window of flies 1-3 per genotype: FFBS (M=200)
  posterior-count Xn, Yn under fitted P and under the feasible P' maximizing and
  minimizing Delta; log-likelihoods under each (transformed pi0).

## Block statement
Narrowest block: one fitted window. Genotype summaries pool windows over flies and days.
