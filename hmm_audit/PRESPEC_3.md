# PRESPEC 3 — is the N312 vs N213 comparison identified? (before any result)
Date 2026-10-19. Same data, window acceptance and grid as PRESPEC 1/2; all 4 genotypes.
Exact H1==H2 windows (S on {1,2}) and exact H1==H3 windows (S on {1,3}), reported separately.

## E (primary). Model-implied route rates over the equivalence set
- Stationary pi of each P (pi' = pi S^-1 for P' = S P S^-1; computed directly from P').
- Jump-chain committors kX = P12/(P12+P13), kY = P13/(P12+P13).
- r312 = pi3 * P31 * kX,  r213 = pi2 * P21 * kY   (expected routes per frame).
- Contrast K = (r312 - r213)/(r312 + r213) in [-1, 1]; also log ratio L = log(r312/r213).
- Per window: range of K over the feasible set. "Invariant" if width < 0.05; "sign determined"
  if the range excludes 0. Per genotype: fractions; window-mean identified set
  [mean of minima, mean of maxima]. Genotype difference identified only if intervals do not overlap.
- Same for r312 and r213 separately (relative range max/min).
## Spot check (counts, not rates)
- First accepted exact-H1==H2 window of flies 1-3 per genotype: FFBS (M=200) posterior mean
  N312 and N213 under fitted P and under feasible P' maximizing and minimizing K; logL each.
