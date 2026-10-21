# PRESPEC 15 — what dissipation do the data REQUIRE? (step 1) and model-free transition signatures (E)
Date 2026-10-21, before running. User's question: posit non-equilibrium (cycle flux != 0) and look for the
hidden routes. Likelihood is flat along the equivalence family, so the posit cannot select a route; the
well-posed version is an extremal one: the MINIMUM dissipation consistent with the data.

## Step 1 — scored 3-state model, exactly degenerate windows (all 4 genotypes)
Windows and equivalence family exactly as PRESPEC 1/2 (H1==H2: S on {1,2}; H1==H3: S on {1,3}),
grid 201 x 201 on [-3,4]^2 plus identity, then Nelder-Mead refinement of (a,b) from the 5 lowest-EPR grid
points (infeasible = any P' entry < 0, penalised). For every feasible P':
- pi = stationary distribution; J_ab = pi_a P_ab.
- EPR (nats/hour) = 3600 * sum_{a<b} (J_ab - J_ba) ln(J_ab/J_ba); a pair with both fluxes 0 adds 0; a pair
  with exactly one zero flux makes EPR infinite (excluded from the finite maximum, reported as a count).
- Cycle flux J = pi_3 P_32 - pi_2 P_23 (pipeline's J_P), per hour.
Per window: EPR at the fitted P; EPR_min, finite EPR_max; J range; Xn_pred and Yn_pred at the EPR_min point.
Read-out fixed now: "data require no dissipation" in a window iff EPR_min <= 1e-3 nats/hour.
Report per genotype: fraction of windows with EPR_min <= 1e-3; median EPR_min; median fitted EPR; fraction
of windows whose J range contains 0.

## E2/E3 — count data (plan details written as an addendum after step 1, before running E)
E2 time-reversal asymmetry of activity around transitions; E3 windowless irreversibility lower bound of the
coarse-grained count sequence (Roldan & Parrondo). Step 1's outcome does not change E's design.

## Addendum — E2/E3 details (written after step 1 results, before running E; step 1 does not inform E)
Count data, 1 s; sleep onset = first second of a zero-count run >= 300 s; waking = first second after it.
### E2 — time-reversal asymmetry of wake-bout edges
Eligible wake bout = interval between two sleep bouts lasting >= 20 min. For each: EXIT profile = movement
fraction in 1-min bins at seconds +61..+600 after waking; APPROACH profile = the same bins counted backward
from the next sleep onset (seconds -61..-600). Minute 1 is excluded on both sides (a movement is guaranteed
there by construction). Under time-reversibility the two profiles are exchangeable within each bout.
Per fly x phase (phase of the wake bout's start): D = mean over bouts and bins of (EXIT - APPROACH).
Genotype x phase statistic: mean of per-fly D (flies with >= 3 eligible bouts).
Null: 2000 permutations swapping EXIT and APPROACH per bout (independently), recomputing the statistic.
Two-sided p. Also reported: the mean profiles, bin by bin.
### E3 — windowless irreversibility lower bound
Each fly's record is cut at FIXED 4 h clock times into segments (no bout-based cuts). Within a segment,
consecutive runs of still (count 0) and moving (count > 0) seconds are classed:
still <10 s 'S0', 10-299 s 'S1', >=300 s 'S2' (sleep); moving 1-2 s 'M0', >=3 s 'M1'.
Statistic: KL divergence D between the distribution of run triples (r_i, r_i+1, r_i+2) and of their
reversals (r_i+2, r_i+1, r_i), +0.5 pseudocount per cell, pooled over flies within genotype x phase.
Null: each segment kept or time-reversed with probability 1/2, 1000 replicates (exact under stationary
reversibility). p = (1 + #null >= obs)/1001; bias-corrected D* = D - mean(null), floored at 0.
Lower bound in nats per hour = D* x (runs per hour) / 3 is NOT used (overlapping triples); report D* in
nats per triple and the number of triples. Irreversibility detected iff p < 0.05.
