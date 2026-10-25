# PRESPEC 20j — per-fly EPR vs J_cycle (hidden): linear vs quadratic, slopes by genotype (before running)
Date 2026-10-25. User request. Descriptive in origin (group medians were seen in Fig. 6), definitions fixed here.
Points: one per fly x phase (light open, dark/DD filled), 20b fits, theta 0.5. J = count J_cycle (FFBS path mean,
cycles/h, forward +); EPR = rule-of-3 bound, per path then averaged (nats/h). Excluded: degenerate flies, 'no cycle'
fly-phases (|J_model| < 1e-12), fly-phases without a Q/W reference (PRESPEC 20i).
## Models per genotype (light and dark pooled; also each phase separately)
- M1 constant affinity: EPR = b + a |J|      (a = typical affinity, nats/cycle)
- M2 linear response:   EPR = b + c J^2      (c = 1/L, Onsager conductance L)
b is a free intercept: the hidden count estimator has an additive offset of ~+0.01-0.025 nats/h (PRESPEC 20f), so
forcing the fit through the origin would bias both models. Unweighted least squares; both have k = 2, so they are
compared by residual sum of squares (equivalently AIC). Confidence intervals: fly bootstrap (2000, seed 20),
resampling flies with both their phases.
- Power-law exponent: log(EPR - 0.02) = log k + n log|J| on fly-phases with |J| >= 0.05 and EPR >= 0.05;
  n ~ 1 constant affinity, n ~ 2 linear response; bootstrap CI.
- Slopes between genotypes: difference in a (M1) between each genotype and CS, bootstrap CI.
## Caveats fixed now
- J and EPR are computed from the same counts, so their errors are correlated; for a unicyclic count network
  EPR_count ~ J x (sum of log count ratios). The question answerable here is whether the affinity is roughly constant
  across flies or grows with J, not whether EPR and J are related at all.
- Weak-circulation fly-phases sit at the estimator offset; the exponent fit excludes them by the thresholds above.
