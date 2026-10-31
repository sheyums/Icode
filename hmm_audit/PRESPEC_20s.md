# PRESPEC 20s-A — Plan A: the paper's quantities at the observable (S/Q/W) level (before running)
Date 2026-10-31. Plan A of PLANS.md (user-approved): paper 1 is observable-first. This PRESPEC applies EXISTING
procedures (20b, 20c, 20d, 20j, 20l) unchanged to the observed S/Q/W symbols and the observed first-order chain.
No new tests are introduced; the HMM values are carried beside every observable value as a model-based bracket.
Every quantity is per fly x phase; group values are medians over flies with fly-bootstrap 95% intervals (2000, seed
20); group comparisons, where made, are paired or per-fly as before.

## Data
CS, per0LD, fumin, sss (light, dark), per0DD (whole) - the paper's genotypes - plus CS DD (whole; data uploaded
2026-10-30), reported as an additional row. theta 0.5; symbols and exclusions exactly as 20b / 20i: fly-phases
without a Q/W reference excluded; for chain-based quantities, fly-phases whose observed chain has no cycle
(|J_chain| < 1e-12) excluded and counted (as 20j). HMM-degenerate flies are NOT excluded from observable quantities
(degeneracy is a property of the HMM), only from the bracket.

## A1 Observable table (per fly x phase)
- Direction: minute-level observable current c_obs (sqw_ref.c_obs_minutes; 20d definition).
- Count J (cycles/h) and count EPR (nats/h, rule-of-three bound) from the observed pair counts
  (sqw_ref.observable_pair_counts -> j_cycle_count, epr_count_bound). Sign: forward W->Q->S->W positive.
- Exact observed-chain J, affinity A and EPR from the fitted first-order chain (20n mc_fit; model_energetics).
- Bracket: the HMM's count J and count EPR (20b FFBS, j_points / route2) and exact model EPR, same fly-phase.
- Reported: genotype x phase medians and CIs; fraction of flies with c_obs > 0 and with J > 0; per-fly
  HMM/observable ratios (median, IQR).
- C3 (CS switches on at night): c_obs and count J in the four 6-h ZT blocks (sqw_ref.phases '6h'), CS only, with the
  same statistics; plus the paired comparison night (ZT 12-24 blocks) vs day per fly (Wilcoxon), as 20c.

## A2 EPR vs J and Onsager L (20j procedure, observable points)
- Points: count J and count EPR of A1 (one per fly x phase). Models M1 EPR = b + a|J|, M2 EPR = b + c J^2 with free
  intercept, RSS comparison, fly bootstrap; exponent fit with |J| >= 0.05 and EPR >= 0.05, offset 0.02 as in 20j.
- Onsager L by the user's Eq. 25 (EPR = J^2/L through the origin), onsager_L.py procedure: OLS SE, fly bootstrap CI,
  Wald vs CS with BH across the 4 contrasts.
- Also with the exact observed-chain J and EPR as points (no estimator offset), as a second line.
- Reading fixed: the paper reports whatever these give; if the hidden-level genotype differences in L (20j) do not
  hold at the observable level, the paper says so.

## A3 Energy landscape (20l procedure, observed chains)
- (1) E2 - E1 = ln(pi_Q/pi_W), E3 - E1 = ln(pi_Q/pi_S) from the observed chain's stationary distribution, per fly x
  phase; medians and CIs.
- (2) Landscape prediction F_S->W, F_W->S (sqw_ref.landscape) vs counted routing on the observed symbols (route_events
  on each run of non-missing minutes; >= 3 routes per denominator); median (counted - predicted), Spearman.
- (3) Not repeated (the duration curves are already observable, 20k).

## Output
sqw17.py -> sqw_out/s_observable.json, sqw17.log; figures in A5 (separate step). Ledger entry with the commit hash.
