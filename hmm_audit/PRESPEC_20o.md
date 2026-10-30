# PRESPEC 20o — what does the HMM's gain over the plain chain consist of? Memory and drift checks (before running)
Date 2026-10-30. Follows 20n (and its correction): the first-order chain on S/Q/W (MC1) is rejected in 16/20 flies, but
the source of the HMM's gain is unidentified: (1) Q/W misclassification, (2) memory beyond one minute, (3) drift of the
rates within a 12-h phase. Genotypes CS (25) and fumin (53); symbols, theta 0.5 and missing coding exactly as 20b/20n.

## Models (all fitted per fly to the same symbols; every likelihood is the full-record P(y_1..y_T), missing minutes
## marginalised by the forward algorithm, so all are directly comparable)
| name | structure | phases | k |
| MC1_12 | 1st-order chain (20n) | light/dark | 2 + 2*6 = 14 |
| HMM_12 | 3-state HMM (20b fit, not refitted) | light/dark | 20 |
| MC2_12 | 2nd-order chain: P(y_t given y_t-2, y_t-1) | light/dark | 8 + 2*9*2 = 44 |
| MC1_6  | 1st-order chain | four 6-h ZT blocks | 2 + 4*6 = 26 |
| HMM_6  | 3-state HMM, shared emissions, 20 starts, seed as 20b | four 6-h blocks | 2 + 4*6 + 6 = 32 |
| MC2_6  | 2nd-order chain | four 6-h blocks | 8 + 4*18 = 80 |
MC2 is fitted exactly as a 9-state chain on pairs (y_t-1, y_t) with deterministic emissions (second element) and
structural zeros for inconsistent pairs; em_step preserves both (asserted). On complete data it must equal the
closed-form count MLE (checked before the run). Phase of a transition = phase of its arrival minute, as everywhere.

## Quantities (per fly)
- dLL of every model against MC1_12; per 1000 min.
- Share of the HMM's gain matched by one more minute of memory:  R_mem = (LL_MC2_12 - LL_MC1_12) / (LL_HMM_12 - LL_MC1_12)
- Share of the HMM's gain that disappears once 6-h drift is allowed:  R_drift = 1 - (LL_HMM_6 - LL_MC1_6) / (LL_HMM_12 - LL_MC1_12)
- Drift alone: LL_MC1_6 - LL_MC1_12.
- Ranking of all six models by BIC (n = minutes; primary, as in 20n) and by AIC (secondary, named here in advance):
  fraction of flies in which each model is best.
- Medians and IQR per genotype; R_* reported only for flies with LL_HMM_12 - LL_MC1_12 > 3 nats (denominator not ~0).

## Predictive checks (20 simulations per fly per model, same mask, same statistics as 20n: switches/h, mean Q/W/S
## runs, fraction of 1-min Q and W runs, c_obs; per phase light/dark of the arrival minute)
For MC2_12 and MC1_6 (MC1_12 and HMM_12 from 20n, same seeds-independent statistic). Key statistic: the 1-min run
fractions that neither MC1 nor HMM reproduced (17-52% of fly-phases outside the central 95%).

## Reading fixed now
- MEMORY accounts for the HMM's gain if median R_mem >= 1 in both genotypes, and MC2_12 brings the 1-min run fractions
  inside the simulated 95% in clearly more fly-phases than MC1 and the HMM did.
- DRIFT accounts for it if median R_drift >= 0.5 (at least half the gain disappears with 6-h rates).
- Whatever gain remains (HMM_6 over MC1_6 with R_mem < 1) is CONSISTENT WITH misclassification but not proof of it.
- Limitation stated now: a 2nd-order chain also captures isolated one-minute misclassifications (W Q W patterns), so
  R_mem >= 1 cannot by itself separate short genuine memory from misclassification. What it can show: if MC2 beats the
  HMM by a margin (AIC), the data contain second-order structure that a 3-state noisy Markov chain does not produce.
- Narrowest block: one fly's light or dark (or 6-h block) minutes pooled over its days.
