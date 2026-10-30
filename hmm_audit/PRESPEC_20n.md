# PRESPEC 20n — is the hidden layer needed? Markov chain on S/Q/W symbols vs the 3-state HMM (before running)
Date 2026-10-30. User question: "you have 3 observable states and 3 hidden states, what's the big deal?"
Genotypes: CS (25 flies) and fumin (53 flies), as a start. theta 0.5. Symbols, phases and HMM fits exactly as 20b
(fit2_0.5.jsonl, K = 3; unknown-type wake minutes coded missing, as in 20b).

## Models
- **HMM**: the 20b K = 3 fit (shared emissions B, light/dark transition matrices, 20 starts). Not refitted.
- **MC**: a first-order Markov chain on the symbols themselves = the same model with B fixed to the identity
  (Q -> H1, W -> H2, S -> H3). Fitted by maximum likelihood with the same forward algorithm, so missing minutes
  are marginalised exactly. Iterate em_step from the transition-count estimate until the log-likelihood gain is
  < 1e-10 (relative); B stays the identity under em_step, which is asserted.
- Free parameters: HMM k = 2 + 2*6 + 6 = 20; MC k = 2 + 2*6 = 14. The MC is nested in the HMM (B = I).

## 1. Likelihood comparison, per fly (all 78 flies)
- dLL = LL_HMM - LL_MC (nats), and per 1000 minutes.
- dBIC = BIC_MC - BIC_HMM with n = number of minutes (primary) and n = symbol changes (sensitivity).
- Report median and range per genotype, and the fraction of flies with dBIC > 10.

## 2. Parametric bootstrap LRT (subset; the boundary B = I makes chi-square invalid)
- Flies: 10 per genotype drawn with numpy default_rng(2026).choice(n_flies, 10, replace=False), drawn before
  any result. Cost was measured first (~90 s per 20-start HMM fit): 20 flies x 10 replicates ~ 75 min on 4 cores.
- For each replicate: simulate symbols from the fly's fitted MC (same length, same phases, same missing mask),
  refit MC (as above) and HMM (em_fit, K = 3, 20 starts, same settings as 20b), LR* = 2(LL_HMM - LL_MC).
  The replicate goes through every step the observed statistic went through.
- Per fly: p = (1 + #{LR* >= LR_obs}) / 11. Across flies: the number with LR_obs above all 10 replicates,
  compared with Binomial(20, 1/11).

## 3. Posterior-predictive checks (all 78 flies; 20 simulated records per fly per model, same mask)
Statistics per fly x phase, the same function applied to real and simulated symbols:
  (a) number of symbol changes per hour; (b) mean Q-run and W-run length (min); (c) fraction of Q runs and of W
  runs lasting exactly 1 min (single-minute flicker); (d) observable current c_obs (c_obs_minutes);
  (e) mean S-run length.
For each model and statistic: fraction of fly-phases whose real value lies outside the central 95% of the 20
simulations (expected ~5% if adequate, larger with 20 simulations because of their own spread). The comparison
of interest is MC vs HMM on the same statistic.

## 4. Thermodynamics from the MC (all 78 flies)
Exact stationary J, EPR and affinity per phase from the MC matrices (model_energetics, states mapped as above),
paired against the HMM's (20b). Paired Wilcoxon per genotype x phase; fraction of flies with infinite EPR.

## Reading fixed now
- The hidden layer is REQUIRED by the data if dBIC > 10 in >= 90% of flies in both genotypes and, in the
  bootstrap subset, LR_obs exceeds all replicates in clearly more flies than Binomial(20, 1/11) allows (p < 0.01).
- If the MC passes the predictive checks as well as the HMM, the HMM is a refinement, not a necessity, for those
  statistics, and the paper should say so.
- Expected from the noise argument (not guaranteed, since the MC is a separate fit and not a coarse-graining of the
  HMM): MC J and EPR BELOW the HMM's, because single-minute Q/W misclassifications enter the MC as reversible
  transitions; and MC underpredicts single-minute Q/W flicker relative to longer runs (geometric runs). No
  prediction for the sign of J.
- Narrowest block: one fly's light (ZT 0-12) or dark (ZT 12-24) minutes pooled over its days (shared emissions).
