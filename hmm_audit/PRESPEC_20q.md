# PRESPEC 20q — is the sleep-wake cycle more precise than its dissipation allows? (TUR; before running)
Date 2026-10-31. Physics: for any Markov jump process in steady state, every current obeys
sigma >= 2 J^2 / D (thermodynamic uncertainty relation; Barato & Seifert 2015 PRL 114:158101; used on experimental
fluctuations by Li, Horowitz, Gingrich & Fakhri 2019 Nat Commun 10:1666). An observed current whose precision implies
a bound B = 2J^2/D ABOVE the observed chain's entropy production would require dissipation the observed chain does not
have. If not, B/sigma is the cycle's thermodynamic efficiency of precision (<= 1).
Per the user's point (2026-10-31): EVERY observable-vs-hidden comparison is made per fly x phase, and group statements
come only from paired tests on per-fly quantities.

## Data, blocks
- PRIMARY: per0DD and CS DD (steady state assumed; whole record, one matrix, as 20b/20p). CS DD has a free-running
  rhythm: it inflates D and makes B smaller, so it biases AGAINST finding B > sigma (conservative) and makes the
  efficiency an underestimate. Stated now; no split (3-h blocks within 12-h halves would leave too few blocks).
- DESCRIPTIVE ONLY: LD genotypes (CS, per0LD, fumin, sss), per light/dark phase. Periodic driving can legitimately
  break the steady-state TUR (Barato & Seifert 2018 PRX 8:031045), so no TUR verdict is drawn in LD.
- Symbols, theta 0.5, missing coding, Q/W-reference exclusion (20p amendment B) as before; degenerate-H1 flies excluded
  from HMM quantities only.

## The current and its estimate (per fly x phase)
- Increments on transitions between consecutive non-missing minutes, attributed to the arrival minute's phase:
  +1 for W->Q, Q->S, S->W (the forward cycle A->Q->S->A), -1 for the reverses, 0 otherwise. N(t) = cumulative sum.
- J = total increments / analysed hours. D(tau) = sample variance of the increments summed over consecutive,
  non-overlapping blocks of tau minutes lying in one phase with no missing minute, divided by 2 tau (hours).
  tau = 3 h PRIMARY; 2 h and 6 h reported as sensitivity. B_obs = 2 J^2 / D, in k_B per hour.
## Exact model quantities (per fly x phase, from that phase's matrix at its stationary state)
- sigma_MC: observed first-order chain (20n fit), exact. B_MC: exact TUR bound of the SAME current under that chain,
  from the tilted matrix M(s)_ij = P_ij exp(s d_ij): J = lambda'(0), D = lambda''(0)/2, lambda = log Perron root
  (central differences, h = 1e-3), per minute -> per hour.
- sigma_HMM (20b fit, exact; CS DD fit from 20p) and B_HMM: the same current on the observed symbols of the HMM, via
  the 9-state chain on (hidden, symbol) pairs with transitions P_zz' B_z'y'.
- Sanity checks before the run: Q = B/sigma <= 1 for 1000 random 3-state chains (TUR must hold), and B_obs -> B_MC on a
  long simulated record (T = 1e6 min).
## Nulls (same estimator, same mask, same blocks)
- 200 records per fly from its MC1; 100 from its HMM.
## Statistics and reading fixed now (Holm across the 2 primary DD blocks)
- S1 PRECISION BEYOND THE OBSERVED CHAIN: per fly p = (1 + #{B* >= B_obs}) / 201; group Wilcoxon of
  log(B_obs / median B*) against 0. Supported in a block if Holm p < 0.05 with median > 0 AND the fraction of flies with
  p < 0.05 exceeds Binomial(n, 0.05) (p < 0.05). Meaning: the observed current is more regular than the observed chain
  produces.
- S2 TUR DEMANDS HIDDEN DISSIPATION: bias-corrected bound B_c = B_obs * B_MC / median(B*) (the chain calibrates the
  estimator's bias; uncorrected B_obs also reported). Per fly log(B_c / sigma_MC); group Wilcoxon against 0. Supported in
  a block if Holm p < 0.05 with median > 0. Also: fraction of flies with B_c > sigma_MC, and with sigma_HMM >= B_c
  (infinite sigma_HMM counts as >=).
- S3 EFFICIENCY (descriptive): per fly Q_MC = B_MC / sigma_MC (exact), Q_obs = B_c / sigma_MC, Q_HMM = B_c / sigma_HMM;
  medians and IQR.
- Expected, stated now: weakly driven cycles usually sit far below the TUR bound (Q << 1), so S2 is unlikely; S3 is then
  the result. Any outcome is reported.
- Narrowest block: one fly's whole DD record (LD: one fly's light or dark minutes pooled over days).
