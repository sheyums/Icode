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

## Amendment C (2026-10-31, before any data run; found by the pre-run sanity check)
1. FORMULA: the bound is B = J^2 / D, not 2 J^2 / D. Barato-Seifert: eps^2 Sigma >= 2 with eps^2 = Var/<X>^2,
   Var = 2 D t, <X> = J t, Sigma = sigma t  =>  sigma >= J^2 / D. Every B above and below is J^2 / D.
2. DISCRETE TIME: the TUR is proven for continuous-time jump processes. Measured with the corrected formula: over 1000
   random 3-state chains, B/sigma <= 1.000 for the generator P - I, but up to 1.074 for the 1-min discrete chain. On
   the flies' fitted observed chains (42 fly-phases, every 7th fly) the discrete/continuous inflation of B is median
   1.0008, max 1.0051. Rule added: per fly, kappa = B_MC(discrete) / B_MC(generator P - I), exact; S2 uses
   B_c / kappa (deflated) against sigma_MC. B_MC, B_HMM reported in both versions.
3. DISCLOSED: that check computed the exact Q_MC = B_MC / sigma_MC (an S3 quantity) for those 42 fly-phases: median 0.98,
   max 1.00. The stated expectation "Q << 1" is withdrawn: near equilibrium a single cycle nearly saturates the TUR, so
   the observed chain sits at the bound. No data estimate (B_obs) had been computed. Tests and decision rules unchanged.

## Amendment D (2026-10-31, AFTER the run; the report crashed before printing any result)
Fly-phases whose fitted observed chain has no cycle (sigma_MC <= 1e-9 or its continuous-time bound B_MC_ct <= 1e-12:
some transition absent in both directions, so J = sigma = 0) have no defined TUR ratio and are EXCLUDED and counted.
Fly-phases where B_obs cannot be estimated (fewer than 4 blocks, or D <= 0) are excluded and counted. Mechanical:
undefined quantities, decided before any ratio was seen (the report crashed on a 0/0 in kappa).

## Amendment E (2026-10-31, AFTER the results; user request "use rule of 3 for 0 cycles"; reported BESIDE amendment D)
The fly-phases excluded under D are of three kinds (classified from the phase's directed transition counts):
(i) one-way edge(s), cycle present, sigma_MC infinite; (ii) all edges two-way but counts exactly symmetric, so the net
cycle count is exactly 0 and the fitted chain is in detailed balance (sigma = B = 0: TUR trivially satisfied, no ratio);
(iii) no cycle: an edge absent in both directions, or a state absent in that phase.
Rule of three is applied to (i) only, as in sqw_ref.epr_count_bound: an unseen direction gets b* = min(3, k) counts,
k = the observed reverse count; P_r3 = row-normalised counts of that phase; sigma_MC, B_MC, B_MC_ct from P_r3; 200 new
null records simulated with P_r3 in that phase (the fly's MC1 matrix in its other phase). (ii) and (iii) stay excluded and
are reported by kind. Direction: b* gives a 95% LOWER bound on sigma, so it inflates B/sigma (toward S2 'violation';
an upper bound on efficiency). Primary DD blocks contain no (i) fly-phases, so the primary verdicts cannot change.
The ledger label "observed chain without a cycle" for all D exclusions was wrong; corrected in the ledger.
