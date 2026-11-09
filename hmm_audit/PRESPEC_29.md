# PRESPEC 29 — do sss and Sh^mns contain distinct forward- and backward-running flies?

Date 2026-11-07. This is an exception to the freeze, approved by Sheyum. Written and committed BEFORE any per-day or
test-retest quantity has been computed. Already known: whole-record per-fly J and its group tests (C2, PRESPEC 23), and
that one fly's direction is below its own noise for days (PRESPEC 25). Motivation: the per-fly J of sss and Sh^mns
looks like the old per0 "two groups" (a Fisher test on signs, later traced to H1 degeneracy). Splitting by sign would
manufacture groups from a weak common bias plus noise, so it is NOT used.

## Quantities
- Observed J per day d and phase p: (N_SW - N_WS) / 12 h over consecutive scored minute pairs arriving in that
  day-phase. Run-compressed symbols are not needed: an S -> W change is one pair.
- Day split (fixed): A = days 0 and 2, B = days 1 and 3 (whatever exists; 3-day flies A = {0, 2}, B = {1}; 2-day
  flies A = {0}, B = {1}). J_A and J_B are the means over each half's days.
- Hidden J per day-phase: the same count on each of 200 FFBS posterior paths from the fly's K = 3 labelled HMM,
  averaged over paths. Existing fits (fit2_0.5) for CS, per0 LD, per0 DD, fumin and sss; degenerate fits excluded;
  no fit exists for CS DD. Sh^mns is fitted now with the 20b procedure (EM, 20 starts, phase-specific transitions,
  shared emissions, label3b), on its PRESPEC 23 records (death rule, up to 4 days).

## Tests
T1 TEST-RETEST: Spearman rho(J_A, J_B) across flies with a scored phase in both halves; values rounded to 10 decimals
   before ranking. Null: 10,000 permutations of the A-B pairing between flies. One-sided (rho > 0 = flies keep their
   direction).
T2 EXCESS SPREAD (observed level only): the SD across flies of the whole-record per-phase count J, against 500
   simulated cohorts. Each fly is simulated with its own length and missing mask from the GENOTYPE-POOLED phase chain
   (counts pooled over flies, row-normalised). One-sided p = (1 + #{SD_sim >= SD_obs}) / 501.
PRIMARY (Holm over 4): T1 and T2 for sss light and Sh^mns light.
REFERENCE (uncorrected p, descriptive): T1/T2 for CS L/D, CS DD, per0 LD L/D, per0 DD, fumin L/D, sss D, Sh^mns D.
HIDDEN (secondary, Holm over 2): T1 on hidden J for sss light and Sh^mns light; other blocks descriptive.

## Decision rule (fixed now)
- "Distinct fly-to-fly differences in daytime direction" in a genotype: T1 AND T2 Holm < 0.05. Then, descriptively,
  a 1- vs 2-component Gaussian mixture (BIC) on the per-fly whole-record J asks whether there are two groups or a
  continuum.
- Otherwise: "the spread of directions is what one shared process gives at this record length". For sss, with 2-day
  records, a null is "not resolved", not "no subgroups".
- No other day split, phase, threshold or genotype is added after looking.

## Output
`sqw29.py` -> `sqw_out/s29_hetero.json`, `figures/fig29_retest.png`. Ledger entry with commit hashes.

## Amendment A (2026-11-09; written and committed BEFORE running): corrected hidden-level test-retest
The original hidden T1 is not interpretable: each fly's HMM was fitted to its whole record, so both halves' posterior
paths share one fitted transition matrix (shared estimation). Corrected design:
- For each fly, two masked records: y_A keeps days {0, 2} and sets the other days' minutes to missing (-1); y_B keeps
  days {1, 3}. EM marginalises missing minutes, so no transition is invented across the gaps.
- Each half is fitted separately with the 20b procedure (K = 3, 20 starts, phase-specific transitions, shared
  emissions, label3b). 200 FFBS paths are drawn from each half's own fit and record; J per day-phase is counted on the
  half's own days and averaged over paths. The halves share nothing but the fly.
- A fly enters the hidden T1 for a block only if neither half's fit is degenerate and both halves have the phase.
- T1 as before (Spearman, values rounded to 10 decimals, 10,000 permutations, one-sided). PRIMARY hidden: sss light
  and Sh^mns light, Holm over 2. Reference: CS L/D, per0 LD L/D, per0 DD, fumin L/D, sss D, Sh^mns D (descriptive).
- sss halves are single days, so its fits may be unstable. The number of degenerate or failed fits is reported per
  genotype.
- The decision rule of the observed level is unchanged. The hidden level is a check: a positive hidden T1 with null
  observed T1 and T2 is reported as "possible hidden-level consistency", not as subgroups.
