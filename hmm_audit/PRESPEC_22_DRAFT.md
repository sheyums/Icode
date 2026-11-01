# DROPPED (2026-11-01, Sheyum's decision): follow-up B will not be run. Kept as a record of the draft only.

# PRESPEC 22 — DRAFT, NOT APPROVED: does the bout process need time-irreversible hidden dynamics? (follow-up B)
Date 2026-11-01. Drafted by cClaudeII at Sheyum's request, for review by cClaude and approval by Sheyum. Nothing runs
until it is approved and committed as PRESPEC_22.md. Number 22 to be confirmed free on the shared branch.
Trigger: PRESPEC 21 fixed B as a conditional follow-up, eligible only if 21 detected a shape asymmetry; 21 did (per0 DD,
CS DD). The limits below are cClaude's, from the PRESPEC 21 review, fixed before 21 ran.

## Question
20u and 20v amendment F found that the best time-REVERSIBLE aggregated HMM of each fly's 1-s activity (K <= 6) does not
reproduce the bout process: it makes 22.5% (per0 DD) and 40% (CS DD) too many wake bouts, while getting total sleep
right. Does the same model class WITHOUT the reversibility constraint reproduce it? A pass where the reversible model
failed, with BIC support, says that within K <= 6 the bout process needs time-irreversible hidden dynamics.
It compares likelihood and adequacy only. No hidden EPR, routing or rate is reported as a result: aggregated models are
identified only up to equivalence classes (Fredkin & Rice 1986; the audit's own entries; C4 lesson).

## Data, blocks
- per0 DD (25 flies; 20u's retained traces, first traceLen seconds) PRIMARY: stationary, undriven.
- CS DD (25 flies; whole 4.000-d files) SECONDARY: the circadian drive is not excluded as a source of irreversibility.
Observation: 1 if count > 0 else 0, per second (as 20u).

## Model and fit (identical to 20u except the constraint)
- Ka active-emitting + Ki inactive-emitting hidden states, deterministic emissions, orders (Ka, Ki) in {1,2,3}^2.
- Transition matrix P: a GENERAL row-stochastic matrix. M-step: P_ij = N_ij / sum_j N_ij from the expected counts, with
  no symmetry.
- E-step: standard scaled per-second forward-backward (numba), with deterministic emissions. Validated BEFORE any data
  run against a brute-force forward-backward on a short synthetic sequence (logL to 1e-10, counts to 1e-8). 20u's
  eigenbasis kernel is not reused, because a non-reversible block can have complex eigenvalues.
- Starts: the SAME 10 starts per order as 20u (sqw-planB random_start, seed default_rng(3000 + fly) for per0 DD and
  default_rng(3000 + 100 + fly) for CS DD), drawn in the same order. Stop at relative gain < 1e-9 or 2000 iterations.
- Order selection: BIC with n = the fly's number of runs, k = K(K-1) free parameters (against K(K+1)/2 - 1 for the
  reversible model). The same BIC is reported for the reversible fits already stored (20u; 20v-F CS DD).

## Adequacy gate (identical to 20u stage 1b)
20 simulated traces per fly from its selected unconstrained model, seeds 4000 + fly (per0 DD) and 4000 + 100 + fly
(CS DD), scored by the same Python extractor (Stage 0-validated in 20u). Per-fly ratios simulated/observed of the
wake-bout count and total sleep. PASS if both medians are in [0.85, 1.15]. Death criterion as 20u amendment A.

## Readouts
- (1) Adequacy pass or fail of the unconstrained class, per block.
- (2) Per fly: dBIC = BIC(best unconstrained) - BIC(best reversible); group median and fraction of flies with dBIC < 0.
- (3) Per fly, descriptive: at the reversible fit's selected order, LR = 2 (logL_unc - logL_rev), reported with
  df = (K-1)(K-2)/2. No chi-square p is used for a decision.
- (4) Descriptive: selected orders, and the per-fly wake-bout ratio for both classes side by side.

## Decision rule and wording (fixed now)
Per block:
- Unconstrained PASSES the gate AND dBIC < 0 in more than half of the flies: "Within aggregated hidden Markov models
  of up to 6 hidden states, the observed bout process requires time-irreversible hidden dynamics: the best reversible
  model fails to reproduce the wake-bout count (ratio ...), the best unconstrained model reproduces it (ratio ...) and
  is preferred by BIC in x of n flies." per0 DD: stationary and undriven. CS DD adds "the circadian drive is not
  excluded as its source".
- Unconstrained PASSES but dBIC < 0 in at most half of the flies: "Adequacy is gained only with extra parameters that
  BIC does not support; no attribution to irreversibility." (INCONCLUSIVE)
- Unconstrained FAILS: "No model in either class (K <= 6) reproduces the bout process." There is no further
  enlargement of K and no other class without a new PRESPEC.
Never: "proves broken detailed balance", "the hidden EPR is ...", or any hidden rate or route as a finding.

## Time box and run logistics
6 h of fitting per block on 4 cores, or "not completed"; orders and starts are not reduced after the fact. Per-fly
checkpoints; a job stopped by the harness's 2-h limit is resumed from its checkpoints (as in 20u and 20v-F) and the
stops are recorded. One run, one verdict per block. Ledger entry "PRESPEC 22" in RESULTS_LEDGER_planB.md.

## Open points for the reviewers
1. Should a pass need BOTH blocks, or is per0 DD alone enough for the headline? Proposed: per0 DD alone carries it;
   CS DD is secondary.
2. Is "dBIC < 0 in more than half of the flies" the right BIC criterion, or should it be the group median dBIC < 0?
   (They usually coincide.)
3. Seeds reuse 20u's starts so that the two classes start from the same points. A reversible start is a feasible
   start for the unconstrained EM, so the unconstrained logL at a matched order can only be >= the reversible one if
   both converge. Readout (3) checks that.
