# PRESPEC T1 — DRAFT (2026-11-12) — SUPERSEDED: adopted into PRESPEC.md as entry T1 (with decisions)

**Status: DRAFT for the user's review. Not a pre-specification until it is moved into `PRESPEC.md` as a
dated entry and committed. Nothing in it has been run.** Open choices are marked **[DECIDE]**.

## Question
Does the 10-s behavioural HMM (D1-B structure: 4 activity symbols, K states, emissions fixed in time) need
transition rates that vary with clock time beyond the light/dark step? If so, in what form? And do the
D1 candidate states survive once clock time is modelled (a state that disappears was rate modulation, not a
state)?

## Data
- Pilot, channels 17-32, PRIMARY set (P1-D; week 1: 12 flies), SENSITIVITY set (8 never-artifact flies)
  reported beside it.
- Weeks 1-2 (as D1).
- Same preprocessing as D1-B: spike filter, P1-C hub repair, handling/death exclusions.
- Preconditions:
  - the D1-B gate result is known;
  - handling periods are settled;
  - week 2 is available;
  - the T1 synthetic gate below has PASSED.

## Time window
Fly x weeks 1-2, fly-days ZT0-ZT0. Each week is also fitted separately, as an age check: descriptive only,
no age model in T1.

## Model: input-output HMM
Emissions: categorical over the 4 symbols, shared across all times (a state means the same thing at every
hour; PRESPEC 14 showed what drifting emissions do to state identity).

Transitions at window t, from state i:

    P(s_t = j | s_{t-1} = i, z_t) = exp(b_ij . z_t) / sum_l exp(b_il . z_t),   with b_ii = 0 (self-transition
    is the reference)

The covariates `z_t` are a function of ZT only, so the model is periodic with period 24 h. The four nested
models:

| model | z_t | parameters per off-diagonal transition |
|---|---|---|
| M0 | intercept, light/dark step (= D1-B) | 2 |
| M1 | M0 + 1st harmonic: cos, sin(2 pi ZT/24) | 4 |
| M2 | M1 + 2nd harmonic: cos, sin(4 pi ZT/24) | 6 |
| M3 | M2 + 3rd and 4th harmonics | 10 |

- **Why these terms:**
  - The step keeps the lights-on and lights-off changes sharp.
  - Harmonics 1-2 give the bimodal morning/evening shape with a midday siesta.
  - Harmonic 4 (6-h period) is the modulation planted in D1's G2 control. It is included so that T1 can
    say whether that kind of modulation is absorbed.
- **[DECIDE]** Whether to add lights-on and lights-off startle indicators (the first 30 min after each
  switch) as M4.

## Estimation
- Generalized EM, with the E-step as in `hmm_cat`.
- **M-step for `b`:**
  - The covariates depend on ZT only, so the expected transition counts are aggregated by
    window-of-day: 8640 x K x K sufficient statistics.
  - For each origin state, a weighted multinomial logistic regression is fitted by L-BFGS.
  - It carries a ridge penalty `lambda = 1e-2` on all `b` except the intercepts **[DECIDE value]**. This
    keeps rare transitions from diverging.
- 10 restarts. Each starts from the M0 fit plus a small random perturbation of the extra coefficients
  (seeds `SeedSequence([5001 + r, fly, K, model, fold])`).
- Maximum 500 iterations, tolerance 1e-4 nats. Best training penalized logL kept. All restart logLs reported.

## Model selection (per fly; flies are the unit)
- **K:** fixed at the fly's D1-B K* for the clock comparison.
- **Clock model:** M* = argmax of held-out logL per window, two-fold by day (odd/even). BIC on all days is
  reported beside it.
- **Across flies:** "clock dependence beyond the light/dark step" is supported if held-out logL of the best
  of M1-M3 exceeds M0 in a one-sided exact sign test across flies at p < 0.05. That needs at least 10 of 12
  flies in the primary set and 7 of 8 in the sensitivity set. Both sets are reported.
- **Description of the clock effect:** for each transition, the fitted rate profile across ZT with a fly
  bootstrap over days (20 resamples) to give pointwise bands. Phase of peak and peak-to-trough ratio per
  fly, summarised across flies. Descriptive only.

## Proposed new D1 criterion (to be entered as D1-C if this draft is adopted)
**C6 survives the clock model.**
- Repeat D1's K selection with M* transitions.
- The candidate must still be present at the new K* (same TYPE, Jaccard >= 0.5 with its windows).
- The candidate must also pass C1 with the bootstrap null simulated from the fly's own M* model at K* - 1.
- A candidate that fails C6 is reported as **"clock-explained"**.

## Synthetic gate (must PASS before any data run; seeds 6001+, 20 synthetic flies per scenario, 14 days)
The truth is built in the 0.4-s simulator (`sim_positions.py`) and passed through the full pipeline,
including the spike filter.

| scenario | truth | PASS iff |
|---|---|---|
| T-G1 recovery | 3 states (SLEEP, QUIET, TRAVEL as D1-A); log sleep-exit rate = step + 1st + 2nd harmonics, peak-to-trough ratio 3, evening peak ZT 12.5 | M1-M3 chosen over M0 in >= 16/20; fitted sleep-exit peak within +-2 h of truth in >= 16/20; median fitted/true peak-to-trough ratio in [0.67, 1.5] |
| T-G2 false positive | the same 3 states with the light/dark step only | M1-M3 chosen over M0 in <= 3/20 |
| T-G3 spurious states | 3 states, per-fly heterogeneity (as D1 G2), clock modulation as T-G1 | under M*, a spurious candidate passes C1-C3 in <= 1/20; the same count under M0 is reported beside it |

If any scenario fails, the failure is logged and the model is changed in a new dated entry.

## Decision rule
- Clock dependence beyond the step: the across-fly sign test above, in each set.
- Candidate states: promoted only if they also pass C6. Promotion still only creates a hypothesis for the
  repeat experiment.

## Matched null
- **Clock comparison:** held-out likelihood. The null model's extra coefficients are zero, so no simulation
  is needed for the selection. T-G2 measures the false-selection rate at our sample size.
- **C6:** the fly's own M* fit at K* - 1, simulated with the same length, missing pattern and covariate
  sequence.

## What T1 does NOT do (later entries)
- **Age (T2):** a covariate or per-week fits. Survivor bias needs alignment by age AND by days before death,
  so T2 waits for the weeks 2-9 data and the death times.
- **Homeostatic covariates (T3):** time awake so far, sleep in the preceding 12 h. These are
  history-dependent, so the model is no longer a fixed-input IOHMM; it needs its own gate, because the
  covariate is computed from the same sequence the model explains.
- **Clock versus light:** this pilot is LD only, so "clock dependence" here means light plus clock. Separating
  the two needs a DD segment, which the repeat experiment should include.

## Compute estimate
- Per fly: 4 models x 2 folds x 10 restarts at K*, plus the C6 K sweep under M*. This is about 2x D1's cost.
- Gate: 60 synthetic flies, about 8-10 h on 4 cores.
- Resumable, as in the D1 gate.
