# Design review — sleep/wake routing through H1 (2026-10-21)

Purpose: take stock before choosing the analysis design, instead of running more tests.
Every number below is in `RESULTS_LEDGER.md` with its plan (`PRESPEC_*.md`) and code.

## 1. The scientific question

Do flies move between sleep and wake through an intermediate hidden state (H1), and is that
routing asymmetric between sleep→wake and wake→sleep? In the paper's terms:
Xn = (N312 − N32)/(N312 + N32) and Yn = (N213 − N23)/(N213 + N23), and the thermodynamic
quantities built on the same transition matrix (cycle flux, EPR, cycle affinity).

## 2. What is established

### 2a. The scored 3-state HMM (the paper's current model)

- **Emissions are deterministic (0/1).** H1 emits exactly like H2 in 68% (CS), 57% (per0), 99%
  (fumin), 99% (sss) of windows, and like H3 in most of the rest.
- **In those windows, routing is not identified.** Transition matrices with identical likelihood
  (checked to 6 decimals on real sequences) span Xn and Yn over about [−1, +1]. The undetermined
  share is 99–100% of windows in every genotype (PRESPEC 1, 2).
  - The Xn − Yn asymmetry is not identified either: its sign is fixed in 0% of windows.
  - Cycle flux, EPR and cycle affinity are not identified either. An equally good fit has P(2→3) ≈ 10⁻²⁰.
- **Viterbi adds a large bias on top.** It suppresses indirect routes, e.g. fumin N312 = 155
  (Viterbi) vs 1312 (posterior). FFBS removes the bias but not the non-identifiability.
- **The windows were extended to finish a bout,** which length-selects the last bout of each
  window. This produced apparent irreversibility; it was removed by dropping edge bouts (8b).

### 2b. Scored bout durations (model-free and pooled)

- **Wake has a robust brief component** (~5–17 s, 4–8% of wake bouts) in CS, per0 and sss.
  It exceeds between-window heterogeneity (p ≤ 0.006). fumin is the exception (0.6%).
- **Sleep: a "second kind" of bout within one fly at one time of day only in CS.**
  - Fisher p = 0.005 at 4 h blocks, 0.013 at 12 h blocks, with edge bouts removed.
  - Its mean is ~320 s, i.e. ~20 s above the 300 s threshold. This is consistent with brief rest
    crossing the 5-minute criterion.
  - per0 and fumin lean towards a rising hazard instead.
- **Irreversibility of consecutive bout pairs:** after removing window-edge bias, only fumin
  passes (whole record, dark), with a tiny lower bound (≤ 0.012 nats/cycle). Post hoc.

### 2c. Per-second movement counts (independent of scoring)

- **CS quiets down in the minutes before sleep onset,** in light (p = 0.001) and dark
  (p = 3 × 10⁻⁵). Model-free.
- **Inactivity durations of all lengths form a ≥5-component mixture.**
  - In CS light, 18% of 5-min+ inactivity belongs to a ~99 s rest component.
  - In fumin, all of it does.
- **A flat count HMM needs 4 states, not 3** (held-out days, all bins).
  - 3 states predicts no better than 2 in CS.
  - The 4 states form 2 slow groups (sleep-like, wake-like) × fast movement states.
  - No third long-lived (≥ 60 s) group appears in CS or fumin.
- **The 3-state count HMM's H1 is a brief pause (CS, nearly identical to H3) or a low-activity
  state (fumin).** Neither is a sleep↔wake regulator.

### 2d. Hierarchical (slow × fast) count HMM, synthetic validation

- **Per fly-bin, routing is not estimable even with the true model:** replicate SD 0.6–0.7.
- **Pooling 25 flies × 4 days recovers routing precisely,** with SD 0.06 (0.02–0.07 with a day
  trend) and every fit at the truth's likelihood.
- **Fly-to-fly heterogeneity creates a FALSE transitional state** under full pooling (V2, V6 fail).
- **Behaviour drifts across days.**
  - CS: sleep bouts shorten about 1 min per day.
  - fumin: large day-1 effect.

## 3. Design lessons

1. **Identifiability needs distinct emissions.** A state that looks like its neighbours cannot be
   routed through, whatever the optimiser does. This is the main obstacle.
2. **Data per fit:** slow-state routing needs many slow transitions. One fly's 4 h × 4 days has about 100.
3. **Pooling buys precision but imports heterogeneity,** which mimics extra states.
4. **Non-stationarity mimics slow states:** time of day within 12 h, and drift across days.
   Keep blocks short, or model the drift.
5. **Cut windows at fixed clock times, never at bout boundaries.**
6. **Every model needs a negative control** (data without the effect). Two of ours failed it and
   would otherwise have produced a spurious H1.

## 4. Candidate designs

| | Design | Answers routing? | Main risk | Cost |
|---|---|---|---|---|
| A | Partial-pooling hierarchical count HMM (built; shared routing, FREE fly-specific texture, day factors) | not as built: its unit test FAILED (wrong routing; T and W lost their distinct signature) | free per-fly emissions erase the shared signature that defines T | would need redesign → B |
| B | Fully Bayesian random-effects HMM (custom Gibbs/FFBS) | yes, with honest per-fly uncertainty | development time; mixing | weeks |
| C | Richer observables: per-beam position (multibeam) or video | yes, and per fly, if H1 gets a distinct signature | data availability | depends on data |
| D | Sleep-deprivation recovery (H1 emits differently) as an identifying contrast | in those windows | transfer to baseline needs a shared-parameter assumption | moderate |
| E | Model-free transition signatures on counts: aligned activity, conditional bout statistics, irreversibility tests without windows | "is there a transitional regime, and is it asymmetric", without routing parameters | not a route count | hours |
| F | Reframe the paper's quantities: identified sets for Xn and Yn; EPR lower bounds from observables | reports what the data determine | weaker-looking numbers | days |

## 5. Recommendation (for discussion)

1. **Primary evidence:** E (model-free, on counts), extended with a windowless irreversibility test.
   It is cheap, robust, and directly about "is there an asymmetric transitional regime".
2. **Model-based estimate:** A failed its unit test as built. The workable version is B-like: per-fly
   deviations SHRUNK toward shared values (random effects), validated with V1 + V2-type tests.
3. **If per-beam or video data exist:** C is the best lever. It turns routing from inferred to
   nearly observed.
4. **For the paper:** the current 3-state routing and energetics numbers need F (identified sets
   or bounds) or replacement by A/C-based estimates. Mapping the 3-state model onto the count
   model is feasible once A is validated: compare co-occupancy of the scored H1/H2/H3 with the
   count model's coarse S/T/W on the same flies.
