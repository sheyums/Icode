# beam_hmm — status of every finding so far (2026-11-19)

Every result in this project, labelled and grouped by genotype (and condition). Numbers and commit hashes are in
`RESULTS_LEDGER.md`; this page only classifies them. When this page and the ledger disagree, the ledger is right.

**Labels**
- **PRELIMINARY** — a pre-specified test (PRESPEC committed before the run) on the pilot. Supported or not, it is a
  hypothesis-generating result: one run, ~11 flies, thresholds partly set after viewing week 1 (P1 disclosure). It
  becomes a finding only when the repeat experiment confirms it under the same rules.
- **DESCRIPTIVE** — computed after looking, not pre-specified as a test, or reported "descriptive only" by its own
  PRESPEC. No inferential claim.
- **METHOD** — synthetic gates, tool validations and data-quality facts. Not about fly biology.

**Genotype is never compared across rows here.** The CS pilot and the per0 data differ in genotype, light condition
(LD vs DD), cohort, lab session and analysis window at once; any CS-vs-per0 difference below is confounded in all of
these and is labelled DESCRIPTIVE.

---

## 1. Canton-S (CS), males, LD 12:12 — the 2026 pilot (one run, channels 17-32; weeks 2-3 = test set)

| Finding | Label | Result | Ledger |
|---|---|---|---|
| A1: by day, flies rest at the tube ends more when asleep than in short pauses | PRELIMINARY | supported (11/12, Holm p 0.0067); not in the dark | 2026-11-14 A1 |
| A1 on week 1 | DESCRIPTIVE | same pattern (viewed before A1 was written) | 2026-11-14 |
| A2: pauses become more frequent before sleep | PRELIMINARY | not supported (week 1's strong pattern did not replicate) | 2026-11-14 A2-A4 |
| A3: edge flicker differs day vs night (behaviour vs detector) | PRELIMINARY | not supported (consistent with a detector effect) | 2026-11-14 A2-A4 |
| A3 secondary: 5-min sleep measure robust to flicker (1-5 min/day) | DESCRIPTIVE | | 2026-11-14 A2-A4 |
| A4: slower walking at night | PRELIMINARY | supported (0/10 positive, Holm p 0.015); path length not supported | 2026-11-14 A2-A4 |
| A5: sleep, activity, end-resting, walking speed change with age (weeks 2-5) | PRELIMINARY | none survives Holm (3 nominal, same direction) | 2026-11-14 A5 |
| A5 weekly means, death-aligned changes | DESCRIPTIVE | survivor-biased | 2026-11-14 A5 |
| A2b: activity shape around sleep boundaries asymmetric (primary 300-s setting) | PRELIMINARY | "no shape asymmetry" | 2026-11-14 A2b |
| A2b secondary 120-s setting: asymmetric, slow wake-up | DESCRIPTIVE (secondary) | direction opposite to per0 DD | 2026-11-14 A2b |
| A2b robustness (lab conversion, no switch exclusion) | DESCRIPTIVE (post hoc) | same direction | 2026-11-14 |
| T1-C1: shared time-of-day modulation of sleep onsets / wake-ups beyond the light step | PRELIMINARY | supported (11 flies, rotation p 0.001); LD: light and clock not separable | 2026-11-14 T1-C1 |
| T1-C1 shapes (siesta onset peak ZT ~5, wake-up peak ZT ~11) | DESCRIPTIVE | | 2026-11-14 T1-C1 |
| O1-B: T1-C1 per week 1-6 (night-onset peak flattens weeks 4-5) | DESCRIPTIVE | survivor bias, n varies | 2026-11-19 O1-B |
| O1-A (a): wake-bout hazard rises with time awake | PRELIMINARY | not supported (RiseNullP 0.84 / 0.97); hazard falls | 2026-11-19 O1-A |
| O1-A: no model in the closed library fits the POOLED wake bouts | PRELIMINARY (pre-specified fit) | all rejected; excess 2-s bouts | 2026-11-19 O1-A |
| Per-fly median bout lengths (4-fold spread between flies; light vs dark within flies) | DESCRIPTIVE (post hoc) | | 2026-11-19 correction |
| O1-A (b): POOLED sleep bouts = three time scales (~2 min, ~22 min, ~65 min) | DESCRIPTIVE | pooled over flies and phases: may be heterogeneity (O1-C tests within flies) | 2026-11-19 O1-A + correction |
| O1-A secondary: light vs dark bout hazards | DESCRIPTIVE | censoring caveat | 2026-11-19 O1-A |
| O1-C: more than one time scale WITHIN a fly and light phase | PRELIMINARY | supported for wake (20/20 units) and sleep (14/21: light 11/11, dark only 3/10) | 2026-11-19 O1-C |
| O1-C: dark-phase sleep mostly one time scale per fly (fly-specific rate); wake short time scale varies seconds-minutes across flies; wake often needs 3 (6/9 large units) | DESCRIPTIVE | 30/41 p at the B = 499 floor; 2 wake units with a 1-s component at the rate ceiling | 2026-11-19 O1-C |
| Q1: Q at the start of a wake bout predicts its remaining duration | PRELIMINARY | not supported (p 0.13; L = 600 s: 0.18) | 2026-11-19 Q1 |
| Q1 direction (Q-start bouts shorter in 9/11 flies) | DESCRIPTIVE | | 2026-11-19 Q1 |
| Sleep 726-1229 min/day; E0 / P1 week-1 exploration | DESCRIPTIVE | | 2026-11-12, 2026-11-14 |
| Position-0 occupancy per fly-week; long stays at 0 before death (2 flies) | DESCRIPTIVE | orientation unknown in the pilot | 2026-11-16 |

**Data-quality facts from this run (METHOD, not biology):** jump artifacts concentrate on a hub beam; P1-C dead-hub
repair; 14 fly-weeks of live flies hidden behind a blocked beam 16; the lab conversion loses up to 261 min/day of sleep
where hubs exist; lab sleep scoring = beam_hmm scoring second for second; LoadingCheck catches 14/14 hidden flies and
~half of the blocked beams in 2 h (1 h the same); clock drift between computer and incubator; mixed tube orientation.

## 2. per0, DD — the earlier dataset (25 flies; analysed in this repo's FINDINGS.md, not by beam_hmm)

| Finding | Label | Source |
|---|---|---|
| Sleep bouts: hyperexponential K = 3 (LRT p < 0.001; third component ~2% of bouts) | DESCRIPTIVE (single dataset; not re-tested on held-out flies) | FINDINGS.md |
| Wake-bout hazard falls, rises (peak ~1000 s), falls; hyper-Erlang with a series branch | DESCRIPTIVE (family chosen after the library failed; fit p optimistic, FINDINGS.md says so) | FINDINGS.md |
| Wind-down into sleep / abrupt waking (planB_21, earlier project) | as reported there; beam_hmm only reproduced the computation exactly | ledger 2026-11-14 repro |

## 3. CS, DD — earlier project

Mentioned in the earlier project's description of the wind-down profile; NOT analysed by beam_hmm. No result here.

## 4. Across genotypes (CS LD pilot vs per0 DD) — DESCRIPTIVE, confounded

- Activity around sleep boundaries: per0 DD winds down into sleep and wakes abruptly; CS LD (secondary setting) wakes
  slowly. Opposite direction.
- Wake-bout hazard: per0 DD rises after ~2 min; CS LD falls throughout.
- Sleep bouts: three time scales in both; CS LD has larger short and long shares.
Genotype, LD/DD, cohort and window all differ. The repeat experiment's DD segment (and an optional per0 group run
alongside) is what can separate them (DESIGN_NOTE.md).

## 5. Synthetic gates and tools (METHOD; no genotype)

| Gate / check | Result |
|---|---|
| R1 hub repair (P1-B / P1-C) | FAIL / PASS |
| D1-B, D1-D state discovery (10-s windows / sample level) | FAIL, FAIL |
| T1 per-fly clock (IOHMM) | FAIL (power) |
| T1-B cohort clock test, hidden states, K fixed | PASS |
| T1-C / T1-C1 cohort clock test, observed level | FAIL / PASS |
| B1 direction-aware AR-HMM state discovery | FAIL (G1 5/20, G2 spurious 19/20) |
| Q1 random-label calibration | PASS (KS p 0.66) |
| MATLAB tools (LoadingCheck, RawToActivity1s, BeamFingerprint) vs Python twins | identical on all checks run |
