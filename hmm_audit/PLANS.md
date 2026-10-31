# Plan A and Plan B (written 2026-10-31)

Written at the end of the S/Q/W audit (PRESPEC 19-20q). Every result referred to here is in
`RESULTS_LEDGER.md` under the PRESPEC number given. This file is the brief for any session that picks the
work up; it does not replace the ledger.

## Where we stand, in five lines
1. Wake splits into quiet (Q) and active (W) wake; with sleep (S) the fly cycles W -> Q -> S -> W with a
   direction, in most genotypes and phases. These are observations at the S/Q/W level.
2. The 3-state HMM's hidden states are almost the observed symbols (emissions 80-100% diagonal). Its gain over
   a plain Markov chain on S/Q/W is real but small and short-range (20n, 20o), and BIC prefers the plain chain
   in ~75% of flies.
3. Model-free probes of what lies beneath: extra time-irreversibility beyond the transitions only in per0 DD
   (and fumin dark, with a caveat) (20p); no thermodynamic-uncertainty evidence for hidden dissipation anywhere
   (20q). The HMM's larger EPR is therefore an upper, model-based bracket, not the measured value.
4. Slow modulation of the rates (6-h drift, circadian) is real and separate from the HMM's gain (20o), and it
   shows as extra noise in the cycle current in CS DD and LD but not per0 DD (20q, post hoc observation).
5. Therefore: paper 1 is observable-first (Plan A); the regulators live in the slow modulation of the rates,
   not in fast hidden states (Plan B).

---------------------------------------------------------------------------------------------------------------

## PLAN A — paper 1: a directed sleep-wake cycle and its thermodynamics, measured at the observable level

Arc: (1) the observation (wake splits) -> (2) the cycle and its direction -> (3) physical quantities (J,
affinity, EPR, Onsager L, landscape) -> (4) mutants dissect it -> (5) memory and hazards -> (6) what lies
beneath the observed cycle (hidden layer thin; per0 DD anomaly; TUR saturation).

### Frozen claim list (each claim, the level it rests on, the analysis that supports it, status)
| # | Claim | Level | Support (ledger) | Status |
|---|---|---|---|---|
| C1 | Wake contains two states, quiet (Q) and active (W), defined per fly by theta x r_f | observed | 20b theta sensitivity; 20m (30 s) | done |
| C2 | Sleep-wake transitions circulate W->Q->S->W: forward in CS night, per0 (LD dark, DD), fumin, CS DD; weak in CS light and per0LD light; ~0 or reversed in sss | observed (minute-level c_obs and count J, 20s-A); HMM agrees in sign in 81-100% of flies | 20b, 20d, 20e, 20m, 20s-A | DONE. NOTE per0LD: transit c (20b) ~0 but minute-level c_obs +0.04/+0.10 |
| C3 | In CS the circulation switches on at night | observed, 6-h blocks | 20c (HMM); 20s-A observable: night - day c_obs +0.12 (80% of flies, p 0.0004) | DONE |
| C4 | Cycle flux J, affinity A and EPR per genotype x phase: new quantities for sleep-wake | observed (count-based, rule of three; exact observed-chain values); HMM value as upper bracket (J x1.0-1.2, EPR x1.2-3.6 per fly) | 20s-A | DONE |
| C5 | EPR grows faster than linearly with J (CS, per0LD, sss; ambiguous per0DD, fumin); Onsager L lower in sss (p < 0.001) and higher in per0 DD (p 0.026) than in CS by fly-bootstrap difference; per0 LD, fumin, CS DD not distinguishable from CS (the hidden-level per0LD > CS does NOT hold) | observed | 20s-A + addendum 2 (20j hidden for comparison) | DONE, claim REVISED |
| C6 | Energy landscape: E3 - E1 negative in CS/per0 (and CS DD), positive in fumin/sss | observed chain | 20s-A (20l HMM agrees) | DONE |
| C7 | Memory: how a fly leaves sleep depends on how long it slept (CS night) | observed (within-fly null), confirmed hidden (400-sim null) | 20k, 20m | done |
| C8 | Wake-bout hazard rises after a trough in per0 DD (bin-sensitive); CS DD to be tested | observed bouts, 1 s | FINDINGS; root ledger 2026-09-15; 20r | IN PROGRESS (lClaude). 20r gate failed: the stored per0 DD bouts (3989/3987) cannot be regenerated (file lost; current and pre-09-22 pipelines both give 4017 wake / 4022 sleep on the same cut files). Route: one documented extractor (validated on golden G1), per0 DD rise test REPLICATED on re-extracted bouts, then CS DD, 20r-B, 20r-C with the same extractor |
| C9 | Beneath the observed cycle: hidden layer thin; extra irreversibility only in per0 DD; cycle sits at the TUR bound | model-free probes | 20n, 20o, 20p, 20q | done |

Supplementary (done, keep): 20e flips, 20f estimator calibration, 20g round trips, 20h dwell times, 20i
missing-reference sensitivity, 20m resolution, 20o drift.

NOT claimed in paper 1: hidden states are the regulators; the HMM's EPR is the dissipation; anything from the old
per-window S/W HMM (old Fig 3, which never left the lab).

### To do (finite; no new tests)
- A0  PRESPEC 20s-A, committed before running: the 20b/20c/20d/20j/20l procedures applied unchanged to the
      observed S/Q/W chain and counts; per fly x phase; HMM values carried beside as a bracket.
- A1  Observable table per fly x phase: c_obs, count J (cycles/h), count EPR with rule of three, exact
      observed-chain J, A, sigma; same for 6-h blocks (C3). Genotype medians with fly-bootstrap CIs.
- A2  EPR vs J (linear vs quadratic, as 20j) and Onsager L (onsager_L.py procedure) at the observable level.
- A3  Energy landscape from the observed chains (sqw_ref.landscape).
- A4  20r (lClaude, MATLAB): gate, CS DD hazard rise, renewal control (root PRESPEC.md, 2026-10-31).
- A5  Figures: redo at the observable level Figs 6/6b, 11, 12, 14/15; keep 7, 8, 9, 13d, 16, 17; combine 18 + 19
      into one 'beneath the surface' figure. Genotype order rule for final figures.
- A6  Methods (tracked docx): observable-first framing; HMM as robustness and upper bracket; validation items for
      20n-20q.
STOPPING RULE: a new analysis runs only if a referee asks for it or a sentence in the draft cannot be written
without it. Anything else goes to Plan B or 'future work'.

---------------------------------------------------------------------------------------------------------------

## PLAN B — paper 2 (or merged): sleep-wake regulators as drives on the observed cycle

Question: which transitions of the observed S/Q/W cycle do the circadian clock and the sleep homeostat act on,
and what do they cost?

Model: the observed S/Q/W chain with time-dependent rates
  k_ij(t) = k_ij^0 * exp[ a_ij cos(2 pi t / tau + phi_ij) + b_ij H(t) ],
tau = the fly's free-running period (24 h in LD; per0 DD: no clock term), H(t) = Process S built from the fly's
own sleep history (rises in wake, decays in sleep, fitted time constants). Few parameters per rate; regulators
pinned by observables rather than free hidden states (why this can work where PRESPEC 12-13 did not).

Readouts:
- a_ij (which transitions the clock modulates), b_ij (which the homeostat modulates); genotype contrasts.
- Thermodynamics with time-dependent rates: housekeeping vs non-adiabatic (driving) entropy production
  (Esposito & Van den Broeck 2010; Hatano & Sasa 2001); A(t), J(t), L along circadian phase and H;
  TUR for time-dependent driving (Koyuk & Seifert 2019) - makes the LD and CS DD blocks testable.
- Lower bound on the clock's own dissipation from the precision of the free-running rhythm (TUR on the circadian
  cycle count: sigma >= 2 tau / s^2 per hour, s = period jitter), per0 DD as negative control.
- Sleep deprivation: dose-response (sleep lost per fly -> change in recovery-hour rates, ZT 0-6), controls as a
  time-trend check (user's within-fly design primary).
Controls: per0 DD (backcrossed to CS; ask: number of generations) vs CS DD; LD genotypes.
Data still needed: deprivation recordings; longer CS DD records (>= 7 days) for the clock-precision bound.

### Feasibility checks (step 3; each is a yes/no, time-boxed; synthetic first)
| # | Check | Pass criterion (fix before running) |
|---|---|---|
| F1 | Synthetic recovery of a_ij, b_ij at our data size (4-day DD, 2-4-day LD flies) | the modulated row is identified in >= 80% of synthetic flies; b within a stated tolerance |
| F2 | Separability of clock and homeostat in real flies: correlation between cos(2 pi t / tau) and H(t) | |r| < 0.8 in most flies (else the two cannot be told apart from behaviour) |
| F3 | Informativeness of the clock-precision bound: period jitter s measurable with the CS DD records at hand | enough cycles per fly, or a pooled estimate with a stated CI |
| F4 | Deprivation data: format, n, spread of sleep lost across flies | enough spread for a dose-response |
| F5 | Synthetic recovery of the housekeeping / driving split of entropy production | recovered within a stated tolerance |
| F6 | Real-data pilot, CS DD vs per0 DD only: covariate model vs constant chain (bootstrap); does the 20q extra variance vanish once modulation is modelled? | pre-registered before the pilot |

Risk / reward: higher risk (identifiability, data length), higher reward (which transitions each regulator acts
on; a measured cost of being driven; a bound on the clock's cost from behaviour).

---------------------------------------------------------------------------------------------------------------

## Step 4 - decision: one merged paper, or two
Merge if F1, F2 and F6 pass and Plan B adds at least one clean, figure-level result (e.g. the clock acts on
specific transitions, with per0 DD null) without reopening Plan A's claims. Otherwise paper 1 = Plan A with a short
outlook; Plan B becomes paper 2.
