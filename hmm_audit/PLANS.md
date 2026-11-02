# Plan A and Plan B (written 2026-10-31; updated 2026-11-01: ONE merged paper)

**ANALYSIS FREEZE 2026-11-01** (ledger): no new PRESPEC unless a referee asks or the controlled deprivation experiment is run.

**Status 2026-11-01.** Plan A and Plan B are merged into ONE paper (Sheyum). The final claim list is C1-C13 below;
the paper outline is in 'The merged paper' at the end. Plan B's per-fly regulator model failed feasibility (20t F1,
F5); its group-level tests 1-4 have run (planB-feasibility, RESULTS_LEDGER_planB.md); test 5 is a Discussion
paragraph. Open: a controlled deprivation experiment (C11). All pre-registered analyses are complete.

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
   (20q). The HMM's EPR is therefore a model-based estimate, not a bound and not the measured value.
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
| C1 | Splitting wake at theta x r_f into quiet (Q) and active (W) minutes reveals a directed cycle, robust to theta (0.33/0.5/0.75) and to 30-s resolution. NOT claimed: that wake has two discrete states (no bimodality test) | observed | 20b theta sensitivity; 20m (30 s) | DONE, REWORDED 2026-11-01 (20t critique) |
| C2 | Sleep-wake transitions circulate W->Q->S->W: forward in CS night, per0 (LD dark, DD), fumin, CS DD; weak in CS light and per0LD light. HIGHLIGHT FOR RESULTS (2026-11-02): sleepless REVERSES the cycle in the light phase - flies fall asleep from active wake and wake into quiet wake. Observed chain J mean -0.038 cycles/h over 50 flies; 13 have no cycle (|J| < 1e-12, set to 0); of the 37 with a cycle, J < 0 in 70%, Wilcoxon p 0.002 (corrected 2026-11-02 after lClaude found the earlier 64% / p 0.004 counted round-off signs of no-cycle flies); HMM J -0.062, p 7e-5 (check); model-free current (20d, pre-registered) mean -0.093, p 0.007; holds at theta 0.33/0.5/0.75. At night sss has no net direction (chain J ~0, p 0.59). Presentation: per-fly MEDIANS hide it (many sss flies have no qualifying transitions, value exactly 0), so the main figure shows per-fly chain J with the HMM beside it (fig_sss_reversal) | observed (minute-level c_obs and count J, 20s-A); HMM agrees in sign in 81-100% of flies | 20b, 20d, 20e, 20m, 20s-A | DONE. NOTE per0LD: transit c (20b) ~0 but minute-level c_obs +0.04/+0.10 |
| C3 | In CS the circulation strengthens at night (mainly nocturnal); it is already forward late in the light (ZT 6-12 c_obs +0.043, J +0.125) | observed, 6-h blocks | 20c (HMM); 20s-A: night - day c_obs +0.12 (80% of flies, p 0.0004) | DONE, REWORDED (was 'switches on'); see C10 for the around-the-clock test |
| C4 | Cycle flux J, affinity A and EPR per genotype x phase: new quantities for sleep-wake | observed (count-based, rule of three; exact observed-chain values, chain-EPR medians with Inf ranked highest); model-based quantities (chain A, chain EPR, barriers, HMM model EPR) use the RULE OF THREE for one-way edges since 2026-11-02 (finite 95% lower bounds; no 'infinite'; ledger amendment); the HMM's values are a MODEL-BASED ESTIMATE, not a bound: typically larger (median x1.2-3.6 per fly) but below the observed chain in CS, fumin and sss light | 20s-A + correction 2026-11-01 | DONE, REWORDED ('upper bracket' removed) |
| C5 | EPR grows faster than linearly with J (clear in CS, per0LD, sss count points; ambiguous per0DD, fumin) - the near-equilibrium expectation for a single cycle, a consistency result. Onsager L lower in sss than CS (fly-bootstrap difference p < 0.001, BH q 0.003); per0 DD higher than CS only nominally (p 0.026, BH q 0.065, and DD vs LD); per0 LD, fumin, CS DD not distinguishable from CS. Like-for-like per0 DD vs CS DD not tested | observed | 20s-A + addenda | DONE, REWORDED |
| C6 | Energy landscape: E3 - E1 negative in CS/per0 (and CS DD), positive in fumin/sss | observed chain | 20s-A (20l HMM agrees) | DONE |
| C7 | Memory: how a fly leaves sleep depends on how long it slept (CS night) | observed (within-fly null), confirmed hidden (400-sim null) | 20k, 20m | done |
| C8 | per0 DD wake-bout hazard rises after a trough (supported but bin-sensitive: RiseNullP 0.001/0.009/0.059/0.086 at 12/16/20/24 bins; reproducible from raw files and documented code). Not produced by burst-pause scoring (renewal null makes LESS of a rise than the falling-hazard null; by rule inconclusive, p_ren 0.015/0.015/0.050/0.050). CS DD shows no rise, pooled or in either 12-h half (p ~ 1, RiseRatio < 1), with power to see a per0-sized rise at 12 bins of 1.00 at CS DD's pooled n (3309), 0.83 and 0.77 at the half sizes (1773, 1536); the per0 DD series branch (>= 2 stages) collapses to exponentials in CS DD | observed bouts, 1 s | root ledger, branch 20r-lclaude (45abec1) | DONE; all 20r analyses run and recorded (20r-lclaude head 32fe85b); CS DD 20r-B and 20r-C not run, per their rules. Not reportable: 'the rise needs a broken clock' (one cohort each); 'beyond renewal' (rule not met). 2026-11-01, test 1 (20u): the best time-reversible aggregated HMM (<= 6 hidden states per fly) gave 22.5% more wake bouts than observed (gate 15%; total sleep right, ratio 1.015) -> NULL INADEQUATE, stage 2 not run; C8 unchanged; Kmax not raised (stopping rule) |
| C10 | The circadian clock modulates the cycle: a 24-h cosine of the pooled 2-h cycle-flux profile, against a null that rotates each fly's record by a random offset (keeps each fly's rhythm and noise, destroys shared clock phase), detects clock-locked modulation in CS LD (A 0.077 cycles/h, peak ZT16, Holm 0.015) and CS DD (A 0.063, peak 23 h clock time from DD-day start, Holm 0.030), not in per0 LD (Holm 0.059, near miss, not a result) and not in per0 DD (A 0.006, p 0.96, clean negative control). Calibration 0.045 at alpha 0.05. LD = light plus clock; CS DD amplitude is a lower estimate (free-running smear). c_obs agrees in sign and phase | observed, pooled 2-h clock bins | 20w-B (planB-feasibility 803a3c4; replaces 20w's insensitive range summary, which stays recorded) | DONE |
| C11 | Sleep deprivation (9 h, ZT15-24 of day 3; 26 CS males): NO change of c_obs, J, EPR, the six transitions or wake activity in recovery ZT0-6 beyond the day-to-day drift of the same flies (placebo = day 2 vs days 0-1). ZT0-6 sleep rises through the whole recording (+28 min on day 2 before deprivation, p 0.007; +26..+58 on days 5-8), so the rebound is not separable from the drift (day 4 beyond placebo +11.6 min, p 0.094). The 20x transition changes (W->S x1.49, W->Q x1.17, S->W x0.75) do NOT survive the post hoc placebo arm | observed | 20x + post hoc P1-P2 (planB-feasibility c98ffb4, 57d0483) | DONE as a NULL (supplement + one Results line). Homeostat claims wait for a controlled experiment (parallel undeprived controls, ideally staggered; Sheyum to decide ~2026-11-08) |
| C12 | REVISED 2026-11-01: a short circuit-level control, not a headline. Silencing reduces sleep (dfb ~19% at night only; vnc ~40% day and night) through fewer entries into sleep and earlier exits; E_S - E_Q follows from occupancy (bookkeeping, not an independent finding); no detectable change of the cycle's drive in a background where the cycle is weak. 20z Part A withdrawn (freeze). Original entry: Circuit silencing (constitutive Kir2.1; empty / dfb / vnc split-GAL4; common non-CS background; males; LD; 19/25/34 flies): vnc>Kir sleeps ~40% less in light and dark, dfb>Kir ~19% less at night only; NEITHER changes c_obs, cycle flux or EPR detectably (12 tests, all Holm >= 0.66; largest J dark dfb +0.10 [+0.04, +0.17], MW p 0.055). The changes sit in the sleep transitions (Q->S down, S->W and S->Q up) and in E_S - E_Q (raised: +0.4 to +0.6 nats dfb, +1.1 vnc), i.e. a shallower, less stable sleep state, not a change of the cycle's drive (secondary, 'consistent with'). Exit-route-vs-sleep-length dependence present in all groups, unchanged by silencing (20z C). Sleep-substate test on bout durations pending (20z Part A, lClaude) | observed | 20y (ledger cf777b1), 20z Part C (7626ecc) | 20y DONE; 20z A PENDING. Compare ONLY with empty>Kir; constitutive; split-GAL4 vs full-GAL4 control |
| C9* | (C9, 2026-11-01 addition) Waiting times: Q and W dwells depend on the direction of passage beyond the observed chain in per0 DD (E_WT +0.083 nats/h, Holm 0.00056) and CS DD (+0.078, Holm 0.00018), comparable to the chain's own EPR. The 1-s -> symbol scoring is EXACTLY time-symmetric (amendment F1: 0 of 275,040 minutes differ), so the 5-min criterion cannot create the dependence from reversible dynamics. Amendment F F2/F3 (planB-feasibility 4d27f50): NO verdict in either block, both NULL INADEQUATE (reversible K <= 6 models give 22.5% (per0 DD) and 40% (CS DD) excess wake bouts). FINAL wording: '... beyond what the observed chain allows. The scoring is exactly time-symmetric (verified), so the sleep criterion cannot produce this from reversible dynamics; a reversible 1-s model adequate to the bout process could not be built within six hidden states, so no scoring-matched verdict exists.' NOT added: 'together with C13, irreversibility at the seconds scale' (C13 does not go beyond the chain). LD descriptive: CS dark +0.14, fumin L/D positive, sss ~0 | observed, 1-min dwells | 20v + amendment F (planB-feasibility b4467b0, e93b0d4) | DONE (F: no verdict, null inadequate in both blocks) |
| C13 | REVISED 2026-11-01 (Sheyum's question): this is a 1-s, model-free CONFIRMATION of the cycle's irreversibility (C2/C4), NOT new irreversibility beyond it. PRESPEC 21's null is a time-reversible process, which the observed minute chain (J > 0, EPR > 0) already rejects; centring removed only the mean level difference, and the forward cycle itself (Q before sleep, W after waking, with Q more likely the closer to onset) predicts a ramp before onset and a flat return after waking. What C13 adds: (i) the asymmetry is present in raw 1-s counts without the minute grid, theta or r_f; (ii) the timescales, descriptive: wind-down over ~5 min, recovery within ~20 s. 'Beyond the observed chain' rests on C9* (20v, chain null), not on C13. Wording: 'falling asleep is not waking played backwards, already in the raw 1-s activity'. Original entry follows. | Seconds-scale irreversibility: within wake bouts (>= 600 s), the 1-s activity trajectory out of sleep is not the time-reverse of the trajectory into sleep beyond their level difference (shape statistic, level removed; null = per-bout sign flips, robust to per-day block flips and a group fly-day null). per0 DD: per-bout z median +0.58 (Holm 0.014), 32% of flies p < 0.05 (binom 2e-05), day-block Holm 0.005, group p 0.0005. CS DD: z +1.40 (Holm 4e-05), day-block Holm 8e-06, group 0.0005. Secondary window (120 s, 5-s bins) agrees; count-intensity shape agrees; level difference (C2 counterpart) +0.04. Synthetic reversible traces (20u fit) gave no positive. Wording: per0 DD 'time-irreversibility of the activity process at the seconds scale, under stationarity and the verified time-symmetric scoring'; CS DD adds 'circadian drive not excluded'. Never 'proves broken detailed balance'. Direction (DESCRIPTIVE, fig_planB21_profiles_primary, planB-feasibility): activity winds down GRADUALLY over the last ~5 min before sleep onset (per0 DD ~0.14 -> ~0.08 active fraction in the last 30-60 s) and returns ABRUPTLY on waking (back to ~0.13-0.15 within ~20 s, then flat); CS DD the same ordering, largest within ~60 s. Paper sentence: 'flies wind down gradually into sleep but wake up abruptly'. LD descriptive: CS, per0LD, fumin positive; sss weak, level reversed (as its C2) | observed, 1-s counts, model-free | PRESPEC 21 + amendment B (planB-feasibility a44ba47, 7429217) | DONE (figure done). Supports C9*: the 20v waiting-time excess now has a model-free seconds-scale counterpart |
| C9 | Beneath the observed cycle: the hidden layer is thin (20n, 20o); extra irreversibility from 3-min memory in per0 DD (also without sleep) and in fumin dark (sleep-criterion caveat) (20p); no TUR evidence of hidden dissipation: the observed current is at most as precise as the observed chain allows, reaching that level only in per0 DD and CS dark (the chain itself sits at its bound by construction near equilibrium) (20q) | model-free probes | 20n, 20o, 20p, 20q | DONE, REWORDED |

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
- A6  Methods (tracked docx): observable-first framing; HMM as robustness check and model-based estimate; validation
      items for 20n-20q. DONE (V7). 2026-11-01: wake-bout hazard and controls (20r, 20u), waiting-time
      irreversibility (20v + F), clock-face test (20w-B), deprivation (20x) added; old deprivation HMM section replaced.
STOPPING RULE: a new analysis runs only if a referee asks for it or a sentence in the draft cannot be written
without it. Anything else goes to Plan B or 'future work'.

---------------------------------------------------------------------------------------------------------------

## PLAN B — merged into ONE paper (Sheyum 2026-11-01). Per-fly Plan B failed feasibility (20t F1, F5 on planB-feasibility). Group-level tests, all run 2026-11-01:
| test | PRESPEC | outcome | claim |
|---|---|---|---|
| 1 reversible-null control of the per0 DD rise | 20u | NULL INADEQUATE (22.5% excess wake bouts vs 15% gate); no verdict | C8 unchanged |
| 2 waiting-time irreversibility | 20v + amendment F | SUPPORTED per0 DD and CS DD; scoring exactly time-symmetric (F1); reversible-model control: per0 DD no verdict, CS DD pending | C9* |
| 3 clock face | 20w (range: insensitive, recorded) -> 20w-B (cosine) | clock-locked J modulation CS LD, CS DD; per0 DD null | C10 |
| 4 sleep deprivation | 20x + post hoc P1-P2 | NULL; rebound not separable from drift | C11 |
| 5 clock-cost bound | 20t F3 | B = 2 tau/s^2 >= 18.7 k_B/h [13, 30] in CS DD (weak: centroid marker, 3 periods per fly; per0 DD 0.9 = method floor) | Discussion paragraph, not a claim |

## (original Plan B text below): sleep-wake regulators as drives on the observed cycle

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
DECIDED 2026-11-01 (Sheyum): ONE merged paper. The criterion below is kept as written for the record.
Merge if F1, F2 and F6 pass and Plan B adds at least one clean, figure-level result (e.g. the clock acts on
specific transitions, with per0 DD null) without reopening Plan A's claims. Otherwise paper 1 = Plan A with a short
outlook; Plan B becomes paper 2.

---------------------------------------------------------------------------------------------------------------

## The story (Sheyum, 2026-11-02): broken detailed balance in three acts. SUPERSEDES the section order of the outline below
Thesis: *sleep-wake transitions run a directed cycle W->Q->S->W that breaks detailed balance; the clock shapes it but does
not create it, and how much a fly sleeps sets neither its direction nor its strength.* Genotype order as the rule
(CS -> CS DD -> per0 LD -> per0 DD -> fumin -> sss). The observed S/Q/W chain carries every act; the HMM is brought in
ONLY where a result could be a scoring artefact or where the hidden level adds something, otherwise it runs in the SI.

| act | question | content (claims) | HMM brought in (main text) | SI |
|---|---|---|---|---|
| 0 Setup | what is measured | S/Q/W symbols, observed chain, J / A / EPR; HMM introduced with its three reasons (C1) | emissions 80-100% diagonal -> hidden states = observed symbols, which licenses the chain (fig1_hmm_emissions) | K=2/3/4 order selection (K=3 in 76-96%); BIC chain vs HMM |
| 1 Wild type (CS LD) | does the cycle exist, and what shapes it | forward cycle, strong at night, already forward late in the light (C2, C3); clock-locked daily rhythm of J, peak ~ZT16 (C10); A, EPR, superlinear EPR(J), Onsager L (C4, C5); landscape, sleep a deep well (C6); the forward step S->W is carried by LONG sleep bouts (C7) | direction agrees in 81-100% of flies; memory C7 also at the hidden level (CS night, beyond all 400 Markov simulations) | hidden A_eff ~2x (observed EPR a lower bound); landscape obs vs hidden (rho 0.83-0.97) |
| 2 Remove the clock (CS DD, then per0 DD; per0 LD one line/SI) | does the drive need the clock | per0 DD still circulates forward with finite A and EPR, but its daily rhythm is flat (A 0.006, p 0.96) while CS DD keeps one (C10): the clock modulates, does not generate. Beneath the chain, in BOTH CS DD and per0 DD: direction-dependent dwells (C9*), 1-s asymmetry (C13). per0 DD only: wake hazard falls then rises, sequential exit route (C8); the waiting-time signature expected of broken detailed balance (reversible finite Markov sojourns have falling hazards, Keilson 1979), with the sleep-criterion and bin-sensitivity caveats | the hidden layer is thin: its gain over the chain equals one extra minute of memory (C9, 20n/20o) - the irreversibility beneath the chain is NOT a missing hidden state of the K=3 kind | K=4 composition (second sleep state in 16/25 per0 DD flies; lower bound); TUR; 20u null inadequate |
| 3 Dial sleep down (fumin, sss) | does sleep amount set the drive | both invert the sleep well (E_S > E_Q, C6); opposite thermodynamic phenotypes: fumin drives the cycle hardest (chain abs A 1.04, forward in both phases), sss REVERSES it by day, no net direction at night, lowest Onsager L (C2, C5). Sleep amount and the breaking of detailed balance are dissociable | sss reversal at the hidden level (p 7e-5; fig_sss_reversal triangles); L lowest at both levels; inverted well at both levels -> none of these is a Q/W scoring artefact in the mutants, where scoring is most suspect | HMM creates more one-way edges than it rescues (31 vs 6, mostly fumin/sss); routing prediction unreliable in sss |
| Support | - | dFB/VNC silencing: a circuit-level sleep reduction that leaves the drive unchanged, supports act 3 (one sentence; C12); memory not abolished by either line. Deprivation null with placebo arm (one line; C11) | none | C12, C11 sections |

Structural point (decided with the story): the reference for per0 DD is CS DD, not CS LD (light and clock otherwise
confounded). Act 2 panels put CS DD and per0 DD side by side, DESCRIPTIVELY; the like-for-like per0 DD vs CS DD L
comparison was never tested and is not added (freeze).

## The merged paper (outline, section order SUPERSEDED by the story above; REVISED 2026-11-02: HMM introduced early, brought in where a result could be a scoring artefact)
Working title (PNAS): *Fruit flies circulate through sleep and wake with a direction, and the circadian clock modulates
the cycle.* Everything at the observable level, per fly; genotype order rule for figures. Draft HMM text:
DRAFT_HMM_motivation.md.

| section | content | claims | HMM role | figures |
|---|---|---|---|---|
| 1 | Three behavioural states and a model of them: wake splits into Q and W; each fly described by the observed S/Q/W chain AND a 3-state HMM allowing misclassification | C1 | INTRODUCED with three reasons (scoring error, hidden dissipation, internal states of sleep/wake); hidden states = observed symbols (emissions 80-100% diagonal) -> licenses the observed chain | old 1-2 + fig1_hmm_emissions |
| 2 | A directed cycle; SLEEPLESS REVERSES IT BY DAY | C2 | sign agreement 81-100% of flies; sss reversal at both levels | A1 (sss mean marker) + fig_sss_reversal (HMM triangles) |
| 3 | Thermodynamics: J, EPR, Onsager L, affinity per cycle A_eff | C4, C5 | one line + SI table: same L ranking (sss lowest at both levels); hidden A_eff ~2x (observed EPR is a lower bound) -> motivates section 7 | A2 |
| 4 | Energy landscape | C6 | none | A3 |
| 5 | The clock modulates the cycle | C3, C10 | none | C10 |
| 6 | Memory and hazards | C7, C8 | one line: memory also at the hidden level (CS night) | C7 observed, C8 |
| 7 | What lies beneath the observed cycle | C9, C9*, C13 | full test: BIC vs chain, gain = one extra minute of memory, TUR, waiting-time asymmetry (in active wake) | A4 (HMM panels), dwell-direction figure, 1-s profiles |
| 8 | Circuits (dfb/vnc control), deprivation null, Discussion (clock-cost bound; what neural recordings could add to a hidden layer) | C11, C12 | none | C12 (SI or small main panel) |

NOT claimed: hidden states are regulators; the HMM's EPR is the dissipation; the rise needs a broken clock; detailed
balance proven broken; any per-fly regulator parameter; any deprivation effect.
