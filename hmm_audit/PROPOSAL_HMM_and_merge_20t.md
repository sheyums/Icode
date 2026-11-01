# What an HMM on S/Q/W can add to paper 1, and which Plan B results can be merged (proposal, 2026-11-01)
Written by the Plan B feasibility session at the user's request ("think about including an HMM component in Plan A";
"merge Plan A with ideas/results from Plan B"; "over the top"). NOTHING here has been run. Each item says what
it would claim, what supports it already (ledger), the decisive control, and the cost. Every item needs its own PRESPEC
committed before it runs (root CLAUDE.md). Plan A's stopping rule ("a new analysis runs only if a referee asks or a
sentence cannot be written without it") is the user's to override. This document asks for that decision; it does not
take it.

## 1. Why the earlier verdict on the HMM was right, and also too narrow
The audit asked one question of the HMM: are its hidden states noisy versions of the symbols? The answer was "barely"
(emissions 80-100% diagonal; BIC prefers the plain chain in ~75% of flies; 20n, 20o). That verdict is correct. But it
tested only one KIND of hidden structure: a hidden layer BESIDE the symbols (misclassification). The ledger points
to a different kind: hidden structure INSIDE the symbols, i.e. several hidden states that emit the SAME symbol.
That is an aggregated Markov model, the standard model of single ion channels: several closed states, one current level
(Colquhoun & Hawkes 1981 Proc R Soc B 211:205; Fredkin & Rice 1986 J Appl Probab 23:208). The evidence already in hand:
- **Dwell times are not geometric.** 20h: observed S, wake, Q and W run hazards lie outside the Markov model's 95%
  band in 17-100% of bins. Wake runs worst: 56-80% of bins in CS, per0LD and fumin.
- **Wake hazards RISE after onset in several blocks.** 20h, minute level: CS dark 0.061 -> 0.088 per min at 5-6
  bins; per0LD dark 0.065 -> 0.095; fumin light 0.015 -> 0.030; fumin dark 0.009 -> 0.028. At 1 s, the per0 DD
  wake-bout hazard rises after a trough (C8, bin-sensitive) and CS DD's does not.
- **Exit routing remembers bout length** (C7, 20k). This is the signature of hidden sub-states within S.
- **One extra minute of memory beats the hidden-noise layer** (20o: R_mem 1.25-2.0), which is what an aggregated
  model produces.
- **Extra irreversibility beyond pair statistics in per0 DD** (20p), present without S (wake-only blocks).
So the HMM was looked for in the wrong place. The question that remains is not "are the symbols noisy?" but
"how many hidden stages does each symbol hide, and are they traversed in a direction?"

## 2. The flagship idea: a model-free theorem turns the hazard shape into a statement about irreversibility
**Theorem (Keilson 1979, Markov Chain Models: Rarity and Exponentiality, Springer; ch. 5).** For a stationary,
time-reversible (detailed-balance) Markov chain, the sojourn time in any subset of states, entered from stationarity,
has a completely monotone density. That density is a positive mixture of exponentials, so its hazard can only FALL.
Why: the generator restricted to the subset is similar to a symmetric matrix, so the dwell density is
sum_k c_k^2 lambda_k e^(-lambda_k t).
**Consequence.** A dwell hazard that RISES in a block with time-independent rates cannot come from ANY reversible
hidden Markov model, however many hidden states it has. Either the hidden dynamics inside that symbol break detailed
balance (hidden stages traversed in a preferred direction, i.e. dissipation), or the rates are driven in time.
- No information criterion is involved and no particular HMM is fitted. Like the structural argument in root
  CLAUDE.md, it is the stronger kind of claim.
- **per0 DD is the one block with neither clock drive nor light** (20q primary steady-state block). A rising hazard
  there therefore means broken detailed balance in hidden wake dynamics. It is independent of, and agrees with, 20p,
  which found extra irreversibility ONLY in per0 DD (and fumin dark), present in wake-only blocks.
- In LD and CS DD, the clock or light can produce a rise, so there it reads as "driven". That contrast IS the Plan B
  question (clock = drive), asked of a shape rather than a fitted rate.
- **Mechanism, and the bridge to Plan B.** A homeostat (Process S, Borbely 1982 Hum Neurobiol 1:195) that builds up
  during wake makes the hazard of falling asleep rise with time awake. A chain of m hidden wake stages traversed in one
  direction is the discrete form of exactly that build-up: the "linear chain trick" (MacDonald 1978, Time Lags in
  Biological Models, Springer). This is the hyper-Erlang family already in the bout suite (C8). So the homeostat that
  F1 could not fit as a covariate per fly REAPPEARS as hidden stages with a measurable direction and cost. That is the
  Plan A / Plan B merge, in one object.

**Decisive control (must run first; this is where the idea can die).** A wake bout is not a plain sojourn of the 1-s
process. It ends only when a run of >= 300 s of inactivity begins, and that is decided by the future. The scoring
rule itself (5-min criterion, bridging) could manufacture a rising hazard from reversible dynamics. Test:
fit a REVERSIBLE (detailed-balance-constrained) hidden Markov model to each fly's 1-s activity, simulate 1-s traces,
pass them through the IDENTICAL sleep_bouts / sqw_symbols pipeline, and compute the same hazard and the same RiseNullP
statistic. The theorem's argument stands only where the observed rise exceeds that reversible-null distribution.
- This is the same logic as 20r-B (renewal control, pending).
- CS DD's absent rise and the hyper-Erlang collapse there (C8) are weak evidence that scoring alone does not produce
  rises, but they are not the control.
- Cost: moderate. A 1-s HMM per fly, then simulation through existing code. per0 DD first (25 flies). Time-boxed by
  a fixed number of null replicates set in the PRESPEC.
- Pre-specified reading: per0 DD rise beyond the reversible null => "hidden wake dynamics break detailed balance"
  (model-free). Not beyond it => C8 is a scoring property and must be worded as such. Either outcome improves the paper.

## 3. The HMM component proper: aggregated (sub-state) model with a thermodynamic bracket
Once item 2 says hidden stages are needed:
- **Model.** Hidden states S1..Sm_S, Q, W1..Wm_W. Each hidden state emits ONE symbol, so emissions are fixed, not
  fitted (which removes the near-identity problem 20n met). Hidden transitions are free within the cycle. Stage counts
  come from the hyper-Erlang shape sweep already in the bout suite (C8) rather than from BIC fishing.
- **Fit per genotype, pooled over flies, with per-fly random offsets.** F1 just showed that 4-day records cannot pin
  modulation per fly, while group-level recovery worked (median clock amplitude ratio 0.99-1.06, phase error under 2 h).
  The same lesson applies here.
- **Readouts.**
  - Hidden-stage mean dwells (identified even where the integer stage count is not; root CLAUDE.md step 3).
  - The hidden cycle affinity, and EPR as a model value.
  - A posterior "stage clock" inside each wake bout (how far along the build-up a fly is): an inferred sleep pressure
    per minute, per fly.
  - **sss.** Its observed cycle current is ~0 (C2), yet hidden stages can still dissipate (Martinez, Bisker, Horowitz &
    Parrondo 2019 Nat Commun 10:3542, "broken detailed balance in the absence of observable currents"). Whether
    sleepless flies run a hidden cycle with no visible current is an open, testable question.
- **The entropy-production ladder, per genotype (the figure).**
  - Lower bounds, all model-free: sigma_MC from transitions (20s-A); the triplet bound from block KL (20p); and a
    waiting-time bound, which uses dwell distributions between observed transitions (Skinner & Dunkel 2021 PRL
    127:198101; van der Meer, Ertel & Seifert 2022 PRX 12:031025; Harunari, Dutta, Polettini & Roldan 2022 PRX
    12:041026). The waiting-time bound is the one that sees hidden stages.
  - Upper rung: the aggregated model's EPR.
  - Coarse-graining can only lower the apparent EPR (Esposito 2012 PRE 85:041125; Roldan & Parrondo 2010 PRL
    105:150607), so the true value is bracketed.
  - This replaces the old "HMM EPR as upper bracket" wording, which C4 does not support per fly (see
    PLANA_CRITIQUE_20t.md), with a ladder whose rungs each have a theorem behind them.
- **Synthetic recovery first** (an F1-style check at pooled group size), with pass criteria fixed beforehand.

## 4. Plan B results that CAN be merged after F1 failed
F1 killed per-fly covariate fitting, not everything in Plan B.
- **B1, a circadian "thermodynamic clock face" (group level, descriptive; no per-fly identification).**
  - J(t), A(t), EPR(t) and the landscape E3 - E1 in ZT bins, from transitions pooled over the flies of a genotype,
    with fly-bootstrap bands. Bin width is fixed in the PRESPEC before looking (rule 4); 2 h is proposed.
  - CS LD vs per0 LD vs per0 DD (a flat clock face is the negative control), and CS DD when uploaded.
  - It extends C3 (6-h blocks) to a full phase portrait. F1's DD-CH run shows group-level clock modulation is
    recoverable at this data size.
  - The periodic-driving TUR (Koyuk & Seifert 2019 PRL 122:230601) then makes the LD and CS DD blocks testable, which
    20q could not do.
- **B2, the cost of the clock from behaviour (F3).**
  - Bound B = 2 tau / s^2 from period jitter: Barato & Seifert 2016 PRX 6:041053, "Cost and precision of Brownian
    clocks". A one-number result: the minimum dissipation a clock needs to keep the observed precision. per0 DD is the
    control.
  - Needs the CS DD files; the 4-day records give about 3 intervals per fly, so a pooled estimate. Longer records,
    if you have them (>= 7 d), make it far stronger.
- **B3, sleep deprivation as a step perturbation (F4).** Dose-response of recovery-hour rates against per-fly sleep
  lost. In the hidden-stage picture, deprivation should push flies further along the wake stages, so the prediction is
  quantitative and testable. Needs the recordings.

## 5. Suggested order (each item a yes/no with its own PRESPEC; stop at the first item that fails)
| step | item | needs | what a pass buys | cost |
|---|---|---|---|---|
| 1 | Reversible-null control of the wake-hazard rise (item 2), per0 DD first, then LD dark blocks | data in hand | the model-free headline "hidden wake dynamics break detailed balance" | moderate |
| 2 | Waiting-time EP lower bound per genotype x phase (item 3, ladder) | data in hand | a new lower rung; tests the sss "hidden cycle without current" question | low-moderate |
| 3 | Synthetic recovery, then fit of the aggregated sub-state model, pooled per genotype | data in hand | stage dwells, inferred sleep pressure, model EPR (upper rung) | high |
| 4 | Clock face (B1) | data in hand (+ CS DD) | circadian phase portrait of the cycle's thermodynamics | low |
| 5 | Clock-cost bound (B2) and deprivation (B3) | CS DD; deprivation files | a cost number; a perturbation test of the stage picture | low each |

**The paper this would make.** (i) A directed sleep-wake cycle with measured flux, affinity and dissipation (Plan A,
unchanged). (ii) Its hidden layer, found where it actually is: inside wake, as stages traversed in one direction,
proven model-free from dwell-time shape in the one undriven genotype, and bracketed by a ladder of entropy-production
bounds. (iii) The regulators: homeostatic build-up as those hidden stages, the clock as a periodic drive seen on the
clock face, and the clock's own minimum cost. Each layer is built on the one before.

**Risks to state now.**
- Item 2 can fail on the scoring control. Then C8 must be reworded and the flagship goes.
- Every added analysis adds to the multiplicity the ledger must show.
- Keilson's theorem needs stationarity. per0 DD qualifies (no clock). For LD it supports only "driven or
  nonequilibrium", not which.
- The homeostat is feedback, so it has no clean housekeeping/driving split (20t F5 caveat). The hidden-stage
  formulation avoids that by making the build-up part of the state.
