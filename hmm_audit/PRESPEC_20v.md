# PRESPEC 20v — merged paper, test 2: waiting-time irreversibility beyond the observed chain
Date 2026-11-01. Written before anything runs; the go-ahead to RUN is Sheyum's (not given by this file).
Owner: Plan B feasibility session. Design agreed with the S/Q/W audit session. Results go to RESULTS_LEDGER_planB.md.

## Question
Do the dwell times in Q and W depend on where the fly came from and where it goes next, more than a first-order chain
allows? In a first-order Markov chain the dwell in a state is independent of its entry and exit neighbours, so any
dependence is irreversibility the transition counts miss. It is the waiting-time route to hidden entropy production
(Martinez, Bisker, Horowitz & Parrondo 2019 Nat Commun 10:3542; Skinner & Dunkel 2021 PRL 127:198101;
van der Meer, Ertel & Seifert 2022 PRX 12:031025; Harunari, Dutta, Polettini & Roldan 2022 PRX 12:041026).

## Data, symbols, blocks
1-min symbols exactly as sqw_ref.sqw_symbols (theta 0.5). DD: light_dark=False (one whole-record reference). LD:
light_dark=True. Fly-phases without a Q/W reference are excluded (20p amendment B).
PRIMARY: per0 DD whole record; CS DD whole record if uploaded before the run (Holm across the primary blocks present).
DESCRIPTIVE: CS, per0LD, fumin, sss, light and dark phases, with the 6-h-rate null (below), because within-phase drift
(20o) makes LD records non-stationary.

## Statistic (per fly x phase)
Runs of Q (and separately of W) with a known entry symbol e and exit symbol x (neighbouring minutes non-missing, in
the same phase). Only the mixed pairs carry direction:
- Q dwells: (e, x) = (W, S) vs its reverse (S, W).
- W dwells: (Q, S) vs (S, Q).
S dwells are NOT used (5-min floor); S enters only as an entry or exit event.
Dwell bins fixed now (minutes): 1, 2, 3, 4-5, 6-9, >= 10; pseudocount 0.5 per bin.
sigma_WT = sum over s in {Q, W} and the two mixed pairs of (n_{s,e,x} / hours) x KL( psi_{s,e,x} || psi_{s,x,e} ),
in nats/h, where psi is the binned dwell distribution and hours = minutes in the phase / 60.

## Null and excess (as 20p)
Null A: 200 records per fly-phase simulated from that fly's own mc_fit chain (12-h phases in LD; whole record in DD),
same length, missing minutes reinserted at their observed positions, scored identically.
LD descriptive null: the same from the 6-h-rate chain (mc_fit with phases(T, '6h')).
E_WT = sigma_WT(obs) - median sigma_WT(null); per-fly p = (1 + #{null >= obs}) / 201.

## Decision rule, fixed now
Per primary block: SUPPORTED if Holm p < 0.05 for the group Wilcoxon of E_WT against 0 with median E_WT > 0, AND the
fraction of flies with p < 0.05 exceeds Binomial(n, 0.05) at p < 0.05. Otherwise NOT SUPPORTED.
Wording if supported: "dwell times in quiet / active wake depend on the direction of passage beyond what the observed
chain allows: a waiting-time lower bound on hidden irreversibility of E_WT nats/h." E_WT is an EXCESS over the null,
not an EPR. The bound statement uses sigma_MC + E_WT, as 20p did with the triplet excess.
Wording if not: "no waiting-time evidence of hidden irreversibility at 1-min resolution." The coarse resolution is
stated (Q runs median 1-2 min), so a null result is weak.
sss gets one descriptive line: E_WT where the observed current is ~0 (C2) bears on a "hidden cycle without a
visible current" (Martinez et al. 2019). It is not part of the decision rule.

## Matched null, narrowest block, time box
Same fly, same length, same mask, same estimator. Narrowest block: one fly's whole DD record (LD: one fly's light or
dark minutes pooled over days). Expected run time under 1 h; no extension.

## Amendment F (2026-11-01; approved by Sheyum; brief by cClaude) — the scoring-matched null. Committed BEFORE any part runs.
20v's recorded verdict stands. This adds one check against the stated caveat ("a sleep-scoring origin is not excluded")
and fixes its wording rule in advance.

F0. The argument. If the 1-s activity process is stationary and time-reversible, and the scoring is time-symmetric
(scoring the reversed trace gives the reversed symbol sequence), the scored symbol process is reversible in
distribution. The reversal of a Q run entered from W and left to S is a Q run entered from S and left to W, with the
same length and probability. So psi_{s,e,x} = psi_{s,x,e} exactly, and sigma_WT has expectation equal to the
estimator's positive bias only.
Each scoring step is symmetric on paper:
- the 300-s silence criterion;
- bridging of single active seconds (>= 300 s of silence on both sides);
- the >= 30-of-60-s rule for S;
- Q/W by the active fraction;
- r_f from the middle 300 s of long wake bouts (symmetric up to one second of rounding at mid = (w0+w1)//2).
So the scoring cannot by itself create direction-dependent dwells from reversible dynamics. If F1 and F2 confirm
this, the 20v excess comes from irreversible activity at the sub-minute scale (e.g. winding down before sleep vs
abrupt waking) or from non-stationarity, not from the 5-min rule.

F1. Gate: time symmetry of the scoring IN THE CODE (deterministic, no simulation).
- Inputs: every per0 DD and CS DD retained trace, trimmed to whole days (per0 DD: first traceLen seconds of the 20r
  reference, floored to whole days; CS DD: the whole 4.000-d files).
- Compare y1 = reverse(sqw_symbols(x)) with y2 = sqw_symbols(reverse(x)), both light_dark=False, theta 0.5.
- PASS: y1 = y2 in >= 99.9% of minutes for every fly. Differing minutes are listed per fly; the only expected cause is
  Q/W flips at the threshold from the one-second r_f asymmetry.
- FAIL: name the asymmetric step and how many minutes it moves. F2 still runs and then measures what the asymmetry
  produces. The scoring is not altered.

F2. Null distribution; per0 DD primary, CS DD secondary.
- per0 DD, no new model: 20u stage 2's simulated 1-s traces (each fly's selected reversible aggregated HMM, 199
  replicates, seeds 5000 + b). They are regenerated bit-identically from the same fits and seeds, because 20u does not
  store traces. Fly-replicates excluded by 20u's death criterion are left out and counted.
- Scoring: each replicate with sqw_symbols(light_dark=False, theta=0.5) on the same traceLen (whole days), using the
  replicate's own r_f, exactly as the observed trace is scored.
- Count matching (rule 5): in each replicate, the runs in every (s, e, x) cell are subsampled without replacement to
  the observed fly's count (seed default_rng(25_000 + 100 b + fly)). A cell with fewer runs than observed keeps all
  its runs and is flagged. sigma_WT uses the observed fly's hours.
- Observed side, recomputed on the same traces: sigma_WT on 20u's retained traces (first traceLen s, whole days),
  with the original 20v per0 DD value reported beside it. The matched recomputation is the one that counts.
- CS DD (secondary): reversible HMMs per CS DD fly by 20u's stage 1 procedure, unchanged:
  - orders (Ka, Ki) in {1,2,3}^2, 10 starts, BIC;
  - seeds 3000 + 100 + fly (fits), 4000 + 100 + fly (adequacy), 5000 + 100 + b (199 null replicates);
  - 20u's stage 1b adequacy gate (wake-bout count and total sleep, medians in [0.85, 1.15]), with CS DD's observed
    bouts from the same Python extractor on its whole-day files.
  Adequacy FAIL => CS DD "NULL INADEQUATE", no verdict.
- Added adequacy check for both blocks, reported, not a gate: per fly, simulated / observed numbers of Q, W and S runs.

F3. Statistic and decision rule, fixed now.
- Per fly: E_score = sigma_WT(obs) - median sigma_WT(scored reversible null, count-matched);
  p = (1 + #{null >= obs}) / (N_kept + 1).
- Per block, the 20v rule unchanged: group Wilcoxon of E_score vs 0 with median > 0 at Holm p < 0.05 (Holm over the
  blocks that pass adequacy) AND the fraction of flies with p < 0.05 above Binomial(n, 0.05) at p < 0.05.
- Descriptive: median sigma_WT under the scored reversible null against the median under 20v's chain null (both
  should sit near the estimator bias under F0; any difference is recorded).

F4. Wording, fixed now.
- F1 PASS and F3 SUPPORTED: "Dwell times in quiet and active wake depend on the direction of passage beyond what the
  observed chain allows. Scoring is time-symmetric (verified), and a time-reversible model of each fly's 1-s activity,
  scored identically, does not produce the dependence (E_score ...). The excess therefore reflects time-irreversible
  activity below the 1-min resolution, not the sleep criterion." The scoring caveat is then lifted.
  - per0 DD: stationary and undriven.
  - CS DD: add "the circadian drive is not excluded as its source". A slowly driven process with reversible
    instantaneous rates can also produce direction dependence, which would be a driving cost, not hidden-state
    dissipation.
- F3 NOT SUPPORTED: "A time-reversible 1-s process, scored identically, accounts for the waiting-time excess; 20v does
  not evidence hidden irreversibility." If F1 passed, the source is the reversible model's dwell structure, not
  scoring asymmetry; record which.
- F1 FAIL: add "the scoring is not exactly time-symmetric (step ...); E_score is measured against a null that carries
  the same asymmetry."
- Never: "proves broken detailed balance" (one model class, up to 6 hidden states, was searched).

F5. Order and time box.
- F1 first (minutes).
- per0 DD F2 as soon as 20u's stage 1 fits exist; it does not depend on 20u's verdict.
- CS DD fits afterwards, under 20u's time box: 6 h on 4 cores, or "not completed". Orders and starts are not reduced
  after the fact.
- One run, one verdict. Ledger entry "PRESPEC 20v amendment F".
