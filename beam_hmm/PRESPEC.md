# beam_hmm pre-specifications

One dated entry per analysis, committed BEFORE it runs. Use the template in the root PRESPEC.md.

## Entries

<!-- Append below. Do not edit an entry once its analysis has run. -->

### 2026-11-12  E0 — exploration of week 1 (descriptive only)

Question        What does the week-1 file contain, and what structure is visible at 0.2-0.4 s that 1-min
                S/Q/W symbols hide? Exploration to decide what is possible; NOT a test of anything.
Data            File 1 (week 1, flies 3-4 d old at start), all 64 channels, CS males, LD 12:12.
Time window     Week 1 only (file 1, start = lights-on as stated by the user).
Fixed choices   Format QC: time-step regularity, -1 lines, out-of-range values, constant channels, jumps of
                > 1 beam per sample, silence >= 12 h. Descriptive: per-fly occupancy by beam and ZT hour,
                step-size distributions, pause durations below 5 min, activity profile by ZT.
Decision rule   None: no p-values are computed in E0. Every pattern seen here is POST HOC with respect to
                any later analysis, and becomes a claim only if re-tested under its own PRESPEC on data
                not used here (weeks 2-9, or week-1 flies held out).
Matched null    None (descriptive).
Notes           Week 1 is thereby spent as an exploration set. Anything about food-end inference that is
                derived from week 1 may only be used on later weeks.

### 2026-11-12  P1 — pilot pipeline: preprocessing, inclusion, scoring, S/Q/W descriptives

Status of this dataset (user decision 2026-11-12): most channels could not be read; the experiment will be
repeated. Channels 17-32 of this run are a PILOT: used to build and validate tools and to GENERATE
hypotheses. Nothing from the pilot is a confirmed result. Each hypothesis it yields gets its own PRESPEC and
is tested on the repeat experiment with the rules below unchanged (or changed in a new, dated entry).
Disclosure: the thresholds below were chosen AFTER seeing week 1 (E0). That is acceptable for a pilot and
is exactly why confirmation must come from the repeat.

Question        (a) Which pilot flies and fly-days are usable? (b) What do 0.4-s position and speed look
                like inside S, Q and W minutes scored from the same record?
Data            Channels 17-32, every weekly file supplied (week 1 now; weeks 2-9 when uploaded).
Time window     Fly-day (ZT0-ZT0, file start = lights-on). Every summary is per fly-day, then per fly
                (median over its days in a stated week), then across flies. Flies are the unit.
Fixed choices
  Preprocessing   Drop the all -1 lines (true step 0.4 s). Keep 0 as its own position ("beyond beam 1").
                  Spike filter, one pass: a sample B with equal neighbours A (A-B-A) and |B - A| >= 4 beams
                  is set to A. Nothing else is altered. Beam 15: no action (recorded as rarely/never read).
  Artifact flag   Per fly-week, after the filter: if > 5% of position changes are jumps of >= 4 beams, the
                  fly-week is ARTIFACT-HEAVY: reported in a separate column, excluded from pooled summaries.
  Empty/dead      Channel with < 100 position changes in its first 24 h: excluded as empty or dead at load.
  Handling        Per file, a HANDLING event = the first 10-min window in which >= 50% of the included
                  channels exceed the 99th percentile of their own 10-min change counts (that file). From the
                  start of that window to the end of the file is excluded for every fly. The first 2 h of
                  each file (re-insertion) are excluded. If the user supplies transfer times, they replace
                  this detector and the change is recorded in a new entry.
  Death           Lab rule: dead at the start of a silence >= 12 h (silence measured outside handling
                  periods) after which the fly never again has >= 100 changes in a 24-h window; the 24 h
                  before death are excluded. A >= 12-h silence followed by recovery is reported, not cut.
  Activity        1-s activity: second s is active if any position change occurs at a sample in [s, s+1).
  Scoring         Exactly sqw_ref.py (hmm_audit): 5-min sleep criterion with lab bridging, S if >= 30 s of
                  the minute are in a sleep bout, else Q if active fraction < 0.5 x r_f (fly reference per
                  light/dark block), else W. Fly-days with any excluded period are scored only on whole
                  remaining minutes; a sleep bout touching an excluded period is censored (not counted).
  Descriptives    Per minute: position changes; distance (sum |dpos|, beams); speed while moving
                  (beams per 0.4 s over samples with a change); time at 0 and at the ends (pos <= 1 or
                  >= 15); longest still run; EDGE FLICKER = single-sample 1-beam excursions (A-B-A, |B-A| = 1)
                  and the share of 1-s active seconds that contain nothing but edge flicker.
                  Summarised by S/Q/W and by light/dark.
Decision rule   None for (b): descriptive, no p-values. For (a) the rules above decide; their output is
                the inclusion table, reported in full with reasons.
Matched null    None (descriptive). Any later test proposed from these figures needs its own PRESPEC.
Notes           Edge flicker matters because a fly resting on a beam boundary produces it, and it counts as
                activity and breaks sleep under the 1-s rule. How much 1-s "activity" is edge flicker is the
                first question 0.4-s position can answer that 1-s counts cannot.

### 2026-11-12  D1 — discovery of candidate observable states in 0.4-s position (pilot)

Question        Does 0.4-s position contain discrete, persistent, reproducible behavioural states that the
                S/Q/W symbols merge or miss, as Q was once found inside wake? A candidate is promoted only if it
                passes all five criteria below; promotion makes it a HYPOTHESIS for the repeat experiment, with
                its definition frozen. No pilot result is a confirmed state. "No candidate passes" is a
                possible and reportable outcome.
Data            Pilot channels 17-32, fly-weeks INCLUDED under P1 (artifact-heavy fly-weeks reported separately,
                never pooled). Discovery uses weeks 1-2 (youngest, most flies alive). Weeks 3-9 are used only
                for a descriptive check of whether promoted candidates persist with age.
Preconditions   (i) Handling periods are settled (user's transfer times, or a new dated detector entry tested
                on weeks 2-9; the P1 detector failed on week 1). (ii) Week 2 is available. (iii) The synthetic
                gate below has PASSED. If any is unmet, D1 does not run on data.
Time window     Fly x weeks 1-2, ZT0-ZT0 fly-days, handling and death exclusions per P1. Light (ZT0-12) and
                dark (ZT12-24) are separate transition blocks inside every model.

Fixed choices
  Window        10 s (25 samples). Windows containing any excluded sample are missing (contribute 1 to
                the likelihood), never dropped (the sequence keeps its time base).
  Symbol        Joint categorical over 12 cells = activity class (4) x zone (3):
                activity  STILL  no position change;
                          FLICK  changes, all of them single-sample 1-beam excursions (A-B-A);
                          LOCAL  other changes, path length (sum |dpos|) <= 3 beams;
                          TRAVEL path length > 3 beams.
                zone      BEYOND1 position 0 for >= half the window (detection-limited zone);
                          END     otherwise, modal position <= 2 or >= 14 (beam 15 is dead/rare);
                          MIDDLE  otherwise.
                Positions after the P1 spike filter.
  Model         Per fly (flies are the unit), flat HMM on the 10-s symbols, K hidden states, categorical
                emissions shared across phases, transition matrices separate for light and dark (sqw_ref
                convention: the step INTO window t uses the phase of window t). Baum-Welch, 10 random restarts
                (seeds 1001-1010), max 500 iterations, tolerance 1e-4 nats; best training logL kept; all
                restart logLs reported.
  K range       2..8.
  K selection   Two-fold by day: fit on odd days, score even days, and vice versa. K* = argmax of the mean
                held-out logL per window over the two folds. BIC on all days reported beside it, not used.
  Degeneracy    Two states in one fit with emission Jensen-Shannon distance < 0.05 bits = DEGENERATE fit.
                K* is taken over non-degenerate fits only; degenerate fits are listed.
  Mapping       Every state of the K* fit gets: argmax activity class, argmax zone, and its S/Q/W composition
                (share of its posterior mass falling in S, Q, W minutes scored by P1). TYPE = (activity class,
                zone, dominant S/Q/W label). A type is a CANDIDATE if it is present in >= half of the
                included flies AND either shares its dominant S/Q/W label with another state of the same fit
                (it splits a label, as Q split wake) or has no S/Q/W label with >= 60% of its mass (it
                straddles labels).

Criteria (each evaluated per candidate; all five needed for promotion)
  C1 discrete   Per fly: held-out gain of K* over K*-1 (mean over folds, nats per window) exceeds the 95th
                percentile of its parametric-bootstrap null: 19 replicates simulated from the fly's own K*-1
                fit, same length, same missing pattern, same phase sequence, each passed through the SAME
                pipeline (two-fold fits at K*-1 and K*, 3 restarts, seeds 2001+). Candidate passes if C1 holds
                in >= 3 of its flies (one-sided binomial vs 0.05 at n = 8: P(X >= 3) = 0.0058).
  C2 persistent Median expected dwell 10 s / (1 - p_kk), averaged over light and dark, >= 30 s in >= half of
                the flies carrying the candidate.
  C3 reproduc.  Fit odd days and even days SEPARATELY at K*. Candidate present (same type) in both halves in
                >= half of its flies, and its JS distance between halves is smaller than its JS distance to
                every other state of the same fly.
  C4 not artif. Refit at K* with (a) FLICK merged into STILL and (b) BEYOND1 merged into END. A candidate
                that is no longer found under (a) or (b) is labelled DETECTION-DEFINED and not promoted
                (reported in full).
  C5 useful     Minute level, per fly, fit on odd days and score even days and vice versa: multinomial logit
                for the next S/Q/W symbol given the current symbol and phase, versus the same plus the fraction
                of the current minute in the candidate state (posterior). Gain = held-out logL difference.
                Passes if gain > 0 in a one-sided exact sign test p < 0.05 across its flies (at n = 8: >= 7).

Synthetic gate (must PASS before any data run; seeds 3001+, 20 synthetic flies per scenario, two weeks each,
generated as 0.4-s positions on a 17-position tube, then passed through the FULL pipeline incl. spike filter)
  G1 planted    Truth: 4 states (sleep-still, quiet-local, travel, and BOUNDARY REST = still with 1-beam flicker
                at rate 1/min, dwell 3 min), light/dark-specific rates. PASS iff the planted state is promoted
                (C1-C3 pass; C4 is expected to flag it as detection-defined under (a), which is the correct
                answer and is reported) in >= 16/20 flies, and K* = 4 in >= 14/20.
  G2 negative   Truth: 3 states, no extra state, but (i) per-fly rate heterogeneity (each rate x exp(N(0,0.3)))
                and (ii) a 6-h cosine modulation of the sleep-exit rate WITHIN each 12-h phase (amplitude 0.5)
                that the light/dark model cannot absorb. PASS iff a spurious candidate passes C1 in <= 2/20 flies
                and is promoted (C1-C3) in <= 1/20. This is the failure mode that sank PRESPEC 13 (heterogeneity)
                plus the non-stationarity a stationary-in-phase HMM can mistake for a state.
  If G1 or G2 fails, the failure is logged and the model is changed in a new dated entry before any data run.

Decision rule   Per candidate: promoted / not promoted, with the criterion it failed. Promoted candidates go to
                a new PRESPEC for the repeat experiment with their symbol definition, K and criteria frozen.
Matched null    C1: the fly's own K*-1 fit, simulated at its own length, missing pattern and phase sequence,
                through the same pipeline. G2: a synthetic population with no extra state.
Notes           - Window (10 s), symbol cells, zone cut-offs and all thresholds are fixed here; none may be
                  changed after the first data fit. Changes need a new entry labelled post hoc.
                - Week 1 has already been viewed (E0, P1); the candidates named in chat on 2026-11-12 (end-zone
                  rest, boundary rest, sub-criterion pauses, two kinds of walking, repeated end visits) were
                  suggested by it. D1 tests them only through the generic procedure above, never by a rule
                  written for one of them.
                - Sleep depth cannot be established from position (arousal threshold is the standard); a
                  promoted sleep sub-state is named by its observable signature, not as "deep" or "light".
                - Compute estimate: ~3 h for the K sweep, ~5 h for C1 bootstraps (8 flies); results are saved
                  per fly as they complete, so a restart resumes.

### 2026-11-12  D1-A — implementation details for D1 (written before any D1 code has run)

Purpose: D1 leaves some operational details open. They are fixed here, before the simulator or the HMM
code has been run once. Nothing in D1's criteria or thresholds changes.

1 TYPE (refines D1 "Mapping", because argmax activity is STILL for most quiet states and would merge them):
  still band   P(STILL) >= 0.9 QUIESCENT; 0.5-0.9 MIXED; < 0.5 ACTIVE  (P summed over zones)
  movement     the most probable of FLICK / LOCAL / TRAVEL (summed over zones); "-" if QUIESCENT
  zone         argmax zone (summed over activity classes)
  label        dominant S/Q/W share of the state's posterior mass
  TYPE = (still band, movement, zone, label). Candidate rule as in D1.
2 Sequences: each fly-day (ZT0-ZT0) is a separate sequence with its own initial distribution pi; held-out
  days are scored independently. Held-out logL per window = total / number of non-missing windows.
3 Random starts: seed = SeedSequence([base + r, fly index, K, fold code]); base 1001 (fits), 2001
  (bootstrap), 3001 (synthetic data). Init: emission rows = 0.5 x Dirichlet(1) + 0.5 x pooled symbol
  frequencies; transition rows: 0.9 on the diagonal, off-diagonal 0.1 x Dirichlet(1); pi uniform.
4 C3 matching: the state of the other half with the same TYPE and the smallest JS distance; if none, absent.
5 C4 "still found": some state of the refit has Jaccard overlap >= 0.5 with the candidate's windows
  (windows whose posterior argmax is the candidate).
6 Gate scope: in G1/G2 each synthetic fly is evaluated on its own (C1 per fly, no across-fly binomial);
  C5 is not computed in the gate. G1 "planted state promoted in a fly" = the fitted state with the largest
  Jaccard overlap with the TRUE planted-state windows (a window counts if >= half of it is in the planted
  state) is a candidate and passes C1, C2 and C3 in that fly. G2 "spurious candidate" = any candidate.
7 Synthetic truth (0.4-s hidden chain; positions 0..16; true position 15 reads as the previous reading,
  i.e. beam 15 dead; reflecting ends). Mean dwells (min) light / dark:
     SLEEP 20 / 40   QUIET 3 / 2   TRAVEL 0.67 / 0.5   BOUNDARY 3 / 3
  Jump probabilities (row = from): SLEEP -> Q .6 T .3 B .1 | QUIET -> S .4 T .4 B .2 |
     TRAVEL -> Q .6 S .1 B .3 | BOUNDARY -> S .5 Q .3 T .2   (G2: BOUNDARY removed, rows renormalised)
  Movement per 0.4 s: SLEEP none; QUIET +-1 with p .03; TRAVEL moves with p .5, keeps direction with
  p .85, step 2 with p .1 (else 1); BOUNDARY none, plus single-sample flicker to an adjacent position at
  rate 1/min (p = 0.4/60 per sample). G2 heterogeneity and 6-h modulation exactly as in D1.

### 2026-11-12  D1-A2 — gate evaluation fix, written after a code smoke test and BEFORE the gate runs

A smoke test of the gate code (toy settings: 4 days, K <= 4, 1 restart, 60 iterations, 2 bootstrap
replicates; 1 synthetic fly per scenario; not the gate) showed that D1-A item 6 counts the planted state as
"recovered" by whichever fitted state overlaps it most, even at Jaccard 0.06. The gate report therefore
adds a SECONDARY count: promoted AND Jaccard >= 0.5 (the threshold C4 already uses). The pre-specified
count is still computed and reported first. The same smoke test showed the HMM splitting sleep by ZONE
(quiescent states at END, MIDDLE and BEYOND1); that observation is logged in the ledger and is NOT acted on
here: the gate runs as specified.

### 2026-11-12  D1-B — behaviour separated from place (user decision; replaces D1's symbol and type)

Why: the D1 code smoke test (ledger 2f6d5a1) showed the HMM spending states on WHERE a resting fly sits
(persistent position makes zone look like hidden state), which crowds out behavioural sub-states. User
decision 2026-11-12: discover states from behaviour only; treat place as a separate layer.
The D1 gate run (as specified) was STOPPED at the user's decision before any synthetic fly finished; no
gate result for D1 exists or was seen. D1's gate is therefore neither passed nor failed: it is superseded.

Changes to D1 (everything not named here stays exactly as in D1, D1-A and D1-A2):
1 Symbol       Activity class only: STILL / FLICK / LOCAL / TRAVEL (M = 4), definitions unchanged.
2 TYPE         (still band, movement, S/Q/W label): D1-A item 1 without the zone.
3 Place layer  For every state of the K* fit, reported descriptively: share of its posterior mass in
               BEYOND1 / END / MIDDLE windows (zones as defined in D1), light and dark separately. Any
               claim about place within a state (e.g. "rest at the food end") needs its own PRESPEC.
4 C4           (a) unchanged: refit at K* with FLICK merged into STILL; candidate must still be found
               (Jaccard >= 0.5). (b) replaced: a candidate with >= 50% of its posterior mass in BEYOND1
               windows (where movement is unseen) is DETECTION-LIMITED and not promoted.
5 Gate         G1, G2 and their pass rules unchanged (truth unchanged, D1-A item 7); evaluation per D1-A
               item 6 with the D1-A2 secondary count (Jaccard >= 0.5) reported.

### 2026-11-12  P1-B — include all pilot flies except ch30; hub repair (user request; before any D1 data run)

Why: user request 2026-11-12 to analyse every channel 17-32 except 30 (dead). No D1/D1-B analysis has
touched real data, so the change cannot be result-driven. Disclosure: the repair rule below was designed
after seeing the week-1 hub structure (ledger, E0 follow-up); it is validated on synthetic data (R1) and,
like everything in the pilot, confirmed only on the repeat experiment.

1 Inclusion   Channels 17-29 and 31-32 (15 flies). P1's artifact flag is kept as a LABEL, not an exclusion.
2 Hub repair  Applied to every fly-week by the same rule, after the P1 spike filter:
              hub = the value receiving the most >= 4-beam jumps, IF it receives >= 30% of them AND >= 20% of
              the runs at it that were entered by such a jump end by returning to the pre-jump value.
              Every run at the hub entered by a >= 4-beam jump is replaced by the last value before it
              (the fly is taken not to have moved; if it exits to a different value, that exit is the move).
              A fly-week with no qualifying hub is unchanged. Hub, its two shares and the number of samples
              replaced are reported per fly-week.
3 Primary and sensitivity sets  Every pilot result is reported for ALL 15 flies (primary) and, beside it,
              for the 8 flies never labelled artifact-heavy (sensitivity). A conclusion that holds in one set
              only is reported as such.
4 R1 synthetic repair check (must PASS before repaired data are used; seeds 4001+). 10 synthetic G1 flies x
              7 days; a hub artifact is added to each: hub value uniform on 2..14, blocking episodes start at
              rate 1/min and last Exp(mean 1.5 s), during which the reading is the hub value. PASS iff, after
              spike filter + repair, compared with the same fly WITHOUT the artifact: (a) S/Q/W minute
              agreement >= 95% in every fly; (b) D1-B window-symbol agreement >= 95% in every fly; and (c) on
              the 10 artifact-free flies the repair changes <= 0.1% of samples. Agreement before repair is
              reported beside it.

### 2026-11-12  P1-C — hub beam declared dead (replaces P1-B item 2 after R1 FAILED; ledger)

Rule: hub identification exactly as P1-B (most frequent target of >= 4-beam jumps, share >= 30%, returns
>= 20%). If a fly-week qualifies, the hub beam is treated as DEAD, as the hardware already does for beam 15:
EVERY sample reading the hub value is replaced by the last value before it that is not the hub. A fly
genuinely at the hub position therefore appears at its previous position, and a crossing of the hub appears
as a 2-beam step. Cost, stated now: time genuinely spent at the hub beam and movement onto it are invisible
in that fly-week.
R1 is re-run unchanged (same synthetic flies, seeds, artifact model and pass criteria). If it fails again,
the 7 artifact-labelled flies stay out of pooled analyses (P1 rule) and are reported separately.

### 2026-11-12  P1-D — post-repair artifact rule (written after seeing the week-1 repair QC)

P1-B/C left open what happens to a fly-week that is still artifact-heavy after the repair. Rule: P1's
artifact threshold is re-applied AFTER the repair: a fly-week with > 5% of changes being >= 4-beam jumps
after spike filter + hub repair is excluded from pooled analyses and reported separately. Disclosure: written
after seeing that ch27, ch31 and ch32 stay above it in week 1; it reuses P1's threshold unchanged.
Sets: PRIMARY = fly-weeks passing this rule (week 1: ch18-26, 28, 29 = 12 flies); SENSITIVITY = flies never
labelled artifact-heavy (week 1: ch18, 19, 20, 22, 24, 26, 28, 29 = 8 flies). Every pilot result is reported
in both. Applied per fly-week, so a fly can enter or leave in later weeks.

### 2026-11-12  T1 — clock-dependent transitions (input-output HMM)

Adopted from beam_hmm/PRESPEC_T1_DRAFT.md (commit 3eb0231) with the user's decisions of 2026-11-12:
NO startle terms (the light/dark step only; M4 dropped), ridge lambda = 0.01, and the pilot's across-fly
thresholds accepted (only large, consistent clock effects are detectable). Everything in the draft holds
except where this entry is more specific; where they differ, this entry wins.

Question, data, time window, preconditions, model family M0-M3, estimation, model selection, decision rule,
matched null and exclusions: as in the draft. Fixed here in addition:
1 Covariates  z(ZT) at window-of-day w (ZT = 10 w s): M0 [1, dark]; M1 + cos, sin(2 pi ZT/24);
              M2 + cos, sin(4 pi ZT/24); M3 + cos, sin(6 pi ZT/24), cos, sin(8 pi ZT/24). dark = ZT >= 12.
              Ridge 0.01 on every coefficient except the intercept and the dark step.
2 C6 (new D1 criterion, entered here as D1-C): as in the draft.
3 Gate truth  0.4-s simulator, states SLEEP, QUIET, TRAVEL with D1-A item 7 dwells and jumps (BOUNDARY
              removed, rows renormalised). Clock modulation (T-G1, T-G3) multiplies the SLEEP exit
              probability by exp(f(ZT)), f = a [cos(th) + cos(2 th)], th = 2 pi (ZT - 12.5)/24,
              a = ln 3 / 3.125 = 0.35156 (peak-to-trough of the modulation = 3; main peak ZT 12.5, secondary
              ZT 0.5). T-G2: no modulation. T-G3: per-fly heterogeneity as D1 G2 (each rate x exp(N(0,0.3))).
4 Gate criteria, made exact
  T-G1 peak     the SLEEP state of the fit = the state with the largest Jaccard overlap with true-SLEEP
                windows (>= half the window in SLEEP). Its fitted exit profile 1 - P_ii(ZT) is compared with
                the TRUE total exit profile (light/dark base x modulation, per 10-s window). PASS needs the
                circular distance between their argmax ZTs <= 2 h in >= 16/20 flies, the median ratio of
                fitted to true peak-to-trough in [0.67, 1.5], and M1-M3 chosen over M0 in >= 16/20.
  T-G2          M1-M3 chosen over M0 in <= 3/20.
  T-G3          the D1-B procedure run with M* transitions (K sweep, C1 bootstrap from the M* model at
                K*-1, C2, C3): a candidate passing C1-C3 in <= 1/20. The same procedure under M0 is run on
                the same flies and reported beside it.
  Seeds: data SeedSequence([6001, scenario, fly]); fits SeedSequence([5001 + r, fly, K, model, fold]).
5 Order       The T1 gate starts after the D1-B gate finishes (they would compete for the same 4 cores).

### 2026-11-13  D1-D — sample-level (0.4-s) behavioural HMM; null simulated at the level it is fitted

Why: the D1-B gate FAILED (run on the user's machine; numbers to be logged in the ledger from the pushed
files): G1 chose K* = 4 in only 7/20 flies (most 5-6) although the planted state was promoted in 17/20;
G2 produced a spurious promoted state in 18/20. G1 has neither clock modulation nor heterogeneity, so its
extra states point to MISSPECIFICATION by windowing: a 0.4-s process switching state inside 10-s windows
creates mixed windows that a window-level HMM absorbs with transitional states, while the C1 null
(simulated directly at the window level) has no mixed windows and cannot see the effect. User decision
2026-11-13: fit at the sample level instead. This replaces D1-B's symbol and procedure; D1, D1-A (except
item 3 seeds, which carry over unchanged), D1-A2, D1-B items 3-5 and P1-B/C/D hold unless changed here.

1 Observation  Per 0.4-s step after spike filter + P1-C hub repair: symbol = min(|x_t - x_{t-1}|, 3), i.e.
               M = 4 classes {no change, 1-beam, 2-beam, >= 3-beam}. A purely local function of two
               consecutive readings: no windows, no aggregation. (Beam 15 being dead turns a 14<->16
               crossing into a 2-beam step in data and simulator alike.) Missing samples = -1.
2 Model        Per fly, categorical HMM on the step symbols, transition matrices separate for light and dark
               (step INTO sample t uses the phase of t), each fly-day a separate sequence. K = 2..8.
3 Fitting      Random starts as D1-A item 3. Convergence: tolerance 0.01 nats (absolute, per fit), max 300
               iterations. Observed-data fits (K sweep, both folds): 5 random starts. Full fit at K*,
               K*-1 null, half fits, C4 refit: 5 random starts. Bootstrap fits: 2 starts each, one random and
               one warm start (K*-1 fits start from the generating null; K* fits start from the fly's
               observed K* fit). A warm start can only add a candidate optimum.
               Reason for the reduction from D1's 10 restarts / 1e-4: measured cost (0.10-0.24 s per EM
               iteration on 7 days; ~200 core-hours for the gate at D1 settings).
4 K selection  Two-fold by day, held-out logL per step; K* = argmax over non-degenerate fits (JS < 0.05 bits
               between two emission rows = degenerate). BIC reported beside it for K*.
5 TYPE         still band on P(no change): QUIESCENT >= 0.999; LOW 0.95-0.999; ACTIVE < 0.95. Movement = most
               probable of {1-beam, 2-beam, >= 3-beam} ("-" if QUIESCENT). Label = dominant S/Q/W share of
               the state's posterior mass (each sample takes its minute's P1 label). Candidate rule as D1.
6 Criteria     C1 as D1 but the 19 bootstrap replicates are simulated from the fly's K*-1 SAMPLE-LEVEL fit with
               the same sequence lengths, phases and missing pattern, i.e. at the same level as the fitted
               data. C2: expected dwell 0.4 s / (1 - p_kk), mean of light and dark, >= 30 s. C3 as D1 (halves
               = the two fold fits). C4(a) replaced: positions with single-sample 1-beam excursions removed
               (A-B-A, |B-A| = 1, B set to A), recoded, refit at K*; candidate still found (Jaccard >= 0.5
               over samples). C4(b) as D1-B (>= 50% of posterior mass at position 0 = detection-limited).
               C5 as D1 (not computed in the gate).
7 Gate         G1 and G2 truths, seeds and pass rules exactly as D1 / D1-A item 7 (G1: planted state promoted
               in >= 16/20 and K* = 4 in >= 14/20; G2: spurious C1 pass <= 2/20 and promoted <= 1/20), with the
               D1-A2 secondary count (Jaccard >= 0.5, over samples). Remaining known gap, stated now: the
               spike filter is applied to the data but not re-applied inside bootstrap replicates (they are
               generated as symbols); the number of samples it changes is reported per fly.
8 If the gate fails again, next step (stated in advance, not run without a new entry): a minimum held-out
               gain per added state, calibrated on simulated flies, before any change to the model family.

### 2026-11-13  D1-D2 — gate evaluation corrected (after the D1-B gate results; BEFORE the D1-D gate runs)

The D1-B results (ledger) show that D1-A item 6's G2 rule ("spurious candidate = any candidate") counts TRUE
states as spurious whenever the S/Q/W scoring merges them (true QUIET and TRAVEL are both scored W). A
negative control must count only states that do not exist. Disclosure: this correction was made after seeing
the D1-B gate; it changes how the gate is SCORED, not the model, the data, the criteria or the pass counts.
For D1-D (G1 and G2):
- Each fitted state gets its best Jaccard overlap (over samples, posterior argmax) with each TRUE state.
- A state is SPURIOUS if its best Jaccard with every true state is < 0.5.
- G2 PASS iff a SPURIOUS candidate passes C1 in <= 2/20 flies and is promoted (C1-C3) in <= 1/20.
- G1 unchanged (planted state promoted in >= 16/20, K* = 4 in >= 14/20), plus the same spurious count reported.
- The original D1-A item 6 count is reported beside the corrected one.
The T1 gate's TG3 is already running under the old rule; its TG3 verdict is therefore to be read as
uninformative, and TG3 will be re-scored the same way only if re-run (TG1 and TG2 are unaffected).

### 2026-11-13  T1-A — TG3 deferred (before any T1 gate result exists)

The T1 gate was stopped at the user's decision after ~1.3 h with 0/60 flies finished (no T1 result exists or
was seen). TG3 is deferred: as specified it reruns D1 discovery with the D1-B window model, which FAILED its
own gate, and scores spurious states with the D1-A item 6 rule that D1-D2 corrected. TG3 will be re-specified
on the D1-D sample-level model with D1-D2 scoring, in a new entry, only if D1-D passes its gate.
TG1 and TG2 run exactly as specified in T1 (truth, seeds, criteria unchanged). gate_t1.py gains an option to
run named scenarios only; the T1 gate verdict is then reported for TG1 and TG2 alone.
