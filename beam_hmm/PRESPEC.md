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
