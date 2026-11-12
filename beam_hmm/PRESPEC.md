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
