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
