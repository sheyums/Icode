# PRESPEC 18 — user-requested variants of figures 1-3 (written after seeing PRESPEC 15-17 results)
Date 2026-10-22. Two changes requested by the user; both are POST HOC with respect to the E2/PRESPEC 17
analyses and are reported ALONGSIDE the pre-registered results, never replacing them.
1. Sleep scoring with 1 s bridging (lab convention): on the binary active/inactive trace, every active run of
   exactly 1 s flanked by inactive seconds is relabelled inactive; THEN sleep = inactive run >= 300 s
   (interpretation (a); user may revise). Activity itself is still measured from the raw counts.
2. Window narrowed to +-5 min: seconds 61-300 after waking and before sleep onset (minute 1 excluded as
   before), profiles in 30 s bins (8 bins). Chosen after seeing that the asymmetry concentrates near the
   transition, so significance on this window is DESCRIPTIVE, not confirmatory.
Kept fixed: wake bouts >= 20 min (same eligibility as E2), flies with >= 3 eligible bouts, blocks as
PRESPEC 16, mid-wake reference = middle 5 min of the bout.
Computed for three versions, all reported:
 V0 pre-registered: no bridging, minutes 2-10 (= existing results, recomputed by the same code as a check);
 V1 bridging, minutes 2-10;  V2 bridging, seconds 61-300 (the requested figures).
Per version and block: per-fly Delta (median, % flies > 0, sign test), E2 permutation p (2000, per-bout swaps),
quieting before sleep x and after waking y (mid - mean over the window), mean +- s.e.m. across flies.
Figures 1-3 regenerated for V2 (files *_pm5_bridged.*); the V0 figures are kept.

## Amendment (user, during the run): bridging convention is (b), not (a)
A 1 s active run is bridged ONLY when both flanking inactive runs are >= 300 s (each already a sleep bout);
e.g. 301 s + 1 s + 300 s -> one 602 s sleep bout. Results of the option-(a) run are kept in v18_optionA.log and
NOT used. Under (b) the start and end of every wake bout >= 20 min are unchanged (only 1 s wake bouts between
two sleep bouts disappear), so V1 is expected to equal V0 exactly; this is checked, not assumed.
V1b = bridging (b), minutes 2-10; V2b = bridging (b), seconds 61-300 (the requested figures).
