# PRESPEC 20y — DRAFT, NOT APPROVED: is the wake bout time-asymmetric at 1 s beyond its activity level?
Date 2026-11-01. Drafted by cClaudeII at Sheyum's request for review by Sheyum and cClaude. NOTHING runs until it is
approved and committed as PRESPEC_20y.md. No data from this analysis has been computed.

## Motivation
- 20v: Q and W dwell times depend on the direction of passage beyond the observed first-order chain (~0.08 nats/h,
  per0 DD and CS DD).
- 20v amendment F, F1: the scoring is exactly time-symmetric, so the 5-min rule cannot create that from reversible
  dynamics.
- 20u: reversible hidden models of the 1-s activity over-fragment wake by 22.5%.
- All three point to structure at the seconds scale, which the minute symbols coarse-grain. This test looks at the
  1-s counts directly, with no hidden model and therefore no identifiability problem (Fredkin & Rice 1986).

## What is NOT new, stated before running
A difference in activity LEVEL between the start and the end of a wake bout is essentially C2. c_obs (20d) is P(last
wake minute before onset is Q) - P(first wake minute after waking is Q), and Q is a minute's active fraction below
theta r_f. "Quieter before sleep than after waking" IS the forward cycle W -> Q -> S -> W. The level comparison is
therefore reported only as a consistency check (secondary). The primary readout is the SHAPE of the 1-s trajectory with
the level removed: whether activity winds down gradually into sleep and rises abruptly out of it (or the reverse),
which the minute chain cannot express.

## The symmetry argument (why the null is exact)
Take a wake bout [w0, w1) between two consecutive sleep bouts (sqw_ref.sleep_bouts, bridged, >= 300 s). Its start
window is x_start(tau) = activity at second w0 + tau - 1; its end window is x_end(tau) = activity at second w1 - tau,
tau = 1..L. If the 1-s process is stationary and time-reversible, and bout detection is time-symmetric (F1 showed
sqw_symbols, which contains sleep_bouts and bridging, maps reversed traces to reversed symbols in every minute),
then x_start and x_end are EXCHANGEABLE within each bout. Reversing the trace maps a bout's start onto its end. A
per-bout sign flip of d = x_start - x_end is then an exact null, with sample size matched by construction.

## Data and blocks
- PRIMARY (Holm over 2): per0 DD (25 cut files) and CS DD (25 cut files), whole retained records. per0 DD is
  undriven; CS DD carries the circadian drive, which a positive cannot be separated from.
- DESCRIPTIVE: CS, per0LD, fumin, sss (LD; light switches break stationarity).
- Bouts: wake bouts between consecutive sleep bouts with length >= 2L = 600 s, so the two windows do not overlap
  (the same 600-s threshold as sqw_ref.wake_bout_events). Flies with fewer than 5 such bouts are excluded and counted.

## Readouts (per fly)
L = 300 s, binned into 30 bins of 10 s, tau-bin 1 = the 10 s next to the sleep boundary. Activity = 1 if count > 0
(as the scoring uses).
- PRIMARY, shape: per bout, d_b = x_start bin means - x_end bin means (30 values); centre d_b by its own mean over the
  30 bins (this removes the level difference, i.e. the C2 component). S = sum over bins of (mean over bouts of
  centred d_b)^2.
- Null: 2000 per-bout sign flips (seed default_rng(27_000 + 1000 g + fly)), S* recomputed;
  p_fly = (1 + #{S* >= S}) / 2001; z_fly = (S - mean S*) / sd S*.
- SECONDARY (outside the Holm family):
  - (a) level, mean over bouts and the first 60 s of (x_start - x_end), with the same sign-flip null: the C2
    consistency check;
  - (b) the shape statistic on count intensity (counts per second, not active/inactive);
  - (c) descriptive pooled profiles m_start(tau) and m_end(tau) per genotype, with fly-bootstrap bands (the figure).

## Decision rule (fixed now)
Per primary block, "TIME-ASYMMETRIC SHAPE" if the group Wilcoxon of z_fly against 0 has Holm p < 0.05 with
median z > 0, AND the fraction of flies with p_fly < 0.05 exceeds Binomial(n, 0.05) at p < 0.05. Otherwise
"no shape asymmetry detected".
Wording if detected (per0 DD): "Within wake bouts, the 1-s activity trajectory into sleep is not the time-reverse of the
trajectory out of sleep, beyond their difference in level; under stationarity and the verified time-symmetric
scoring, this is time-irreversibility of the activity process at the seconds scale." The direction (gradual into
sleep, abrupt out, or the reverse) is read from the pooled profiles as description.
CS DD adds "the circadian drive is not excluded as its source". Never "proves broken detailed balance": stationarity
is assumed, not tested within bouts.
Wording if not detected: "No seconds-scale shape asymmetry of wake bouts beyond their level difference was detected."
Then the 20v excess has no 1-s trajectory counterpart in this statistic.

## Narrowest block, matched null, time box
Narrowest block: one fly's whole DD record (the bouts it contains). Matched null: the same bouts, sign-flipped.
Expected run time: minutes. One run, one verdict.

## Open points for the reviewers
1. L = 300 s and 10-s bins are proposed, not derived. Fix them, or set a small pre-declared set reported in full
   (rule 4: not chosen after looking).
2. Centring removes the mean level over 300 s. If the reviewers consider the 61-300 s level difference part of
   "shape", drop the centring and accept that the primary then overlaps C2.
3. Within-bout stationarity: long wake bouts in DD could carry slow drifts (circadian in CS DD) that make the start
   and end windows differ for reasons other than irreversible fast dynamics. The per0 DD result is the clean one.
