# PRESPEC 21 — is the wake bout time-asymmetric at 1 s beyond its activity level? (formerly draft 20y)
Date 2026-11-01. Drafted by cClaudeII; reviewed by cClaude (open points answered, amendments below folded in).
NOTHING runs until Sheyum gives the go-ahead. No data from this analysis has been computed.

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
then reversing the WHOLE record maps every bout's start onto its end: the exact symmetry is a GLOBAL flip of all bouts
together. A per-bout sign flip of d = x_start - x_end additionally needs the bouts to be independent given |d_b|, and
20k/C7 shows that bouts carry memory. The per-bout flip is therefore the primary null only TOGETHER with the
whole-day block null below; neither alone is called exact.
What this needs, and what it does not (cClaude's review):
- It does NOT need stationarity WITHIN a bout. The rise and fall of activity inside a bout is exactly what is measured,
  and a reversible process can have a non-flat profile, as long as start and end mirror each other.
- It DOES need stationarity of the whole record plus reversibility. That holds for per0 DD; for CS DD the circadian
  drive is the stated caveat.
- The weak point is INDEPENDENCE OF BOUTS, which the per-bout flip assumes (20k/C7: bouts carry memory). Hence the
  per-day block sensitivity null below.

## Data and blocks
- PRIMARY (Holm over 2): per0 DD (25 cut files) and CS DD (25 cut files), whole retained records. per0 DD is
  undriven; CS DD carries the circadian drive, which a positive cannot be separated from.
- DESCRIPTIVE: CS, per0LD, fumin, sss (LD; light switches break stationarity).
- Bouts: wake bouts between consecutive sleep bouts with length >= 2L = 600 s, so the two windows do not overlap
  (the same 600-s threshold as sqw_ref.wake_bout_events). Flies with fewer than 5 such bouts are excluded and counted.

## Readouts (per fly)
PRIMARY setting: L = 300 s in 30 bins of 10 s (bouts >= 600 s). SECONDARY setting, also pre-declared and reported in
full, outside the Holm family: L = 120 s in 24 bins of 5 s (bouts >= 240 s, the same no-overlap rule). Bin 1 = the bin
next to the sleep boundary. Activity = 1 if count > 0 (as the scoring uses). Wake bouts used per fly are reported per block
(median, range, and the per-fly list) for both settings, because S grows with the number of bouts, in the null and in
the data alike.
- PRIMARY, shape: per bout, d_b = x_start bin means - x_end bin means (30 values); centre d_b by its own mean over the
  30 bins (this removes the level difference, i.e. the C2 component). S = sum over bins of (mean over bouts of
  centred d_b)^2.
- Null (primary): 2000 per-bout sign flips (seed default_rng(27_000 + 1000 g + fly)), S* recomputed;
  p_fly = (1 + #{S* >= S}) / 2001; z_fly = (S - mean S*) / sd S*.
- Sensitivity null, reported beside: per-DAY block sign flips. All bouts whose midpoint falls in the same 24-h day of
  the record are flipped together. All 2^(number of days) configurations are enumerated exactly (16 for a 4-day
  record), giving p_fly_block and z_fly_block. This guards against dependence between bouts within a day.
- SECONDARY (outside the Holm family):
  - (a) level, mean over bouts and the first 60 s of (x_start - x_end), with the same sign-flip null: the C2
    consistency check;
  - (b) the shape statistic on count intensity (counts per second, not active/inactive);
  - (c) descriptive pooled profiles m_start(tau) and m_end(tau) per genotype, with fly-bootstrap bands (the figure).

## Decision rule (fixed now)
Per primary block, "TIME-ASYMMETRIC SHAPE" if the group Wilcoxon of z_fly against 0 has Holm p < 0.05 with
median z > 0, AND the fraction of flies with p_fly < 0.05 exceeds Binomial(n, 0.05) at p < 0.05. Otherwise
"no shape asymmetry detected". The same group rule is applied to the block null (Holm over the same 2 blocks) and
reported. If the two nulls DISAGREE (one positive, the other not), the per-bout result is NOT reported as clean: the
verdict reads "not robust to dependence between bouts within a day", and no shape asymmetry is claimed. The secondary L = 120 s setting is reported with the same rule and does not
change the verdict.
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

## Resolved review points (cClaude, 2026-11-01; two notes)
0. Renamed from draft 20y because 20y and 20z are taken on the shared branch (the Kir2.1 circuit-silencing test, cf777b1,
   and its follow-up).
1. L and bins: L = 300 s / 10-s bins primary; L = 120 s / 5-s bins secondary; both pre-declared, both reported.
2. Centring: yes (without it the shape statistic duplicates C2).
3. Within-bout stationarity: not required (see the symmetry argument); whole-record stationarity is, so per0 DD is
   the clean block and CS DD carries the circadian caveat. Bout dependence is handled by the whole-day block null
   (one block per 24-h day: 4 blocks for a 4-day record, 3 for a 3-day record); disagreement between the nulls means
   no clean claim.
4. Exclusions: flies with < 5 qualifying bouts are excluded and counted; bouts per fly per block are reported.

## Conditional follow-up B (NOT part of this test; its own PRESPEC, written only if PRESPEC 21 detects a shape asymmetry)
Refit 20u's models WITHOUT the reversibility constraint: same orders, starts, seeds, adequacy gate and time box.
Limits fixed now:
- a pass where the reversible model failed supports only "within K <= 6, the bout process needs time-irreversible
  hidden dynamics" (the class is named);
- BIC is reported beside the gate, so a pass is not simply extra flexibility (more free parameters at the same K);
- if it also fails, the verdict is "no model in either class reproduces the bouts", with no further enlargement;
- B runs only after PRESPEC 21's verdict, so it is not chosen in light of anything else.

## Amendment A (2026-11-01, after Sheyum's go-ahead, BEFORE any data run; found while implementing)
Whole-DAY blocks give 2^4 = 16 sign configurations for a 4-day record (2^3 = 8 for 3 days). The smallest attainable
p_fly_block is then 1/16 = 0.0625 (1/8 for 3 days), so the binomial-excess half of the decision rule (flies with
p < 0.05) could NEVER pass under the block null. Any per-bout positive would "disagree" by construction.
The block is therefore the 12-h HALF-DAY, within the "8 or 4 blocks per fly" cClaude's review allowed: blocks =
floor(bout midpoint / 43200 s). That gives 8 blocks and 256 configurations for a 4-day record (minimum p 1/256), and
6 blocks and 64 configurations for 3 days (1/64). Configurations are enumerated exactly, identity included;
p_fly_block = #{S* >= S} / (number of configurations), z_fly_block from the same set. Nothing else changes.

## Amendment B (2026-11-01, BEFORE any data run; cClaude's fix, relayed by Sheyum) — SUPERSEDES amendment A
cClaude independently found the same flaw (whole-day blocks: p_fly_block >= 1/16, so the binomial arm cannot pass) and
fixed it differently. This rule replaces amendment A; blocks are back to whole 24-h days (floor(midpoint / 86400 s)).
(a) Per-fly block null: all 2^(days) configurations enumerated (16 for 4 days, 8 for 3), z_fly_block from them. The
    block-null rule uses ONLY the group Wilcoxon of z_fly_block against 0 (Holm over the 2 primary blocks, median > 0).
    The binomial arm is DROPPED for this null, because with 8-16 configurations p_fly_block cannot fall below 0.0625.
(b) Group-level block null: across all flies of a genotype, each fly-day is flipped independently (e.g. 25 x 4 = 100
    blocks); 2000 random configurations (seed default_rng(30_000 + genotype index)); statistic = the group sum of the
    flies' S; p = (1 + #{sum S* >= sum S}) / 2001. This is the second robustness check.
DISAGREE now means: the per-bout rule passes AND either (a) (Holm p >= 0.05 or median <= 0) or (b) (p >= 0.05) fails.
A disagreement is reported as "not robust to dependence between bouts within a day", with no claim.
The per-bout rule is unchanged (Wilcoxon of z_fly with Holm, median > 0, binomial excess).
