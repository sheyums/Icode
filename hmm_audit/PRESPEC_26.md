# PRESPEC 26 — does the direction of successive sleep-wake loops persist, alternate, or behave as the chain predicts?

Date 2026-11-03. This is an exception to the analysis freeze, which Sheyum approved. Written and committed BEFORE any
loop-sequence statistic has been computed. Already known: per fly x phase loop RATES and durations (cycle_times.py,
descriptive), and that 4-6-h windows cannot resolve direction (PRESPEC 25). Not yet looked at: the ORDER of loop
directions. Numbering: 23 and 24 are reserved for Sh^mns and the extra silencing lines.

## Question
Each fly completes a sequence of loops, each forward (W->Q->S->W) or backward. Is the direction of one loop
independent of the previous one, given the fly's own chain? Or do same-direction loops cluster (EPISODES: a drive
that switches between forward and backward modes)? Or do they alternate (a RESTORING tendency)? The comparison of
clockless per0 with CS separates intrinsic switching from clock-driven changes of the bias.

## Loops (identical for data and null)
- Loop erasure on the minute symbols, run-compressed (cycle_times.loops). A loop erased with 3 distinct states is a
  completed cycle: forward if it is a rotation of W->Q->S->W, else backward. Loops are taken in order of completion
  over the whole record.
- Consecutive pairs: successive loops in the same uninterrupted stretch of scored minutes. A missing minute (-1,
  sss only) resets the path, and no pair spans it.
- A fly enters a test only if it has >= 10 consecutive pairs.

## Statistic
Per fly: S = the fraction of consecutive loop pairs with the SAME direction. Excess E = S_obs - mean(S_null).

## Null (sample-size matched; every step identical for real and simulated data)
- PRIMARY: a clock-matched first-order chain, one transition matrix per 4-h clock window (6 per day), pooled over
  the fly's days. Counts are taken on pairs of scored minutes, normalized per row; a row with no counts takes the
  fly's 12-h phase row. This keeps clock-time changes of the bias (including day/night), so an excess is NOT the
  known daily rhythm.
- SECONDARY: the fly's 12-h phase chain.
- 500 simulated records per fly. Each is the same length, starts at the first scored state, switches matrices at
  window boundaries, and gets the fly's missing minutes blanked. Loops, pairs and S are computed exactly as for the
  data.
- Descriptive only: a within-fly shuffle of the observed loop order (2000 shuffles). It shows how much the chain
  itself departs from independent loops.

## Tests (flies are the unit)
- PRIMARY: per0 DD and per0 LD. Wilcoxon signed-rank of E against 0, two-sided, Holm over the 2.
- REFERENCE: CS LD and CS DD, the same test, Holm over the 2.
- SECONDARY: per0 DD against CS DD (like for like, no light), Mann-Whitney on E, two-sided. Then per0 LD against
  CS LD.
- DESCRIPTIVE: fumin and sss, plus the distribution of same-direction run lengths against the null, per genotype.

## Decision rule (fixed now)
- EPISODES (switching drive): Holm < 0.05 with median E > 0.
- RESTORING tendency: Holm < 0.05 with median E < 0.
- Otherwise: "successive loop directions are as the observed chain predicts", a constant bias at the loop scale.
- An effect in CS but not in per0 is read only through the secondary per0-vs-CS test, never as a difference of
  significance.
- Interpretation limit, stated with any positive result: an excess over a clock-matched chain is memory beyond the
  first-order minute chain (of the same kind as C7 and C9*). It is not separable from slower day-to-day drift.
- No change of the loop rule, pair rule, window size of the null or genotype set after looking.

## Output
`sqw26.py` -> `sqw_out/s26_loopseq.json`, `figures/fig26_*.png`. Ledger entry with commit hashes.
