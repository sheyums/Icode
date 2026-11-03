# PRESPEC 25 — when a fly measurably cycles within a few hours, which way does it turn?

Date 2026-11-03. This is an exception to the analysis freeze, which Sheyum approved: a new question, not a re-test.
Written and committed BEFORE any window-level count has been computed. Already known: 12-h and pooled results (C2,
C10, 20c). Those pool over days or over flies, so they cannot see backward EPISODES within a fly that are outweighed
by forward ones. Numbering: 23 and 24 are reserved for Sh^mns and the extra silencing lines.

## Question
Within single flies, in single short windows (one window on one day): when the cycle is measurable, is it always
forward, or do flies run it backward for hours at a time more often than a constant forward process would produce?

## Data
- Genotypes. PRIMARY: per0 LD, per0 DD. REFERENCE: CS LD, CS DD. DESCRIPTIVE: sss (LD).
- Symbols as 20s-A: sqw15.seqs, theta 0.5, whole days. Records are 4 days per fly (sss 2), with no missing minutes;
  the script asserts this.

## Windows (fixed now; both sizes reported; no other size, offset or merging)
- 4 h (6 per day) and 6 h (4 per day), contiguous, on single days, starting at ZT0. In DD they start at the record
  start, which is clock time from the first DD day.
- In a window, count N_SW (sleep -> active wake) and N_WS (active wake -> sleep) over consecutive minutes that both
  lie inside the window. Net forward count n = N_SW - N_WS. This is the edge used by the count J.
- MEASURABLE: two-sided exact binomial test of N_SW against N_WS (p = 0.5) gives p < 0.05. This needs at least 6
  such transitions. The rule does not depend on the sign, so if there were no preferred direction, measurable
  windows would split 50:50 between forward and backward.
- A measurable window is FORWARD if n > 0 and BACKWARD if n < 0.

## Null (sample-size matched; every step identical for real and simulated windows)
- PRIMARY NULL, clock-matched. For each fly, a first-order chain for each window POSITION in the day (e.g. ZT4-8),
  pooled over that fly's days: transition counts normalized per row. With no missing minutes this equals mc_fit. A
  row with no counts takes the fly's 12-h phase row. This null keeps each fly's clock-time changes in the cycle and
  removes only the day-to-day variation within the same clock window.
- SECONDARY NULL: the fly's 12-h phase chain (the chain of every other analysis).
- Simulation: 500 replicates per fly and window size. Each replicate is a continuous minute-level record of the same
  length, starting at the observed first state, with the chain switching at window boundaries. Windows, counts and
  the measurable rule are applied exactly as to the data.

## Statistics (flies are the unit; windows are nested in flies)
- PRIMARY, per genotype x window size: B = the total number of measurable BACKWARD windows over the genotype's
  flies. Monte Carlo p = (1 + #{replicate totals >= B}) / 501, one-sided (an excess of backward windows). Holm over
  the 4 primary tests (per0 LD and per0 DD x 4 h and 6 h). The reference genotypes get their own Holm over 4. sss is
  reported with uncorrected p, descriptively.
- Reported alongside: F (measurable forward windows) with its own one-sided excess p; the backward fraction
  B/(B+F) against the null's; the per-fly excess (observed minus null-mean B), with a Wilcoxon test as secondary; the
  number of flies with at least one measurable backward window, against the null.

## Decision rule (fixed now)
- "Backward EPISODES within flies": primary p_B Holm < 0.05 AND the observed backward fraction B/(B+F) above the
  null's median fraction. Both are needed. Day-to-day variation in the cycle's SIZE alone, with a constant direction,
  adds measurable windows in the fly's own direction, so it raises F and leaves B/(B+F) at or below the null's.
- An excess of B with a backward fraction at or below the null is reported as "day-to-day variation of the cycle's
  size", not as reversal.
- No excess of B: "when per0 flies measurably cycle within hours, they turn forward, as often as a constant forward
  process predicts". This is direction stability at the window scale.
- Interpretation limit, stated with any positive result: an excess over the clock-matched null means variation
  BETWEEN DAYS at the same clock time (episodes). It is not separable from slow drift of the fly's state across days.
- sss: forward and backward are defined as for every other genotype. The descriptive question is whether its
  daytime reversal is a stable backward drive or a mix of episodes.

## Figures (descriptive)
1. A per-fly raster: rows are flies, columns are windows in time order, colour = forward / backward / not
   measurable. Per genotype, 6-h and 4-h.
2. Observed counts of measurable forward and backward windows against the null distributions, per genotype x size.

## Output
`sqw25.py` -> `sqw_out/s25_windows.json`, `figures/fig25_*.png`. Ledger entry with commit hashes.
