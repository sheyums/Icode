# PRESPEC 27 — through which state's exits does the clock modulate the cycle?

Date 2026-11-06. This is an exception to the analysis freeze, which Sheyum approved: a cleaner reformulation of the
old Fig. 3A-C question ("which state's transition-probability oscillations are the primary source"). Written and
committed BEFORE any decomposition has been computed. Already known: C10, the pooled 2-h clock-bin cycle flux and its
24-h cosine amplitude (CS LD 0.077, CS DD 0.063 cycles/h, per0 DD 0.006). Not yet computed: any per-row or per-state
counterfactual. Numbering: 23 Sh^mns, 24 the 23E10 silencing, 25/26 used.

## Why not repeat the old analysis
In the old 4-h-window HMM, constraining one hidden state's rows changed the others through the refit (degenerate H1).
In the current HMM the hidden states are the observed symbols (emissions 80-100% diagonal), so a refit would barely
move. The question is therefore asked of the observed chain directly, as a counterfactual on fitted matrices, with no
refitting.

## Data and quantities (identical to C10 / PRESPEC 20w-B)
- Symbols: sqw_ref.sqw_symbols(read_counts(file), theta 0.5, light_dark per genotype), whole days. Files are those
  of C10.
- Per genotype: minute transitions pooled over flies and days by arrival minute of day, summed into twelve 2-h clock
  bins (ZT in LD; clock time from the start of each day in DD). P_b is the row-normalised pooled matrix of bin b.
- Bin flux J_b = the stationary cycle flux of P_b (C10's J_bins, cycles/h). Rhythm amplitude A = the amplitude of the
  24-h cosine fitted to the 12 J_b (C10's cosfit).
- GATE: the full model must reproduce C10's A to within 1e-3 cycles/h for CS LD and CS DD (it is the same quantity).
  If it does not, stop and report.

## Counterfactuals (the decomposition)
For each state i in {Q, W, S}:
- FLATTEN i: replace row i of every P_b by the row of the day-pooled matrix (all 12 bins' counts summed), keeping the
  other rows bin-specific. A_flat(i) = amplitude of the resulting J_b. R_i = 1 - A_flat(i)/A = the share of the rhythm
  REMOVED when state i's exits stop following the clock.
- ONLY i: flatten the two OTHER rows, keep row i bin-specific. O_i = A_only(i)/A = the share of the rhythm that state
  i's exits produce on their own.
R and O need not sum to 1, because effects interact and can cancel in phase. Both are reported, with the peak time of
each counterfactual rhythm.

## Uncertainty
Fly bootstrap: resample flies with replacement (1000 replicates, seed 27), recompute everything from the pooled
resampled counts. 95% percentile intervals for A, R_i, O_i. Replicates where A < 0.01 cycles/h are kept, and their
number is reported.

## Genotypes
- PRIMARY: CS DD (the clock without light).
- SECONDARY: CS LD (clock plus light; light-driven rhythm is not separable here).
- DESCRIPTIVE: per0 LD (light without clock; C10 not detected, Holm 0.059).
- CONTROL: per0 DD (no detectable rhythm, A 0.006). Its decomposition is reported, but NOT interpreted.

## Decision rule (fixed now)
- "The clock acts on the cycle mainly through the exits from state i": R_i's 95% interval lower bound >= 0.5, AND
  O_i's point estimate is the largest of the three.
- If no state meets it: "the rhythm is shared across states". The R and O estimates are reported with intervals and a
  ranking only where the intervals do not overlap.
- No other binning, flattening scheme (single edges, phase shifts) or genotype is added after looking.

## Secondary, descriptive
- The daily rhythm of each of the six exit probabilities itself (24-h cosine amplitude relative to the day-pooled
  value, with fly-bootstrap intervals): which transitions oscillate most. This is the direct analogue of old Fig. 3A.
- The same decomposition for the sleep occupancy rhythm (stationary pi_S of P_b), to see whether the states that carry
  the flux rhythm also carry the sleep rhythm.

## Output
`sqw27.py` -> `sqw_out/s27_decomp.json`, `figures/fig27_*.png`. Ledger entry with commit hashes.
