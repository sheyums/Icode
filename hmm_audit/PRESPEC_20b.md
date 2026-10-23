# PRESPEC 20b — S/Q/W HMM per fly x light/dark (before running; replaces PRESPEC 20's per-window fits)
Date 2026-10-23. Reason: PRESPEC 20's synthetic positive control failed (one 4 h window: BIC chose K=3 in 4%,
routing SD 0.75), and a steps-per-fit check (ledger) put estimable routing at ~2880 minutes per fit. The user chose
per fly x light/dark. This change was made on SYNTHETIC results only; no fly data had been fitted.
Unchanged from PRESPEC 20: section 1 (S/Q/W minute coding, bridging (b), theta 0.5 primary, 0.33/0.75
sensitivity, r_f per fly x phase), symbol codes, state labelling rule, degeneracy rule (TV < 0.05), routes and
their FFBS counting, energetics formulas, comparison with the model-free c and its sign prediction, outputs.

## 1. Model (one fit per fly)
- K hidden states; ONE emission matrix B (K x 3) shared by both phases; TWO transition matrices, P_light used for
  transitions into minutes at ZT 0-12 and P_dark into ZT 12-24; initial distribution free. per0DD: a single P
  (whole record), as in PRESPEC 19.
- The whole record (all full days) is one sequence. EM (Baum-Welch with phase-split transition counts),
  20 random starts, seed 20 + fly index; best log-likelihood kept.
- Explicit assumption: each phase is time-homogeneous. Changes within a phase (ZT drift) and across days are
  averaged. Sections 4 (V3b) and 5 test how much this matters.

## 2. Model order
- K = 2, 3, 4 per fly. k = (K-1) + 2K(K-1) + 2K = 9, 20, 35 (per0DD: (K-1) + K(K-1) + 2K = 7, 14, 23).
- BIC with n = number of minutes (primary) and n_eff = number of symbol changes (sensitivity).
- Flag rule unchanged: a fly is "strongly not 3-state" if K=2 or 4 beats K=3 by > 10 BIC units; reported, not
  dropped. Report fractions of flies best at each K per genotype, and how often the two n rules disagree.
- The parametric-bootstrap LRT of amendment A is DEFERRED (cost: ~7 h at this record length); it is run only if
  > 20% of flies in some genotype are flagged, on up to 5 flagged flies per genotype, 50 replicates.

## 3. Routing and energetics (K = 3)
- Per fly x phase: Xn_pred, Yn_pred from P_light / P_dark; FFBS (200 paths) and Viterbi route counts, each route
  attributed to the phase of its completion minute; Xn, Yn from counts (>= 3 routes per denominator);
  c_model = (Yn - Xn)/2. Flies whose H1 row is degenerate are reported and excluded from routing summaries.
- Energetics per fly x phase from P_light / P_dark (zero transitions -> infinite, counted, never floored).
- Group summaries: mean with 95% fly-bootstrap interval (2000, seed 20).
- Uncertainty: parametric bootstrap on 5 non-degenerate flies per genotype (seed 20), 30 replicates, refit at
  K = 3 with 10 starts; SD and bias of Xn_pred, Yn_pred per phase.

## 4. Synthetic validation (FIRST; fly data untouched until logged)
20 synthetic flies x 5760 minutes per case, shared B as in PRESPEC 20's V1.
- **V1b** positive: light Xn = -0.3, Yn = +0.4; dark Xn = +0.2, Yn = -0.1 (V1's dwell structure).
  Pass: K=3 chosen in >= 70% of flies; |mean bias| < 0.1 for all four quantities; SD < 0.3.
- **V2b** negative: 2-state truth (no H1), different light/dark rates. Pass: K=2 in >= 80% of flies.
- **V3b** within-phase drift: 2-state truth whose sleep-exit and wake-exit rates change linearly across each
  12 h phase by a factor 2 (no H1). Pass: K=2 in >= 80% AND the K=3 fit's H1 degenerate or its routing
  stationary H1 occupancy < 5% in >= 80% of flies. A failure means
  within-phase drift can manufacture H1, and fly results must be read with that caveat.
Failures are logged, and the design is reviewed before fly data.

## 5. Sensitivity and model checks (descriptive)
- Day split: refit K = 3 on the first and the second half of the days separately; report Xn, Yn per half.
- Within-phase adequacy: from the fitted K = 3 model simulate 100 records; for each 4 h window, the observed
  fraction of S, Q, W minutes and number of symbol changes vs the simulated 95% band. Report the fraction of
  windows outside the band by ZT window (a time-homogeneous phase should give ~5% everywhere).
- theta 0.33 and 0.75: K = 3 only.

## 6. Export
`.mat` per genotype as PRESPEC 20 section 9: `R{fly, window}` 1 x 240 symbols; `C{3*fly-2, w}` = 3,
`C{3*fly-1, w}` = P of that window's phase, `C{3*fly, w}` = the shared B. Scratchpad only.
