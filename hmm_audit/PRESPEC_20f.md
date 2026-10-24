# PRESPEC 20f — hidden vs observable entropy production: difference and ratio, with estimator calibration
Date 2026-10-24, before running. Requested by the user after Fig. 6b (post hoc figure; the medians in it were seen).
Data and definitions as PRESPEC 20b (theta 0.5, 12 h phases), counting as sqw_ref.py.
- EPR_h: hidden, rule-of-3 lower bound per FFBS path (200), averaged over paths (Jensen).
- EPR_o: observable, the same bound on run-compressed S/Q/W minute symbols (Q<->H1, W<->H2, S<->H3); transitions
  touching a missing minute dropped.
- Per fly: D = EPR_h - EPR_o (nats/h), r = EPR_o / EPR_h (only if EPR_h > 0).

## Blocks (rule fixed now)
"Driven" = genotype x phase blocks whose 20b count J_cycle 95% interval excludes 0: CS dark, fumin light, fumin dark,
per0LD light, per0LD dark, per0DD, sss light. Excluded (current ~0, ratio is noise): CS light, sss dark (reported
descriptively only).

## Real data (per driven block)
Median D and median r over flies with 95% fly-bootstrap intervals (2000, seed 20); Wilcoxon signed-rank of D vs 0;
fraction of flies with D < 0 (where the estimators violate the true-rate inequality).

## Calibration on synthetic flies (per driven block)
For every non-degenerate fly: 20 records simulated from its fitted model (its length and phases), seed 2020 + fly.
On each record:
- E_model: the fitted model's exact stationary hidden EPR per phase (nats/h), as sqw_ref.model_energetics;
- E_path: rule-of-3 bound on the TRUE simulated hidden path's counts;
- E_ffbs: the estimator used on real data, FFBS (100 paths) under the TRUE parameters (no refit; parameter-estimation
  error is therefore NOT included -- 20b's bootstrap covers that for routing);
- E_obs: rule-of-3 bound on the simulated symbols.
Reported per block: median relative bias of E_ffbs vs E_path and vs E_model; simulated distribution of the group
median r (replicate-wise, flies pooled) and of the fraction with D < 0; real values set beside them.
## Reading (fixed now)
- The real median r is called CONSISTENT with coarse-graining of the fitted dynamics if it lies inside the central 95%
  of the simulated group-median r.
- The hidden-minus-observable difference is reported as a physical quantity ("hidden dissipation") only for blocks
  where (i) Wilcoxon D > 0 at p < 0.05 (Holm over the 7 driven blocks) and (ii) |median relative bias of E_ffbs vs
  E_path| < 0.25. Otherwise descriptive.
