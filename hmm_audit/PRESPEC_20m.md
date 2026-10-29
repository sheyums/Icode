# PRESPEC 20m — time-resolution check: 30-s steps vs 1-min steps (before running)
Date 2026-10-29. User question: did moving from 1-s data to 1-min symbols lose resolution? Genotypes: CS and fumin
(user's choice of CS and fmn/sss; fumin has more flies and fewer phases without a Q/W reference). theta 0.5.
## Symbols on a 30-s grid
Same rules as PRESPEC 20 scaled to the bin: S if >= 15 of the bin's 30 s lie in a (bridged) sleep bout; otherwise Q
if the bin's active fraction < theta r_f, else W; r_f unchanged (computed from the 1-s data, per light phase).
Missing (no r_f) as before. Whole days only; 2880 steps per day; phase = light if step mod 2880 < 1440.
## Model
K = 3 only (model order not re-examined), shared emissions, light/dark transition matrices, EM with 20 starts,
seed as 20b; labelling and degeneracy rule as 20b. FFBS 200 paths.
## Comparisons with the 1-min (20b) results, per fly x phase, paired
1. H1 time scale: mean dwell 1/(1 - P_11) converted to minutes; fraction of decoded H1 runs lasting a single step.
2. Count J_cycle (cycles/h) and EPR bound (nats/h, per path then averaged).
3. Direction: c_pred and FFBS c; fraction of fly-phases with the same sign at both resolutions (non-zero cycles).
4. Memory, R1 in CS dark: hidden D with the matched definition rescaled (H3 >= 10 steps = 5 min, gaps >= 20 steps
   = 10 min); per-fly D at both resolutions, Spearman, and group mean. Descriptive (no new simulation null).
Paired Wilcoxon (30 s vs 1 min) for items 1-2 per genotype x phase; descriptive.
## Reading fixed now
- If faster dynamics were hidden by the 1-min grid: H1 mean dwell (in minutes) decreases, J and EPR increase.
- Results are robust to the step size if routing directions and the CS-dark R1 effect keep their sign.
- The 30-s Q/W call rests on half the data per symbol, so more Q/W flicker (shorter H1/H2 runs) is expected from
  noise alone; an increase in J or EPR is therefore an upper bound on what the finer grid reveals.
