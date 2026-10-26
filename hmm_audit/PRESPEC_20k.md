# PRESPEC 20k — does routing depend on bout duration? (before running)
Date 2026-10-26. User request. Question: is the choice of route (via quiet wake / H1 vs direct) history-dependent,
which a Markov model forbids for exits from sleep?
## Relations (quiet = Q at the observed level, H1 at the hidden level; active = A / H2)
- R1 exit route vs the duration of the sleep bout it ends.
- R2 entry route vs the duration of the wake period that precedes it.
- R3 entry route vs the duration of the sleep bout it begins.
## Data (theta 0.5; 20b fits; phase = phase of the sleep onset for R1, R3, of the wake period's end for R2; per0DD whole)
- Observed: bridged sleep bouts (>= 300 s). Exit classified only if the following wake bout lasts >= 600 s, entry
  only if the preceding wake bout lasts >= 600 s; Q if active fraction in the 61-300 s window < theta r_f
  (r_f of the onset's phase), else A. Durations in minutes.
- Hidden: on posterior paths (200 per fly), H3 runs >= 5 min with a previous and a next H3 run; exit classified if
  the gap to the next H3 run is >= 10 min, entry if the gap since the previous one is >= 10 min; entry/exit =
  neighbouring state (H1 quiet, H2 active); wake duration = the gap.
## Statistic
Per fly x phase: D = mean ln(duration | quiet route) - mean ln(duration | active route); needs >= 2 bouts of each
kind (hidden: per path, then averaged over paths where defined). Group statistic = mean D over flies.
D > 0: quiet routes go with longer durations.
## Nulls
- Observed: permute route labels among the fly's bouts within phase (2000 permutations, seed 20); two-sided p.
- Hidden: 40 records simulated from each fly's fitted model (Markov: R1 null by construction), decoded with the true
  parameters (FFBS, 20 paths), same statistic; p = rank of the real group mean among the 40 simulated group means
  (two-sided). Also reported: the simulated mean D (what a Markov model produces through the pipeline).
Holm over the genotype x phase blocks, per relation and level.
## Descriptive figure
P(quiet route) in fixed duration bins with Wilson 95% intervals: sleep bouts 5-10, 10-20, 20-40, 40-80, > 80 min;
wake periods 10-20, 20-40, 40-80, 80-160, > 160 min; observed and hidden, per genotype x phase.
## Reading (fixed now)
A relation is called history-dependent if the real D differs from its null at Holm-corrected p < 0.05. If so, the
landscape-vs-observed routing plot is examined next.

## Amendment A (2026-10-26, after the first run; power fix, not a change of statistic)
With 40 simulated records the smallest attainable hidden-level p is 2/41 = 0.049, so after Holm over 9 blocks no
hidden-level result could reach significance (min Holm p ~0.44). The first run is kept (sqw10_nsim40.log). The
hidden-level null is rerun with 400 simulated records per fly (same statistic, same decoding), seed offset +1.
Observed-level results unchanged.
