# PRESPEC 20h — dwell-time distributions, observed (S/Q/W) and hidden (H1/H2/H3) (before running)
Date 2026-10-25. User request, for a main figure. Descriptive: no distribution fitting or model selection here (the
user's MATLAB bout library will do that later). theta 0.5, 20b fits (12 h phases; per0DD whole).

## Runs
- Observed: maximal runs of one symbol on the 1-min S/Q/W sequence: S, Q, W, and wake = Q or W (a run of non-S
  minutes). A missing minute (-1) ends a run, which is then censored there.
- Hidden: maximal runs of H1, H2, H3 on each of 200 FFBS paths per fly (non-degenerate flies).
## Phase boundaries (fixed now)
A run is assigned to the phase of its first minute. A run still going at the end of that phase (ZT 12 or ZT 24) is
RIGHT-CENSORED at the boundary, as is a run cut by the end of the record or by a missing minute. Durations are in
minutes (grid resolution 1 min). Kaplan-Meier survival handles the censoring. The uncensored version (whole runs by
start phase) is reported beside it, because it is what a naive analysis would give.
## Summaries per genotype x phase (runs pooled over flies; flies also summarised separately)
- Kaplan-Meier survival S(t) of each run type; median run length per fly (KM).
- Life-table hazard per minute on FIXED bins (edges in min): 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 721;
  h = events / minutes at risk in the bin; 95% interval from the Poisson count of events. Bins with < 10 runs at
  risk are not shown.
- Observed runs, model check (posterior-predictive): 20 records simulated from each fly's fitted model (same length
  and phases), processed identically; band = central 95% over replicates of the pooled group survival and hazard.
  Reported: fraction of hazard bins where the data lie outside the band.
- Hidden runs: KM survival and histogram (log bins as above) per state, pooled over paths and flies (each fly's
  paths weighted 1/200); model mean dwell 1/(1 - P_ii) minutes per fly and phase, and its group mean; the
  geometric survival implied by each fly's P_ii, averaged over flies, overlaid on the decoded survival.
- Stated now: decoded hidden dwell distributions are shaped by the fitted model (close to geometric by construction);
  they show the states' time scales, not evidence for the model.
## Figures
Fig. 9: observed S, wake, Q, W -- survival (log scale) and hazard with the model band, per genotype x phase.
Fig. 10: hidden H1, H2, H3 -- survival/histograms and mean dwell times.
