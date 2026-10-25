# PRESPEC 20i — sensitivity: exclude fly-phases without a Q/W reference (before running)
Date 2026-10-25. Found while answering a user question: 15 fly-phases (fumin 4, sss 11) have fewer than 3 qualifying
wake bouts, so r_f is undefined and every wake minute of that phase is missing; their phase matrix is then not
informed by wake data. Recompute, with those fly-phases EXCLUDED (everything else unchanged, theta 0.5, 20b fits):
c_model (FFBS, per-path mean), count J_cycle, EPR bound (per path), and observed / hidden-matched round-trip
forward - reverse. Same tests as originally (bootstrap interval of the mean; Wilcoxon for round trips).
A result HOLDS if its sign is unchanged and its interval still excludes 0 (or its Wilcoxon p < 0.05) when it did before.
