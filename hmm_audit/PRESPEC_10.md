# PRESPEC 10 — per-second movement counts: model-free transitions, inactivity durations, count HMM
Date 2026-10-20 (UTC). Data: counts/CS (25 files, 345600 s = 4 d) and counts/fumin (53 files,
259199 s = 3 d); 1 s, start ZT0, LD 12:12, integer counts, no missing data. Only format and value
shares were inspected (0: CS 92.8%, fumin 73.3%). Files cannot be mapped to the scored .mat rows;
the user judges alignment unnecessary: this is a stand-alone analysis with its own sleep definition.
Fly = file (sorted filename order). Day d covers [86400(d-1), 86400 d); light = first 43200 s.

## Shared definitions
- Inactivity run = maximal run of zero counts. Sleep bout = inactivity run >= 300 s (count-based,
  no bridging). Sleep onset = first second of a sleep bout; wake onset = first second after it.

## C1 (model-free): activity around transitions
- Movement fraction f = share of seconds with count > 0.
- Before sleep onset: minutes -6..-2 (seconds -360..-61); minute -1 excluded because second -1 is
  a movement by construction. After wake onset: minutes +2..+6 (seconds +61..+360); minute +1
  excluded likewise. Only transitions whose adjacent wake bout (time between sleep bouts) is
  >= 20 min, so the windows lie inside wake.
- Reference: for each such wake bout, the 5-min stretch centred on its midpoint.
- Per fly and LD phase: contrast_pre = f(pre) - f(reference), contrast_post = f(post) - f(reference),
  averaged over that fly's eligible transitions (flies with >= 3 eligible transitions).
- Test per genotype x phase: two-sided Wilcoxon signed-rank across flies, for each contrast.
  4 tests per genotype (pre/post x light/dark); all reported, no correction applied, count stated.

## C2: inactivity-run durations of ALL lengths (no 300 s cut)
- Interior inactivity runs (first and last run of each fly excluded), assigned to phase by start.
- Discrete hyperexponential, xmin = 1 s (n_min = 1), K = 1..5 (pooled.py fit, 30 starts), per
  genotype x phase; K chosen by 5-fold CV by fly (fold = fly index mod 5); report q and means.
- For the CV-chosen K: fraction of runs >= 300 s whose posterior responsibility is assigned to
  components with mean < 300 s ("short-kind" tail crossing the sleep threshold).

## C3: 3-state count HMM (the user's model, with count emissions)
- Observation per second: symbol 0,1,2,3 = count 0,1,2,>=3 (categorical emissions, no shape
  assumed). Time-homogeneous per fly x phase: each day's light (or dark) half is a separate
  sequence sharing one parameter set. hmmlearn CategoricalHMM, EM, n_iter <= 500, tol 0.01,
  8 random restarts (random_state 0..7), start probabilities fitted.
- Labels: H3 = state with largest P(symbol 0), H2 = smallest P(symbol 0), H1 = the third.
- Identifiability checks (primary question: does counting break the degeneracy?):
  (i) emission separation: L1 distance of H1's row to H2's and to H3's (0..2); "degenerate" if < 0.05;
  (ii) restart spread: among restarts with logL within 0.5 nats of the best, the range of
       Xn_pred and Yn_pred (formulas of PRESPEC 1); "route-unidentified" if either range > 0.2.
- Route counts at the best restart: FFBS, 100 draws per sequence, routes do not span sequences;
  N32, N312, N23, N213 posterior means summed over sequences; Xn, Yn from summed means, with 95%
  posterior interval from draw-level ratios.
- Model adequacy (reported, not used for C3's K): held-out last day. Fit K = 2, 3, 4 (2 restarts)
  on days 1..D-1 and report the held-out log-likelihood per second for each K.
## Reporting
Per fly tables (CSV) and per genotype x phase summaries; every check reported for every fly.
