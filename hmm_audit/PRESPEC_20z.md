# PRESPEC 20z — hidden structure inside sleep: does silencing the dfb or vnc cells remove a deep sleep substate?
Date 2026-11-01. Written and committed BEFORE any sleep-bout distribution, hazard or route-vs-duration statistic of the
Kir2.1 flies has been computed. What has been seen: 20y (sleep amount, S/Q/W transitions, landscape; ledger cf777b1).
Sheyum's go-ahead 2026-11-01.

## Why this question, and why it is identifiable
20y: silencing raises E_S - E_Q and makes sleep easier to leave and harder to enter, without a detectable change of the
cycle's drive. The S/Q/W chain treats S as ONE state. Inside a scored sleep bout there is no activity by definition, so
sleep substates have no emission signature; they can only show in TIMING. If S is an aggregate of hidden states, the
left-truncated sleep-bout duration law is a mixture of exponentials, and the number of components bounds the number of
hidden sleep states from below (Fredkin & Rice 1986). Mixture weights (observed q) and time constants are identified;
the aggregated transition matrix is not (PRESPEC 13). Hence the test is on the duration law, not on an HMM.
Literature on depth structure of fly sleep: van Alphen et al. 2013 J Neurosci 33:6917; Tainton-Heap et al. 2021 Curr
Biol 31:578; Wiggin et al. 2020 PNAS 117:10024. All tests are TWO-SIDED (the two line papers are still unread).

## Data
The 20y flies (empty 19, dfb 25, vnc 34; none dead; 8 whole days from ZT0). Same caveats as 20y: constitutive Kir2.1;
enhancerless full-GAL4 control vs split-GAL4 drivers; common non-CS background; compare ONLY with empty>Kir.
Sleep bouts: sqw_ref.sleep_bouts(x, bridge=True) on the whole-day trace (8 x 86400 s); bouts starting at second 0 or
ending at the last second are dropped (the 20u G1 convention, shown identical per fly to the MATLAB 20r extractor).
Durations in seconds; xmin = ceil(300/dt)*dt = 300 s; dt = 1 s. Each bout carries fly id, monitor, and the ZT of its
start. Exported as derived data (kir_sleep_bouts_20z.csv); raw counts are not pushed.

## Part A — the sleep-bout duration law (MATLAB bout suite; run by lClaude; Icode root CLAUDE.md applies)
Block (primary): all sleep bouts of a group pooled over flies and the whole day (bouts assigned by start; only record
edges censor, and those bouts are dropped).
A1 (descriptive) CompareBoutModels(x, 300, SamplingInterval=1, MaxComponents=5) per group, default library. Report in
   CLAUDE.md order: R.Skipped, R.Table, R.Selected and SelectionPath, guards of the winner, q (not w).
   Order comparison matched in n (rule 5): each group subsampled without replacement to the smallest group's bout
   count, 20 subsamples (RandomSeed 2000 + s), the hyperexponential order selected by AICc and by BIC (K = 1..5)
   recorded; report the modal order and its frequency per group. Descriptive, no test.
A2 (PRIMARY) FitHyperexponentialMLE(x, 300, K=2, SamplingInterval=1, discrete) per group. Readouts: tau_slow (the
   larger time constant) and q_slow (its observed weight). Uncertainty: fly bootstrap, B = 500 per group (resample flies
   with replacement, pool their bouts, refit K = 2; RandomSeed 3000 + b, independent per group). Contrasts dfb - empty
   and vnc - empty of ln(tau_slow) and of q_slow: difference of point estimates, 95% percentile CI of the difference of
   bootstrap replicates, p = 2 min(P(diff* <= 0), P(diff* >= 0)), floored at 1/B. Holm over the 4 tests.
   A group whose K = 2 fit fails its guard (GuardOK false or smallest branch q x n < 10) gets "not identified" for
   that group's contrasts; replicates whose guard fails are dropped and counted; > 10% dropped is stated.
   Reading: "the slow (long) sleep component is shorter / smaller after silencing" only if Holm p < 0.05. Wording
   about substates: "consistent with loss of a long-lived sleep state"; never "the deep sleep state is removed" (the
   mixture bounds the number of states from below and does not locate them).
A3 (secondary) A2 on bouts STARTING in ZT 12-24, because dfb's sleep loss is night-only (20y). Holm over 4,
   separately. Caveat stated: long bouts starting late at night run into the light phase (censoring gap, CLAUDE.md).
A4 (descriptive) BoutHazard(x, 300, SamplingInterval=1, NumBins = 12, 16, 20, 24, Plot=true) per group: the hazard
   of waking against time asleep, with Wilson bands. All four bin counts reported. No test.

## Part C — the observable footprint of sleep depth: exit route vs time asleep (Python, sqw_ref; run by cClaude)
sqw_ref.route_duration_observed(x, theta=0.5, light_dark=True) and route_duration_D(records, 'R1_exit_vs_sleep',
phase), exactly as 20k: per fly x phase D = mean ln(sleep length | quiet exit) - mean ln(sleep length | active exit),
>= 2 bouts of each kind.
C1 per group x phase: group mean D against the within-fly label-permutation null (2000 permutations, seed 20), as
   20k. Descriptive here (does the dependence exist in this background?).
C2 (PRIMARY for Part C) dfb vs empty and vnc vs empty, per phase: two-sided Mann-Whitney on per-fly D. Holm over 4.
   Reading: "silencing changes how the exit from sleep depends on time asleep" only if Holm p < 0.05.

## Narrowest blocks
A2: all sleep bouts of a group over 8 days (pooled; the fly bootstrap carries the fly structure). A3: bouts starting
in ZT 12-24. C: one fly's light (or dark) sleep bouts over 8 days.

## Bounds
One run of each part, one verdict each. No other K, no other xmin, no other windows, no re-binning. If Part A cannot be
run in MATLAB (lClaude unavailable), it waits; it is NOT ported to Python for this test.
