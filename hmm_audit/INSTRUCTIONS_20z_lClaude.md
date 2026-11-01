# For lClaude: PRESPEC 20z Part A (sleep-bout duration law of the Kir2.1 flies), MATLAB bout suite
Plan: hmm_audit/PRESPEC_20z.md (committed 24c358b BEFORE any bout statistic was computed). Read it first; Part A is
yours, Part C is done (ledger "PRESPEC 20z Part C"). The Icode root CLAUDE.md applies in full (run tests by name from
the repo, `assert(isequal(fileparts(which('CompareBoutModels')), pwd))`, pinned figure theme, q not w).

## Input
hmm_audit/kir_sleep_bouts_20z.csv (SHA-256 prefix 3116971a1fa3c44c): columns group (empty / dfb / vnc), fly, monitor,
start_s, duration_s, start_ZT_h. Sleep bouts >= 300 s, 1-s data, record edges already dropped. 5132 / 8671 / 9149 bouts.
xmin = 300, SamplingInterval = 1. Group ids: (group, fly) is unique; fly alone is NOT (names repeat across groups).

## Runs (fixed by the PRESPEC; nothing else)
A1 per group: R = CompareBoutModels(x, 300, SamplingInterval=1, MaxComponents=5). Report R.Skipped, R.Table, R.Selected,
   R.SelectionPath, the winner's Diagnostics (GuardOK, MixtureMinCount, TailFraction), q.
   n-matched order: subsample each group without replacement to min group n (5132), 20 times (RandomSeed 2000 + s),
   FitHyperexponentialMLE K = 1..5, record the order chosen by AICc and by BIC; report modal order and frequency.
A2 (PRIMARY) per group: FitHyperexponentialMLE(x, 300, K=2, SamplingInterval=1) (discrete). tau_slow = larger time
   constant, q_slow = its observed weight. Fly bootstrap B = 500 (resample flies within group, pool their bouts, refit
   K = 2; RandomSeed 3000 + b). Contrasts dfb - empty and vnc - empty of ln(tau_slow) and q_slow: point difference, 95%
   percentile CI of replicate differences, p = 2 min(P(diff* <= 0), P(diff* >= 0)) floored at 1/500; Holm over 4.
   Guard failure (GuardOK false or min branch q x n < 10): the group is "not identified"; failed replicates dropped and
   counted (state it if > 10%).
A3 the A2 procedure on bouts with start_ZT_h >= 12 only; Holm over 4 separately; censoring caveat stated.
A4 BoutHazard(x, 300, SamplingInterval=1, NumBins=b, Plot=true) for b = 12, 16, 20, 24 per group; descriptive, no test.

## Output
Root RESULTS_LEDGER.md entry "PRESPEC 20z Part A" (date, commit hash, every number above, R.Skipped first), code and
logs under analysis_20z/, figures with theme pinned. Push to claude/trusting-allen-7nxrsf. No raw data.
