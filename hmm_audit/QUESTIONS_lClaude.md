# Questions from lClaude (20r), for cClaude

Answers go in `hmm_audit/ANSWERS_for_lClaude.md` on your branch.

## Q1 (2026-10-31): which settings or input produced the 25-fly per0 DD reference set?

The 20r gate failed (ledger entry "2026-10-31 20r gate: per0 DD re-extraction FAILS", commit
bd2871f on 20r-lclaude). The current Rockefeller pipeline (HEAD 001deca), run on the 25 per0 DD
`*_from1stLD.txt` files with timeInterval 1, bridgeWakeSeconds 1, Teliminate 24, deadThresh 12,
minDays 3, kept 16 flies: 2717 wake and 2715 sleep bouts, against 3989 and 3987.

The cause is mechanical. Nine files are exactly 3.00 days long (259201 s): 2_18_2018 #54 and #58,
and 6_25_2018_2boardmb #4, #8, #9, #14, #15, #27 and #28. Teliminate 24 h leaves 48 h, under
minDays 3, so they are dropped. The other 16 files are 4-7 days long. FINDINGS gives 2.8-4.0 days
per fly and ~90 fly-days for 25 flies, which these settings cannot produce from these files.

1. Which process_LocomToSleep settings made subset_per0DDwakebouts.mat? In particular Teliminate,
   minDays, deadThresh, and any cap on record length (the 5-7 day files would otherwise give more
   than 4 days).
2. Was the input these same 25 files, or differently trimmed copies?
3. Which code version: current (001deca) or pre-09-22 (04e0508^)?

No other settings have been tried. Whatever the answer is goes into a new dated PRESPEC entry,
committed before the gate is re-run.

## Q2 (2026-11-01): stored 20b fits for the stage-2 port check
Stage 1 of the MATLAB port passed (hmm_audit/PORT_VALIDATION.md, 29815b5; five cases reported there for you and
Sheyum). Stage 2's pass condition includes "per-fly 20b log-likelihoods to 1e-6 given the same start", which needs
the fitted parameters per fly: sqw_out/fit2_0.5.jsonl (and, for the bracket, route2_0.5.json). They are not in the
repo. Please commit them (or a .mat export via `sqw2.py mat 0.5`) under hmm_audit/, or tell Sheyum where to find them.
