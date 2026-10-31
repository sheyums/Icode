# 20r for lClaude (MATLAB): hidden steps between waking and sleep

Prepared 2026-10-31 in the S/Q/W audit session. The pre-specifications are the three entries dated
2026-10-31 at the end of the repo-root `PRESPEC.md` (20r-A, 20r-B, 20r-C). They were committed BEFORE
anything was run; do not edit them. If something in them cannot be done as written, add a NEW dated
entry saying what changed and why, commit it, and only then run.

## Why this exists (two paragraphs, so you do not need the audit history)

The physics question is whether waking-to-sleep passes through sequential hidden stages. The
specific signature is a wake-bout hazard that RISES after a trough: no mixture of memoryless
processes (and no between-fly heterogeneity) can make a hazard rise, while a series of stages can.
FINDINGS.md shows this rise in per0 DD wake bouts (hyper-Erlang series branch); the root ledger
(2026-09-15) has it "supported but bin-sensitive" at 12/16/20/24 bins and not carried by any single
fly. Two things are open: (A) is it there in CS DD (clock intact)? (B) is it in the flies, or made by
scoring bursty movement into wake bouts? (C, exploratory) does it track the per-fly excess
irreversibility found in the S/Q/W audit (PRESPEC 20p: per0 DD carried extra time-asymmetry)?

The S/Q/W audit randomness-parameter idea was dropped on purpose: with the FINDINGS per0 mixture
(81% memoryless 660 s, 8% memoryless 9.5 s, 11% five-stage Erlang 720 s) Var/mean^2 = 1.04, so
heterogeneity hides the stages from r. The hazard is the right statistic.

## Rules that apply (CLAUDE.md; read it and FINDINGS.md first)
- Run from the Icode folder, by name, behind `assert(isequal(fileparts(which('BoutHazard')), pwd))`.
  Check `Icode_sync_stamp.txt` in the analysis folder against `git log -1` before trusting copies there.
  Do not delete or "clean up" copies in the user's analysis folder.
- `SamplingInterval` = 1 everywhere; discrete mode; report q, not w; read `R.Skipped` first.
- Pin figure themes (`theme(fh,'light')` before axes) if you plot.
- Never change the bin set, never widen the window, report every bin count (rule 4).
- Results go into the repo-root `RESULTS_LEDGER.md` in its template, one entry per PRESPEC entry,
  with the commit hash, including null and inconclusive results.

## Step 0 — gate (must pass before ANY CS DD number is looked at)
Re-extract per0 DD wake and sleep bouts from the raw 1-s count files with the SAME code and settings
that produced `subset_per0DDwakebouts.mat` (1-s ironout bridging, bridgeWakeSeconds = 1, dt = 1 s;
find the producing script in the analysis folder; record its name and hash). Confirm:
- 3989 wake bouts and 3987 sleep bouts, identical after sorting to the stored sets (and FINDINGS);
- the 16-bin hazard reproduces the FINDINGS table (e.g. 1.16e-2 at 1.41 s, trough near 109 s).
If either fails: stop, record the failure in the ledger, report to Sheyum.

## Step 1 — 20r-A (CS DD, blind)
Data: the 25 CS DD files (`*_from1stLD.txt`, from CS_DD_activfiles.zip, uploaded 2026-10-30).
Record a SHA-256 per file. Extract wake (xmin 2 s) and sleep (xmin 300 s) bouts with the gated code.
1. Rise test at NumBins {12,16,20,24} exactly as the 2026-09-15 per0 entry (null = K = 2
   hyperexponential fit's SurvivalHandle, 999 replicates, RandomSeed 1). Report trough, peak,
   RiseRatio, null median and 95th, RiseNullP, NullDropped at every count.
2. Per-fly influence: the 2026-09-15 removal design unchanged.
3. Sleep-bout control: the same rise test on CS DD sleep bouts.
4. `CompareBoutModels(x, 2, SamplingInterval=1, MaxComponents=5)` on CS DD wake bouts; report in the
   CLAUDE.md order (Skipped, Table, Selected/SelectionPath, guards, q, ShapeSweep flatness).

## Step 2 — 20r-B (renewal control, per0 DD then CS DD)
Per fly, from the raw trace (same trimming as the extraction): run-length encode into active runs
(count > 0) and inactive runs (count == 0, all lengths). Build a synthetic trace of the fly's record
length by alternating iid draws (with replacement) from its own active and inactive run lengths,
starting in the state of its first second, truncating the last run. Score it with the identical
pipeline (ironout bridging + extraction, xmin 2 s). Pool the 25 synthetic flies; subsample to the
observed n if larger (report if smaller). BoutHazard RiseRatio at the four bin counts. B = 199,
RandomSeed 2000 + b. Pre-run check: one synthetic fly's active and inactive run-length histograms
match the real fly's. Decision rule in the PRESPEC entry. Also report the renewal null's median
RiseRatio beside the 2026-09-15 two-exponential null medians (1.083, 1.136, 1.186, 1.235 for per0).
For CS DD run it only if 20r-A gave "supported" or "supported but bin-sensitive".

## Step 3 — 20r-C (exploratory, only if hyper_erlang is the per0 / CS DD wake model)
Per fly: mean posterior responsibility of the series branch over its wake bouts (pooled fit).
Spearman with `E_nats_per_h` and with `E_wakeonly_nats_per_h` from `hmm_audit/e_perfly_20p.csv`
(match by file name; skip flies with qw_reference False). Report rho, p, n. Label it exploratory.

## What to send back
The four ledger entries (gate, 20r-A, 20r-B, 20r-C) and any figures, pushed to the repo.
