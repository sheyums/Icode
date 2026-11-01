# PRESPEC 20u — merged paper, test 1: is the per0 DD wake-hazard rise produced by reversible dynamics plus scoring?
Date 2026-11-01. Written before anything runs. Owner: the Plan B feasibility session (branch planB-feasibility).
Design agreed with the S/Q/W audit session (messages of 2026-11-01). The go-ahead to RUN is Sheyum's and is NOT given
by this file. Results go to RESULTS_LEDGER_planB.md.

## Question
20r replicated a rise of the per0 DD wake-bout hazard after a trough (RiseNullP 0.001/0.009/0.059/0.086 at 12/16/20/24
bins). 20r-B: a burst-pause renewal null (iid runs, statistically time-reversible) gave p_ren 0.015/0.015/0.050/0.050
(INCONCLUSIVE). Here the null is generalised to CORRELATED but time-reversible 1-s dynamics: does the best-fitting
reversible hidden Markov model of each fly's activity, scored by the identical pipeline, produce the rise?
Motivation, NOT the test: in a stationary reversible chain the sojourn in a set of states has a completely monotone
density, so its hazard can only fall (Keilson 1979, ch. 5; discrete time needs a lazy chain, i.e. a positive spectrum).
A SCORED wake bout is not such a sojourn. It ends at a silence that LATER reaches 300 s, a hitting time of an
irreversible silence counter, so the theorem does not apply to it. The simulation through the scoring IS the test.

## Data, block
per0 DD, 25 flies: the cut 1-s count files uploaded 2026-11-01 (per0DD_activity_files_1sec_res.zip), the 20r input.
Each file name must match analysis_20r/per0DD_bouts_20r_replication.mat `files`, and its length must be >= traceLen.
Narrowest block: one fly's retained trace (first traceLen seconds); the statistic pools the 25 flies (as 20r, 20r-B).
per0 DD only: it is the one stationary, undriven block. LD minute-level rises (20h) are excluded because within-phase
drift (20o) can raise a hazard without irreversibility.

## Stage 0 — gates (run first; if either fails, STOP and report; the port is NOT adjusted to match)
G1 extractor. A Python port of the 20r extractor: lab bridging of single active seconds (bridgeWakeSeconds 1,
  ironout), 5-min sleep criterion, Teliminate 0, deadThresh 12, minDays 3, phase 'entire', bouts touching the first or
  last retained sample dropped. PASS: per fly, sorted wake bouts AND sorted sleep bouts identical to the reference .mat
  (SHA-256 b14ba340...; 4017 wake bouts in total).
  The Rockefeller MATLAB code is NOT in the repository; the port is written from sqw_ref.sleep_bouts and the 20r
  documentation. If G1 fails, the fallback is NOT to tune the port. Instead the synthetic count files are scored in
  MATLAB by extract_20r (as 20r-B did), which needs Sheyum or lClaude to run them locally. That is a logistics
  decision for Sheyum.
G2 statistic. A Python port of BoutHazard RiseRatio (SamplingInterval 1, Spacing "log", MinAtRisk 10, edges from each
  dataset itself). PASS: on the reference wake bouts it reproduces 20r-B's observed RiseRatio 1.416 / 1.402 / 1.385 /
  1.435 at 12/16/20/24 bins to 3 decimals. If G2 fails, BoutHazard is run in MATLAB on the null bouts (same fallback).

## Stage 1 — the reversible null model (per fly)
Observation a_t = 1 if count > 0 else 0, per second (only count > 0 vs 0 enters the scoring; 20r-B note).
Aggregated hidden Markov model: Ka hidden states emit 1 (active), Ki emit 0 (inactive), deterministic emissions.
Hidden 1-s transition matrix P REVERSIBLE by construction: P_ij = F_ij / pi_i with F symmetric, non-negative,
sum F = 1, pi_i = sum_j F_ij. Fitted by EM (forward-backward on the run-length-compressed sequence; M-step = the
reversible MLE fixed point x_ij <- (N_ij + N_ji) / (N_i / x_i + N_j / x_j), Bowman, Beauchamp, Boxer & Pande 2009
J Chem Phys 131:124101; Prinz et al. 2011 J Chem Phys 134:174105). 10 random starts per order, seed 3000 + fly;
stop when the relative log-likelihood gain < 1e-9 or after 2000 iterations.
Orders (Ka, Ki) in {1, 2, 3} x {1, 2, 3}, so Kmax = 6. Why modest: fits are per fly, with a few thousand runs each.
Ki >= 2 is needed for long silences alongside short pauses; Ka >= 2 for bursty activity. More states per fly would fit
noise (PRESPEC 12: replicate SD 0.6-0.7 on routing in per-fly hidden models). Order chosen by BIC with
n = the fly's number of runs. Laziness check: the smallest eigenvalue of the selected P must be > 0; flies where it is
not are listed (no refit).

## Stage 1b — adequacy gate (before the verdict)
From 20 simulated traces per fly (seed 4000 + fly), scored identically: per-fly ratio simulated / observed of
(i) wake-bout count and (ii) total sleep (sum of sleep bouts). PASS: median over flies of each ratio in [0.85, 1.15].
FAIL => verdict "NULL INADEQUATE": the best reversible model does not reproduce the bout process, and the test says
nothing either way. Reported regardless: the per-fly ratios, and 20r-B's renewal excess (~15% more wake bouts).

## Stage 2 — the null distribution
B = 199 pooled datasets; replicate b uses seed 5000 + b for every draw. Each fly is simulated from its selected model
for its traceLen seconds, starting from pi (simulation by hidden jump times, not per second), and scored identically.
Pool the 25 flies; if the pooled wake-bout count exceeds the observed 4017, subsample uniformly without replacement to
4017 (if smaller, keep and report). RiseRatio at 12/16/20/24 bins.
p_rev = (1 + #{RiseRatio* >= observed}) / 200 per bin count.

## Decision rule (20r-B's rule, unchanged) and verdict wording, fixed now
Adequacy failed -> "NULL INADEQUATE" (overrides everything below).
p_rev < 0.05 at >= 3 of 4 counts -> RISE BEYOND REVERSIBLE DYNAMICS. Wording: "The per0 DD wake-bout hazard rises
  more than the best-fitting time-reversible hidden Markov model of each fly's 1-s activity (up to 6 hidden states),
  scored identically, produces (p_rev ...). Within that class, reversible dynamics plus the scoring rule do not account
  for the rise." NOT allowed: "proves broken detailed balance" (only one model class was searched).
p_rev >= 0.05 at >= 3 of 4 counts -> EXPLAINED BY REVERSIBLE DYNAMICS PLUS SCORING. Wording: "A time-reversible
  hidden Markov model of each fly's activity, scored identically, reproduces the rise; the rise does not evidence
  irreversible hidden wake dynamics." C8 is then reworded accordingly.
Otherwise -> INCONCLUSIVE (bin-sensitive); all four p_rev are reported, never one alone.
All four bin counts are reported in every case.

## Matched null, narrowest block
Same flies, same record lengths, same scoring, the observed n by subsampling. Narrowest block: one fly's retained
trace; pooled statistic. Time box: if stage 1 has not finished in 6 h wall time on 4 cores, report "not completed".
Orders and starts are not reduced after the fact.

## Amendment A (2026-11-01, before ANY part of 20u ran; go-ahead for tests 1-3 given by Sheyum, relayed by cClaude)
1. DEATH CRITERION ON SYNTHETIC TRACES (raised by cClaude). The pipeline's deadThresh 12 never acted on the real cut
   files, but a reversible-null fly with long silences could trigger it, and the Python port does not reproduce MATLAB's
   trimming. Assumption, stated: deadThresh 12 means a continuous inactive run of >= 12 h (as Teliminate is in hours).
   Check, before stage 2: on the REAL per0 DD retained traces, no inactive run reaches 12 h (otherwise G1 would already
   differ). Rule for the nulls: a synthetic fly-trace with any inactive run >= 12 h is EXCLUDED from that replicate's
   pool. The total count of excluded fly-replicates is reported. If more than 5% of the 25 x 199 fly-replicates are
   excluded, the verdict carries "death-criterion exclusions > 5%" as a stated limitation (no change of rule).
2. FALLBACK (a) APPROVED (Sheyum): if G1 fails, lClaude commits the five Rockefeller extraction files with SHA-256 to
   the shared branch under analysis_20r/rockefeller_extractor/; that branch is merged, and the port is rewritten line by
   line. The switch is recorded as a dated amendment before use. There is still no tuning: the port must match per fly
   exactly, or stop and report. Option (b) is used only if (a) also fails.
3. G1 implementation fixed now: sqw_ref.sleep_bouts(x[:traceLen], bridge=True) for sleep bouts; sleep bouts starting at
   second 0 or ending at traceLen are dropped. Wake bouts are the gaps between consecutive sleep bouts (all of them,
   including gaps next to a dropped edge sleep bout); leading and trailing wake segments touch the edge and are not bouts.
   G2: BoutHazard's buildEdges / hazard / trough-peak logic ported verbatim, with MATLAB round (half away from zero).
4. Reporting conventions (20s-A correction, shared branch 0212253): medians rank Inf highest where infinities occur;
   "no cycle" is split into its three kinds. These apply to 20v and 20w; 20u has neither.
