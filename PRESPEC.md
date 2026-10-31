# Pre-specifications

One dated entry per analysis, **written and committed before the analysis
runs**. An entry added or edited after seeing a result is not a
pre-specification; start a new entry instead and say what changed.

Rules 2, 4 and 5 in `CLAUDE.md` govern this file. Results go in
`RESULTS_LEDGER.md` and reference the entry here that they answer.

## Template

```
### YYYY-MM-DD  <short name>

Question        What is being asked, in one sentence. A question, not a hope.
Data            Genotype, condition, which animals, which recording.
Time window     The block: ZT range / startZT segment / days. The NARROWEST
                window the numbers will be computed from.
Fixed choices   xmin and dt (protocol facts, stated not chosen); models to be
                fitted; NumBins set for any hazard; bootstrap B; seeds.
Decision rule   What outcome would count as support, and what would count
                against. Written so that both are possible.
Matched null    What the result will be compared against, at the SAME n.
Notes           Anything that would otherwise be reconstructed from memory.
```

## Entries

<!-- Append below. Do not edit an entry once its analysis has run. -->

### (none yet)

The per0 DD sleep and wake bout analyses in `FINDINGS.md` predate this file.
They are recorded there with their caveats stated explicitly -- including that
`hyper_erlang` and `weibull_mix` were added BECAUSE the earlier library failed
on those data, so their goodness-of-fit p-values are optimistic. That is
exactly the situation this file exists to prevent repeating.

### 2026-09-15  per0 DD sleep: pooled K=3 against a per-fly null

```
Question        If each fly's sleep bouts have at most two timescales, does pooling the 25 flies
                reproduce the pooled evidence for a third component (K = 3 against K = 2)?
Data            per0, constant darkness, 25 flies, 1-s bridged scoring. Persleep in
                C:\Users\sheyu\Icode\subset_per0DDsleepbouts.mat (SHA-256 AD6EC0B1850095E1...),
                25 x 249: rows are flies, columns are bouts in time order, NaN padding. 3987 bouts
                >= 300 s; discrete one-exponential logL -31665.2380, matching FINDINGS' -31665.24.
Time window     Per-fly fits: each fly's whole analysed record, so the fly is the narrowest block
                for those. The tested statistic: all 25 flies pooled.
Fixed choices   xmin 300 s, dt 1 s, discrete likelihood throughout. Per-fly model:
                FitHyperexponentialMLE, MaxComponents 2, SamplingInterval 1, RandomSeed 1; the
                fly's K = 2 fit, or its K = 1 fit if the K = 2 fit is flagged Degenerate.
                Statistic: LR = 2 (logL at K = 3 - logL at K = 2), both from FitHyperexponentialMLE
                with MaxComponents 3 on the pooled bouts, with the fitter's default multistart for
                the real data. Replicates: each fly's own n bouts drawn from its own model (component
                by q, then the exact discrete draw n_min + Geometric that HyperexponentialLRT uses),
                pooled, and LR computed with HyperexponentialLRT's replicate multistart budget
                (nStartsBase 4, nStartsPerComponent 4, maxStarts 20). RandomSeed 2 for replicates.
                Staged: run 199 replicates; stop if p < 0.02 or p > 0.12; otherwise continue to 999.
Decision rule   p = (1 + #{replicate LR >= observed LR}) / (B + 1). p >= 0.05: spread between
                two-timescale flies is enough to produce the pooled third component. p < 0.05:
                pooling two-timescale flies does not produce it, so something within flies, or the
                scoring, is needed.
Matched null    Every replicate has the same 25 flies with the same bout counts as the data, each
                drawn from its own fitted model.
Notes           Purpose: decide the unit of analysis for the full datasets. p >= 0.05 means a pooled
                mixture fit says nothing about individual flies, so the full datasets are analysed
                fly by fly or with a hierarchical model; p < 0.05 means pooled structure needs more
                than differences between flies. FINDINGS' order test used a pooled K = 2 null with
                no fly structure; this replaces the null, not the statistic, so the observed LR
                should come out near 16.72. Report the fraction of replicate LRs below zero; above
                5%, raise the replicate multistart budget and rerun before reading p. The bouts are
                bridged and cannot be unbridged from these files: per-fly fits absorb merged bouts
                only as far as two exponentials can, so p < 0.05 would not separate sleep biology
                from merging, which needs the raw traces. Also record per replicate whether the
                K = 3 fit has a slow component holding under 5% of bouts (descriptive only). The
                flies have about 3 or about 4 analysed days (9 and 16 flies), which may mark two
                recording runs; not used here.
```

### 2026-09-15  per0 DD wake: does any single fly carry the hazard rise?

```
Question        Does removing one fly lower the pooled wake hazard's rise by more than removing
                the same number of randomly chosen bouts does?
Data            per0, constant darkness, 25 flies, 1-s bridged scoring (ironout). Perwake in
                C:\Users\sheyu\Icode\subset_per0DDwakebouts.mat (SHA-256 039C7BC2E3C8FED1...),
                25 x 248: rows are flies, columns are bouts in time order, NaN padding at row ends
                only. 3989 bouts, the FINDINGS data (verified identical after sorting to
                wakebouts_per0dd.mat); fly i contributes n_i bouts (81-248).
Time window     Each fly's whole analysed record. Every number is computed on the 24 remaining
                flies pooled, or on 3989 - n_i bouts pooled after random removal.
Fixed choices   xmin 2 s and dt 1 s (protocol facts under 1-s bridging). BoutHazard with NumBins in
                {12, 16, 20, 24}, Spacing "log", MinAtRisk 10, Alpha 0.05; the statistic is
                RiseRatio, with no null replicates and no bootstrap for the removal runs. For each
                fly i: RiseRatio with fly i removed, and RiseRatio after removing n_i bouts drawn
                uniformly without replacement from all 3989, repeated 199 times with RandomSeed
                1000 + i. Where a null is needed: FitHyperexponentialMLE on that run's bouts,
                MaxComponents 2, DistributionType "discrete", SamplingInterval 1, RandomSeed 1, its
                K = 2 SurvivalHandle as NullSurvivalHandle, NullReplicates 999, RandomSeed 1. Code:
                Icode main at 8565082 or later, run from the Icode folder behind
                assert(isequal(fileparts(which('BoutHazard')), pwd)).
Decision rule   Fly i is influential at a bin count if removing it gives a RiseRatio below the
                2.5th percentile of its 199 random removals, and influential overall if that holds
                at 3 or more of the four counts. No influential fly: the rise is not carried by any
                single fly. If some are influential, remove all of them together, refit the
                two-exponential null on what remains, and rerun BoutHazard with it. If RiseRatio is
                then at or below that run's own null median (RiseNullQuantiles(1)) at 3 or more
                counts, the rise is carried by a few flies; otherwise some flies matter but the rise
                does not depend on them alone.
Matched null    Random removal of the same number of bouts as the fly (199 draws per fly), on the
                same bins with the same settings.
Notes           Purpose: decide whether wake hazards in the full datasets can be read pooled or must
                be examined fly by fly. The bin-count set is the one proposed as the standing choice
                for every new dataset; it is not blind here, since FINDINGS already has RiseNullP at
                16 and 20 bins on these bouts. By chance about 1 fly in 40 falls below the 2.5th
                percentile at any one count; requiring 3 of 4 counts only partly guards against
                that, because the counts are correlated. The last step compares with the reduced
                data's own null median rather than asking for p >= 0.05, so that the loss of power
                from removing bouts cannot by itself produce the verdict "carried by a few flies".
                Pooling flies cannot create a rise, and nothing here can say whether the rise is made
                by scoring bursty movement into wake bouts; that needs raw activity traces.
```

### 2026-09-15  per0 DD wake: the hazard rise at four fixed bin counts

```
Question        Does the pooled per0 DD wake hazard rise after its trough by more than a strictly
                falling (two-exponential) hazard produces, at each of four bin counts fixed before
                the run?
Data            per0, constant darkness, 25 flies, 1-s bridged scoring. Perwake in
                C:\Users\sheyu\Icode\subset_per0DDwakebouts.mat (SHA-256 039C7BC2E3C8FED1...),
                3989 bouts, the FINDINGS data.
Time window     The whole analysed record of each fly (2.8-4.0 days), all 25 flies pooled. That
                pooled set is the narrowest block for every number here.
Fixed choices   xmin 2 s, dt 1 s. BoutHazard with NumBins in {12, 16, 20, 24}, Spacing "log",
                MinAtRisk 10, Alpha 0.05, Bootstrap 0, NullReplicates 999, RandomSeed 1. Null
                model: FitHyperexponentialMLE on the same bouts, MaxComponents 2, discrete,
                SamplingInterval 1, RandomSeed 1; its K = 2 SurvivalHandle as NullSurvivalHandle.
                Code: Icode main, run from the Icode folder behind the which/pwd assert.
Decision rule   RiseNullP < 0.05 is a rise at that bin count. All four counts: supported. Some but
                not all: supported but bin-sensitive, the wording decided 2026-09-14. None: not
                supported. Any count with NullDropped > 0 is reported as not computable rather
                than given a verdict.
Matched null    999 datasets of n = 3989 drawn from the fitted two-exponential model, each put
                through the whole estimator on its own bin edges.
Notes           Reinstated at Sheyum's request on 2026-09-15, having been dropped earlier the same
                day as a FINDINGS replication; the question of whether the rise is real made the
                few minutes worth spending. NOT blind at 16 and 20 bins, where FINDINGS already
                records 0.010 and 0.068 at RandomSeed 1; 12 and 24 are new. Reporting every count
                is the point of the fixed set. The null model is fitted once and used as if known.
                This tests the rise against a falling-hazard null only: it cannot say whether the
                rise is in the flies or is produced by scoring bursty movement into wake bouts,
                which needs raw traces.
```

### 2026-10-31  20r-A  CS DD wake: the hazard rise at the four fixed bin counts (blind)

```
Question        Does the pooled CS DD wake hazard rise after its trough by more than a strictly
                falling (two-exponential) hazard produces, at each of the four standing bin counts?
                (Is the per0 DD rise clock-independent?)
Data            Canton-S, constant darkness, 25 flies, 4 days each: the 25 files
                *_from1stLD.txt uploaded by Sheyum on 2026-10-30 (CS_DD_activfiles.zip; recorded in
                the same runs as per0 DD, e.g. 12_29_2017 multibeam). Checked DD in the S/Q/W audit:
                no activity step at the old light switches (-0.014 +- 0.008 vs CS LD +0.078).
                1-s counts, 1-s bridged scoring (ironout), extracted with the SAME code and settings
                that produced subset_per0DDwakebouts.mat (gate below).
Time window     The whole analysed record of each fly, all 25 pooled. No subjective-phase split:
                CS DD has a free-running rhythm, but mixing day and night bouts can only make a
                pooled hazard fall faster, never rise, so pooling is conservative for a rise.
Fixed choices   xmin 2 s, dt 1 s (wake; xmin = (floor(bridgeWakeSeconds/dt)+1)*dt with 1-s
                bridging). BoutHazard with NumBins in {12, 16, 20, 24}, Spacing "log", MinAtRisk 10,
                Alpha 0.05, Bootstrap 0, NullReplicates 999, RandomSeed 1. Null model:
                FitHyperexponentialMLE on the same bouts, MaxComponents 2, discrete,
                SamplingInterval 1, RandomSeed 1; its K = 2 SurvivalHandle as NullSurvivalHandle.
                Influence: the 2026-09-15 per-fly removal design unchanged (199 random removals of
                n_i bouts per fly, RandomSeed 1000 + fly). Control: the same rise test on CS DD
                SLEEP bouts (xmin 300 s, dt 1 s), where per0 DD showed no hump.
                Model library: CompareBoutModels on the CS DD wake bouts, default library,
                SamplingInterval 1, MaxComponents 5, read per CLAUDE.md steps 1-3.
Decision rule   As 2026-09-15: RiseNullP < 0.05 at all four counts: supported; some but not all:
                supported but bin-sensitive; none: not supported; NullDropped > 0: not computable.
                Influence verdict by the 2026-09-15 rule. Sleep control read the same way.
Matched null    999 datasets of the observed n from the fitted two-exponential model, each through
                the whole estimator on its own bin edges.
Notes           GATE before any CS DD number: re-extract per0 DD from its raw count files with the
                same code and confirm 3989 wake and 3987 sleep bouts identical (after sorting) to
                subset_per0DDwakebouts.mat / FINDINGS. If the gate fails, stop and report.
                BLIND for CS DD at all four counts. hyper_erlang is now PRE-SPECIFIED for CS DD, so
                its fit p-value there is not optimistic (unlike per0 DD, FINDINGS).
                Prepared in the S/Q/W audit session (hmm_audit/INSTRUCTIONS_20r_lClaude.md).
```

### 2026-10-31  20r-B  per0 DD and CS DD wake: is the rise made by burst-pause renewal and scoring?

```
Question        Does the wake hazard rise survive a null in which each fly's raw 1-s trace is
                replaced by an alternating renewal of its OWN active and inactive runs (same marginal
                run lengths, no correlation between successive runs), scored by the identical
                pipeline? This is the open question of the 2026-09-15 entries ("in the flies, or
                produced by scoring bursty movement into wake bouts; needs raw traces").
Data            per0 DD (25 flies, the 2026-09-15 data) and CS DD (20r-A data), raw 1-s counts.
Time window     Each fly's whole analysed record (the same trimming as the extraction); per fly
                for the null construction, pooled over flies for the statistic.
Fixed choices   Per fly: active runs = maximal runs of count > 0; inactive runs = maximal runs of
                count == 0, ALL lengths (including >= 300 s). Synthetic trace: start in the state of
                the fly's first second; alternate active/inactive runs drawn iid with replacement
                from that fly's own empirical run lengths; stop at the fly's record length (truncate
                the last run). Each synthetic trace goes through ironout 1-s bridging and the same
                bout extraction (xmin 2 s wake, dt 1 s); pool the 25 synthetic flies; if the pooled
                synthetic wake-bout count exceeds the observed n, subsample uniformly without
                replacement to n (if smaller, keep and report). Statistic: BoutHazard RiseRatio at
                NumBins {12, 16, 20, 24}, Spacing "log", MinAtRisk 10, no null replicates.
                B = 199 synthetic pooled datasets per genotype, RandomSeed 2000 + b.
Decision rule   Per bin count p_ren = (1 + #{RiseRatio* >= observed}) / 200. Rise BEYOND renewal if
                p_ren < 0.05 at >= 3 of 4 counts; EXPLAINED by renewal and scoring if p_ren >= 0.05 at
                >= 3 of 4 counts; otherwise inconclusive (bin-sensitive). Applied to per0 DD and CS DD
                separately; for CS DD only if 20r-A gives "supported" or "supported but bin-sensitive".
Matched null    The renewal construction above: same flies, same record lengths, same run-length
                marginals, same scoring, observed n (by subsampling where possible).
Notes           Reading: a rise beyond renewal needs correlations between successive runs (a
                progression through wake), which is what sequential hidden stages would produce. A
                rise explained by renewal means run-length statistics plus the 5-min scoring rule are
                enough; it does not then evidence hidden stages. Also report, descriptively, the
                renewal null's median RiseRatio next to the two-exponential null medians of the
                2026-09-15 entry (does renewal itself make a rise?). Check before the run: one
                synthetic per0 fly's run-length histograms reproduce the real fly's.
```

### 2026-10-31  20r-C  EXPLORATORY: does the hazard-rise branch track the 20p excess irreversibility per fly?

```
Question        Across flies, is the share of a fly's wake bouts attributed to the series (Erlang)
                branch associated with that fly's 20p excess irreversibility E?
Data            per0 DD and CS DD wake bouts (20r-A/B data); E per fly from
                hmm_audit/e_perfly_20p.csv (matched by file name; flies without a Q/W reference
                have no E and are excluded).
Time window     Each fly's whole record.
Fixed choices   The pooled hyper_erlang fit selected for that genotype (per0 DD: FINDINGS shapes
                [1 1 5]; CS DD: 20r-A's fit if hyper_erlang is selected or identified, else not run).
                Per fly: mean posterior responsibility of the series branch over its wake bouts.
                Spearman rank correlation with E, and separately with E_wakeonly.
Decision rule   None. Exploratory and hypothesis-generating: report rho, its p-value, n. Not to be
                cited as evidence unless replicated by a pre-specified test on new flies.
Matched null    Not applicable (exploratory).
Notes           Motivated AFTER seeing 20p (per0 DD carried the excess) and FINDINGS (per0 DD carries
                the series branch): a post hoc link, so it is labelled exploratory here. ~25 flies
                per genotype and noisy per-fly E: low power, and a null result means little.
```

### 2026-10-31  20r gate amendment: stored bout files not found

```
Question        Unchanged from 20r-A: the GATE only is replaced. Does the extraction code, re-run on
                the raw per0 DD files, reproduce the per0 DD bouts that every 20r comparison builds on?
Why             subset_per0DDwakebouts.mat (Perwake, SHA-256 039C7BC2E3C8FED1...) and
                wakebouts_per0dd.mat (w_bouts) were not found: absent from the Icode working tree and
                its git history, and a recursive search of C:\Users\sheyu (the only drive) for
                subset_per0DDwakebouts.mat, wakebouts_per0dd.mat and *per0*wake*.mat found nothing.
                Sheyum, 2026-10-31: treat them as lost. The 20r-A gate ("identical after sorting to
                subset_per0DDwakebouts.mat") therefore cannot be run as written.
Data            per0 DD: the 25 files *_from1stLD.txt in
                C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data\per0 dd
                (the 28 other .txt files in that folder are not used). SHA-256 recorded per file.
Code            The CURRENT extraction code in C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab
                scripts (git repo, HEAD 001deca), chosen by Sheyum on 2026-10-31 over the pre-09-22
                state: process_LocomToSleep.m 7aec365d5b79c859, sleeptrace.m 0deb702eed5cbfcc,
                ironout.m f22f91c262f8c370, segment_sleep_phases.m a627f698159bc97a,
                extractBoutsForModelCompetition.m 64f6f15985d04e5e (SHA-256 prefixes). No new extractor
                is written. Settings: timeInterval 1, bridgeWakeSeconds 1, Teliminate 24,
                deadThresh 12, minDays 3 (defaults, confirmed by Sheyum); phase 'entire'.
                Fitting/hazard code: Icode, run by name from the Icode folder behind
                assert(isequal(fileparts(which('BoutHazard')), pwd)). New 20r scripts live in
                Icode/analysis_20r/.
Gate            ALL of the following, else STOP and report which item failed:
                (1) 25 flies retained; 3989 wake bouts (xmin 2 s) and 3987 sleep bouts (xmin 300 s).
                    Exact.
                (2) The 16-bin hazard table in FINDINGS ("The hazard with bands"): bin centres to the
                    printed precision, and at-risk and event counts in every one of the 15 kept bins.
                    Counts exact.
                (3) The 2026-09-15 two-exponential null fit on the wake bouts
                    (FitHyperexponentialMLE, MaxComponents 2, discrete, SamplingInterval 1,
                    RandomSeed 1): tau 8.4 and 661.9 s, observed weights 0.079 / 0.921 (to the printed
                    precision), logL -29328.2460 to within 1e-4.
                (4) RiseRatio 1.405 / 1.388 / 1.376 / 1.430 at 12 / 16 / 20 / 24 bins (BoutHazard,
                    Spacing "log", MinAtRisk 10, Alpha 0.05), to the 3 decimals recorded.
                Also reported, not gated: per-fly bout counts, fly-days, and the trough/peak bins.
Output          The re-extracted bouts (wake and sleep, per fly, NaN-padded, with the file list) are
                saved to Icode/analysis_20r/per0DD_bouts_20r.mat; its SHA-256 goes in the ledger so
                the reference set cannot be lost again.
Decision rule   Pass: proceed to 20r-A exactly as written, using the same code for CS DD. Fail: no
                CS DD number is looked at; the failure is recorded in RESULTS_LEDGER.md.
Matched null    Not applicable (reproduction check).
Notes           Gate (2)-(4) are recorded statistics rather than the bouts themselves, so they check
                the bout SET only through what those statistics see; (1) and (2) together pin the
                counts in every log bin. Everything else in 20r-A, 20r-B and 20r-C is unchanged.
```

### 2026-10-31  20r gate amendment 2: the pre-cut per0 DD files, Teliminate 0

```
Question        Unchanged. The gate of "20r gate amendment: stored bout files not found" (same items
                (1)-(4), same tolerances) FAILED on the raw *_from1stLD.txt files with Teliminate 24
                (ledger, bd2871f): 16 of 25 flies kept, because nine 3.000-day files fell under
                minDays after the 24 h trim. This entry changes ONLY the input files and Teliminate.
Why             Sheyum and cClaude, 2026-10-31: the per0 DD set behind FINDINGS and all S/Q/W results
                is per0DD activity files_1sec res.zip (C:\Users\sheyu\Box\HMM paper\different
                genotypes\per0DD\flies\; SHA-256 778599dbcff6a94adc19e804b1b186e824f375833957b5a47bba40f0d74dce5f),
                25 files *_from1stLD.txtcutdata_4.8.2025.txt, already cut, each starting at subjective
                day. "Don't cut the first 24 hours" (Sheyum): Teliminate 0.
The cut         Measured, not assumed (analysis_20r/align_cut_20r.m, report
                analysis_20r/per0DD_cut_alignment.txt): every cut file is an exact contiguous piece of
                its raw *_from1stLD.txt, one placement only, at offset 0. The cut keeps the first
                min(N, 345600) samples: 4-7 day files are capped at 4.000 d (345600 s), the nine
                3.000-day files (259201 s) are kept whole. 16 x 4.000 d + 9 x 3.000 d = 91.0 fly-days.
Data            per0 DD: the 25 cutdata files from that zip, unchanged.
                CS DD (20r-A/B): CS_DD_activfiles.zip (SHA-256
                33495ce480fd4ededc020169288aafaeb08b3d9b7a18cdf03dbc1eeec3bb5c33), 25 raw
                *_from1stLD.txt files of 345601 s each, cut by THE SAME RULE: the first
                min(N, 345600) samples, i.e. the final sample dropped. Sheyum: both genotypes are cut
                the same way, starting at subjective day.
Fixed choices   process_LocomToSleep (current Rockefeller code, hashes as in the first amendment):
                timeInterval 1, bridgeWakeSeconds 1, Teliminate 0, deadThresh 12, minDays 3,
                phase 'entire'. 20r-B builds its synthetic traces from the same cut traces.
Gate            Items (1)-(4) and tolerances exactly as in the first amendment. Fail: stop, record,
                report; no CS DD number is looked at.
Output          analysis_20r/per0DD_bouts_20r_T0.mat (wake and sleep, per fly), SHA-256 in the ledger.
Matched null    Not applicable (reproduction check).
Notes           The cut was identified AFTER the first gate failed, from file provenance supplied by
                Sheyum, not by trying settings against the gate. It was not tuned: one rule, read off
                the alignment, applied unchanged to both genotypes.
```

### 2026-10-31  20r gate amendment 3: pre-09-22 pipeline, two pre-listed configurations, stopping rule

```
Question        Unchanged. Amendment 2 (pre-cut files, Teliminate 0, current code) FAILED narrowly
                (ledger 2e09b9e): 25 flies but 4017 wake / 4022 sleep against 3989 / 3987. 4017 equals
                the Python extraction of the documented rule, so the stored set most likely differs in
                edge handling. This entry tries the code that existed when the reference set was made.
Code            The Rockefeller repo at 04e0508^ (the state before "Sleep/wake pipeline: one coding
                convention", 2026-09-22): process_LocomToSleep 31d341d, ironout b54ade9,
                segment_sleep_phases f6c81f5, sleeptrace 2ca51d9, extractBoutsForModelCompetition
                22214ac, exported with git archive into a scratch folder and placed on the path INSTEAD
                of the current folder: addpath(..., '-begin'), which() asserted to resolve into that
                folder for every pipeline file listed by matlab.codetools.requiredFilesAndProducts,
                their which -all output recorded in the ledger, and the folder removed from the path
                when the run ends so later steps use the current code. outputCoding left at its
                default.
                Fitting/hazard code: Icode, by name, behind the which/pwd assert.
Data            per0 DD: the 25 cutdata files (zip SHA-256 778599db...), unchanged.
Configurations  Pre-listed, run IN THIS ORDER, nothing else:
                  C1  timeInterval 1, bridgeWakeSeconds 1, Teliminate 0, deadThresh 12, minDays 2
                  C2  the same with minDays 2.8
Stopping rule   Stop at the first configuration that EXACTLY matches 3989 wake / 3987 sleep and the
                FINDINGS 16-bin hazard table (at-risk and events in every kept bin). That
                configuration then goes through gate items (3)-(4) as in the first amendment; if it
                passes, 20r proceeds with it for both genotypes. If neither configuration matches,
                NO further settings are tried: switch to the replication route ("20r amendment: stored
                per0 DD bouts lost", to be written as its own entry before anything runs): one
                documented extractor for both genotypes, the 2026-09-15 per0 DD rise test repeated on
                the re-extracted bouts as a REPLICATION, then CS DD, 20r-B and 20r-C with that
                extractor.
Diagnostic      For every configuration, descriptive only, no hazard statistics beyond the gate items:
                per-fly wake and sleep bout counts; their min/max against the recorded 81-248 wake and
                80-249 sleep per fly (2026-09-15 entries); and per fly whether the first and last
                bouts of the record are kept (first kept bout starts at sample 1; last kept bout ends
                at the final sample).
Matched null    Not applicable (reproduction check).
Notes           minDays cannot remove a fly here unless the old code trims differently (all 25 have
                >= 3.000 d with Teliminate 0); it is listed because the recorded per-fly records were
                2.8-4.0 d. C2 runs only if C1 does not match (stopping rule).
```

### 2026-10-31  20r amendment: stored per0 DD bouts lost (replication route)

```
Why             The stopping rule of gate amendment 3 was reached (ledger 0d007e8): no pre-listed
                configuration of the current or pre-09-22 pipeline reproduces the stored 3989 wake /
                3987 sleep set, whose files are lost. The gate is therefore REPLACED by a replication.
                Text from Sheyum and cClaude, 2026-10-31.
Extractor       ONE extractor for both genotypes: the current Rockefeller pipeline (HEAD 001deca;
                process_LocomToSleep 7aec365d5b79c859, sleeptrace 0deb702eed5cbfcc, ironout
                f22f91c262f8c370, segment_sleep_phases a627f698159bc97a, extractBoutsForModelCompetition
                64f6f15985d04e5e), timeInterval 1, bridgeWakeSeconds 1, Teliminate 0, deadThresh 12,
                minDays 3, phase 'entire'. It implements the lab rule (sleep = >= 300 s of zero
                counts; a single active second bridged only when the inactive runs on both sides are
                each >= 300 s; wake bout = gap between consecutive sleep bouts; xmin 2 s wake, 300 s
                sleep, dt 1 s). Validated on hmm_audit/golden/G1_sleep_bouts.json (ledger bf40518):
                bridging agrees in every case; bouts touching the first or last sample are DROPPED as
                censored, the edge convention DESIGNATED BY SHEYUM on 2026-10-31 ("Drop edge
                bouts"), where the Python reference keeps them.
Data            per0 DD: the 25 cutdata files (zip 778599db...). Extracted: 4017 wake / 4022 sleep,
                90.45 fly-days covered by bouts (recorded set: 3989 / 3987, ~90).
                CS DD: CS_DD_activfiles.zip (33495ce4...), each file cut to its first min(N, 345600)
                samples (the cut measured on per0 DD, amendment 2), Teliminate 0. SHA-256 of every
                cut file recorded.
Steps           (1) Save the per0 DD bouts (wake and sleep, per fly) to
                    analysis_20r/per0DD_bouts_20r_replication.mat, SHA-256 in the ledger.
                (2) REPLICATION of "2026-09-15 per0 DD wake: the hazard rise at four fixed bin
                    counts", unchanged: null = FitHyperexponentialMLE MaxComponents 2, discrete,
                    SamplingInterval 1, RandomSeed 1, its K = 2 SurvivalHandle; BoutHazard NumBins
                    {12,16,20,24}, Spacing "log", MinAtRisk 10, Alpha 0.05, Bootstrap 0,
                    NullReplicates 999, RandomSeed 1. Same decision rule (all four: supported; some:
                    supported but bin-sensitive; none: not supported; NullDropped > 0: not
                    computable). Reported whatever it gives, next to the 2026-09-15 numbers.
                (3) Only then CS DD: 20r-A exactly as written (rise test, per-fly influence, sleep
                    control, CompareBoutModels), with this extractor.
                (4) 20r-B and 20r-C exactly as written, with this extractor. 20r-B's synthetic traces
                    go through the same pipeline (process_LocomToSleep, same settings). 20r-C's per0
                    fit is hyper_erlang with Shapes [1 1 5] (the FINDINGS shapes) refitted on the
                    re-extracted per0 bouts; responsibility of the series branch for a bout of k
                    bins = w_3 [S_3(k-1) - S_3(k)] / sum_j w_j [S_j(k-1) - S_j(k)], the discrete
                    likelihood's own convention.
                (5) 20r-A's per-fly influence and 20r-B also run on per0 DD's re-extracted bouts
                    where 20r-A/B call for per0 (influence: 20r-A asks only CS; not added).
Decision rules  Unchanged from 20r-A, 20r-B, 20r-C and the 2026-09-15 entry.
Matched null    As in those entries, at the re-extracted n.
Notes           FINDINGS' per0 DD numbers were computed on the lost set; the replication says whether
                they carry to the re-extracted set (28 more wake bouts). 20r-A's "hyper_erlang is
                pre-specified for CS DD" and its blind status are unchanged: no CS DD number has been
                looked at.
```

### 2026-10-31  20r-A follow-ups: power at CS DD's n, and the CS DD rise within 12-h clock halves

```
Question        (a) At CS DD's sample size, how often would the rise test detect the per0 DD rise?
                (b) Is a CS DD wake-hazard rise present within either 12-h half of the record clock,
                    masked in the pooled test by the free-running rhythm?
                Asked by cClaude after 20r-A (not supported), to tell "no rise" from "not detectable".
                Written after the 20r-A result was read; the questions follow from it, the settings
                below are fixed before either runs.
Data            (a) the per0 DD reference bouts, analysis_20r/per0DD_bouts_20r_replication.mat
                    (SHA-256 b14ba340...), 25 flies, 4017 wake bouts.
                (b) the CS DD bouts, analysis_20r/CSDD_bouts_20r.mat (SHA-256 80402183...), 3309 wake
                    bouts, with their onset bins (meta.OnsetBins, 1 = the first second of the cut
                    record, which starts at subjective day).
Time window     (a) each subsample, pooled. (b) each half, pooled over the 25 flies and 4 days.
Fixed choices   (a) n_CS = 3309. B = 199 subsamples, subsample b drawn with rng(3000 + b, 'twister'):
                    flies in random order (randperm), whole flies added while the running total stays
                    <= 3309, then the remainder taken as a uniform random subset, without replacement,
                    of the next fly's bouts. Each subsample through rise_test_20r unchanged EXCEPT
                    NullReplicates 199 (Sheyum, 2026-10-31: ~4 h instead of ~20 h; p in steps of
                    0.005): null = K = 2 hyperexponential fitted to that subsample (RandomSeed 1),
                    NumBins {12,16,20,24}, Spacing "log", MinAtRisk 10, BoutHazard RandomSeed 1.
                (b) A wake bout is in the first half if mod(onset - 1, 86400) < 43200 (hours 0-12 of
                    each record day), else the second. Each half through rise_test_20r unchanged
                    (NullReplicates 999). The 20p amendment-A split; the free-running period is not
                    estimated, so the halves approximate subjective day and night to ~1-2 h over 4 d.
Decision rule   (a) Descriptive: per bin count, power = fraction of the 199 with RiseNullP < 0.05; also
                    the fraction whose four p-values give "supported" and "supported or supported but
                    bin-sensitive". Reading: if power is high at 12 and 16 bins (where per0 detects
                    the rise), CS DD's null is informative; if low, it is not.
                (b) Each half read by the 2026-09-15 rule. Either half "supported" or "supported but
                    bin-sensitive": a rise is present within a half and the pooled CS DD verdict may
                    reflect masking. Neither: no rise within a half either. Two halves tested, no
                    correction applied; say so when reporting.
Matched null    (a) the per0 rise test's own null, at n = 3309, 199 replicates per subsample.
                (b) per half, 999 datasets of that half's n from the half's own K = 2 fit.
Notes           (b) inherits the censoring gap (CLAUDE.md step 0): bouts are assigned by START and not
                clipped, so a long bout starting late in a half carries time from the other half, and
                those long bouts sit in the tail where a rise would be read. The direction of the
                bias is not established; it must be stated with any (b) result.
                (a) subsamples share flies and bouts, so the 199 are not independent datasets; the
                power is an estimate for these flies, not for a new cohort.
```
