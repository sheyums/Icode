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
