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

### 2026-11-05  per0 DD: is there a time arrow in the sleep and wake bout sequences?

```
Question        In per0 flies in constant darkness, are successive bouts of the same type
                (sleep-to-sleep, wake-to-wake) asymmetric in time -- i.e. is the sequence
                irreversible -- beyond what a time-reversible process with the same bouts and
                the same memory produces?
Data            per0, constant darkness, raw per-second beam counts in
                C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data\per0 dd
                (53 .txt files; read only). Each file is one fly read from its FIRST sample, which
                is CT0 (Sheyum, 2026-11-05). 48 FLIES, chosen by Sheyum 2026-11-05: the 25
                "_from1stLD" files, plus the 23 plain files that have no "_from1stLD" version.
                Five plain files are EXCLUDED as the same recordings as their cut versions
                (checked: each cut file is an exact contiguous stretch of its plain twin, from
                sample 327432): 12_29_2017_multibeam_pruned2_{3,5,7,10,12}.txt. All 48 kept
                files are distinct recordings (content hashes all differ; no file's first day
                occurs inside any other). The kept files are copied, unchanged, into one staging
                folder that process_LocomToSleep reads with selectAll=true. SHA-256 prefixes
                (first 8 bytes) of the raw files:
                  from1stLD (suffix "_from1stLD.txt" omitted below), days:
                  12_29_2017_multibeam_pruned2_10.txt  4  41CC721B2D8DFD9D
                  12_29_2017_multibeam_pruned2_12.txt  4  C26D406A3AFDA43E
                  12_29_2017_multibeam_pruned2_3.txt   4  7F61EF0F4941A65D
                  12_29_2017_multibeam_pruned2_5.txt   4  8F820D66DCEE7BDD
                  12_29_2017_multibeam_pruned2_7.txt   4  59096AD2AEDC79B5
                  1_26_2018_multibeam_pruned_32.txt    7  912A907366266D2D
                  1_8_2018_multibeam_pruned_4.txt      4  4351F21E46301B75
                  1_8_2018_multibeam_pruned_5.txt      4  EC8EFFEEFD47A254
                  1_8_2018_multibeam_pruned_7.txt      4  4148B341F9B2DC6F
                  2_18_2018_multibeam_pruned2_54.txt   3  D6E41245FF3D58B1
                  2_18_2018_multibeam_pruned2_58.txt   3  97477867ACB393F0
                  2_6_2018_multibeam_pruned_14.txt     6  1AB8505CA0A9C512
                  2_6_2018_multibeam_pruned_7.txt      6  7F9A8D322962862A
                  5_28_2020_4b_wc_pruned_11.txt        7  628CBFB2E9BF446F
                  5_28_2020_4b_wc_pruned_2.txt         7  B1093D178E7D1E31
                  6_11_2020_4b_rc2_pruned_15.txt       5  FE56C92151D0A19B
                  6_25_2018_2boardmb_pruned_14.txt     3  22E5AE9EAC6DE5FC
                  6_25_2018_2boardmb_pruned_15.txt     3  4B064FF6B6CD5942
                  6_25_2018_2boardmb_pruned_27.txt     3  099EC89590525589
                  6_25_2018_2boardmb_pruned_28.txt     3  C6836B6E75DC6F40
                  6_25_2018_2boardmb_pruned_4.txt      3  AE48ADC8810EFC99
                  6_25_2018_2boardmb_pruned_8.txt      3  24E2B49AF98BE4F6
                  6_25_2018_2boardmb_pruned_9.txt      3  41166C90119451DD
                  6_25_2020_4b_rc_pruned_19.txt        4  015D0186F1EF96C1
                  6_25_2020_4b_rc_pruned_20.txt        4  86418C5C2A37923B
                  plain, days:
                  12_23_2017_multibeam_pruned2_{3,4,5,8,9,10,12,13,15,16}.txt  5.96 each
                    078CE757DBB9A9B6 5761334946DF1E0E 33E71866E76DC050 B71942395BC4E672
                    F0BC5F2FC89B81CB 9E37D0F888029ADF 72DFB11B69634349 26D1B1F3E6CD1147
                    116146277FA5A0F6 D9819B286C4A3E6E   (same order as the numbers)
                  12_29_2017_multibeam_pruned2_{4,6,8,11,13,15}.txt  7.90 each
                    CB19611F8835E93B B73B178D5780F383 630279263E3CC4F2 B7EDB264E6EC1F4F
                    D40B71FB9E31411D 2CB511618A124C2E
                  1_17_2018_multibeam_pruned_{18,...,24}.txt  6.48 each
                    A8B758F92EEED474 1FAA9F080288A015 CA00C730BAFE6691 452DF25E69EBB897
                    2A1B73F0F6B5DD9B 986DA5F3360AF6C3 D624866850220F6C
Time window     Each fly's whole retained record, one block per fly: the pipeline's dead-fly trim
                (deadThresh 12 h) and minDays 3 apply, nothing else is cut. A difference is never
                formed across flies. The whole record of one fly is the narrowest block every
                per-fly number is computed in; the group numbers pool these.
Fixed choices   Scoring: process_LocomToSleep(timeInterval=1, segmentMode=1, Teliminate=0,
                startZT=0, dayLength=24, deadThresh=12, minDays=3, bridgeWakeSeconds=1,
                selectAll=true) on the staged folder; then extractBoutsForModelCompetition(S,
                "entire", "sleep") and (S, "entire", "wake"). Sleep bouts >= 300 s, wake bouts
                >= 2 s (protocol facts). Before any statistic: assert meta.OnsetBins strictly
                increases along every row (bouts in time order); if not, stop and report.
                A fly enters a bout type's analysis if it has >= 40 bouts of that type.
                Test: BoutIrreversibility(bouts, Lags=[1 2 4 8], NumSurrogates=1000, Null="both",
                Seed=1, MinPairs=20, Plot=true), sleep and wake run separately. Code: Icode
                f46b912 or later, run from the Icode folder behind
                assert(isequal(fileparts(which('BoutIrreversibility')), pwd)).
Decision rule   PRIMARY (4 tests): the group UpFraction and the group Bowley skewness at lag 1,
                each against the IAAFT null, for sleep and for wake. Holm correction across the
                four at family-wise 0.05. Any primary rejection: a time arrow is supported for
                that bout type, and its SIGN (direction of the asymmetry) is reported, not
                predicted -- the test is two-sided. No rejection: no time arrow detected at lag 1
                at this sample size; that is a finding, not a cue to change lags or blocks.
                SECONDARY, descriptive, no correction: lags 2, 4, 8; the number of flies
                individually significant (p < 0.05, IAAFT) for each statistic; memory Rho
                against the shuffle null; moment skewness (no p-value, by design). The shuffle
                p-values for the skew statistics are reported but NOT used for any verdict:
                they over-reject once memory is present (measured, BoutIrreversibility header).
Matched null    IAAFT surrogates of each fly's own sequence (same bout values, approximately the
                same autocorrelation, time-reversible), 1000 per fly; the group null is the mean
                over the same flies of their r-th surrogates. Same n by construction.
Notes           Purpose: a model-free check of the "noise skew is irreversibility" reading
                (memory: project_noise_skew_irreversibility). The skew Sheyum reported is in
                bout-space residuals of an AR(1) on PIT-transformed bouts, which are skewed even
                for a reversible process; this test does not depend on that model. Calibration
                was measured only on simulated AR(1)-type processes (n = 150): real flies may
                have longer memory or regime switches. Power at lag 1 against the simulated
                irreversible process: UpFraction 0.875, Bowley 0.285 per fly -- so a Bowley null
                with an UpFraction rejection is not a contradiction. Nothing handles censoring;
                segment_sleep_phases already drops bouts touching the record ends. Blocks finer
                than the whole record (days, subjective day/night) are NOT part of this entry;
                if wanted, they need their own entry before they are run.
                Known open point, recorded rather than resolved: the five excluded plain
                12_29_2017 files are 7.90-day records whose cut versions START 327432 s (about
                91 h) into them, and six other 12_29_2017 plain files of the same length are
                used WHOLE here, per Sheyum's instruction that the per0 DD files start at CT0.
                If any plain file in fact opens with non-DD days, those days are inside its
                block.
```
