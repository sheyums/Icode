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

### 2026-11-05  CS LD: is there a time arrow in successive bouts beyond the day/night profile?

```
Question        In wild-type (Canton-S) flies in 12:12 light:dark, are successive bouts of the
                same type asymmetric in time BEYOND what each fly's own time-of-day profile
                produces -- i.e. does bout order carry an arrow that the light cycle and clock
                do not explain?
Data            CS, LD 12:12, raw per-second beam counts in
                C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data\CS ld
                (read only). All 40 files, every one a "_from1stLD" file; each is one fly read
                from its FIRST sample, which is ZT0 (Sheyum, 2026-11-05). Checked before writing
                this entry: 40 distinct recordings (content hashes all differ; no file's first
                day occurs inside another), and in all 8 experiments the largest 15-min activity
                bins fall within 30 min after ZT0 and ZT12 -- the lights-on and lights-off
                startles -- consistent with a ZT0 start. Files are copied unchanged into one
                staging folder. Name ("_from1stLD.txt" omitted), days, SHA-256 prefix:
                  10_1_2020_4b_rc2_pruned_55.txt         4.00  065FFA25C770DF19
                  10_5_2017_multibeam_pruned2_11.txt     5.00  F47A9793561E1DFF
                  10_5_2017_multibeam_pruned2_20.txt     5.00  4E4A8F876533B76D
                  10_5_2017_multibeam_pruned2_21.txt     5.00  804C786C4955817E
                  10_5_2017_multibeam_pruned2_23.txt     5.00  DAACD5A8F11F12CD
                  10_5_2017_multibeam_pruned2_24.txt     5.00  0B2A51496EE80E99
                  10_5_2017_multibeam_pruned2_27.txt     5.00  09900A9B1BFE23A4
                  10_5_2017_multibeam_pruned2_29.txt     5.00  7BA4308969B6D1D5
                  10_5_2017_multibeam_pruned2_30.txt     5.00  8FBE79BA2DAAAB59
                  10_5_2017_multibeam_pruned2_32.txt     5.00  0C506A9B5BFF6280
                  10_5_2017_multibeam_pruned2_4.txt      5.00  37D7F3F078549456
                  10_5_2017_multibeam_pruned2_9.txt      5.00  220646F92DB64F97
                  6_30_2017_multibeam_pruned10.txt       5.00  A6673E46C5985731
                  6_30_2017_multibeam_pruned11.txt       5.00  6DBF1878C7B5C325
                  6_30_2017_multibeam_pruned13.txt       5.00  6B7BE34FDD45F526
                  6_30_2017_multibeam_pruned6.txt        5.00  F07B9D48460878EC
                  6_30_2017_multibeam_pruned7.txt        5.00  9C6C972A81D48117
                  6_30_2017_multibeam_pruned9.txt        5.00  34BBAE255F5E29EA
                  6_9_2017_mb_pruned10.txt               5.00  280938BDDB3ACD3B
                  6_9_2017_mb_pruned15.txt               5.00  F7F5DB9302410443
                  6_9_2017_mb_pruned6.txt                5.00  EA8F9F4D22D40AEE
                  6_9_2017_mb_pruned8.txt                5.00  AF186CA3643F9B49
                  7_10_2017_multibeam_pruned15.txt       5.00  FD705DE7D0F9A62E
                  7_10_2017_multibeam_pruned5.txt        5.00  BA99998D8D5F7A49
                  7_10_2017_multibeam_pruned7.txt        5.00  69D8079519E9E2DE
                  7_10_2017_multibeam_pruned9.txt        5.00  F206A18D19E54F63
                  7_19_2018_2boardmb_pruned_24.txt       4.00  88EEC651D3F68100
                  7_19_2018_2boardmb_pruned_27.txt       4.00  0F881EB8C5FE3EE1
                  7_19_2018_2boardmb_pruned_29.txt       4.00  8EC51C265D191D1F
                  7_19_2018_2boardmb_pruned_31.txt       4.00  57FF785A981253A2
                  8_25_2017_multibeam_pruned25.txt       5.00  6982E8E396F3E25E
                  8_25_2017_multibeam_pruned26.txt       5.00  954313FD0BCCBF28
                  8_25_2017_multibeam_pruned28.txt       5.00  AB22ED4D6A31C341
                  8_25_2017_multibeam_pruned31.txt       5.00  E097376D7717C3B7
                  8_25_2017_multibeam_pruned32.txt       5.00  478FB83C5FE3B5CE
                  9_29_2017_multibeam_pruned2_10.txt     5.00  BE763A5F7F663D18
                  9_29_2017_multibeam_pruned2_6.txt      5.00  B8778313A57A812E
                  9_29_2017_multibeam_pruned2_7.txt      5.00  DEED7E8FB1F3D3F0
                  9_29_2017_multibeam_pruned2_8.txt      5.00  B65DC8429F44835F
                  9_29_2017_multibeam_pruned2_9.txt      5.00  74EB8F2CC1ED1813
Time window     Each fly's whole retained record as ONE block (4-5 days); pairs DO span lights-on
                and lights-off, because the null keeps every position's time of day and so
                carries the transitions. Strata are 1-h ZT bins of bout onset. The fly's whole
                record is the narrowest block of every per-fly number; the count pools flies.
Fixed choices   Scoring: process_LocomToSleep(timeInterval=1, segmentMode=1, Teliminate=0,
                startZT=0, dayLength=24, daytime=12, nighttime=12, deadThresh=12, minDays=3,
                bridgeWakeSeconds=1, selectAll=true) on the staging folder; then
                extractBoutsForModelCompetition(S, "entire", "sleep") and (S, "entire", "wake").
                Sleep >= 300 s, wake >= 2 s. Assert OnsetBins strictly increase along every row.
                A fly enters a bout type's analysis if it has >= 40 bouts of that type.
                Strata: bout onset ZT = mod((OnsetBins - 1)*1/3600 + 0, 24) hours, stratum =
                floor(ZT), i.e. 24 one-hour bins, built separately for sleep and wake.
                Test: BoutIrreversibility(bouts, Strata=strata, Lags=[1 2 4 8],
                NumSurrogates=1000, Null="both", Seed=1, MinPairs=20), sleep and wake run
                separately. Code: Icode at the commit adding the strata null or later, run from the
                Icode folder behind assert(isequal(fileparts(which('BoutIrreversibility')), pwd)).
Decision rule   PRIMARY (2 tests): for sleep and for wake, the number of flies whose per-fly
                UpFraction p at lag 1 against the STRATA null is < 0.05, against
                Binomial(NFlies, 0.05), upper tail (UpFraction_CountP_Strata). Holm across the two
                at family-wise 0.05. A rejection: order beyond the day/night profile carries an
                arrow in more flies than chance; report how many flies, and in which direction
                (UpFraction above or below its strata null). No rejection: no arrow beyond the
                profile detected by this per-fly count. That is a finding, not a cue to change
                strata, lags or blocks.
                SECONDARY, descriptive, no correction: lags 2, 4, 8 (same count); per-fly Bowley
                against the strata null (number of flies only, no count p -- its per-fly rate
                under the null was measured at about 0.068); memory beyond the profile (per-fly
                Rho against the strata null); the IAAFT and shuffle results, reported as what a
                clock-blind test says, NOT used for any verdict (measured to be fooled by the
                profile: per-fly UpFraction rejection 0.070-0.135 on reversible LD-like flies).
Matched null    Per fly, 1000 profile-matched surrogates: the fly's own bouts permuted only among
                positions in the same 1-h ZT bin, across all its days. Same n, same bouts, same
                time-of-day profile. The count's null is Binomial(NFlies, 0.05).
Notes           Why per-fly counts and not a group mean: the strata null's group mean was
                measured to over-reject (0.15-0.18 on 10-fly LD-like groups) because permuting
                within ZT bins cuts the link between a bout's length and the next bout's onset
                time; the count was measured calibrated for UpFraction (0.030 on 40-fly groups,
                2.06 flies significant on average against 2.00 expected). Sheyum chose per-fly
                tests over building a timeline-regenerating null (2026-11-05).
                POWER: the count is far weaker than the per0 DD group-mean test. Per-fly power
                was 0.535 against a strong simulated arrow; an effect the size of per0 wake (one
                percentage point of increases, 2 of 48 flies individually significant) would be
                invisible here. A null result here therefore does NOT contradict per0 DD, and the
                two results are not directly comparable. For comparison only, the per0 DD per-fly
                counts at lag 1 were sleep 4/48 and wake 2/48 (IAAFT null, not strata).
                ASSUMPTION for review: Teliminate=0 (no entrainment trim), reading from ZT0 as
                for per0 DD; the 2026-08-30 CS work used Teliminate=24. Calibration was simulated
                for sleep bouts only, with one profile shape; wake bouts and other profiles were
                not simulated. Nothing handles censoring.
```

### 2026-11-05  per0 DD and CS LD: can x(t+1) = A + B*x(t) + eps generate the bouts, and what is eps?

```
Question        For each fly and bout type, (1) does a linear AR(1) in BOUT SPACE (seconds, no
                transform) with noise independent of the previous bout describe successive bouts,
                and (2) what is the distribution of its noise eps? This is the question the
                collaborator's PIT/AR(1) scheme set out to answer, asked directly.
Data            (a) per0 DD: the same 48 flies, files, staging and scoring as the 2026-11-05 entry
                "per0 DD: is there a time arrow ..." (committed 2df6bef), from the first sample.
                (b) CS LD: the same 40 flies, files, staging and scoring as the 2026-11-05 entry
                "CS LD: is there a time arrow ... beyond the day/night profile?" (58ac4d7), from ZT0.
Time window     per0 DD: each fly's whole record is one block. CS LD: the two circadian phases are
                fitted SEPARATELY (Sheyum, 2026-11-05): a bout belongs to the phase of its onset
                (light ZT 0-12, dark ZT 12-24); every 12-h phase of every day is its own block, so
                no pair spans lights-on or lights-off; one fit per fly per phase uses all that
                phase's blocks. The narrowest block of a number is the fly (per0) or the fly's
                light or dark phase (CS).
Fixed choices   Bouts: as in the two entries above (bridged 1 s; sleep >= 300 s, wake >= 2 s; a fly
                enters a bout type if it has >= 40 bouts of that type; OnsetBins order asserted).
                CS onset ZT = mod(OnsetBins - 1, 86400)/3600; day = floor((OnsetBins - 1)/86400).
                BoutAR1Noise(bouts, Xmin = 300 (sleep) or 2 (wake), NumSurrogates = 1000, Seed = 1,
                MinPairs = 20, ReturnResiduals = true); per0 with the shuffle null; CS with
                Strata = floor(onset ZT) (1-h bins) within each phase. Then FitResidualFamilies on
                each fly's (per phase's) residuals: shifted exponential, Normal, Laplace,
                asymmetric Laplace, shifted lognormal, shifted Weibull, GenNorm K = 1, GH family.
                Code: Icode 328116c, run from the Icode folder behind the which/pwd assert (the six
                noise fitters resolve to Icode's copies there).
Decision rule   PRIMARY, question (1): the count of flies whose Hetero p (Spearman of |eps| with
                x(t), against the null) is < 0.05, against Binomial(NFlies, 0.05). Tests: per0
                sleep, per0 wake (Holm across these two); CS sleep-light, sleep-dark, wake-light,
                wake-dark (Holm across these four). A rejection: the noise's spread depends on the
                previous bout in more flies than chance, so additive iid noise is wrong there.
                No rejection: the linear model's independence assumption is not contradicted.
                Question (2), DESCRIPTIVE with a pre-set wording rule: per fly, the BIC winner among
                the eight families and each family's DeltaBIC from the best. If the shifted
                exponential wins in more than half the flies of a set, eps is reported as "the bout
                distribution shifted by a constant" for that set; otherwise the family winning most
                often is named, with its count, and the median DeltaBIC of the shifted exponential
                is reported. Calibration behind this: the shifted exponential won 100/100 when true
                (test_FitResidualFamilies), so a flexible family winning is not chance flexibility.
                SECONDARY, descriptive, no correction: B per fly (with B_P; per0 group mean against
                the shuffle null; CS counts only); ShapeKS count (fires whenever there is memory, so
                read as "the dynamics change eps"); BelowHalfB; regeneration -- the fraction of
                generated bouts below the floor (300 s sleep, 2 s wake; never clipped), the
                smallest generated bout, and GenKS against GenKS_Ref; residual SD, skewness and
                Bowley skewness against their no-dynamics values.
Matched null    Per fly (per phase), 1000 surrogates of its own bouts: within-block shuffles (per0)
                or permutations within 1-h ZT bins across the phase's days (CS). Same n and bouts.
                Count null Binomial(NFlies, 0.05).
Notes           Measured (BoutAR1Noise header): Hetero rejects 0.045 on iid bouts, 0.085 on EAR(1)
                (linear, flat spread) and 0.695 on a transformed-Gaussian AR(1) -- the collaborator's
                model class. In LD the day/night profile alone produced B ~ 0.18 and shuffle-null
                rejection 0.465; with strata 0.030. That LD calibration was simulated on whole
                records, NOT per phase; per-phase blocks were not simulated. CS phases have about
                half the bouts of the whole record (and lose one pair per block), so power per fit
                is lower than in per0. OLS on exponential-like bouts is dominated by the long bouts.
                eps are treated as iid by the family fits although they overlap in time.
                Expectation stated before the run, not a decision rule: per0 memory is very weak
                (lag-1 rank correlation +0.03 sleep, -0.06 wake), so B ~ 0 and eps ~ the bout
                distribution minus a constant is the anticipated answer for per0; a negative B
                alone breaks the 2-s wake floor on regeneration (measured), so per0 wake
                regeneration may put bouts below 2 s.
```

### 2026-11-06  Four homeostat analyses: shared data and preliminaries (applies to entries A-D)

```
Datasets        All raw per-second beam counts under
                C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data
                (read only). 140 recordings, all distinct (content hashes; no file's first day
                occurs inside another; checked across the four sets on 2026-11-06).
                per0 DD  48 flies: as in the 2026-11-05 per0 DD entries (2df6bef), from sample 1.
                CS LD    40 flies: as in the 2026-11-05 CS LD entry (58ac4d7), from ZT0.
                CS DD    25 flies, folder "CS DD", the 25 "_from1stLD" files (4.00 days each). The
                         zip in that folder holds the same flies (Sheyum) and is not used. Pooled
                         activity peaks 10.25-11.75 h after the record start, consistent with a CT0
                         start and an evening peak. START PHASE: CT0 (confirmed by Sheyum,
                         2026-11-06).
                per0 LD  27 flies, folder "per0 ld", plain files, 4.74-8.71 days. Sheyum's
                         impression: records start at or near ZT0. Checked 2026-11-06 from the
                         lights-on/off startles (per0 has no clock, so its activity peaks in LD
                         mark the light transitions; per0ld_phase_check_2026-11-06.png in the
                         analysis folder): 9_21_2020 (11 flies) and 8_18_2017 (1) show
                         transitions within about 0.5-1.5 h of the record start -- consistent
                         with a ZT0 start; 8_25_2017 (3), 9_22_2017 (3) and 9_29_2017 (4) show
                         stripes at about 4 h and 16 h every day -- transitions about 4 h off a
                         ZT0 start; 9_1_2017 (5) is unclear (one stripe near 14 h). (The
                         2026-08-30 work estimated 9_21_2020 at ZT 11.5 from a different copy of
                         that experiment, Per0_LD_all.) START PHASE, decided 2026-11-06 (Sheyum:
                         estimate from the startles): EstimatePhaseStamp(folder,
                         ClockState="arrhythmic"), analysis folder copy, gives startZT
                           8_18_2017 11.0 (1 fly; = the schedule-validated 2026-08-30 value),
                           9_21_2020 11.5 (11; = the schedule-validated 2026-08-30 value; partner
                                     rise 6.5 h away, anchor margin 10.9),
                           8_25_2017 20.0 (3; partner rise 12.0 h away, no flag),
                           9_22_2017 20.5 (3; no flag),
                           9_29_2017  8.0 (4; no flag),
                           9_1_2017   7.5 (5; FLAG partner rise 9.0 h away) -> DROPPED as not
                                     reliably phased (no schedule).
                         None starts at ZT0. per0 is more active in the dark (EstimatePhaseStamp
                         header: the sleep-in-dark measure is inverted for arrhythmic flies),
                         which is why 9_21_2020 is busy in its first 12 h. Teliminate per
                         experiment = hours to its first ZT0: 13, 12.5, 4, 3.5, 16 respectively,
                         startZT as above, so every retained record starts at ZT0. 22 flies. Each
                         experiment is processed in its own staging folder (startZT is one value
                         per call). Every retained day must show its lights-on and lights-off
                         rises within 1 h of the stated ZT (checked at run time; a failing
                         experiment is dropped and named in the ledger).
Scoring         process_LocomToSleep(timeInterval=1, segmentMode=1, startZT and Teliminate as
                above (both 0 for the three ZT0/CT0 sets), dayLength=24, daytime=12,
                nighttime=12, deadThresh=12, minDays=3, bridgeWakeSeconds=1); sleep and wake bouts
                with OnsetBins from extractBoutsForModelCompetition; sleep >= 300 s, wake >= 2 s.
                A fly enters if it has >= 40 bouts of each type.
Alternation     Sleep and wake bouts are merged in onset order per fly. Two bouts are ADJACENT
                only if the later starts exactly where the earlier ends (onset = previous onset +
                previous duration, in bins); non-adjacent pairs are never formed, and their
                number is reported. Which state starts each record is read from the onsets.
Time of day     "Profile sets": CS LD, per0 LD (ZT) and CS DD (CT, assuming a 24-h free run).
                The CS free-running period is NOT estimated per fly, so over 4 days its rhythm
                can slide against the CT bins (about 1 h per 0.25 h of period error): a stated
                limitation. per0 DD has no profile.
Code            New functions, each with a test suite and a calibration on simulated flies,
                written, committed and run BEFORE any data analysis: BoutCrossType (A),
                SleepAmountACF (C), ProcessSProxy (D). B uses existing code. If a calibration
                fails (null rejection above 0.095 at nominal 0.05 over 200 replicates), the
                affected test is NOT run as a verdict: a dated amendment below the entry,
                committed before the run, demotes it to descriptive.
```

### 2026-11-06  A. Cross-type homeostat test: does a sleep bout depend on the wake bout before it more than on the one after it?

```
Question        Per fly, is a sleep bout's length correlated with the preceding wake bout's length
                (homeostat: longer wake -> more pressure -> longer sleep), and is that correlation
                larger than with the FOLLOWING wake bout (the time arrow a homeostat predicts;
                reversing time would swap the two)?
Data            All four datasets (preliminaries). Whole record, one block per fly.
Statistics      Per fly, over adjacent triples (wake w_pre, sleep s, wake w_post):
                r_pre = Spearman(s, w_pre); r_post = Spearman(s, w_post); Delta = r_pre - r_post.
                Mirror for wake bouts: q_pre = Spearman(w, s_pre), q_post = Spearman(w, s_post),
                DeltaW = q_pre - q_post.
Nulls           per0 DD: circular shift of the fly's wake sequence against its sleep sequence, all
                shifts k with 2 <= k <= n-2 (keeps each sequence's own memory, breaks the pairing).
                Profile sets: 1000 permutations of the fly's wake bouts among wake positions in
                the same 1-h ZT/CT bin of onset, across days, sleep bouts fixed.
Decision rule   PRIMARY, per0 DD: group mean of r_pre (one-sided, > 0) and group mean of Delta
                (one-sided, > 0) against the shift null's group means; Holm across the two.
                PRIMARY, profile sets (only if their calibration passes): count of flies with
                per-fly one-sided p < 0.05 for r_pre, and for Delta, against Binomial(n, 0.05);
                Holm across the six (3 sets x 2).
                Reading: r_pre > 0 AND Delta > 0 supports a homeostat. r_pre > 0 with Delta ~ 0 is
                a symmetric co-variation of sleep and wake (a shared slow state), not a homeostat.
                SECONDARY: r_post, q_pre, q_post, DeltaW; per-fly values; lag-2 versions (s with
                the wake two bouts back).
Calibration     Before the run, 200 simulated flies each: (i) alternating independent exponential
                bouts (every test ~0.05); (ii) a shared slow state, no homeostat (sleep and wake
                means modulated by one latent AR(1)): r_pre > 0 expected, Delta must stay
                calibrated; (iii) LD-like flies with a ZT profile and no coupling: the profile
                null must stay <= 0.095 (the timeline-coupling bias measured for the strata null
                in BoutIrreversibility may reappear: permuting a wake bout cuts its link with the
                next sleep's onset time); (iv) a homeostat (Process-S-modulated rates) for power.
Matched null    Shifts or permutations of each fly's own bouts, same n; count null
                Binomial(n, 0.05).
Notes           The fly-labelled per0 DD subset in C:\Users\sheyu\Icode had an unknown first state
                for 13 flies; this entry uses the raw files, where onsets resolve it.
```

### 2026-11-06  B. Cox-process check: do memory and state-dependent noise go together, as a slowly drifting rate predicts?

```
Question        If each bout is exponential with a mean set by a slowly varying hidden state (a
                doubly stochastic, Cox-type process), flies with positive memory B should show
                noise that SPREADS after long bouts (Hetero > 0), and negative B should go with
                shrinking spread. Across flies, are B and Hetero positively related beyond what
                independent bouts give at the same n?
Data            The per-fly B and Hetero of the 2026-11-06 ledger entry (625dc08): per0 DD sleep
                and wake (primary); CS LD sleep/wake x light/dark (secondary).
Statistic       Per set, Spearman correlation across flies of B with Hetero.
Null            1000 simulated sets, each fly matched to its real n (pairs), xmin and mean bout,
                bouts independent exponential; B and Hetero by the same formulas (checked to equal
                BoutAR1Noise's on the real data to 1e-12). This measures how far the two
                estimators are tied MECHANICALLY, with no biology.
Decision rule   PRIMARY: per0 sleep and per0 wake, one-sided (> 0) Monte Carlo p against that
                null; Holm across the two. Support: the relation exceeds the mechanical one.
                SECONDARY: the four CS LD sets; and the correlation a Cox process would give --
                simulated flies whose log-mean follows a latent AR(1) (phi 0.5), its SD set so the
                set's median B is reproduced.
Notes           NOT BLIND. The sign concordance that motivated this (per0 sleep: 9 of 11 flies with
                significant B positive, 7 of 8 with significant Hetero growing; wake: 7 of 9
                negative, 4 of 6 shrinking) was seen in the 2026-11-06 results before this entry.
                The test uses ALL flies and a new null, but a positive result confirms a pattern
                already noticed and must be reported that way.
```

### 2026-11-06  C. Sleep amount in windows: compensation beyond independent bouts?

```
Question        In per0 flies in constant darkness, are the minutes of sleep in one window
                followed, hours later, by deviations the other way (compensation, the hallmark of
                a homeostat), more than the same bouts in random order produce?
Data            per0 DD only (48 flies): the one dataset with neither a clock nor light. The
                profile sets need a null that keeps a free-running or entrained rhythm together
                with its timeline, which does not exist yet; they are NOT part of this entry.
Statistics      Per fly, from the first bout onset to the last bout end: sleep minutes per window,
                30-min windows (PRIMARY) and 60-min (secondary), both reported. Autocorrelation at
                lags 0.5-12 h. Deviation curve = observed ACF minus the null mean ACF.
Null            Renewal surrogates, 1000 per fly: sleep bouts permuted among the fly's sleep
                positions and wake bouts among its wake positions (alternation and first state
                kept), the trace rebuilt end to end and binned identically. Keeps both bout
                distributions; removes all sequence dependence. Binning's own positive ACF (a long
                bout fills several windows) is in the null.
Decision rule   PRIMARY: the minimum over lags 1-12 h of the group-mean deviation curve, against
                the same minimum on the surrogates' group means (a max-statistic null, which
                corrects for scanning lags); one-sided (more negative than the null). Support:
                compensation beyond independent bouts, at the lag of the minimum. SECONDARY: the
                maximum (persistence), the same way; 60-min windows; per-fly curves; the count of
                flies whose minimum is below their own null 5th percentile.
Calibration     Before the run, 200 simulated per0-like flies each: independent alternating
                exponential bouts (<= 0.095), and a homeostat (Process-S-modulated rates) for power
                and to see where the minimum appears.
Notes           Records differ in length (3-7.9 days), so long lags rest on fewer windows in short
                records; the group mean weights flies equally.
```

### 2026-11-06  D. A sleep-pressure proxy: does accumulated wake predict the next sleep bout?

```
Question        In per0 DD, does a Borbely-type pressure proxy S built from each fly's own history
                predict how long the next sleep bout lasts (and the next wake bout, inversely),
                beyond independent bouts -- and on what timescale?
Data            per0 DD only (48 flies), for the reason given in C.
Proxy           S in [0, 1], computed exactly bout by bout: during wake S <- 1 - (1 - S)*exp(-dt/tw);
                during sleep S <- S*exp(-dt/ts); S = 0.5 at the record start; bouts starting in
                the first 12 h are excluded from the statistics. Grid fixed now: tw and ts each in
                {0.5, 1, 2, 4, 8, 16} h (36 pairs).
Statistics      Per fly and grid pair: rho_S = Spearman(S at sleep onset, that sleep bout's
                length); rho_W = Spearman(S at wake onset, that wake bout's length).
Null            The renewal surrogates of C (1000 per fly), S recomputed on each.
Decision rule   PRIMARY: the group mean of rho_S maximised over the 36 grid pairs, against the
                same maximum on the surrogates' group means (max-statistic null, corrects the grid
                scan); one-sided (> 0). SECONDARY: rho_W minimised the same way (< 0); the best
                (tw, ts), descriptive only (a grid peak is not an estimate with an interval);
                whether S adds to the immediately preceding wake bout (rho_S against A's r_pre,
                per fly).
Calibration     Before the run: independent alternating bouts (<= 0.095); a homeostat with known
                tw and ts (power, and whether the grid maximum lands near the truth); a shared slow
                state without a homeostat (does S pick it up? reported).
Notes           With short tw, S at sleep onset is mostly the preceding wake bout, so D overlaps A
                at that end of the grid; long time constants carry several bouts of history, which
                A cannot see. A positive result says the proxy predicts sleep; it does not
                identify the circuit that integrates pressure.
```

### 2026-11-06  AMENDMENT to entries A, C, D (calibration results; written and committed BEFORE any data run)

```
Basis           The calibrations the preliminaries entry required, run on SimulateSleepWake flies
                (independent bouts; a shared slow state with no homeostat; a homeostat with Tw 2 h,
                Ts 1 h, gain 2; an LD-like profile) with the code committed alongside this
                amendment (test_BoutCrossType 9/9, test_SleepAmountACF 5/6, test_ProcessSProxy 6/6).
                Rejection rates at nominal 0.05 (rule: above 0.095 = fails):
                A per0 DD (shift), 48-fly groups x 100: r_pre 0.050, Delta 0.040 -> calibrated.
                  Homeostat power (50 groups): r_pre 1.000, Delta 0.980. Shared slow state, per
                  fly (200): r_pre 0.375, Delta 0.020.
                A profile sets (strata), 40-fly count test x 100: r_pre 0.060 -> calibrated;
                  Delta 0.120 -> FAILS; DeltaW 0.150 and r_pre2 0.120 (secondary) also above 0.095.
                C independent bouts, 48-fly groups x 100: compensation (MinP) 0.100 -> FAILS;
                  persistence (MaxP) 0.030. Homeostat: 1.000 (minimum at 1.0 h). Shared slow state:
                  0.580.
                D independent bouts, 48-fly groups x 60: rho_S 0.017, rho_W 0.017 -> calibrated
                  (conservative). Homeostat: 1.000 and 1.000. Shared slow state: rho_S 1.000,
                  rho_W 0.367. Best grid pair: homeostat (2, 1) h in 9/30 groups; shared slow
                  state (0.5, 0.5) h in 10/30.
Changes         A, per0 DD: unchanged; both primaries stand. The reading rule stands and is now
                  measured to work: a shared slow state gives r_pre > 0 but not Delta > 0.
                A, profile sets: the Delta count is DEMOTED to descriptive. The r_pre count stays
                  primary, Holm across the three sets (not six). A positive r_pre there can come
                  from a shared slow state as well as a homeostat, so it will be reported as
                  coupling between a wake bout and the next sleep, NOT as evidence of a homeostat.
                C: the whole entry is DEMOTED to descriptive (compensation over-rejects on
                  independent bouts, and a shared slow state also produces it 58% of the time).
                  Curves, extremes and p-values are reported with these calibration numbers and
                  no verdict.
                D: the primary test stands (calibrated), but its READING changes: a positive
                  result says a slow variable links successive bouts -- a homeostat OR a shared
                  slow state -- not a homeostat. The best (tw, ts) is reported descriptively,
                  noting that the two models placed it differently in simulation.
                B: unaffected (its null is its own simulation).
Consequence     The only pre-specified test that is both calibrated and specific to a homeostat is
                A's Delta on per0 DD.
```

### 2026-11-06  E. Is cumulative sleep mean-reverting? Variance-time against a renewal null

```
Question        In per0 flies in constant darkness, does the running sleep surplus (time asleep
                minus the fly's average sleep fraction times elapsed time) wander LESS over long
                stretches than the same bouts in random order allow -- i.e. are sleep deficits and
                surpluses repaid (mean reversion, the regulation of sleep AMOUNT a homeostat
                implies) -- and on what timescale?
Data            per0 DD only (48 flies), as in entries C and D (preliminaries entry, fe307b0); the
                rhythmic sets need a null that keeps a free-running or entrained rhythm, which does
                not exist, and a 24-h rhythm dominates long-scale variance.
Statistic       SleepVarianceTime (Icode e327aab or later): per fly and scale L in {1, 2, 4, 8, 16,
                24} h, V(L) = variance over start times (15-min grid, overlapping windows) of the
                sleep surplus change over L; logR(L) = log(V(L) / mean renewal-surrogate V(L)).
                Per fly, Slope = least-squares slope of logR on log2(L). Group: mean Slope.
Null            1000 RenewalSurrogate draws per fly (Seed 1, substream = fly index: the SAME
                surrogates as entries C and D).
Decision rule   PRIMARY: group mean Slope, one-sided, more negative than the surrogates' group
                means -> mean reversion that grows with timescale (SUPPORTED if p < 0.05). The
                opposite tail (persistence) is reported as the second, two-sided reading of the
                same statistic, not as a separate test. SECONDARY, descriptive: logR at each
                scale (where the bend happens); per-fly slopes and p; the runs test per bout type
                (more runs than the surrogates = reversion at the bout level).
Calibration     MEASURED before this entry was committed (test_SleepVarianceTime 6/6, Icode
                e327aab; SimulateSleepWake 48-fly groups, 199 surrogates; one-sided at 0.05):
                  independent bouts (100 groups): reversion 0.040, persistence 0.020 -> calibrated;
                    median group Slope -0.031 (log ratio drifts below 0 at long scales even under
                    the null: a log-of-variance bias the null shares).
                  homeostat, Tw 2 h, Ts 1 h, gain 2 (50): reversion 1.000; median Slope -0.235,
                    logR -0.10 at 1 h down to -1.16 at 24 h.
                  shared slow state, no homeostat (50): reversion 0.140; median Slope -0.038 (about
                    the independent-bout value); persistence 0.020 -- it does NOT bend up.
                  runs test: shared slow state 1.000 for sleep and wake, homeostat 0.800 / 1.000:
                    not specific, so the runs test stays descriptive.
                Reading rule fixed by this: a supported reversion is reported as regulation of
                sleep AMOUNT, with its Slope set against the two simulations (homeostat about
                -0.24, shared slow state about -0.04); it is called homeostat-like only if the
                Slope is below -0.10, a value the shared slow state did not reach in median
                (and is reported with that caveat, since only one shared-state structure was
                simulated).
Notes           Records are 3-7.9 days, so the 16- and 24-h scales rest on few independent
                stretches in short records; the group mean weights flies equally. A renewal null
                already contains the bouts' own variance-time behaviour, so the test is about
                ORDER only. The descriptive dip of entry C at 1-1.5 h is the short-lag view of the
                same question; E integrates over lags instead of scanning them.
```

### 2026-11-09  F. Collaborator's latent-AR residuals: does the model reproduce the noise shape, and does memory change it?

```
Question        In Shreenath's latent-AR(1) (PIT) model of CS LD sleep bouts, (1) do the observed
                bout-space residuals have the skewness the fitted model produces, and (2) does memory
                (fitted phi) change residual skewness compared with phi = 0? (His reading: "memory
                reduces the skewness".)
Data            His export, read only: C:\Users\sheyu\Box\Data from Shreenath\
                memory_vs_no_memory_residuals_raw_export\{day,night}\{residuals,sim_residuals}.csv
                and notes.json (dated 2026-11-06; sizes day 964162 / 432764116 B, night 729130 /
                326537424 B). CS LD sleep bouts, K = 1; day 56 flies, night 55. Residual =
                observed minus predicted duration (s); prediction = empirical inverse CDF of
                Phi(refitted intercept + phi * z), a conditional median. Simulations: 200 replicates
                per fly per version, "fitted" (his phi) and "phi0" (iid latent, same empirical
                marginal), paired draws, each scored by refitting exactly like the data. His
                extraction and scoring are used as delivered (bridging and record start not stated
                in the export). 24 of his 56 flies are in our 40-fly CS LD set; no cross-matching here.
Time window     His day and night phases separately; each fly's residuals within a phase.
Statistics      Per fly, phase and sequence (observed, or one replicate of a version): moment
                skewness and Bowley skewness (Parzen mid-quantiles) of residual_seconds, and SD.
Decision rule   PRIMARY (4 tests, Holm across them):
                (1) model adequacy, per phase: per fly, two-sided Monte Carlo p of the observed
                    moment skewness among the fly's 200 "fitted" replicate values; count of flies
                    with p < 0.05 against Binomial(n, 0.05). Exactly calibrated if the model is right
                    (replicates are exchangeable with the data). Rejection: the noise's shape is NOT
                    what the fitted model produces.
                (2) memory effect, per phase: per fly, mean moment skewness over "fitted" replicates
                    minus over "phi0" replicates; across flies, two-sided sign-flip permutation test
                    of the mean difference (10000 flips, seed 1). Rejection: memory changes the
                    noise's skewness, in the direction of the sign.
                SECONDARY, descriptive: the same with Bowley skewness and SD; observed vs phi0;
                per-fly phi; effect sizes (mean difference and its ratio to the phi0 skewness).
Matched null    (1) the fly's own fitted-model simulations, same n and scoring; (2) the paired phi0
                simulations of the same fly.
Notes           Expected from our bout-space results (not a decision rule): phi is tiny (|phi| <
                ~0.1), so memory should change skewness by a few percent at most. A change in the
                best-fitting family (GH vs shifted exponential) is NOT tested here: family choice is
                not a skewness measure.
```

### 2026-11-09  G. Collaborator's latent residuals: are they Gaussian, as his model assumes?

```
Question        In Shreenath's latent-AR(1) model of CS LD sleep bouts, do the observed latent
                residuals (latent coordinate minus its AR(1) prediction) depart from the Gaussian
                innovations the model assumes -- in skewness or tail weight -- beyond what his own
                fitted-model simulations, scored identically, produce?
Data            The export of entry F (same files; column latent_residual), day 56 / night 55 flies.
Statistics      Per fly and phase: moment skewness and excess kurtosis of latent_residual, for the
                observed sequence and for each of the 200 "fitted" replicates.
Decision rule   PRIMARY (4 tests, Holm): per phase and statistic, the count of flies whose observed
                value is outside its 200 fitted replicates (two-sided Monte Carlo p < 0.05), against
                Binomial(n, 0.05). Rejection: the latent noise is not Gaussian in that respect, i.e.
                the dependence between successive bouts is not the Gaussian (copula) kind.
                SECONDARY, descriptive: mean skewness, excess kurtosis and SD (observed vs fitted vs
                phi0); SD against sqrt(1 - phi^2).
Matched null    Each fly's own fitted-model simulations (Gaussian innovations), same n and scoring.
Notes           The randomized PIT makes each fly's latent MARGINAL Gaussian by construction, so a
                departure can arise only from the dependence structure; with phi ~ 0.03 (entry F) the
                expected outcome is no departure, and a null result confirms the construction rather
                than adding biology. Skewness of latent residuals is not by itself a time arrow.
```

### 2026-11-09  H. Collaborator's model: do extreme sleep bouts cluster beyond what Gaussian dependence allows?

```
Question        In CS LD sleep bouts, is a very long (or very short) bout followed by another very
                long (very short) bout more often than Shreenath's fitted Gaussian latent AR(1)
                produces? A Gaussian copula has no tail dependence, so excess clustering of extremes
                would be the one dependence feature his model cannot carry (candidate: a shared slow
                state); no excess means a non-Gaussian dependence model has nothing to find here.
Data            The export of entries F and G (same files), day 56 / night 55 flies; column
                observed_duration_seconds, rows of version "observed" (replicate -1) and "fitted" /
                "phi0" (replicates 0-199), ordered by block_index then bout_index_within_block.
Statistic       Per fly, phase and tail: thresholds from the fly's own bouts in that sequence (upper:
                duration >= the ceil(0.9 n)-th smallest; lower: <= the ceil(0.1 n)-th smallest).
                Pairs = consecutive bouts within one block (bout indices differing by 1).
                C = P(next bout extreme | current bout extreme) / (fraction of bouts extreme);
                C = 1 means no clustering. Identical computation on every replicate (its own
                thresholds). NaN if no extreme bout starts a pair.
Decision rule   PRIMARY (4 tests: day/night x upper/lower, Holm): group mean C over flies (flies with
                finite C) against the group means on each fitted replicate index (200 values, each a
                full 56- or 55-fly group with matched bout counts); two-sided Monte Carlo p.
                Rejection with observed ABOVE the replicates = extremes cluster beyond the Gaussian
                model (tail dependence or a slow state); BELOW = extremes avoid each other.
                SECONDARY: per fly, two-sided Monte Carlo p against its 200 fitted replicates, count
                p < 0.05 against Binomial(n, 0.05) (low power: ~11 extreme bouts per fly);
                descriptive: the same group means under phi0 replicates and under 1000 within-fly
                permutations of the observed durations (block structure kept); the same statistics
                at the 20%/80% tails, reported only, never substituted for the primary.
Matched null    Each fly's own fitted-model simulations: same marginal, same block lengths, same n.
Narrowest block One day or night phase segment of one fly (pairs never cross blocks).
Notes           The group-level test is primary because per-fly counts of extreme pairs are ~1.
                The 10% tail is fixed here and will not be changed after the run. A C above 1 under
                the fitted null is expected from phi alone and is already in the null.
```

### 2026-11-09  I. Is the clustering of long sleep bouts (entry H) time of day, a slow state, or memory?

```
Question        Entry H: in Shreenath's CS LD export, long sleep bouts follow long sleep bouts far
                more often than his stationary Gaussian latent AR(1) allows; short bouts do not. A
                drifting mean bout length predicts exactly that asymmetry. Does the clustering
                survive a null that keeps each fly's time-of-day profile of bout lengths? And does
                it appear in per0 DD, which has no clock and no light cycle?
Data            Our own bouts, as staged and scored for the homeostat entries (preliminaries entry
                2026-11-06; flies saved in homeo_results.mat): CS LD 40 flies from ZT0, per0 DD 48
                flies from sample 1. 1-s scoring, sleep >= 300 s. Sleep bouts only, in onset order.
Blocks          CS LD: phase segments, day = onset ZT 0-12, night = ZT 12-24 of the same 24-h day
                (as in H). per0 DD: 12-h segments from the record start, all pooled as one phase.
                Pairs = two sleep bouts consecutive among the fly's sleep bouts within one segment
                (narrowest block: one 12-h segment of one fly).
Statistic       As entry H: per fly and phase, thresholds from that fly's sleep bouts in that phase
                (upper: >= the ceil(0.9 n)-th smallest; lower: <= the ceil(0.1 n)-th smallest);
                C = P(next extreme | current extreme) / fraction extreme; group mean over flies.
Nulls           Permutations of the sleep-bout durations within fly and phase, positions fixed,
                1000 each (Seed 1, substream = fly index):
                N1 ZT-STRATIFIED (time of day kept): only among bouts whose onsets fall in the same
                   1-h ZT bin, pooled across days (each bin draws on bouts from every day of the
                   record, so the permutation is not degenerate; it breaks day-to-day differences,
                   slow states within a day and bout-to-bout memory, and keeps the ZT profile).
                N2 WITHIN-SEGMENT (each day's level kept): only within the same 12-h segment; keeps
                   day-to-day differences, breaks the profile within a phase.
                N3 PLAIN: within fly and phase; reproduces the H comparison in our data.
                For per0 DD, N1 uses 1-h bins of time since record start modulo 24 h.
Decision rule   PRIMARY (3 tests, Holm), upper 10% tail, two-sided Monte Carlo p of the group mean
                against the null group means:
                  1 CS LD day vs N1;  2 CS LD night vs N1;  3 per0 DD vs N3.
                Reading: CS LD clustering that VANISHES under N1 (and is present under N3) = time of
                day; clustering that SURVIVES N1 = a slow state or memory, not the ZT profile. Under
                N2 (descriptive) survival = within-day slow state or memory; vanishing = differences
                between days. per0 DD clustering beyond N3 = clustering without a clock.
                SECONDARY, descriptive: lower 10% tails; every tail x null combination; N1 at 2-h
                and 3-h bins (pre-specified here, reported alongside, never substituted for 1 h);
                the fraction of consecutive pairs whose two bouts share a 1-h ZT bin on the same day.
Calibration     Permutation within strata is exact under exchangeability within strata, so no
                simulation calibration is needed for size; small-sample bias of C (~0.9 under no
                clustering, entry H) is present equally in data and nulls.
Notes           1-h bins: about 9 bouts per bin per fly per phase pooled over days; wider bins keep
                less of the profile. Non-adjacent bouts (gaps) do not break pairs; their count is
                reported. Our CS LD flies overlap 24 of his 56; scoring differs (ours bridged 1 s).
```

### 2026-11-09  J. CS LD daytime long-bout clustering: slow state within the day, or bout-to-bout memory?

```
Question        Entry I: CS LD daytime clustering of long sleep bouts survives keeping each fly's
                time-of-day profile (pooled across days) and survives keeping each day's level.
                Does it survive keeping BOTH -- each fly's bout lengths at that hour on that day --
                and does it decay with bout lag the way memory would, or stay flat like a state?
Data            As entry I (homeo_results.mat; CS LD 40 flies from ZT0, per0 DD 48), sleep bouts
                only, same phase segments, thresholds and statistic C (top 10% primary).
Nulls           N4 DAY-BY-ZT: permute sleep durations only among bouts on the same day whose onsets
                fall in the same ZT bin (bin width 2 h primary). Keeps the time-of-day profile, each
                day's level and any drift slower than the bin; breaks order within the bin.
                Bin width fixed from occupancy, counted BEFORE this entry (no results seen): day x
                1-h bins leave 16% of daytime bouts alone in their bin and only 48% in bins of >= 3
                (a near-degenerate null that would reproduce any clustering); 2 h: 3% alone, 88% in
                bins >= 3; 3 h: 1%, 95%. 1 h and 3 h are reported alongside, never substituted.
                N1 as entry I (ZT 1-h bins pooled across days) for the lag profile.
Lag profile     C_k for pairs k sleep bouts apart within a segment, k = 1..4 (lag 4 ~ 2 h by day).
                Excess_k = observed C_k - mean N1 C_k. Contrast Dlag = Excess_1 - mean(Excess_2..4),
                its null from the N1 permutations scored the same way.
Decision rule   PRIMARY (2 tests, Holm), CS LD day, upper 10%, two-sided Monte Carlo p, 1000
                permutations (Seed 2, substream = fly):
                  1 group mean C vs N4 (2-h).  SURVIVES -> dependence on a sub-2-h scale (memory, or
                    a state faster than ~2 h); VANISHES -> a slow state within the day (hours).
                  2 group mean Dlag vs N1.  POSITIVE -> the excess is concentrated at lag 1 and
                    decays (memory-like); NOT -> flat over lags 1-4 (state-like).
                SECONDARY, descriptive: the same for CS LD night and per0 DD (per0: same 12-h
                segment and same 2-h bin of time since start); lower tails; occupancy per set.
Notes           Memory and a fast state are not separable by any shuffle; a positive Dlag with
                survival under N4 is the closest this design comes to memory. The two primary
                readings may disagree; both are reported.
```

### 2026-11-09  K. Is the clock-free slow variation of sleep-bout length a trend or a fluctuating state?

```
Question        Entries I and J: in per0 DD, long sleep bouts cluster because the typical bout
                length differs between 12-h stretches of the record. Is that a TREND (bouts
                steadily lengthen or shorten over the record, e.g. habituation, food, age; or a
                first-day effect common to flies) or a FLUCTUATING STATE (stretches of consolidated
                and fragmented sleep that come and go)?
Data            As entries I and J (homeo_results.mat): per0 DD 48 flies (2.97-7.89 days, median
                5.91; median 193 sleep bouts), sleep bouts only, y = log(duration in s), t = onset
                in days from the record start, 12-h segments from the record start. CS LD 40 flies
                descriptive only (day-phase and night-phase series separately; segments = phase
                segments).
Statistics      TREND: per fly, least-squares slope b of y on t. Group: mean b; null = random sign
                flips of the per-fly slopes (10000; flies independent), two-sided. Effect reported
                as exp(mean b) - 1 per day.
                FLUCTUATION: per fly, residuals of y after (i) the fly's own linear trend in t and
                (ii) the common time course (group mean of the per-fly segment means, by segment
                index); V = variance over segments (>= 3 bouts) of the segment means of those
                residuals. Null: durations permuted within fly (positions fixed), then the SAME
                detrending (own trend refitted, common course recomputed from the permuted group),
                1000 permutations (Seed 3, substream = fly). Per fly logRatio = log(V / mean null
                V); group mean logRatio against the null group means, one-sided (greater).
Decision rule   PRIMARY (2 tests, per0 DD, Holm): 1 trend; 2 fluctuation.
                Fluctuation SUPPORTED -> a state that comes and goes beyond any per-fly linear drift
                and any time course shared by the flies (incl. a first-day effect). Trend only ->
                the clock-free "state" is a drift. Neither -> the clustering of I/J came from
                something else (reported as such).
Calibration     Run in the same script BEFORE the data, on simulated flies that keep each real
                per0 fly's bout onsets and draw durations 300 s + Exp(fly mean - 300 s): (a)
                stationary; (b) common linear trend +10%/day; (c) fly-specific trends, slopes
                N(0, 0.1)/day; (d) common first-day step x1.5; (e) a fluctuating state, independent
                level per 12-h segment, SD 0.3 in log mean (power). 100 groups x 199 permutations.
                Rejection rate of the fluctuation test above 0.095 in any of (a)-(d) demotes it to
                descriptive (stated in the ledger).
SECONDARY       Descriptive: excess between-segment variance without any detrending vs after
                detrending (share removed by trends); long-bout clustering C (entry H/I) on the
                detrended durations vs their plain within-fly shuffle; dropping the first 24 h;
                4-h segments; per0 split by record length (< 4 d vs >= 5 d); the same statistics
                for CS LD day and night series; lag-1 correlation of successive segment levels.
Narrowest block One 12-h segment of one fly (fluctuation); the whole record of one fly (trend).
```

### 2026-11-09  L. per0 DD after the first day: does long-bout clustering survive, and on what timescale?

```
Question        Entry K: per0 DD half-day variation in sleep-bout length is mostly the first day
                (acclimation) plus a ~2%/day decline. After the first 24 h, do long sleep bouts
                still cluster -- and if so, is it slow drift between 12-h stretches, a clock-free
                state lasting hours (within 12-h stretches but slower than ~2 h), or something
                faster than 2 h (memory or a fast state)?
Data            per0 DD 48 flies (homeo_results.mat, as I-K), sleep bouts with onset >= 24 h from
                the record start (cutoff fixed here; 48 h reported alongside, never substituted).
                Segments = 12-h stretches from the record start; pairs = consecutive sleep bouts
                within one segment (narrowest block: one 12-h stretch of one fly).
Statistic       C, upper 10% (entry H/I): thresholds from the fly's retained bouts. Lower 10%
                descriptive. Lag profile C_k, k = 1..4.
Nulls           1000 within-fly permutations of the durations, positions fixed (Seed 4, substream =
                fly): P plain (all retained bouts); S12 within the same 12-h stretch; W2 within
                the same 2-h window (floor(onset / 2 h)).
Decision rule   PRIMARY (3 tests, Holm), two-sided Monte Carlo p of the group mean C:
                  1 vs P   -- any clustering after day 1 (not robust to trends; descriptive of
                             the data, not of mechanism);
                  2 vs S12 -- clustering beyond differences between 12-h stretches;
                  3 vs W2  -- clustering beyond drift slower than ~2 h.
                Reading: 1 not rejected -> per0's clustering was acclimation; the clock-free
                claim is dropped. 1 and 2 rejected, 3 not -> a clock-free state lasting hours.
                3 rejected -> dependence faster than ~2 h (memory or a fast state; the lag profile
                then reads which). 1 rejected, 2 not -> slow drift between stretches only.
                SECONDARY: lag profile excess over S12 and Dlag = excess_1 - mean(excess_2..4)
                (as J); lower tail; the 48-h cutoff; per-fly counts.
Calibration     In the same script BEFORE the data, simulated flies on each real fly's retained
                onsets, durations 300 s + Exp(m(t)), m from the fly's mean: (a) stationary; (b)
                common trend -2%/day; (c) acclimation fading past 24 h, m x (1 + 0.5 exp(-t/0.75 d));
                (d) fly-specific trends N(0, 0.1)/day; (e) POWER: a state, independent level per
                3-h block, SD 0.3 in log mean. 100 groups x 199 permutations. P must hold size in
                (a); S12 and W2 in (a)-(d); rejection above 0.095 where size is required demotes
                that test to descriptive. Power in (e) reported.
```

### 2026-11-09  M. Sleep-wake and wake-sleep pairs: what drives the coupling, and is there a before/after asymmetry after long wake?

```
Question        (Q1) Adjacent sleep and wake bouts are positively coupled (homeostat entry A: per0
                DD r_pre +0.104, r_post +0.115). Is that coupling the same slow layer that makes long
                sleep bouts cluster (entries I-L), or a bout-to-bout dependence that survives the
                shuffle ladder? (Q2) A homeostat should act mainly after LONG wake bouts: is a long
                wake bout followed by a long sleep bout more often than it is preceded by one? (A
                before/after asymmetry is also a time-reversal asymmetry of the two-state sequence.)
Data            homeo_results.mat (as I-L): per0 DD 48 flies, bouts with onset >= 24 h from the
                record start; CS LD 40 flies, whole records from ZT0. Sleep >= 300 s, wake >= 2 s,
                1-s scoring, 1-s bridging. Blocks: per0 12-h stretches from the record start; CS LD
                phase segments (day ZT 0-12, night ZT 12-24 of one day). A PAIR = two bouts adjacent
                in time (onset = previous onset + duration), both in the same block. pre = (wake,
                following sleep); post = (sleep, following wake). Narrowest block: one 12-h segment
                of one fly.
Statistics      Q1: per fly, Pearson correlation of within-fly ranks (Spearman) for pre and post
                pairs; coupling R = (r_pre + r_post) / 2; group mean R.
                Q2: per fly, long = top 10% of the fly's retained bouts of that type (>= the
                ceil(0.9 n)-th smallest). Cpre = P(sleep after a long wake is long) / P(sleep long);
                Cpost = P(sleep before a long wake is long) / P(sleep long); Dtail = Cpre - Cpost
                (positive = homeostat direction). Group mean Dtail.
Nulls           Q1: sleep durations permuted among sleep bouts and wake durations among wake bouts,
                independently, within fly, positions fixed, 1000 times (Seed 5, substream = fly):
                P plain; S12 same block; W2 same 2-h window (floor(onset / 2 h)); for CS LD also
                Z same 1-h ZT bin pooled over days (descriptive).
                Q2: random sign flips of per-fly Dtail (10000): under time reversal each fly's
                Dtail is symmetric about 0.
Decision rule   PRIMARY (5 tests, Holm), two-sided:
                  1 per0 R vs P;  2 per0 R vs W2;  3 per0 Dtail (sign flip);
                  4 CS LD R vs P; 5 CS LD R vs W2.
                Reading Q1: rejected vs P but not W2 -> the coupling is the slow layer (hours-scale
                drift shared by sleep and wake; positive = a consolidation "tempo" state, negative
                = an arousal trade-off); rejected vs W2 -> a bout-to-bout coupling, the first
                bout-level dependence found. Q2: Dtail > 0 -> homeostat-direction asymmetry after
                long wake; Dtail ~ 0 -> no bout-level homeostat signature even in the tail.
                CS LD Dtail is DESCRIPTIVE only: a time-of-day profile makes the sequence locally
                non-stationary and can produce Dtail without any homeostat (cf. the CS LD time-arrow
                entry), and no profile-respecting null is calibrated for it here.
                SECONDARY: r_pre, r_post and Delta = r_pre - r_post (sign flip); S12 and Z rungs;
                the reverse tail (long sleep -> next vs previous wake); per-fly counts.
Calibration     In the same script BEFORE the data: simulated flies keep each real fly's bout
                positions, types and blocks; durations resampled from the fly's own sleep / wake
                durations: (a) independent; (b) shared slow state, one factor exp(N(0, 0.3^2)) per
                2-h window multiplying sleep and wake; (c) sleep durations declining 2%/day; (d)
                POWER, homeostat: a sleep bout x (preceding wake / median wake)^0.15. 100 groups;
                199 permutations, 2000 sign flips. Size required: R vs P in (a); R vs W2 in (a)-(c);
                Dtail in (a)-(c). Rejection above 0.095 where required demotes that test to
                descriptive.
```

### 2026-11-10  Amendment to K, L, M: calibration rerun with an independent permutation seed per simulated group

```
Why             In the calibration sections of entries K, L and M, every simulated group passed the
                SAME permutation / sign-flip seed to the test (K: 7; L: 9; M: 11), so the 100 groups
                reused identical shuffles and were not independent. Rejection rates were therefore
                not valid size or power estimates (M showed 0.00 almost everywhere). The data
                sections are unaffected (one dataset, its own permutations) and are NOT rerun.
Change          Only the test seed: K seed 40000 + 100 k + g; L 30000 + 100 k + g; M 20000 + 100 sc
                + g (same seed for the M sign flips). Everything else exactly as pre-specified:
                same scenarios, 100 groups, 199 permutations, 2000 sign flips, simulation seeds,
                required-size sets and the 0.095 demotion threshold.
Decision rule   Apply each entry's original rule to the new rates. A test whose required size now
                exceeds 0.095 is demoted to descriptive and its verdict is re-read; a test that was
                demoted (K fluctuation) is restored only if all its required sizes hold. Power
                statements (L scenario e; M scenarios b, d) are replaced by the new rates. Results
                go to the ledger as correction entries; earlier entries are not edited.
```

### 2026-11-10  Second amendment to M: calibration counter bug

```
Why             The M calibration line `rej(sc,:) = rej(sc,:) + [p...] < 0.05` lacks parentheses;
                MATLAB applies + before <, so the count never accumulates and every reported rate
                is 0.00 or 0.01. The original M calibration AND the fresh-seed rerun begun under
                the first amendment are both invalid. K and L count correctly.
Change          `+ ([p...] < 0.05)`, plus the independent seeds of the first amendment. per0 DD and
                CS LD calibrated in two parallel processes; nothing else changes. Same decision
                rule as the first amendment.
```

### 2026-11-10  N. tim0 DD: does the first 24 h carry an LD / acclimation effect (as in per0 DD)?

```
Question        per0 DD half-day differences in sleep-bout length were mostly the first day (entry
                K), so per0 analyses start at 24 h. Does tim0 DD's first day differ from the later
                days in the same way? If yes, all tim0 analyses start at 24 h.
Data            tim0 DD from Zichao (zip of 25 per-fly 1-s activity files, 96 h each, "_from1stLD"
                cut: time 0 = start of the first subjective day). Pipeline exactly as per0 DD
                (run_homeo.m): process_LocomToSleep timeInterval 1, segmentMode 1, startZT 0,
                Teliminate 0, deadThresh 12, minDays 3, bridgeWakeSeconds 1; sleep >= 300 s;
                flies kept with >= 40 sleep and >= 40 wake bouts. Reference: per0 DD 48 flies
                (homeo_results.mat), same statistics on its days 1-4 only, so the windows match.
Statistics      Per fly, per 24-h day d (1-4 from record start), per half H1 = 0-12 h, H2 = 12-24 h:
                L = mean log sleep-bout duration (bouts by onset); A = sleep minutes.
                Level effect   dL = L(day 1) - mean L(days 2-4); dA likewise.
                LD carry-over  contrast c = value(H2) - value(H1) within a day;
                               dcL = cL(day 1) - mean cL(days 2-4); dcA likewise.
Nulls           Random sign flips across flies (10000, Seed 1), two-sided, of the per-fly values.
Decision rule   PRIMARY tim0: the four tests dL, dA, dcL, dcA, Holm. DROP the first 24 h if any is
                Holm < 0.05. If none is but a day-1 effect matches per0's in sign and is >= half
                its size, report it as underpowered and say so (25 flies) -- the cut is then the
                user's call. per0 values are a reference only (descriptive).
                Also shown: the day-by-day profiles (L, A, both halves) per genotype.
```
