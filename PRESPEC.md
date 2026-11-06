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
