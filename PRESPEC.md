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
