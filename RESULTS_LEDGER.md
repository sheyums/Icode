# Results ledger

**Append-only.** One entry per result, **including null and inconclusive
results**. Do not delete or rewrite an entry; add a correcting entry that
references it.

The point of recording nulls is that selection cannot be audited from the
successes alone. A ledger with only positive findings gives no way to know
whether a p-value of 0.03 was the first analysis run or the twentieth.

Rule 3 in `CLAUDE.md` governs this file. Each entry names the `PRESPEC.md`
entry it answers.

## Template

```
### YYYY-MM-DD  <short name>

Commit          Hash of the code that produced it (git log -1 --format=%H).
Prespec         Which PRESPEC.md entry this answers.
Time window     The narrowest block the number was computed from.
n               Sample size actually used, after truncation at xmin.
Result          The number, with its uncertainty. State the units.
Matched null    The same statistic under a sample-size-matched null.
Verdict         Support / against / inconclusive, per the pre-specified rule.
Not reportable  Anything that must NOT be quoted from this run, and why.
```

## Entries

<!-- Append below. Never edit an entry above. -->

### (none yet)

Earlier per0 DD results are in `FINDINGS.md`, written before this ledger
existed. They are not retrofitted here: back-filling a ledger from memory
would give it the appearance of a contemporaneous record without the
substance, which is worse than an honestly empty file.

### 2026-09-15  per0 DD wake: does any single fly carry the hazard rise?

```
Commit          35ee73a358c741f6a1b7209f8b7f35a58585d7d1 (Icode main; code unchanged since 8565082)
Script          wake_loo_rise.m, session scratch, not in the repo
Prespec         2026-09-15 "per0 DD wake: does any single fly carry the hazard rise?"
Time window     The whole analysed record of each fly (2.8-4.0 days, constant darkness). Every
                number was computed on 24 flies pooled, or on 3989 - n_i bouts after random removal.
n               3989 wake bouts from 25 flies; xmin 2 s, dt 1 s, 1-s bridged scoring. 3888 bouts
                after the one influential fly was removed.
Result          Full data RiseRatio 1.405, 1.388, 1.376, 1.430 at 12, 16, 20 and 24 bins.
                One fly, row 14 (n = 101), fell below the 2.5th percentile of its own random
                removals at all four counts; no other fly fell below at more than two of them.
                Removing fly 14 left RiseRatio 1.351, 1.326, 1.323, 1.370 against the null medians
                of that same run, 1.090, 1.139, 1.188, 1.240, so above the median at every count.
                RiseNullP there was 0.003, 0.025, 0.137 and 0.187, with NullDropped 0 throughout.
                Removing fly 10 (n = 239) instead RAISED the ratio to 1.684-1.936, so that fly
                damps the rise rather than causing it.
Matched null    For influence: 199 removals of the same number of bouts as the fly, drawn uniformly
                from all 3989 (RandomSeed 1000 + fly), on the same bins. For the p-values after
                removal: 999 draws from the two-exponential fit to the remaining bouts.
Verdict         Per the pre-specified rule: some flies matter, but the rise does not depend on them
                alone, and no single fly carries it. For the design question, a pooled wake hazard
                in the fuller datasets can be read as such, with per-fly influence checked rather
                than assumed.
Not reportable  The p-values after removing fly 14 are not a test of the rise in the full data, and
                they straddle 0.05 across bin counts (0.003 to 0.187), so no single one should be
                quoted. Nothing here says whether the rise lives in the flies or is produced by
                scoring bursty movement into wake bouts.
```

### 2026-09-15  per0 DD sleep: pooled K=3 against a per-fly null

```
Commit          35ee73a358c741f6a1b7209f8b7f35a58585d7d1 (Icode main; code unchanged since 8565082)
Script          sleep_pooled_null.m, session scratch, not in the repo
Prespec         2026-09-15 "per0 DD sleep: pooled K=3 against a per-fly null"
Time window     Per-fly fits: the whole analysed record of each fly (2.8-4.0 days, constant
                darkness). The tested statistic: all 25 flies pooled.
n               3987 sleep bouts from 25 flies; xmin 300 s, dt 1 s, 1-s bridged scoring. Per fly
                80 to 249 bouts, median 159.
Result          Observed logL -31626.8558 at K = 2 and -31618.4952 at K = 3, so LR = 16.72,
                reproducing FINDINGS. The observed K = 3 fit is the FINDINGS one: tau 54.6, 1040.5
                and 3141.2 s with observed weights 0.0496, 0.9298 and 0.0206. Fitted one fly at a
                time, only 10 of the 25 flies supported two timescales; the other 15 came back
                degenerate at K = 2 and entered the null with one. Pooling 25 such simulated flies
                gave LR well ABOVE the observed value: median 26.68, 95th percentile 48.14, maximum
                58.94 over B = 199, with no negative LR and no failed fit. p = 0.86. The slowest of
                the three fitted components held under 5 percent of bouts in 11 percent of the
                replicates.
Matched null    199 pooled datasets, each with the same 25 flies at the same bout counts, every fly
                drawn from its own fitted model. The staged rule stopped at 199 because p > 0.12.
Verdict         The pre-specified reading is supported: differences between flies are enough to
                produce the pooled third component, and more than enough. A pooled mixture fit on
                these data says nothing about how many timescales an individual fly has, so the
                fuller datasets are to be analysed fly by fly or with a hierarchical model.
Not reportable  Not evidence that no fly has a third timescale: the test asks what pooling can
                manufacture, not what individuals do. The null is if anything too generous, since it
                simulates from unshrunk per-fly fits whose estimation noise widens the spread
                between flies by an unknown amount; a hierarchical, shrunk null would be the sharper
                comparison. The bouts are bridged, so merging sits inside both the data and the
                per-fly fits.
```

### 2026-09-15  Correction to "per0 DD sleep: pooled K=3 against a per-fly null"

```
Commit          35ee73a358c741f6a1b7209f8b7f35a58585d7d1 (no new computation; wording only)
Corrects        2026-09-15 "per0 DD sleep: pooled K=3 against a per-fly null", recorded in 84a5187
What was wrong  That entry said a pooled mixture fit "says nothing about how many timescales an
                individual fly has". Too strong, as Sheyum pointed out.
Correct reading A pooled fit does constrain the population: K > 1 means the flies are not all one
                exponential at a common rate. What it cannot do is attribute the extra timescales,
                because differences between flies and structure within a fly produce the same pooled
                distribution. On these data the third component needs nothing beyond differences
                between flies (p = 0.86), so that component cannot be read as a third timescale in
                any fly. One direction does carry: mixing cannot make a hazard rise, so a
                non-monotone pooled hazard does imply non-monotone hazards in some flies. That is
                what the wake entry rests on.
Separate        The per-fly evidence is a different analysis from the pooled fit. Fitted one at a
evidence        time, 10 of the 25 flies supported two timescales. That count is the fitter
                degeneracy flag at about 160 bouts a fly, not a test, and at that bout count an
                unresolvable second timescale and a genuine rare component look alike, as the 219-fly
                LD work measured.
Unchanged       The numbers and the verdict on the pre-specified question stand: p = 0.86, and the
                unit of analysis for the fuller datasets is the fly or a hierarchical model.
```

### 2026-09-15  per0 DD wake: the hazard rise at four fixed bin counts

```
Commit          172315b (Icode main; analysis code unchanged since 8565082)
Script          wake_rise_fixed_bins.m, session scratch, not in the repo
Prespec         2026-09-15 "per0 DD wake: the hazard rise at four fixed bin counts"
Time window     The whole analysed record of each of the 25 flies (2.8-4.0 days, constant
                darkness), pooled. That pooled set is the narrowest block here.
n               3989 wake bouts; xmin 2 s, dt 1 s, 1-s bridged scoring.
Result          Null K = 2 fit: tau 8.4 and 661.9 s, observed weights 0.079 and 0.921,
                logL -29328.2460.
                bins 12: trough 119.8 s, peak 1094.0 s, RiseRatio 1.405, null median 1.083,
                         95th 1.210, RiseNullP 0.001
                bins 16: trough 109.3 s, peak  997.7 s, RiseRatio 1.388, null median 1.136,
                         95th 1.285, RiseNullP 0.010
                bins 20: trough 161.0 s, peak  944.3 s, RiseRatio 1.376, null median 1.186,
                         95th 1.390, RiseNullP 0.068
                bins 24: trough 144.1 s, peak  910.0 s, RiseRatio 1.430, null median 1.235,
                         95th 1.512, RiseNullP 0.101
                NullDropped 0 at every count; one bin dropped and one saturated at each, as
                expected. The values at 16 and 20 bins reproduce FINDINGS exactly.
                The observed ratio barely moves across bin counts (1.376 to 1.430). What moves is
                the NULL: its median climbs from 1.083 to 1.235 and its 95th percentile from 1.210
                to 1.512, because more bins put fewer bouts in each and a spurious peak-over-trough
                ratio grows.
Matched null    999 datasets of n = 3989 drawn from that two-exponential fit, each put through the
                whole estimator on its own bin edges.
Verdict         Supported but bin-sensitive, the wording fixed on 2026-09-14: a rise at 12 and 16
                bins, none at 20 and 24. Every count reported, not the smallest p.
Not reportable  Do not quote 0.001 on its own. The bin sensitivity here is a loss of power rather
                than a wobble in the data, since the observed ratio is stable while the null widens;
                that reading is an observation from these four counts, not a separate test. The test
                is against a falling-hazard null only: it says nothing about whether the rise is in
                the flies or is produced by scoring bursty movement into wake bouts.
```

### 2026-11-05  per0 DD: time arrow in successive sleep and wake bouts

```
Commit          f46b912 (Icode branch claude/shifted-lognormal-mle-wc03ob: BoutIrreversibility)
                with the analysis folder's scoring code at its commit 23240e7
                (process_LocomToSleep, extractBoutsForModelCompetition); run from the Icode
                folder behind the which/pwd assert, MATLAB R2026a, 2026-11-05.
Script          run_per0dd.m, session scratch, not in the repo
Prespec         2026-11-05 "per0 DD: is there a time arrow in the sleep and wake bout
                sequences?" (committed 2df6bef before the run)
Time window     Each fly's whole retained record, one block per fly (3-7.9 days, constant
                darkness, read from the first sample); group numbers pool the 48 blocks.
n               48 of 48 flies passed the pipeline and the >= 40-bout rule for both types.
                Sleep: 83-426 bouts per fly (median 193), 9790 lag-1 pairs; wake: 84-425
                (median 193), 9775 pairs. Ties (zero differences) 0.1%. Bouts verified in time
                order (OnsetBins strictly increasing in every row).
Result          Group means at lag 1 (IAAFT p; shuffle p for memory):
                  sleep  UpFraction 0.5070 (p 0.020)  Bowley -0.0075 (p 0.55)
                         Moment -0.073 (descriptive)  Rho +0.028 (shuffle p 0.002)
                  wake   UpFraction 0.4903 (p 0.004)  Bowley +0.0178 (p 0.14)
                         Moment +0.052 (descriptive)  Rho -0.060 (shuffle p 0.002)
                Lags 2, 4, 8, UpFraction: sleep 0.5050 / 0.4998 / 0.4977 (p 0.12 / 0.95 / 0.48);
                wake 0.4996 / 0.4985 / 0.5007 (p 0.92 / 0.62 / 0.82). Bowley p >= 0.45 at every
                lag in both. Memory: sleep Rho 0.024 / 0.007 / 0.018 (p 0.008 / 0.26 / 0.036),
                wake 0.021 / 0.004 / 0.010 (p 0.014 / 0.38 / 0.19).
                Flies individually p < 0.05 at lag 1 (IAAFT): UpFraction sleep 4 of 48 (3 above
                0.5), wake 2 (both below); Bowley 1 and 1. Memory (shuffle) sleep 9 (7 positive),
                wake 8 (all 8 NEGATIVE). Per-fly UpFraction at lag 1 spans 0.460-0.555 (sleep),
                0.449-0.524 (wake).
Matched null    1000 IAAFT surrogates per fly (same bouts, about the same autocorrelation,
                reversible); group null = mean of the flies' r-th surrogates. Group null 95%
                band for UpFraction at lag 1 is about 0.494-0.506 in both types.
Verdict         PRIMARY, Holm across four: wake UpFraction p 0.004 -> 0.016, SUPPORTED: a time
                arrow in successive wake bouts at lag 1, direction = slightly more decreases
                than increases (49.0% increases). Sleep UpFraction 0.020 -> 0.060, sleep Bowley
                0.55, wake Bowley 0.28: not supported. So: a small but detectable arrow in wake,
                none detected in sleep at the family-wise level.
Not reportable  Do not quote sleep UpFraction p = 0.020 as significant; it fails the Holm
                correction. Do not quote the shuffle p-values for the skew statistics (they
                over-reject with memory, as measured). The effect is ONE percentage point of
                increases, carried by the group, not by individual flies (2 of 48 individually),
                so it is not a per-fly trait at this n. It concerns the ORDER of successive
                bouts; it says nothing about the GH skew of the collaborator's bout-space
                residuals, whose null is different (simulation from the fitted model). IAAFT
                calibration was measured only for positive memory (latent AR phi 0.5); wake
                memory here is NEGATIVE at lag 1 (-0.06), which was not simulated. The six
                12_29_2017 plain files were used whole (see the PRESPEC's open point).
```

### 2026-11-05  CS LD: time arrow beyond the day/night profile

```
Commit          4426e08 (Icode branch claude/shifted-lognormal-mle-wc03ob: BoutIrreversibility with
                the strata null) with the analysis folder's scoring code (process_LocomToSleep,
                extractBoutsForModelCompetition; unchanged since its commit 23240e7); run from the
                Icode folder behind the which/pwd assert, MATLAB R2026a, 2026-11-05.
Script          run_csld.m, session scratch, not in the repo
Prespec         2026-11-05 "CS LD: is there a time arrow in successive bouts beyond the day/night
                profile?" (committed 58ac4d7 before the run)
Time window     Each fly's whole record (4-5 days, 12:12 LD, from ZT0; DataStartZT = 0 checked),
                one block per fly; strata = 1-h ZT bins of onset. The count pools the 40 flies.
n               40 of 40 flies passed the pipeline and the >= 40-bout rule for both types. Sleep:
                98-234 bouts per fly (median 178); wake 97-233 (median 177). Bouts in time order
                (checked). 95-100% of bouts per fly sit in a stratum with at least one other bout.
Result          PRIMARY, flies with per-fly UpFraction p < 0.05 vs their own strata null, lag 1:
                  sleep 1 of 40 (count p 0.87), wake 0 of 40 (count p 1.00). 2 expected by chance.
                Lags 2, 4, 8: sleep 3, 2, 2; wake 1, 2, 2 (count p >= 0.32 everywhere).
                Median per-fly UpFraction minus its own strata-null mean, lag 1: sleep +0.0072,
                wake +0.0021.
                Clock-blind nulls, group mean (NOT for a verdict; shown to be fooled by the
                profile): sleep UpFraction 0.5107 at lag 1, IAAFT p 0.004, shuffle p 0.002; wake
                0.4996, IAAFT p 0.88.
                Memory: group Rho at lag 1 sleep 0.152, wake 0.094 -- larger than per0 DD (0.028,
                -0.060) -- but only 2 (sleep) and 1 (wake) flies keep significant memory against
                the strata null; median Rho minus strata-null mean +0.02 (sleep), +0.04 (wake).
                Bowley vs strata, flies significant at lag 1: sleep 3, wake 3 (descriptive; its
                per-fly null rate was measured at about 0.068).
Matched null    Per fly, 1000 surrogates permuting its own bouts within 1-h ZT bins across its
                days (same bouts, same time-of-day profile); count against Binomial(40, 0.05).
Verdict         NOT SUPPORTED for either bout type (Holm 1.00 and 1.00): no time arrow beyond the
                day/night profile detected by the per-fly count.
Not reportable  Do not quote the IAAFT/shuffle sleep p = 0.004 / 0.002 as a time arrow: on
                LD-like simulated flies the profile alone made those nulls reject 0.07-0.135 per
                fly, and their group mean is exactly where that bias accumulates. Do not read this
                null as contradicting per0 DD wake: the count is far weaker (per-fly power 0.535
                against a strong simulated arrow; an effect of per0's size, 2/48 flies, is
                invisible to it), so the two results are not comparable. "Memory is mostly the
                profile" is descriptive (no pre-specified test). Teliminate = 0: the first,
                possibly less entrained, day is included.
```

### 2026-11-06  per0 DD and CS LD: the bout-space AR(1) x(t+1) = A + B*x(t) + eps and its noise

```
Commit          328116c (Icode branch claude/shifted-lognormal-mle-wc03ob: BoutAR1Noise,
                FitResidualFamilies; the six noise fitters at their audited versions) with the
                analysis folder's scoring code (unchanged since its commit 23240e7); run from the
                Icode folder behind which/pwd asserts for all seven functions, MATLAB R2026a,
                2026-11-05/06.
Script          run_ar1.m, session scratch, not in the repo
Prespec         2026-11-05 "per0 DD and CS LD: can x(t+1) = A + B*x(t) + eps generate the bouts,
                and what is eps?" (committed 14adef6 before the run)
Time window     per0 DD: each fly's whole record (one block). CS LD: each fly's light phase and
                dark phase fitted separately, every 12-h phase of every day its own block.
n               per0: 48 flies, median 192 pairs per fly per bout type. CS: 40 flies; median
                pairs per fly sleep light 94, sleep dark 70, wake light 95.5, wake dark 68.
                Sleep xmin 300 s, wake 2 s; DataStartZT 0 and time order checked.
Result          PRIMARY (Hetero = Spearman of |eps| with x(t); flies with p < 0.05 vs null):
                  per0 sleep 8/48 (7 with spread GROWING with x(t)), count p 0.0024, Holm 0.0049
                  per0 wake  6/48 (4 shrinking, 2 growing),           count p 0.0317, Holm 0.0317
                  CS sleep light 2/40, sleep dark 3/40, wake light 3/40, wake dark 2/40;
                  count p 0.32-0.60, Holm 1.00 for all four.
                eps FAMILY, BIC winners per set (pre-set wording rule):
                  per0 sleep: shifted Weibull 18, shifted exp 13, asym Laplace 7, shifted
                    lognormal 7, GH 3 -> named "shifted Weibull (18/48)"; shifted-exp median
                    DeltaBIC 13.8.
                  per0 wake: shifted exp 18, shifted Weibull 11, asym Laplace 7, shifted
                    lognormal 7, GH 5 -> most frequent shifted exp (18/48, not a majority);
                    median DeltaBIC 4.9.
                  CS sleep light: asym Laplace 18, shifted exp 14, shifted Weibull 4, GH 3,
                    shifted lognormal 1 -> "asymmetric Laplace (18/40)"; DeltaBIC 7.2.
                  CS sleep dark: shifted exp 21/40 (majority) -> "the bout distribution
                    shifted by a constant"; DeltaBIC median 0.0.
                  CS wake light: GH family 21, shifted lognormal 11, shifted exp 5, asym
                    Laplace 3 -> "GH family (21/40)"; shifted-exp DeltaBIC 43.0.
                  CS wake dark: shifted exp 17, shifted lognormal 11, GH 6, asym Laplace 4 ->
                    most frequent shifted exp (17/40); DeltaBIC 2.2.
                SECONDARY: B median (IQR): per0 sleep 0.037 (-0.015, 0.113), group mean 0.041,
                  p 0.002, 11/48 flies (9 positive); per0 wake -0.035 (-0.106, 0.016), mean
                  -0.042, p 0.002, 9/48 (7 negative); CS sleep light 0.080, 8/40 (count p
                  0.0007); sleep dark 0.114, 4/40 (0.14); wake light 0.122, 8/40 (0.0007); wake
                  dark -0.041, 4/40 (0.14).
                  eps against the no-dynamics eps (same bouts, random order), medians: skewness
                  per0 sleep 1.29 vs 1.31, per0 wake 2.03 vs 2.06, CS 1.46-3.20 vs 1.54-3.24;
                  SD within 1-2% in every set (panel F). ShapeKS significant flies: per0 11 and
                  6, CS 9, 5, 3, 1.
                  Regeneration, median fraction of generated bouts below the floor: per0 sleep
                  0.034, wake 0.050; CS 0.052-0.094 (max 0.27); generated vs observed KS 0.072-
                  0.146 against a sampling reference of 0.058-0.092.
Matched null    1000 surrogates per fly (per phase): within-block shuffles (per0), permutations
                within 1-h ZT bins across the phase's days (CS). Count null Binomial(n, 0.05).
                Family comparison calibrated: shifted exponential wins 100/100 when true.
Verdict         Question (1): the additive-noise assumption is CONTRADICTED in per0 for both bout
                types -- more flies than chance have noise whose spread depends on the previous
                bout (sleep mostly growing, wake mostly shrinking) -- and NOT contradicted in any
                CS LD phase. Question (2), by the pre-set rule: per0 sleep shifted Weibull, per0
                wake shifted exponential (plurality), CS sleep light asymmetric Laplace, CS sleep
                dark the bout distribution shifted by a constant, CS wake light GH family, CS
                wake dark shifted exponential (plurality). In EVERY set, eps has nearly the same
                skewness and width as the residuals of the same bouts in random order.
Not reportable  Do not present the winning family as a property of the dynamics: panels E-F show
                eps has the bouts' own shape, so the family describes each set's bout
                distribution as seen through a near-flat line. Do not read GH winning CS wake
                light as support for the collaborator's GH noise without that caveat. Per-fly
                BIC winners at n = 68-192 are noisy; quote counts, not a single fly's family.
                Hetero significant in 8/48 and 6/48 flies says the additive model fails for SOME
                per0 flies, not which model is right. The CS LD calibration was simulated on
                whole records, not per phase. Regeneration numbers are descriptive.
```

### 2026-11-06  Four homeostat analyses (A cross-type, B Cox check, C window ACF, D Process S proxy)

```
Commit          141b620 (Icode branch claude/shifted-lognormal-mle-wc03ob: BoutCrossType,
                SleepAmountACF, ProcessSProxy, RenewalSurrogate; B used 328116c's BoutAR1Noise
                results of ledger 625dc08) with the analysis folder's scoring code (unchanged since
                23240e7); run from the Icode folder behind which/pwd asserts, MATLAB R2026a.
Script          homeo_run/run_homeo.m (A, C, D), homeo_run/run_B.m (B), stage_check.m (per0 LD
                phase check); session scratch, not in the repo
Prespec         2026-11-06 preliminaries and entries A-D (fe307b0) as AMENDED by the 2026-11-06
                calibration amendment (141b620), both committed before any data run.
Time window     Each fly's whole record, one block per fly (per0 DD from sample 1; CS LD from ZT0;
                CS DD from CT0; per0 LD trimmed to its first ZT0). C and D: per0 DD only.
n               per0 DD 48, CS LD 40, CS DD 25 flies; every record strictly alternating with no
                non-adjacent link. per0 LD: the run-time phase check (largest rise within 1 h of
                each transition in the top 10% of that day's rises, fixed before looking) passed
                every day only for 8_25_2017; 8_18 (lights-on 5/6 days), 9_21 (5/8), 9_22 (5/6),
                9_29 (lights-on 4/6, lights-off 5/6) were DROPPED by the committed rule -> per0 LD
                = 3 flies, uninformative. Lights-off was found on 31 of 32 days; the misses are
                the weak lights-on startle of arrhythmic flies.
Result          A per0 DD (shift null), group means: r_pre +0.104 (p 0.001), r_post +0.115
                (0.002), Delta -0.011 (0.74), r_pre2 -0.007 (0.67), q_pre +0.113 (0.001),
                q_post +0.103 (0.002), DeltaW +0.010 (0.25). Flies p < 0.05: r_pre 24/48,
                r_post 20, Delta 5, q_pre 24.
                A profile sets (strata null), flies with r_pre p < 0.05: CS LD 7/40 (count p
                0.0034), CS DD 14/25 (1.6e-12), per0 LD 0/3. Descriptive: CS DD Delta 5/25,
                r_pre2 mean -0.097 (0/25 positive-significant); CS LD q_pre 12/40.
                B (Spearman across flies of B with Hetero; null = independent exponential flies
                matched in n, xmin, mean bout): per0 sleep +0.714 vs null median +0.754 (95th
                0.853), p 0.69; per0 wake +0.630 vs +0.757, p 0.94. Secondary CS LD: sleep light
                +0.714 (p 0.64), sleep dark +0.724 (0.53), wake light +0.932 (0.001), wake dark
                +0.915 (0.002).
                C per0 DD (DESCRIPTIVE): group deviation 30-min min -0.035 at 1.5 h (p 0.005),
                max +0.045 at 10.5 h; 60-min min -0.069 at 1 h (p 0.001), max +0.072 at 4 h.
                D per0 DD: max group rho_S +0.092 at (Tw 0.5, Ts 0.5) h, p 0.001; min group rho_W
                -0.140 at (0.5, 0.5) h, p 0.001; both decline monotonically with longer time
                constants. Per fly, rho_S at (0.5, 0.5) vs A's r_pre: Spearman +0.825.
Matched null    A: all circular shifts (per fly) / 1000 random shifts (group) of each fly's wake
                sequence; strata: 1000 within-1-h-bin permutations; count Binomial(n, 0.05). B:
                1000 matched simulated sets. C, D: 1000 renewal surrogates per fly (identical in C
                and D).
Verdict         A per0 DD: r_pre SUPPORTED (Holm 0.002): a longer wake bout precedes a longer sleep
                bout. Delta NOT SUPPORTED (Holm 0.74): no before/after asymmetry -- by the
                pre-set and simulation-checked reading rule, a symmetric coupling of adjacent
                bouts (the shared-slow-state pattern), NOT the homeostat signature. A profile
                sets: coupling SUPPORTED in CS LD (Holm 0.0068) and CS DD (4.8e-12), read as
                coupling, not homeostat (amendment); per0 LD uninformative (3 flies).
                B: NOT SUPPORTED for per0 (Holm 1.00, 1.00): the B-Hetero concordance is what
                the estimators produce on independent bouts.
                C: descriptive only (no verdict, per the amendment).
                D: SUPPORTED as amended (p 0.001): a short-memory variable links successive bouts.
                Its optimum at the shortest time constants, and its per-fly agreement with r_pre
                (+0.825), say it is the adjacent-bout coupling of A, not multi-hour pressure.
Not reportable  Do not call the per0 DD coupling homeostatic: Delta is flat, and coupling reaches
                only the immediate neighbour (r_pre2 ~ 0), which a multi-hour pressure would not
                do. Do not quote C's p-values as tests (over-rejecting, non-specific). Do not
                present the earlier "B and spread go together, so a slowly drifting rate" reading
                (session discussion, not a ledger result): B refutes it for per0. per0 LD carries
                no result. The shared-slow-state model was simulated with one structure (a
                per-cycle AR(1)); a symmetric coupling confined to adjacent bouts is a different
                pattern from that model's (whose coupling also reaches further), so "shared slow
                state" names the reading rule's category, not a fitted mechanism.
```

### 2026-11-06  E. Mean reversion of cumulative sleep (variance-time), per0 DD

```
Commit          e327aab (SleepVarianceTime, RenewalSurrogate) on the bouts of the 2026-11-06 homeostat
                run (same scoring, 48 flies, contiguous, strictly alternating); MATLAB R2026a.
Script          homeo_run/run_E.m, session scratch, not in the repo
Prespec         2026-11-06 entry E (committed 06b38b9 before the run)
Time window     Each fly's whole record (3-7.9 days), one block; scales 1-24 h within it.
n               48 flies; 1000 renewal surrogates per fly (the same as entries C and D).
Result          Group mean Slope -0.0089; one-sided p (reversion) 0.903; persistence tail 0.098.
                Per-fly slope median +0.002 (IQR -0.084 to +0.141); 5/48 flies p < 0.05 for
                reversion; 10/48 slopes below -0.10.
                Mean logR at 1, 2, 4, 8, 16, 24 h: -0.122, -0.209, -0.280, -0.265, -0.197, -0.183.
                For comparison (calibration medians): independent bouts -0.002 ... -0.160 (the
                log-variance bias); homeostat -0.10 ... -1.16; shared slow state -0.05 ... -0.24.
                Runs (descriptive): mean Z sleep +0.18 (p 0.18), wake +0.36 (p 0.010).
Matched null    Renewal surrogates of each fly's own bouts; group null from the surrogate index.
Verdict         NOT SUPPORTED: no mean reversion that GROWS with timescale (Slope ~ 0, p 0.90);
                by the pre-set rule nothing here is homeostat-like.
Not reportable  The logR profile (descriptive, not a test): the surplus varies LESS than the renewal
                null already at 1-8 h (-0.12 to -0.28, against the null's own bias of -0.002 to
                -0.05 there), then the gap narrows towards 24 h. That shape -- reduced variability
                at a few hours that does not keep growing -- is what the Slope cannot see, and it
                must not be reported as a test result: it was looked at after the primary failed.
                It may be the short-range reversion of entries C (dip at 1-1.5 h) and of the wake
                runs; a scale-specific test would need its own pre-specified entry. The runs test
                is non-specific (a shared slow state triggered it every time).
```
