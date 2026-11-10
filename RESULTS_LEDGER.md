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

### 2026-11-09  F. Collaborator's latent-AR residuals: model adequacy and the memory effect on skewness

```
Commit          50beb18 (PRESPEC only; analysis is a self-contained MATLAB script on his export, no
                Icode function used); MATLAB R2026a.
Script          F_run/run_F.m, session scratch, not in the repo
Prespec         2026-11-09 entry F (committed 50beb18 before the run)
Time window     His CS LD day and night phases separately; each fly's sleep-bout residuals in a phase.
n               Day 56 flies, 6215 observed residuals, 2,486,000 simulated; night 55 flies, 4623 and
                1,849,200. 200 replicates per fly per version, all successful. phi median 0.036
                (IQR -0.000 to 0.079) by day, 0.032 (-0.018 to 0.070) by night.
Result          (1) Adequacy, flies whose observed moment skewness lies outside their 200 fitted-model
                replicates (p < 0.05): day 0/56, night 0/55 (count p 1.00).
                (2) Memory, mean moment skewness fitted vs phi0: day 2.035 vs 2.048 (difference
                -0.013, -0.6%; 39/56 flies lower; sign-flip p 0.0004); night 1.839 vs 1.850 (-0.011,
                -0.6%; 37/55; p 0.0012).
                Descriptive: observed skewness mean 2.18 (day), 1.99 (night) -- close to the
                exponential's 2. Bowley obs/fitted/phi0 day 0.304/0.315/0.310, night
                0.210/0.220/0.217 (memory moves Bowley UP by 0.003-0.006). SD obs/fitted/phi0 day
                1167/1145/1153 s, night 1757/1718/1730 s (memory lowers SD by 7-12 s, under 1%).
Matched null    (1) each fly's fitted-model simulations, scored identically; (2) paired phi0 simulations.
Verdict         (1) NOT REJECTED: the residuals' skewness is what his fitted model produces, in every
                fly. (2) SUPPORTED in direction (Holm 0.0016 day, 0.0036 night): memory lowers the
                moment skewness, but by 0.6%; the paired design detects a negligible effect.
Not reportable  Do not write "memory reduces the skewness" without its size (0.6%, and Bowley moves
                the other way). The GH-versus-shifted-exponential family change he observed cannot be
                attributed to memory on this evidence: a 0.6% skewness change does not move a family
                choice; the phi0 residuals share the same shape. The residuals are, to within 1%,
                the empirical bout distribution re-centred -- the result of the bout-space analysis
                (ledger 625dc08), reproduced in his framework. His scoring (bridging, record start)
                is not stated in the export; 24 of his 56 flies overlap our CS LD set.
```

### 2026-11-09  G. Collaborator's latent residuals: Gaussian as assumed

```
Commit          2559473 (PRESPEC only; self-contained MATLAB script on his export); MATLAB R2026a.
Script          F_run/run_G.m, session scratch, not in the repo
Prespec         2026-11-09 entry G (committed 2559473 before the run)
Time window     His CS LD day and night phases; each fly's latent residuals within a phase.
n               Day 56 flies, night 55; 200 fitted-model replicates per fly.
Result          Flies outside their fitted-model range (p < 0.05): skewness day 4/56 (count p 0.31),
                night 1/55 (0.94); excess kurtosis day 2/56 (0.78), night 1/55 (0.94). Holm 1.00 for
                all four. Means obs/fitted/phi0: skewness day -0.018/-0.001/-0.002, night
                -0.045/+0.002/+0.001; excess kurtosis day +0.059/-0.021/-0.030, night
                +0.030/-0.040/-0.043; SD day 0.984/0.986/0.994, night 0.975/0.985/0.992 (sqrt(1 - phi^2)
                ~ 0.998).
Matched null    Each fly's own fitted-model simulations with Gaussian innovations, scored identically.
Verdict         NOT REJECTED: the latent residuals are Gaussian as far as these flies can show; the
                Gaussian latent AR(1) describes the dependence between successive sleep bouts.
Not reportable  This is largely guaranteed by construction (the randomized PIT makes each latent
                marginal Gaussian, and phi ~ 0.03), so it confirms the model rather than adding
                biology. It does not test time reversal, and says nothing about wake bouts or other
                genotypes (not in the export).
```

### 2026-11-09  H. Long sleep bouts cluster beyond the collaborator's Gaussian latent AR(1); short ones do not

```
Commit          61a040d (PRESPEC only; self-contained MATLAB script on his export); MATLAB R2026a.
Script          F_run/run_H.m (+ fig_H.m), session scratch, not in the repo
Figure          extreme_bout_clustering_2026-11-09.png (analysis folder)
Prespec         2026-11-09 entry H (committed 61a040d before the run)
Time window     His CS LD day and night phases; pairs of consecutive sleep bouts within one phase
                segment of one fly (narrowest block).
n               Day 56 flies (median 113 bouts), night 55 (median 85); 200 fitted and 200 phi0
                replicates per fly; 1000 within-fly permutations.
Result          Group mean C (1 = no clustering), observed | fitted [2.5-97.5%] | phi0 | shuffle:
                day upper 10%   1.880 | 1.155 [0.921-1.398] | 0.912 | 0.913   p 0.0100 (MC floor)
                day lower 10%   0.964 | 1.169 [0.922-1.365] | 0.919 | 0.918   p 0.0995
                night upper 10% 1.402 | 1.059 [0.790-1.318] | 0.881 | 0.886   p 0.0100 (MC floor)
                night lower 10% 0.993 | 1.088 [0.860-1.308] | 0.901 | 0.893   p 0.4478
                PRIMARY Holm: day upper 0.040, night upper 0.040, day lower 0.199, night lower 0.448.
                Descriptive 20% tails: upper day 1.365 vs 1.111 (p 0.010), night 1.463 vs 1.056
                (0.010); lower day 0.989 vs 1.114 (0.040), night 0.970 vs 1.058 (0.159).
                Per fly (secondary, low power): upper 4/56 day (Binomial p 0.31), 4/55 night (0.30);
                lower 0/56, 0/55; observed above its fitted mean in 42/56 and 33/55 flies (upper),
                24/56 and 23/55 (lower). Small-sample baseline without memory is ~0.9, not 1.
Matched null    Each fly's own fitted-model simulations: same marginal, block lengths and n.
Verdict         REJECTED for long bouts in both phases (observed beyond all 200 fitted groups): a
                long sleep bout is followed by another long one ~1.9x (day) / ~1.4x (night) as often
                as chance, against ~1.1x in his model. NOT REJECTED for short bouts, which sit at the
                no-memory level, BELOW his model. The bout-to-bout dependence is asymmetric --
                concentrated in long bouts -- whereas a Gaussian copula spreads its (tiny) phi
                symmetrically over both tails. This qualifies ledger G: the Gaussian latent AR(1) fits
                the residual moments but not the dependence of extremes.
Not reportable  Cannot yet separate tail dependence from a slow drift within a phase segment (e.g.
                long bouts concentrated at one time of night/day): both give clustered long bouts, and
                onsets are not in the export. Not a time-arrow test; CS LD sleep only.
```

### 2026-11-09  I. Long-bout clustering: time of day at night; beyond time of day by day; between 12-h stretches in per0 DD

```
Commit          d66e6df (PRESPEC only; script uses no repo function, data from homeo_results.mat of
                the preliminaries entry); MATLAB R2026a.
Script          I_run/run_I.m (+ fig_I.m), session scratch, not in the repo
Figure          long_bout_clustering_source_2026-11-09.png (analysis folder)
Prespec         2026-11-09 entry I (committed d66e6df before the run)
Time window     Pairs of consecutive sleep bouts within one 12-h segment of one fly (CS LD: day ZT
                0-12 or night ZT 12-24 of one day; per0 DD: 12-h stretches from record start).
n               CS LD 40 flies, per0 DD 48; 1000 permutations per null and fly; gaps per fly median 0.
Result          Upper 10%, group mean C: observed | N1 ZT 1-h | N2 within-segment | N3 plain (p):
                CS LD day    1.511 | 1.104 (0.016) | 1.183 (0.022) | 0.915 (0.002)
                CS LD night  1.629 | 1.550 (0.663) | 0.797 (0.002) | 0.888 (0.002)
                per0 DD      1.340 | 0.930 (0.002) | 1.290 (0.597) | 0.954 (0.002)
                PRIMARY Holm: CS LD day vs N1 0.032; CS LD night vs N1 0.663; per0 DD vs N3 0.006.
                N1 at 2 h / 3 h (pre-specified sensitivity): day 1.086 (0.004) / 1.095 (0.002); night
                1.499 (0.51) / 1.439 (0.27); per0 0.931 / 0.944 (0.002 both) -- same reading as 1 h.
                Lower 10%: no clustering under any null (all p >= 0.19): day 0.817, night 1.079,
                per0 1.041. Consecutive pairs sharing a 1-h bin on the same day: 0.53 day, 0.39
                night, 0.44 per0.
Verdict         CS LD NIGHT: the clustering is time of day -- fully reproduced by keeping each fly's
                ZT profile (pooled across days), not by keeping each night's level. CS LD DAY:
                survives the ZT profile and survives each day's level -> a slow state within the
                day or bout-to-bout memory (not separable with these nulls). per0 DD: clustering
                without a clock (Holm 0.006), fully reproduced by keeping each 12-h stretch's
                level -> slow variation between half-days (a fly-specific state or a trend over
                the record). Short bouts never cluster. Entry H (his export) is reproduced in our
                flies under the plain shuffle.
Not reportable  Day: a combined ZT-by-day null was not pre-specified; whether the day excess is
                within-day slow drift or memory needs its own entry. per0 DD: half-day differences
                cannot be split into trend (e.g. first day) vs fluctuating state here. Sleep only.
```

### 2026-11-09  J. Daytime long-bout clustering is a slow state within the day, not detectable memory

```
Commit          efe3256 (PRESPEC only; script uses no repo function; data as entry I); MATLAB R2026a.
Script          J_run/run_J.m (+ fig_J.m), session scratch, not in the repo
Figure          long_bout_state_vs_memory_2026-11-09.png (analysis folder)
Prespec         2026-11-09 entry J (committed efe3256 before the run)
Time window     Pairs of sleep bouts 1-4 apart within one 12-h segment of one fly; N4 permutes only
                within one 2-h ZT window of one day of one fly (narrowest block).
n               CS LD 40 flies, per0 DD 48; 1000 permutations per null and fly (Seed 2).
Result          PRIMARY, CS LD day, upper 10%:
                  1 observed C 1.511 vs N4 day x 2-h 1.347: p 0.106, Holm 0.212 (1 h 1.427 p 0.29;
                    3 h 1.426 p 0.48). NOT REJECTED.
                  2 lag profile, excess over N1 at lags 1-4: 0.407 0.305 0.101 0.283; Dlag 0.178,
                    N1 95% -0.331..0.308, p 0.262, Holm 0.262. NOT REJECTED (no decay after lag 1).
                Descriptive: per0 DD observed 1.340 vs N4 2-h 1.309 (p 0.63), excess 0.412 0.367
                0.111 0.332, Dlag p 0.19. CS LD night observed 1.629 vs N4 2-h 1.262 (p 0.006; 1 h
                1.414 p 0.016; 3 h 1.201 p 0.012) although N1 (ZT pooled over days) reproduced it in
                entry I (1.550, p 0.66); night lag-2 C 0.952 vs N1 1.476. Lower tails: no clustering
                (all p >= 0.07). Occupancy, day x 2-h bins: bouts alone 0.03 / 0.08 / 0.04, in bins
                >= 3: 0.88 / 0.72 / 0.83 (day / night / per0).
Verdict         CS LD day: keeping each day's 2-h windows brings the null from 1.10 (time of day
                only) to 1.35 of the observed 1.51, and the remainder is not significant; the excess
                does not fall off after one bout. The daytime clustering is a slow state on the
                scale of hours within each day, with no detectable bout-to-bout memory. per0 DD:
                the same, without a clock. CS LD night: UNRESOLVED -- the ZT-pooled null explains
                it and the day x ZT null does not; within-day permutation of small bins samples
                without replacement and can push C down, so N4 may be biased low at night (fewer
                bouts per bin; the night N2 of entry I also fell below the plain shuffle).
Not reportable  Absence of memory is "not detected", with N4 95% ranges about +-0.15 around 1.35.
                A state faster than ~2 h cannot be told from memory by any shuffle. The night
                discrepancy needs a calibrated null (simulated flies with a known ZT profile) before
                any night claim.
```

### 2026-11-09  K. per0 DD half-day variation in sleep-bout length: mostly the first day; fluctuation test failed calibration

```
Commit          04cbce3 (PRESPEC only; script uses no repo function; data as entries I-J); MATLAB R2026a.
Script          K_run/run_K.m (+ fig_K.m), session scratch, not in the repo
Figure          per0dd_trend_vs_state_2026-11-09.png (analysis folder)
Prespec         2026-11-09 entry K (committed 04cbce3 before the run)
Time window     Trend: the whole record of one fly. Fluctuation: 12-h segments of one fly.
n               per0 DD 48 flies (CS LD 40 descriptive); 1000 permutations (Seed 3); trend 10000
                sign flips.
Calibration     100 simulated 48-fly groups on the real per0 onsets, 199 permutations. Fluctuation
                test rejection: (a) stationary 0.01, (b) common trend 0.04, (c) fly trends 0.03,
                (d) common first-day x1.5 0.26 -> FAILS the 0.095 limit, (e) real 12-h state 1.00.
                Trend test: 0.07 / 1.00 / 0.08 / 1.00 / 0.08. Per the prespec, the FLUCTUATION test
                is DEMOTED TO DESCRIPTIVE. Trend test passes.
Result          TREND (primary): mean slope -0.0205/day (-2.0%/day), 18/48 flies positive, p 0.032,
                Holm 0.032. SUPPORTED, small: per0 sleep bouts shorten ~2%/day on average.
                FLUCTUATION (descriptive): group logRatio 0.271 (null 95% upper 0.011), p 0.001,
                36/48 flies > 0; without detrending 0.419; trends remove 27% of the excess variance.
                First 24 h dropped: logRatio -0.041, p 0.048 (null centre is negative), trend -2.5%/day
                p 0.15. 4-h segments: -0.000, p 0.19. Records < 4 d (19 flies) 0.027 p 0.065; >= 5 d
                (28) 0.412 p 0.001. Lag-1 correlation of successive detrended 12-h levels -0.174 vs
                null -0.225 (p 0.10). Long-bout clustering C on detrended durations 1.174 vs shuffle
                0.948 (p 0.020; raw 1.340). Group levels (figure): highest in the second 12 h
                (~1.2x the fly's average), then a slow decline.
                CS LD (descriptive): day trend -6.3%/day (p 0.0001), no day-to-day fluctuation after
                detrending (logRatio -0.42, p 0.64; raw p 0.006); night trend -2.0%/day (p 0.18), no
                fluctuation (p 0.67). C on detrended durations day 1.351 vs 0.900 (p 0.004), night
                1.345 vs 0.886 (p 0.002).
Verdict         The per0 DD half-day differences that reproduced long-bout clustering in entry I
                are mostly the first day of the record (the excess collapses from 0.27 to near the
                null when the first 24 h are dropped, and the calibration shows a first-day effect
                alone fakes this test), plus a small decline. Some clustering remains after
                detrending (1.17 vs 0.95), and the >= 5-day records still show excess variance, so
                a weaker recurring state is NOT excluded -- but it is not established. In CS LD the
                day-to-day differences are a trend (day bouts shorten 6%/day), and the within-day
                clustering of entries I-J is untouched by detrending.
Not reportable  "Clock-free slow state" as a finding: deflated to "mostly a first-day / acclimation
                effect" for per0 DD until a first-day-robust test is pre-specified. Whether the
                first-day effect is acclimation to the tube or to DD cannot be told from these data.
```

### 2026-11-09  L. per0 DD after the first day: long bouts still cluster, explained by differences between 12-h stretches

```
Commit          13ebef5 (PRESPEC only; script uses no repo function; data as entries I-K); MATLAB R2026a.
Script          L_run/run_L.m (+ fig_L.m), session scratch, not in the repo
Figure          per0dd_after_day1_2026-11-09.png (analysis folder)
Prespec         2026-11-09 entry L (committed 13ebef5 before the run)
Time window     Sleep bouts from 24 h after the record start; pairs within one 12-h stretch of one
                fly (narrowest block).
n               per0 DD 48 flies, median 172 retained sleep bouts; 1000 permutations (Seed 4).
Calibration     100 simulated groups x 199 permutations, rejection P / S12 / W2: (a) stationary
                0.02/0.01/0.01; (b) trend -2%/d 0.02/0.03/0.05; (c) acclimation tau 18 h
                0.06/0.08/0.06; (d) fly trends 0.12/0.01/0.03 (P not required to hold here);
                (e) 3-h state, POWER 0.80/0.60/0.17. All required sizes hold; no test demoted.
Result          Upper 10%, from 24 h: observed C 1.356 | P 0.938 (p 0.002, Holm 0.006) | S12 1.211
                (p 0.184, Holm 0.184) | W2 1.203 (p 0.060, Holm 0.120).
                Lag profile excess over S12: 0.145 -0.071 -0.127 -0.098; Dlag 0.244 (S12 95%
                -0.270..0.276), p 0.106. Per fly: outside S12 range 1/48, W2 2/48. Lower 10%: 1.090,
                no test rejected (p >= 0.18).
                From 48 h (descriptive): upper 1.264 | P 0.912 (p 0.008) | S12 1.117 (p 0.23) | W2
                1.184 (p 0.39); lower 1.231 | P 0.904 (p 0.014) | S12 0.994 (p 0.058).
Verdict         Reading "1 rejected, 2 not": after the first day per0 long bouts still cluster
                (Holm 0.006), and keeping each 12-h stretch's bouts together reproduces it -- slow
                differences between 12-h stretches (fly-specific drift, or slow variation on >= 12 h;
                not separable here), not a clock-free state lasting hours within a stretch, and no
                detectable memory. Power caveat: a 3-h state of SD 0.3 would be caught by the S12
                test only 60% of the time, so a moderate hours-scale state is not excluded.
Not reportable  "Clock-free hours-scale state" -- not supported. The CS LD daytime within-day drift
                (entry J) therefore has no clock-free counterpart in per0 at this power.
```

### 2026-11-09  M. Neighbouring sleep and wake bouts: a coupling finer than 2 h, and long wake is followed (not preceded) by long sleep in per0 DD

```
Commit          6cdeca7 (PRESPEC only; script uses no repo function; data as entries I-L); MATLAB R2026a.
Script          M_run/run_M.m (+ fig_M.m, diag_M.m), session scratch, not in the repo
Figure          sleep_wake_coupling_2026-11-09.png (analysis folder)
Prespec         2026-11-09 entry M (committed 6cdeca7 before the run)
Time window     per0 DD bouts from 24 h after the record start; CS LD whole records from ZT0. Pairs =
                time-adjacent bouts within one 12-h block of one fly (narrowest block).
n               per0 DD 48 flies, median 169 wake->sleep and 164.5 sleep->wake pairs per fly; CS LD 40
                flies, median 170 / 174.5. 1000 permutations (Seed 5); sign flips 10000.
Calibration     As pre-specified (100 groups x 199 permutations, 2000 flips), rejection R vs P / R vs
                W2 / Dtail: per0 (a)-(d) all 0.00; CS LD (a) 0.01/0.01/0.00, (b)-(d) 0.00. All sizes
                hold; no test demoted. FLAW FOUND AFTER THE RUN: the permutation seed was the same in
                every simulated group, so the 100 groups were not independent and the near-zero
                rates are not a size estimate. Unplanned check (diag_M.m), 150 per0 groups under (a)
                with fresh seeds, 99 permutations: R vs P rejection 0.033 (null z SD 0.975), Dtail
                0.020 -> both tests are calibrated, Dtail slightly conservative. POWER: (b) a shared
                2-h state of SD 0.3 and (d) a homeostat of power 0.15 were essentially never
                detected -- both are far weaker than the effects seen in the data.
Result          per0 DD: coupling R 0.104 | P -0.000 (p 0.002) | S12 0.020 (p 0.002) | W2 0.047
                (p 0.002). r_pre 0.113, r_post 0.096, Delta 0.017 (sign-flip p 0.53).
                Dtail 0.537 [Cpre 1.85, Cpost 1.32], sign-flip p 0.007; flies > 0: 27/48.
                Reverse tail (long sleep: wake after minus before) -0.44, p 0.028 (same direction:
                long wake sits BEFORE long sleep; shares the long-long pairs, not independent).
                CS LD: R 0.037 | P -0.000 (p 0.002) | S12 0.012 (p 0.004) | W2 -0.002 (p 0.002) | Z
                -0.022 (p 0.002). r_pre 0.044, r_post 0.030, Delta p 0.46. Dtail 0.264 [1.52, 1.26]
                p 0.33, 18/40 (descriptive); reverse tail +0.30 p 0.078.
                PRIMARY Holm: all five 0.010 (per0 R vs P, per0 R vs W2, per0 Dtail, CS LD R vs P,
                CS LD R vs W2).
Verdict         Q1, reading "rejected vs W2": neighbouring sleep and wake bouts are coupled beyond
                anything that shuffling within the same 2-h window keeps -- in per0 about half the
                coupling (0.047 of 0.104) is the slow layer, the rest is at finer than 2 h. It is
                symmetric (wake->sleep = sleep->wake), so it looks like a shared fast "tempo"
                (consolidated vs fragmented stretches), not a one-way homeostat; the test cannot
                tell a bout-to-bout dependence from a state lasting well under 2 h. In CS LD the
                coupling is a third the size but likewise exceeds the 2-h and clock-hour shuffles.
                Q2: in per0 DD after day 1, a long wake is followed by a long sleep 1.85x chance and
                preceded by one 1.32x chance (Holm 0.010) -- a homeostat-direction asymmetry in the
                tail, and a time-reversal asymmetry of the clock-free sleep/wake sequence. It is
                carried by a subset of flies (27/48 positive). CS LD shows the same direction, not
                significant (descriptive).
Not reportable  "Bout-to-bout memory" as such (a sub-2-h state is not excluded). A bulk homeostat
                (the symmetric correlation is not directional). A per-fly homeostat trait (27/48).
```

### 2026-11-10  Correction to K, L, M: calibrations rerun with independent seeds (and M's counter fixed)

```
Commit          Amendments 64d1afe (independent seeds) and c3c11b6 (M counter); scripts
                recal/calK.m, calL.m, calM1.m, calM2.m (session scratch). Data sections not rerun;
                every data p-value in K, L and M stands.
Defects         (1) K, L, M: every simulated group used the same permutation seed, so the 100 groups
                were not independent. (2) M only: `rej + [p] < 0.05` without parentheses never
                accumulated, so every M rate was 0.00 or 0.01. The original M calibration was void;
                the "fresh-seed check" in entry M (0.033 / 0.020) came from a separate, correct
                script and stands.
K               Fluctuation / trend rejection: (a) 0.01/0.08 (b) 0.03/1.00 (c) 0.02/0.06 (d) 0.25/1.00
                (e) 1.00/0.06. Fluctuation test size still fails in (d) -> stays DESCRIPTIVE.
                K unchanged.
L               Rejection P / S12 / W2: (a) 0.03/0.02/0.01 (b) 0.02/0.03/0.05 (c) 0.05/0.10/0.03
                (d) 0.14/0.02/0.04 (e) POWER 0.78/0.61/0.12. S12 size 0.10 in (c) exceeds 0.095 ->
                the S12 test is DEMOTED to descriptive (was primary, p 0.184, not rejected). The
                verdict stands: a non-rejection by a slightly liberal test is no weaker. Power for
                a 3-h state: S12 0.61 (was 0.60), W2 0.12.
M               Rejection R vs P / R vs W2 / Dtail.
                per0: (a) 0.01/0.05/0.03 (b) 1.00/0.04/0.04 (c) 0.05/0.01/0.03 (d) 1.00/1.00/1.00.
                CS LD: (a) 0.03/0.04/0.09 (b) 1.00/0.04/0.05 (c) 0.03/0.08/0.05 (d) 1.00/1.00/0.99.
                All required sizes hold; no test demoted; Holm results of entry M stand.
                CORRECTED CLAIM: entry M said the 2-h state and the homeostat (power 0.15) were
                "essentially never detected" -- wrong (counter bug). The whole-record shuffle
                catches the 2-h state every time, and all three tests catch the homeostat every
                time. What holds: a shared 2-h state does NOT pass the 2-h-window test (0.04) or
                the long-wake test (0.04-0.05), so neither result can come from a slow state.
                New point: a homeostat ALSO passes the 2-h-window test (1.00), so "finer than 2 h"
                alone does not tell a homeostat from a fast tempo; the direction symmetry of the
                data (r_pre 0.113 vs r_post 0.096, p 0.53) is what argues against a bulk
                homeostat, since scenario (d) acts on wake -> sleep only. CS LD Dtail size 0.09 in
                (a) is borderline; it was descriptive already.
```
