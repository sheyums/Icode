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
