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
