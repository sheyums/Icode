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
