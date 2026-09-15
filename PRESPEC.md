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
