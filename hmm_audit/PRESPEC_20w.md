# PRESPEC 20w — merged paper, test 3: the circadian "clock face" of the cycle's thermodynamics (descriptive)
Date 2026-11-01. Written before anything runs; the go-ahead to RUN is Sheyum's (not given by this file).
Owner: Plan B feasibility session. Results go to RESULTS_LEDGER_planB.md. DESCRIPTIVE: no test, no p-values, no
per-fly fits. It replaces Plan B's per-fly covariate model, which failed F1, by group-level description, where F1
showed recovery is possible (DD-CH: median clock amplitude ratio 0.99-1.06, phase error under 2 h).

## Data and bins (fixed now; never re-binned)
1-min symbols exactly as sqw_ref.sqw_symbols (theta 0.5; LD light_dark=True; DD light_dark=False).
Genotypes: CS LD, per0 LD, per0 DD (the flat negative control), CS DD if uploaded before the run. fumin and sss are not
included (not in the agreed design).
Bins: 12 bins of 2 h of clock time, bin k = minutes [120k, 120k + 120) mod 1440. In LD the bins are aligned to ZT0,
so no bin straddles a light switch. Each transition is attributed to the bin of its ARRIVAL minute (sqw_ref
convention).
CS DD alignment, fixed now: CLOCK TIME from the start of the cut file, with no per-fly phase estimate. With a
free-running period tau != 24 h, a fly's phase drifts by (24 - tau) h per day, at most about 2 h over 4 days for
tau in [23.5, 24.5]. This smears the profile; the smearing is stated, not corrected. Per-fly circadian phase would
need a fitted period per fly, which is the per-fly fitting F1 argued against.

## Quantities per genotype x bin
Pooled transition counts over all flies and days in the bin, giving the pooled first-order chain P_k (row-normalised
counts; no pseudocount; a one-way edge is reported as such, never floored). From P_k: J (cycles/h), affinity A,
EPR via sqw_ref.model_energetics, the landscape E3 - E1 via sqw_ref.landscape, and c_obs (the minute-level cycle
current). Bands: 2000 fly-bootstrap replicates (flies resampled with replacement, seed 2026), 95% percentile band.
The pooled-chain EPR is not the mean of the flies' EPRs (Jensen); it is labelled "pooled-chain EPR".

## Summary numbers (descriptive, fixed now)
Per genotype: profile range of J and of EPR (max minus min over the 12 bins) with its bootstrap band, and the bin of
the maximum. Reading planned: per0 DD's range is the noise-level reference. "Rhythmic" is said only where a
genotype's range band lies entirely above per0 DD's range band. In LD the profile mixes clock and light (masking).
per0 LD vs CS LD is the only contrast that speaks to the clock there, and it is descriptive.

## Narrowest block, matched null, time box
Narrowest block: one 2-h clock bin pooled over flies and days (stated with every number). The sample-size-matched
reference is per0 DD at the same bin width and comparable fly-days. Expected run time: minutes. Figure: 12-bin
circular "clock face" of J and EPR per genotype with bands.
