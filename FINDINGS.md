# Findings

Results and methodological decisions, kept with the code so a later session
picks them up. Conventions and traps live in `CLAUDE.md`.

## per0, DD, 3987 sleep bouts (xmin = 300 s, dt = 1 s)

Selected model: **hyperexponential, K = 3**

| | tau (s) | tau (min) | q, observed | surviving truncation | w, untruncated |
| --- | --- | --- | --- | --- | --- |
| 1 | 54.6 | 0.91 | 4.96% | 0.417% | 0.904 |
| 2 | 1041 | 17.4 | 92.98% | 75.0% | 0.094 |
| 3 | 3141 | 52.4 | 2.06% | 90.9% | 0.002 |

- AICc 63247.01, BIC 63278.44, Akaike weight 0.998, k = 5
- Goodness of fit: G = 39.103 on 40 bins, chi-square p in [0.2513, 0.4653]
  (Chernoff-Lehmann bounds). Passes; the bounds do not straddle alpha, so no
  bootstrap was spent.
- Mixture order settled by parametric bootstrap LRT: **LR = 16.72, p < 0.001**
  (B = 999; null 95th percentile 4.71, 18% of the null at zero exactly).

### Why the order needed a test rather than a criterion

AICc separated K=3 from K=2 by **12.71** while BIC separated them by **0.14** --
posterior odds near 52:48. Both are correct: the third component is real but
small (2% of observed bouts), and BIC's penalty of `log(3987) = 8.3` per
parameter very nearly cancels a genuine 8.36-nat improvement. Neither criterion
answers "is the component real"; the LRT does, once its null is simulated
rather than assumed to be chi-square.

**Reporting**: a third exponential component is supported (bootstrap LRT
LR = 16.7, p < 0.001, B = 999) though it accounts for only ~2% of bouts, which
is why BIC is nearly indifferent between K = 2 and K = 3.

### Caution on the fast component

`tau_1 = 54.6 s` sits 5.5x **below** the 300 s threshold, so only 0.417% of
that component's mass survives truncation and its untruncated weight carries a
**240x back-transform gain**. `w_1 = 0.904` does not license "90% of bouts are
short" -- that is extrapolation into the region the protocol excluded by
construction. The defensible numbers are the `q` column.

### Other rows, for context

`gamma` and `powerlaw_cutoff` were bit-identical (-31646.75), as were `erlang`
and `hyperexp K=1` (-31665.24) -- the nesting identities holding on real data.
`powerlaw` lost by 1483 AICc: these bouts are not scale-free. `pearson3` was
excluded, its location having run up against xmin at 299.852; note its AICc
would have placed it third, above gamma, so the gate mattered. `K=4` and `K=5`
were excluded as collapsed or non-identified, `K=4` with logL identical to K=3.

## per0, DD, 3989 wake bouts (xmin = 2 s, dt = 1 s)

The hazard observation came first and stands on its own; the model below was
found afterwards, and reproduces it.

**The empirical hazard has two turning points**: 1.2e-2 at 2 s, falling to
1.1e-3 by 70 s, rising to 1.65e-3 near 500 s, then falling again past 1000 s.
Nothing in the original library can produce that shape, and the exclusion is by
construction rather than by evidence:

- A hyperexponential's hazard is **strictly decreasing at every order** — it is
  a sum of decreasing exponentials — so no K will do. The symptom was exactly
  that: best was K=2 with G = 77.0 on 40 bins (G/df = 2.14, p ~ 1e-4), and
  K=3/4/5 bought 0.33 nats between them.
- Every other family in the original library has a **monotone** hazard, and
  `exp_weibull` allows at most **one** turning point.
- **Mixing cannot rescue it.** `h_mix(t) = sum_i w_i(t) h_i(t)` with weights
  shifting toward the longer-lived components, so pooling flies can only make a
  hazard fall faster, never rise. Between-fly heterogeneity is not the
  explanation.
- A **PH(2)** hazard is monotone (O'Cinneide), so the smallest phase-type that
  can hump has `sum(m_j) >= 3`.

Sleep bouts show no such hump, so the asymmetry is wake-specific.

`hyper_erlang` and `weibull_mix` were added for this.

### The model: hyper-Erlang, three branches, one of them in series

Selected model: **`hyper_erlang`, shapes [1 1 5]**, k = 5, logL = -29309.219.

| branch | q, observed | mean | kind |
| --- | --- | --- | --- |
| 1 | 80.7% (SE 2.1%) | 660 s | memoryless, rate 1.515e-3 |
| 2 | 7.9% (SE 0.6%) | 9.5 s | memoryless, rate 0.1053 |
| 3 | 11.4% (SE 1.9%) | 720 s | **5 sequential stages**, rate 6.949e-3 |

- Goodness of fit: **G = 38.92 on 40 bins, df [34, 39], p in [0.258, 0.474]**.
  Passes on both Chernoff-Lehmann bounds. Akaike weight 0.99886.
- `weibull_mix` second at dAICc 13.55, itself ambiguous on the fit test
  (G = 50.69, p in [0.033, 0.099]).
- **The best hyperexponential is rejected**: K=2 at dAICc 34.04 with
  G = 76.96, p ~ 1e-4. K=3 no better. Exactly as the hazard argument requires.
- **And the exponential branches are exhausted, not merely outscored**: K=2 to
  K=3 buys 0.32 nats for two parameters, K=3 to K=4 another 0.30 and lands on
  the rate ceiling (excluded, tau = 1 s).
- The step-down did NOT fire: J=3 was identified, and the row carries the bare
  `hyper_erlang` name.

**The cleanest form of the evidence is at MATCHED complexity.** `hyperexp K=3`
and `hyper_erlang` both have k = 5:

    hyperexp K=3   three memoryless branches          logL -29327.9229
    hyper_erlang   two memoryless, one 5-stage series logL -29309.2191
                                                           ----------
                                                            18.70 nats

The penalties cancel, so dAICc = 37.41 is exactly twice that, and no argument
about how k is counted or whether the integer shape should be charged can touch
it. Spending the same five parameters on a SERIES branch rather than a third
exponential buys 18.7 nats.

**Three nesting identities check out on this data**, through different code
paths: `hyperexp K=3` equals `hyper_erlang`'s ShapeSweep at m=1 (-29327.9229
both, since shapes [1 1 1] IS the 3-component mixture), `powerlaw_cutoff`
equals `gamma` (-29450.2139), and `erlang` at shape 1 equals `hyperexp K=1`
(-29572.4554). Implementation checks, not coincidences.

**The runner-up is untested, not rejected.** `weibull_mix` has G = 50.69 with
the Chernoff-Lehmann bounds straddling alpha, which is exactly the case
`GoFBootstrap="auto"` exists to resolve -- but the walk-down stops at the first
row that passes, and `hyper_erlang` passed first. Claiming `weibull_mix` is
excluded requires running that bootstrap deliberately.

**The guard reports health here, not just absence of failure.**
`MixtureMinCount` = 315 of 3989 bouts in the smallest branch, sixty times the
floor -- worth quoting given that the phantom-component bug below concerned a
branch holding 2.5e-18 of its data while looking respectable.

**The mechanism, and why it needs a series branch.** Branches 1 and 3 have
nearly the same mean -- 660 s against 720 s -- and opposite hazard shapes. The
9.5 s branch makes the hazard high at 2 s and then depletes; the Erlang branch's
RISING hazard lifts it through the middle; that branch then dies off faster than
the slow exponential (asymptotic rate 6.9e-3 against 1.5e-3), so the hazard falls
back to the slow branch's constant. Two turning points, from three branches, none
of which individually has one.

**The fit reproduces the hazard it was never shown.** Computing h(t) from the
fitted parameters alone gives a trough of 1.31e-3 at t = 126 s and a peak of
1.74e-3 at t = 896 s, falling to 1.52e-3 by 4000 s -- against the 1.2e-2 at 2 s,
1.1e-3 near 70 s, 1.65e-3 near 500 s read off the data. Nothing in the
likelihood targeted the hazard, so this is a check rather than a restatement.

### What must NOT be reported

**Not "five stages".** The shape sweep is flat where it matters:

    m:      1          2          3          4          5          6
    logL:  -29327.92  -29312.11  -29310.06  -29309.29  -29309.22  -29309.48
    gain:       +15.81     +2.05      +0.77      +0.075     -0.26

A series branch is decisive -- m=1 to m=2 buys **15.8 nats**. Past m=3 the
profile is flat, m=5 beats m=4 by 0.075 nats, and `H.k` charges nothing for the
integer so no criterion penalised taking the largest. The defensible claim is
**a series branch of at least two, and probably three or more, stages**.

**The p-value is optimistic.** Both non-monotone families were added AFTER the
original library failed on these same bouts, so p in [0.26, 0.47] is a fit test
on a family chosen having seen this data. The hazard argument is what stands
without qualification; the G test supports it and does not establish it.
Re-testing on held-out flies is the fix, and has not been done.

**Truncation is nearly harmless here, unlike the sleep bouts.**
`TailFraction` = 0.989, so q = [0.807, 0.079, 0.114] against untruncated
w = [0.793, 0.095, 0.111]. The 240x back-transform caution that governs the
sleep result does not bite at xmin = 2 s -- but it is a property of THIS
threshold, not a general reprieve.

### The hazard with bands: `BoutHazard` on these bouts, measured

MATLAB R2026a, `wakebouts_per0dd.mat` (`w_bouts`, n = 3989, xmin = 2 s,
dt = 1 s). `Bootstrap=999`, `RandomSeed=1`, `MinAtRisk=10`, the `hyper_erlang`
fit as `SurvivalHandle`. Rates per second, 95% Wilson bands.

At 16 bins (15 kept; the one dropped is the final bin, saturated by construction):

    center (s)  at risk  events   hazard     lower      upper      model
       1.41      3989      46    1.16e-2   8.70e-3   1.55e-2   9.94e-3
       2.45      3943      34    8.66e-3   6.20e-3   1.21e-2   9.15e-3
       3.87      3909      69    8.91e-3   7.04e-3   1.13e-2   8.11e-3
       6.71      3840      85    5.60e-3   4.53e-3   6.92e-3   6.40e-3
      12.0       3755     112    4.33e-3   3.60e-3   5.20e-3   4.29e-3
      21.2       3643      98    2.27e-3   1.86e-3   2.77e-3   2.50e-3
      36.7       3545     125    1.80e-3   1.51e-3   2.14e-3   1.57e-3
      63.1       3420     160    1.37e-3   1.17e-3   1.60e-3   1.33e-3
     109.3       3260     237    1.24e-3   1.09e-3   1.41e-3   1.31e-3   <- trough
     190.1       3023     400    1.33e-3   1.20e-3   1.46e-3   1.32e-3
     330.8       2623     593    1.39e-3   1.28e-3   1.50e-3   1.42e-3
     574.5       2030     844    1.67e-3   1.56e-3   1.79e-3   1.63e-3
     997.7       1186     731    1.72e-3   1.59e-3   1.85e-3   1.71e-3   <- peak
    1733          455     351    1.52e-3   1.35e-3   1.70e-3   1.57e-3
    3012          104      97    1.60e-3   1.20e-3   2.02e-3   1.52e-3

The fitted model's binned hazard is inside the band in every kept bin at both
bin counts.

| | 16 bins | 20 bins |
| --- | --- | --- |
| trough bin | 109 s, 1.237e-3 [1.089, 1.405]e-3 | 161 s, 1.258e-3 [1.116, 1.418]e-3 |
| peak bin | 998 s, 1.717e-3 [1.591, 1.849]e-3 | 944 s, 1.732e-3 [1.596, 1.876]e-3 |
| `RiseRatio` | 1.388 | 1.377 |
| `RiseCI` (descriptive only) | [1.279, 1.866] | [1.317, 1.824] |
| `RiseDisjoint` | 1 | 1 |
| `RiseNullP` (null = hyperexp K=2 fit, 999 draws) | 0.058 | 0.156 |
| `NonMonotone` (= `RiseNullP` < 0.05) | 0 | 0 |

Before 4f49911 the trough could land on the last kept bin: at 20 bins it did
(2850-4435 s, 35 at risk, 30 events, 1.228e-3 [0.773, 1.748]e-3), giving
`RiseRatio` 1.000 and `NonMonotone` 0 on the same data that read 1.388 at 16
bins. 4f49911 excludes the final bin from the trough search; the table above is
after that fix.

**The criteria, measured on data with NO rise.** 50 datasets of n = 3989
simulated from the hyperexponential K=2 fit to these bouts (strictly decreasing
hazard), each run through `BoutHazard` with the settings above, statistic as of
bbccb5a (before the trough fix):

| claims a rise | 16 bins | 20 bins |
| --- | --- | --- |
| `NonMonotone` (then `RiseDisjoint` OR `RiseCI(1)` > 1) | 10/50 | 8/50 |
| `RiseCI(1)` > 1 | 9/50 | 5/50 |
| `RiseDisjoint` | 3/50 | 3/50 |
| null `RiseRatio` median / 95th pct / max | 1.117 / 1.288 / 1.331 | 1.112 / 1.378 / 1.526 |

That is what retired `RiseCI` as a test in 4f49911. `RiseDisjoint` has not been
re-measured under the new trough rule.

**What `RiseNullP` measures, and why it reads high.** Two defects in its null,
separated by a 2x2 run (999 datasets from the same K=2 fit, observed ratio as
above, p = (1 + #null >= observed)/(M + 1)):

| null datasets drawn by | binned on | 16 bins | 20 bins |
| --- | --- | --- | --- |
| engine convention P(T=k) = S(k-1)-S(k) | each dataset's own bins | **0.019** | **0.057** |
| `BoutHazard`'s null sampler | each dataset's own bins | 0.019 | 0.054 |
| engine convention | the observed data's fixed bins and keep mask | 0.047 | 0.185 |
| `BoutHazard`'s null sampler | fixed observed bins (= `RiseNullP`) | 0.048 | 0.173 |

1. **The null sampler is one grid step early** -- it starts its grid at xmin
   and forces S = 1 there, but the engine's truncated survival is 1 at
   xmin - dt. E[T] 610.518 s against 611.508 s; P(T=2) = 0.019666, exactly the
   engine's P(2) + P(3); every later value equals the engine's at one step
   higher to 2.7e-5. **It does not move the p-value** (rows 1 and 2).
2. **Binning null replicates on the observed data's edges and keep mask does**
   -- p rises from 0.019 to 0.048 at 16 bins and from 0.057 to 0.173 at 20
   (rows 1 and 3). The observed statistic went through its own bin and keep
   selection; the replicates did not, so the two are not exchangeable.
   `RiseNullP` as shipped (0.058, 0.156) agrees with the fixed-bin rows to
   within Monte Carlo error.

Row 1 is the exchangeable version: **p = 0.019 at 16 bins, 0.057 at 20.**
Caveat on every row: the null is the K=2 model FITTED to these bouts, used as
if known; real use would refit it per dataset, which this does not capture.

Both defects are addressed in 2371913 (own bins per replicate, sampler grid
shift). The `RiseNullP` values in the table above (0.058, 0.156) are from
4f49911, BEFORE that fix. **At 2371913 on these bouts** (same settings,
NullReplicates=999, RandomSeed=1):

| | 16 bins | 20 bins |
| --- | --- | --- |
| `RiseNullP` | **0.010** | **0.068** |
| null ratio median / 95th pct | 1.136 / 1.285 | 1.186 / 1.390 |
| row 1 above (independent 999 datasets) | p 0.019, 1.138 / 1.303 | p 0.057, 1.184 / 1.392 |
| `NonMonotone` | 1 | 0 |

The fixed test agrees with row 1 to within Monte Carlo error (differences 0.009
and 0.011, about 1.5 and 1.1 standard errors of a difference between two
999-draw estimates), and the null distributions match. Trough, peak,
`RiseRatio` and `RiseDisjoint` are unchanged from the table above.

**The bin count decides the side of 0.05.** Same data, same null, same
statistic: p = 0.010 at 16 bins, 0.068 at 20. Choosing the bin count after
seeing these would be selecting the analysis on its outcome.

**Decided (user, 2026-09-14): report the rise as supported but bin-sensitive,
across a set of bin counts fixed IN ADVANCE, every count reported rather than
the minimum.** The alternatives were one count pre-specified by a rule, or a
statistic that does not bin. **Not yet run**: the set has not been fixed and
the user is looking at other datasets first, so the rule applies there too
(judgement call 8). Until it is run, quote the two measured values above as
exactly that -- two bin counts -- not as a range. The structural argument above
-- no hyperexponential can hump, and at k = 5 the series branch buys 18.70
nats -- does not rest on this p-value.

**`weibull_mix` fit-test bootstrap: still running.** Observed G = 50.69 on 40
bins, df [34, 39], chi-square p in [0.0327, 0.0995]. Staged rule: 199
replicates; stop if p < 0.02 or > 0.12, else pool 200-replicate batches under
new seeds until the 95% Monte Carlo interval clears 0.05 (cap 4999).
First batch (seed 2, finished 10:31): **B = 199 valid, p = 0.080, MC SE 0.019,
95% [0.042, 0.118]** -- inside the undecided range, so pooling continues at
~5.5 s per replicate. Not a result until the interval clears 0.05.

## Mixture weights: the guard has to see what the data see

A mixture component can hold a large share of the **untruncated** weight and
none of the observations. On synthetic humped data at `xmin = 100`, a
3-component hyper-Erlang returned

    w     = [0.055, 0.809, 0.136]     rates = [0.019, 0.566, 0.0074]

The middle component carries the largest weight of the three and has mean
1.8 s, so `S(100) = 3e-25` and its observed share is `1.5e-24` — **zero of 1500
bouts**. A guard testing `min(w)*n` read 82.3 and passed it; the phantom
component then bought a stage count the data never supported.

Two consequences worth carrying into any mixture result:

1. **Report and judge on `q`**, the observed weight `~ w_j S_j(xmin)`. This is
   the same caution as the 240x back-transform on the sleep bouts' fast
   component, in a form that also corrupts model *selection* rather than only
   interpretation.
2. **The same truncated law has several parametrizations and only some look
   degenerate.** The fit above and one with `w = [0.288, 1.3e-10, 0.712]` have
   identical log-likelihoods to ten digits; they are the same distribution
   above `xmin`. Where the optimizer stops is not evidence about
   identifiability.

Correctly sized, the model was strongly preferred on that data: `hyper_erlang`
at **J=2** gave m=4, logL = -6905.149, AICc = 13816.32 with k=3, against the
best hyperexponential's AICc = 14263.31 — a gap of **447**. The J=3 fits had
the *same* log-likelihood as J=2; they were J=2 plus a branch holding nothing,
charged two extra parameters.

### Stage counts are weakly identified — report the branch mean

Data generated with an Erlang(3) branch are routinely fitted with m=2, in two
independent implementations, with the survival still within 0.015 of truth. The
**branch mean** `m/lambda` is what the data determine; the integer stage count
is not, and no information criterion here charges for having swept over it.
Quote `ShapeSweep` and how flat it is, not the selected integer alone.

## Pearson III: a rejected ridge that outscores every legitimate model

On test 25's synthetic humped data (n = 1500, xmin = 100, MATLAB), Pearson III
returned

    logL = -6828.18    shape 162058   scale 0.738   location -119120
    GuardOK = 0        LocationGap = 1192.2

against the best legitimate model on the SAME data, hyper-Erlang at J=2, with
logL = -6841.39. **The rejected fit scores 13 nats higher than the winner.**

It is rejected by the span-ratio arm of `pearson3Guard`, not the `xmin` arm:
the location has run 119,000 s BELOW the data, which spans a few thousand, so
`(t - location)` is nearly constant over every observation, the density
flattens to `const * exp(-t/scale)`, and shape and location stop being
separately identifiable. What the optimizer reports is a point on a ridge, not
a maximum -- the gamma's normal limit, a huge shape with the location pushed
far out. `pearson3Guard`'s header already documented this at
`location = -5.34e5` on real bouts, where "that fit beat the true model and won
a comparison outright". This is the same failure on different data.

Two things follow, and they pull in opposite directions.

**The row cannot win, and it cost more than the rest of the library
combined.** 516 s for THREE starts, so roughly 70 minutes for the engine's 26,
because every start walks that ridge to its evaluation budget. Pearson III is
therefore **no longer in the default library** -- pass `Models={'pearson3',
...}` to fit it.

**But the ridge is a real hazard, and the default table no longer shows
it.** Anyone who disables the guard, or who reads the AICc column without
checking `Degenerate`, selects a fit whose parameters mean nothing. A
three-parameter family with a free location will do this to left-truncated
data whenever the data are close to exponential over their observed range;
the guard is the only thing standing between that ridge and a published
model. If a reviewer asks why Pearson III is absent, the answer is that it
was excluded for being unidentifiable here, not for fitting badly -- it
"fits" better than anything else.

## Judgement calls, open to revision

These are choices, not results. Worth revisiting before publication.

1. **G-test verdict uses `pLower`**, the conservative Chernoff-Lehmann bound
   (df = B-1-k) rather than the lax one (df = B-1). A model called acceptable
   has passed the harder reading.
2. **`GoFBootstrap="auto"`** spends a bootstrap only when the two df bounds
   straddle alpha, or when the lax p lands in 0.01 to 0.2. Elsewhere both
   bounds already agree.
3. **`OrderLRTThreshold = 2`** on *either* criterion triggers the order ladder.
   Two is Burnham & Anderson's "no meaningful difference" band.
4. **The order ladder does not control a family-wise error rate.** Each rung is
   valid at its own alpha, but climbing several and stopping on the first
   acceptance inflates the overall type-I rate. Quote the rungs individually;
   do not present the endpoint as carrying one alpha-level guarantee.
5. **`HyperErlangComponents` stays at 3**, decided by the user on the
   condition that J=3's suitability is ASSESSED per dataset rather than
   assumed. It is: J=3 is attempted on every dataset and refused only when
   that dataset's own fits come back unidentified, so nothing anywhere
   encodes "J=3 does not work". On two-regime data it is refused and J=2
   wins by 447 AICc (Octave) and 514 (MATLAB), on different draws.

   First decided as RELY ON THE REFUSAL, then revised by the user the next
   morning (2026-09-14) to a **STEP-DOWN ON REFUSAL**, because relying on the
   refusal meant an unsupported J=3 yielded no hyper-Erlang row at all --
   only the `R.Skipped` entry -- rather than the J=2 fit that works and wins
   comfortably, and the reader had to rerun by hand to get it.

   Now: when J is refused, `CompareBoutModels` refits at J-1 down to 2 and
   the first identified order enters the table as `hyper_erlang J=j`. The
   refusal stays in `R.Skipped`, its Reason naming that row.
   `HyperErlangStepDown=false` restores the bare refusal. Still no J sweep.

   **Cost, measured rather than asserted.** An earlier draft of this entry
   said "paid for only when J fails", which is true but invites the wrong
   inference. What J fails ON is the point: data with fewer regimes than
   branches. Ordinary gamma(0.7, 500) at n = 300-500 SUPPORTS J=3, so the
   step-down never fires there and costs nothing (verified two ways: an
   Octave default-library call returns the plain `hyper_erlang` row with
   `R.Skipped` empty, and a MATLAB A/B with the option on and off differs by
   run-to-run noise). Where it does fire -- two-regime data -- it costs about
   1.3 s, a second sweep at the lower order. And the row it produces is the
   right one: the stepped-down `hyper_erlang J=2` has logL -6841.391,
   identical to an explicit J=2 run on the same data.

   What this choice carries, stated so it is not rediscovered:
   - **J is now chosen after seeing the data**, and no information criterion
     charges for it -- the same caution as the stage count. It is mild: the
     step fires only on a refusal, never to improve AICc, and stops at the
     first identified order.
   - **It is not replayed in a bootstrap.** A GoF replicate refits at the
     order the data landed on, not by re-running the step-down.
   - **The order is in the row name** because a J=2 fit reported as
     `hyper_erlang` would be read as the 3-branch model the data refused.
6. **`erlangCDFint` switches to a direct tail sum below `F = 1e-6`.** The
   threshold is a speed/accuracy trade, not a derived constant: `1 - S` has
   relative error about `eps/F`, so 1e-6 keeps that below ~1e-10 while leaving
   the slow branch rare. Lower it if a fitted range ever needs bin
   probabilities finer than that.
7. **`B = 999`** for the LRT. A test has a decision boundary the estimate must
   resolve: the Monte Carlo error is `sqrt(p(1-p)/B)`, so a true p of 0.05 has
   a 95% interval of [0.020, 0.080] at B=199 against [0.036, 0.064] at B=999
   (Davison & Hinkley 1997, sec. 4.2). Even 999 is not tight at the boundary --
   go higher rather than report a borderline p.

8. **A hazard rise is reported across a PRE-SPECIFIED set of `NumBins`,
   every count, never the minimum p.** `BoutHazard`'s `RiseNullP` depends on
   binning: on the per0 DD wake bouts it was 0.010 at 16 bins and 0.068 at 20.
   Fix the set before looking at any p-value (e.g. 12/16/20/24), report all of
   them, and call the rise bin-sensitive if they straddle alpha. Decided by the
   user 2026-09-14; not yet applied to any dataset.

## One thing that does not transfer between implementations

The bootstrap null absorbs the *estimation procedure*, not just the model. On
the same test data MATLAB put **41%** of null LRs at zero where Octave put
**16%**, because MATLAB's multistart collapses the extra component more
readily. That is correct methodology -- replicates get the same treatment as
the data -- but it means a p-value is a property of model *plus* procedure.
Quote it with B and the multistart settings attached, and do not expect exact
agreement with another implementation.

## References

- Chernoff & Lehmann (1954) *Ann Math Statist* 25:579-586 — df bounds for a
  chi-square statistic with estimated parameters
- Davison & Hinkley (1997) *Bootstrap Methods and their Application*, sec. 4.2
- Day (1969) *Biometrika* 56:463-474 — unbounded mixture likelihood
- Hartigan (1985) — mixture LRT non-regularity
- Hurvich & Tsai (1989) *Biometrika* 76:297-307 — AICc
- McLachlan (1987) *Appl Statist* 36:318-324 — bootstrap LRT for mixture order
- McLachlan & Peel (2000) *Finite Mixture Models*, ch. 6
- Mudholkar & Srivastava (1993) — exponentiated Weibull hazard regimes
- O'Cinneide (1990) *Stochastic Models* 6:1-57 — PH(2) hazards are monotone
- Tijms (1994) *Stochastic Models: An Algorithmic Approach* — hyper-Erlangs are
  dense in the distributions on [0, inf)
- Thummler, Buchholz & Telek (2006) *IEEE Trans Dependable Secure Comput*
  3:245-258 — EM fitting for hyper-Erlang (not used here; noted as the
  principled alternative to the shape sweep)
