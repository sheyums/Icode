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

## per0, DD, 3989 wake bouts (xmin = 2 s, dt = 1 s) — NO MODEL YET

Recorded because the observation is solid even though the model is not.

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

`hyper_erlang` and `weibull_mix` were added for this. **Both were added AFTER
the original library failed on these same bouts**, so any G-test p-value on a
winner among them is optimistic — the family was chosen having seen the data.
The hazard argument above is the defensible claim; a passing fit-test on the
new families is not, and should be stated as such or re-tested on held-out
flies.

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
