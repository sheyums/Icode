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
5. **`B = 999`** for the LRT. A test has a decision boundary the estimate must
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
