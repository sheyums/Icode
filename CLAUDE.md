# Icode — left-truncated bout-duration fitting

MATLAB (R2025a+) maximum-likelihood fitters for *Drosophila* sleep and wake
bout durations, plus a model-comparison pipeline. Everything is
**left-truncated**: bouts below a protocol threshold (`xmin`, e.g. 300 s for
the 5-minute sleep criterion) do not exist in the data, so every likelihood
conditions on `T >= xmin`.

## Rules for every analysis

These are working rules, not suggestions. They exist because the failures
recorded in `FINDINGS.md` were not coding errors -- the code was right and the
inference was not.

1. **Read `FINDINGS.md` before analysing data or changing a fitter.** It holds
   the results, the judgement calls, and the things that must NOT be reported.
2. **Write `PRESPEC.md` before running any new analysis, and commit it first.**
   A choice recorded after the fact is not a pre-specification, and a commit is
   what makes the timing checkable by someone who was not there.
3. **Append every result to `RESULTS_LEDGER.md`, including null results**, with
   the date and the commit hash of the code that produced it. The ledger is
   append-only. A file holding only the analyses that worked cannot tell you
   how many were run.
4. **Never widen a block or re-bin after seeing a result.** A **block** is a
   time window -- a contiguous stretch of recording, such as a ZT range, a
   `startZT` segment, or a set of days. If the effect is not there in the
   window you specified, that is the finding.
5. **When reporting a result, state the narrowest block it was computed within,
   and what it looks like under a sample-size-matched null.** The narrowest
   block is the smallest window the reported number was actually computed from
   -- not the widest one that contains it. Matching the null's sample size
   matters because almost every statistic here moves with `n`.

Rule 4 is why `BoutHazard`'s bin count is fixed in advance and every value
reported (FINDINGS judgement call 8: `RiseNullP` was 0.010 at 16 bins and 0.068
at 20 on the same wake bouts -- picking after looking would have let either
conclusion be written up). Rule 5 is why the hyper-Erlang p-value is called
optimistic in FINDINGS rather than quoted plainly: the family was chosen after
the earlier library failed on those data.

## Layout

Two unrelated kinds of fitter live here, and they fit different physical
quantities. The engine and its `Fit*MLE` wrappers fit **bout durations** --
left-truncated positive times. `shiftlognormal_MLE`, `GenNorm_MLE`, `multigaussian_MLE` and
`GeneralizedHyperbolic_MLE` fit **noise distributions** -- the recorded signal
itself, real-valued and untruncated. A log-likelihood from one side is
meaningless against one from the other, and not for the usual reason: it is
not that the criteria penalise differently or the dominating measure differs,
it is that the two are fitted to DIFFERENT DATA. Never let a noise fitter into
a bout comparison, whatever its AIC.

### Bout-duration suite

| File | Role |
| --- | --- |
| `FitTruncatedDiscreteMLE.m` | **Shared engine.** One truncated likelihood, nine families in a `switch` registry. Adding a family costs ~15 lines: `ParamNames`, `cdf`, `sf`, `logpdf`, `unpack`, `valid`, `starts`; a mixture also needs `report`/`nReport` and a `guard`. |
| `FitHyperexponentialMLE.m` | Mixture of exponentials, K=1..N, parametrized in **observed** weights `q`. |
| `FitExponentiatedWeibullMLE.m` | Exponentiated Weibull; `FixAlpha=true` nests the plain Weibull. |
| `FitHyperErlangMLE.m` | Mixture of Erlangs, shapes `[1..1 m]` swept over m. The exponential-native route to a **non-monotone hazard**. The sweep is warm-started (`WarmStart`, `SweepStarts`) via the engine's `StartZ`; cold starts are reduced, never removed, so a warm start can only add a candidate optimum. Refuses rather than return a fit whose guard rejected every shape. |
| `FitWeibullMixtureMLE.m` | Two-component Weibull mixture; the other non-monotone-hazard family. |
| `FitGammaMLE` … `FitBetaMLE` (8 files) | Thin wrappers over the engine. |
| `CompareBoutModels.m` | Fits the whole library, ranks by AICc/BIC, G-tests the winner, plots it. `pearson3` and `beta` are **opt-in**, not part of a default run. |
| `HyperexponentialLRT.m` | Parametric bootstrap LRT for mixture order. |
| `BoutHazard.m` | Life-table hazard with Wilson bands. Does NOT rank models — it decides whether a whole FAMILY can work, by asking whether the hazard rises. |
| `test_*.m` (6 files) | 129 tests. Run each by name from this directory **in MATLAB**. See Testing for the Octave caveat. |

**The bout-duration library is CLOSED.** The families `CompareBoutModels`
fits are the only ones to be used for bout distributions -- the user's
decision, 2026-09-14 local. By default: `hyperexponential` (K=1..N),
`exp_weibull`, `weibull`, `gamma`, `powerlaw_cutoff`, `powerlaw`, `erlang`,
`chisquared`, `weibull_mix`, `hyper_erlang`; `beta` and `pearson3` opt-in by
name. The engine's registry is cheap to extend (see its row above) but that is
a statement about the CODE, not an invitation: do not add a family, and do not
port one in from elsewhere, unless asked for it directly.

The analysis folder holds about 37 further `*_MLE` files from roughly ten years
of separate work, and **most of them DO fit durations** -- `ExpoMLE`,
`GammaMLE`, `LogNormMLE`, the `*ExpoMLE` mixtures, `PlMLE` and others. (An
earlier version of this file said they fit noise or other quantities. That was
wrong; only `shiftlognormal_MLE`, `GenNorm_MLE`, `multigaussian_MLE`,
`GeneralizedHyperbolic_MLE`, `AsymLaplace_MLE` and `shiftweibull_MLE`, in the
section below, are noise fitters here; all six are called by the analysis
folder's `SelectBestModel`, `GenNorm_MLE`, `multigaussian_MLE`,
`AsymLaplace_MLE` and `shiftweibull_MLE` joined Icode from that folder in Nov
2026, and all six are audited.) They are not candidates for three reasons
that have nothing to do with what quantity they fit: the user decided the
library is closed, they are **continuous** where `CompareBoutModels` is
discrete by default, and most condition on **`[xmin, xmax]` with `xmax`
defaulting to `max(data)`** -- double truncation whose normalizer depends on
the largest observation, which is a different likelihood, not a different
parametrization of ours.

If a dataset cannot be described by this library, that is a FINDING to report
-- as the per0 wake bouts were, where every monotone family failed and the
hazard said why -- not a cue to widen the library.

### Noise-distribution fitters — NOT bout fitters

These fit the noise in the recorded signal, not durations. Nothing in
`CompareBoutModels`, `BoutHazard` or the engine applies to them: no `xmin`, no
`SamplingInterval`, no truncation, no hazard argument, and no shared AIC table.

| File | Role |
| --- | --- |
| `shiftlognormal_MLE.m` | Shifted-lognormal noise fitter. Positive support with a free shift, **untruncated**. **Audited Nov 2026** against an independent dense profile likelihood (branch `claude/shifted-lognormal-mle-wc03ob`). Found and fixed: (i) NOT equivariant under offsets -- `x + 1e8` moved the shift by 2.3% and LogLike by 0.049 (`fminbnd`'s relative tolerance on the raw shift, ~1.5 at 1e8); the search now runs over `log(min(x) - shift)` and the shifted data are formed as `(x - min(x)) + gap`, equivariant to ~1e-8; (ii) the upper bound `min(x) - 0.001*range(x)` EXCLUDED the true shift for heavy right tails (n=2000, true shift 3: sigma 1.5 clipped 20/20, median error 0.45) -- it now scales by `median(x - min(x))` (0/20, error 0.006); (iii) fits on a bound stopped just inside it (fminbnd never evaluates its ends) -- the bounds are now candidates; (iv) new `SmallSample` warning below n=15, where the unbounded spike at `min(x)` wins (134/200 fits at n=5, 9/200 at n=10, 0/200 from n=15). Data with no right skew clip to the lower bound (warned). **SEs and CIs added (Nov 2026):** `SE` from the analytic observed information (matches an independent numerical Hessian to 1e-6), `CI` by profile likelihood (recommended), `CIWald` on constraint-respecting scales, `CIOpen` for ends that reach a bound. Measured coverage of nominal 95%: 0.94-0.97 for all three parameters at n=200 and 1000; at n=50 profile 0.93, Wald mu 0.90. Base MATLAB only (`fminbnd`). |
| `GenNorm_MLE.m` | Generalized Normal (generalized Gaussian) MIXTURE noise fitter, K components of `[mu, sigma, alpha, weight]`, real-valued support. Added to Icode Nov 2026 from the analysis folder and **audited** there against an independent closed-form K=1 profile and an independent multistart. Found and fixed: the `AIC` field held AICc (now plain AIC, plus an `AICc` field -- a BREAKING change for code that read `AIC`); tiny samples returned garbage (K=2 at n=2..7 gave finite AICc of -25 to -97; now K needs n >= 4K+1); NOT equivariant (an absolute 1e-6 sigma floor and raw-unit fitting: `1e-6*x` moved the K=2 LogLike by 10.1, `x + 1e8` by 36; now fitted in z-units); mixture fits stuck in local optima (an independent 40-start search beat the old fitter by 11.8, 5.97 and 1.09 log-likelihood units in 3 of 4 K=2 samples; random starts and a polish added -- see its test 12 and LIMITATIONS); alpha <= ~1.5 stalled at the cusps `mu` = data point (up to 0.08 short at K=1); `kmeans` used and advanced the caller's global random stream (now private, restored); `alpha` sat on its bounds 0.5/10 silently (now `AlphaAtBound`); failures were silent (now `Converged`/`NotConverged`). **SEs and CIs added (Nov 2026):** `SE` from the observed information, except `mu` of a component with alpha < 1.5, which uses OPG (the cusp makes the Hessian meaningless: Hessian-Wald coverage of mu was 0.04 at alpha 0.7, 0.33 at 1, 0.88 at 1.3); Wald `CILower`/`CIUpper` by default, profile likelihood for mu/sigma/alpha with the `Profile` mode (about 6x slower). Measured coverage of nominal 95%: Wald 0.93-0.97 at alpha >= 1.3 and for K=2, but mu only 0.80 at alpha 0.7 and 0.91 at 1 (a `CuspSE` warning says so); profile 0.92-0.98 throughout, including alpha 0.7. Overlapping mixture components were not measured. The `Select` LRT for the number of components uses a chi-square reference that is NOT valid for mixtures -- a heuristic; prefer BIC. Needs **Optimization** (`fmincon`) and **Statistics** (`kmeans`, `iqr`, `quantile`, `chi2cdf`). |
| `multigaussian_MLE.m` | Gaussian MIXTURE noise fitter, K components of `[mu, sigma, weight]`, real-valued support. Added to Icode Nov 2026 from the analysis folder, **audited and rewritten** against an independent closed-form K=1 MLE and an independent Nelder-Mead multistart. Old-code defects measured: K=1 used the n-1 SD (not the MLE); K=2 on three populations was 84.7 log-likelihood units short of the optimum, K=3 cases 2.5-5.7; `1e-6*x` moved the LogLike by 191 (absolute sigma floor, raw units); `AIC` held AICc; n=2..7 at K=2 gave AIC -56 to Inf; constant data gave NaN; an unseeded `randn` start used the global stream. New engine: z-units, k-means++/quantile/grid/random/SPLIT seeds (so LL never falls as K grows), SQUAREM-EM screening on <= 5000 points ranked on the full data, `fmincon` with the analytic gradient on the best K+2, Newton polish. **Spurious components** (< 5 effective points, or sigma on the 1e-3*IQR floor -- chance clumps that the unbounded mixture likelihood rewards) are not returned unless nothing else exists (`SmallComponent`, `EffectiveN`). Faster than the old code in every timed case with LL equal or higher (n=50000, K=4: 26.5 s vs 83.7 s). SEs (observed information), Wald CIs, and profile CIs with `Profile`; measured Wald coverage 0.93-0.99 for separated components but **0.81-0.88 for overlapping ones** (profile 0.88-0.97). For over-fitted K an independent search still finds optima 3-7 units higher built on 5-8-point clumps; `Converged` is often false there. Needs **Optimization** (`fmincon`) only. |
| `AsymLaplace_MLE.m` | Asymmetric Laplace noise fitter, `[mu, sigma, kappa]` (Kotz-Kozubowski-Podgorski parameterization; kappa = 1 is the symmetric Laplace), real-valued support. Added to Icode Nov 2026 from the analysis folder (baseline `059f10f` = the working copy) and **audited** against a brute-force profile between the data points and Nelder-Mead profiles. The exact O(n) solver was already global; fixed: `AIC` held AICc (now plain AIC plus `AICc`); the bootstrap consumed the global random stream; an absolute 1e-12 degeneracy filter (scale-dependent); constant data returned `[mu NaN NaN]` (now an error); one-sided data silently gave a near-degenerate fit (now `OneSided` warning, `SEValid` false). New: closed-form Fisher-information SEs (the observed Hessian is meaningless at the cusp), exact profile-likelihood CIs (default), Wald, optional `BOOTSTRAP` on a private stream, `NOPROFILE`. **Measured coverage** (400 samples per row): profile 0.92-0.96 at n >= 100, 0.81-0.92 at n = 30; Wald for mu only 0.84-0.91 even at n = 500. Header has the full table. No toolbox needed except `prctile` in `BOOTSTRAP` mode. |
| `shiftweibull_MLE.m` | Three-parameter (shifted) Weibull noise fitter, `[scale, shape, shift]`, shape floored at 1.001 (below 1 the likelihood is unbounded). Added Nov 2026 from the analysis folder (baseline `059f10f`) and **audited** against an independent 40-start Nelder-Mead and an independent profile-over-gap fit: the old Nelder-Mead estimates were right in all 8 cases probed. Fixed: the bootstrap ended with `warning('on','all')`, switching ON every warning the caller had turned off, and consumed the global random stream; `AIC` held AICc; no flag when the shape sat on its floor. New: exact profile over the gap (global, deterministic), observed-information SEs, profile CIs (default), Wald, flags `ShapeAtFloor`, `ShiftAtLowerLimit` (Gumbel limit), `NonRegular` (shape <= 2), optional `BOOTSTRAP` -- which covers the shift only 0.62 at shape 1.5, n=100 (the old bootstrap-only version: 0.74), so do not use it for the shift. **Measured profile coverage** 0.94-0.97 at n >= 50 for shape 3-5, 0.89-0.96 at shape 1.5. Much faster than the old version (0.2-1.4 s per fit versus 5-8 s). |
| `GeneralizedHyperbolic_MLE.m` | GH / NIG / hyperbolic noise fitter, five physical parameters (`mu`, `lambda`, `alpha`, `beta`, `delta`); `AUTO` takes the better of NIG and hyperbolic, then chooses between it and GH by LRT and BIC. **Real-valued support and untruncated** -- the whole real line, so it does not even share a support with a duration law. **Audited and fixed in `dd7f613`** against ground truth that shared no code with the estimator (Bessel-free density by quadrature, independent Nelder-Mead, exact equivariance identities): open profile bounds now report 0 or Inf instead of overflow edges, the fit is scale-equivariant, and the header's special cases are corrected (hyperbolic is `lambda = 1`; the variance-gamma limit is `delta -> 0` with `lambda > 0`). Needs **Optimization Toolbox** (`fmincon`) and **Statistics Toolbox** (`chi2cdf`, `chi2inv`) -- unlike the bout suite. Three things survive the fix and are in its KNOWN LIMITATIONS. (i) **Exact VG (`delta = 0`) is outside the family** -- `delta` has a floor -- but VG-DISTRIBUTED DATA fit perfectly well: a true VG sample (`lambda = 1`, n=600) fitted with an interior `delta = 0.354`, and the GH density matches VG to 3.6e-5, 4e-9 and 4e-13 at `delta` of 1e-2, 1e-4 and 1e-6. VG-like data CAN drive `delta` to the floor -- at n=1000 the floored fit sat 0.022 nats from the unconstrained optimum -- and `delta_at_bound` then flags it, with an open lower bound. That openness is a MEASURED flatness of the profile down to the floor, not a consequence of construction. (ii) Finite bounds near the Normal ridge can still be ridge artifacts when `NoGainOverNormal` warns. (iii) An `alpha`/`beta` upper bound orders of magnitude above the estimate **can** mark where the nuisance fits stopped converging rather than a true crossing: at two `beta` bounds (1.2e6, 1.3e7) an independent profile was still below threshold, so those are too narrow and should be read as effectively open, while at the `alpha` bounds (1.2e6, 1.3e7, 1.16e9) the independent checker itself failed, leaving them unverified in either direction. Costs about **2.3x** more than before -- median GH fit at n=1000 was 171 s against 75 s, both under concurrent load. **Second audit (Nov 2026, branch `claude/shifted-lognormal-mle-wc03ob`, 8 findings):** `AUTO` now prefers a VALID fit (if both are invalid it applies the usual LRT/BIC rule and warns `BothModelsInvalid`) and computes profile CIs for the selected model only; GH reuses the polished NIG fit as its warm start. If a profile nuisance fit finds a lower NLL than the point estimate (by more than max(1e-3, Bessel noise)), the fit is re-optimized and the CIs recomputed, up to 2 times (`profile_refits`), else `ProfileImprovedOptimum` warns. New `student_t_limit` flag (`lambda < 0` and pushing `psi` 20 lower costs < 1e-3 nats) -- MEASURED to fire rarely, even on t(3), t(1.5) and Cauchy samples at n=1000: 11 of 12 fits sat at an interior gamma ~ 1/max\|x\|, and the one flagged fit (t(3), GH) cleared the tolerance by 0.15 millinats, so it is a rule check, not a t-data detector. Requires **n >= 10** (`TooFewObservations`), warns `SmallSample` below 50. Validity is judged on the log-density (new `LogFit` output), shared by the NLL and `Fit`; the Bessel-argument cap is now n-dependent, `min(1e10, 1e-2/(n*eps))`. On NIG data n=600 (suite 12C), `AUTO` took 17.3 s against 81.1 s for stand-alone GH (64.4 s) plus NIG (16.8 s). **Third revision (2026-11-04, same branch):** new `HYP` model (hyperbolic, `lambda = 1`; alias `HYPERBOLIC`), which `AUTO` now also considers (`LRT.Restricted`, `LogLike_HYP`, `pValue_NIG`, `pValue_HYP`); GH warm-starts from BOTH nested fits and keeps the better nested point, so its LogLike can never fall below theirs. **SEs added**: `SE`, `Cov`, `SEValid`, `CIWald` (log scale for alpha, delta), from a numerical observed information (`information`; central differences + Richardson on the exact NLL) -- `hessian` is still fmincon's quasi-Newton Lagrangian Hessian and must not be used for SEs. At the delta floor all SEs are NaN: conditional-on-delta Wald intervals covered lambda/alpha/delta 0 of 12 there. `NoProfile` mode skips the profiles (10-30x faster, Wald only, `CIMethod`). Profile crossings made by the Bessel-argument guard (alpha/beta upper bounds ~8e8 at n=500) are now reported open. **Measured coverage** (100 samples each, nominal 95%): NIG n=500 profile 0.95-0.96, Wald 0.94-0.98; HYP n=500 profile 0.90 (alpha) -0.97, Wald 0.96-0.97; GH n=1000 (`lambda = -1.5`) profile 0.95-0.98, Wald 0.79-0.80 for lambda and alpha -- prefer profile for GH. 12 of those 100 GH samples had their GLOBAL optimum in the variance-gamma corner (delta floor, lambda 1.9-5.4; confirmed by refitting from the truth). No fit in the 300 ended below the log-likelihood at the true parameters. A seed-paired rerun of 76 of those samples (HYP 50, GH 26) with the final code (profile steps capped at 1 z-unit) changed no fit; the only interval changes were two guard-made alpha bounds (5.9e8, 8.6e8) now Inf, and root shifts of 1e-3 relative or less. **Fourth revision (2026-11-05, same branch): skewness outputs.** `Skewness` (model moment skewness from the GIG cumulants; NIG closed form `3*beta/(alpha*sqrt(delta*gamma))`), `SkewnessSE` and `SkewnessCIWald` (delta method through the z-unit information; NaN when `SEValid` is false), `SkewnessBowleyModel` (quartile skewness of the fitted density, by quadrature and `fzero`) and `SkewnessBowleyData` (quartile skewness of the data with Parzen MID-quantiles, so integer-minute ties do not bias it: plain quartiles gave 0.333 on rounded data with true 0.067, mid-quantiles 0.056). `Skewness` is NaN at `student_t_limit`. Checked against a symbolic reference, the NIG closed form and a brute-force 4e6-point CDF (1e-8 to 1e-10). **Measured** (mu 0, alpha 2, beta 0.8, delta 1): `SkewnessCIWald` coverage NIG n=150 0.92, n=500 0.965, HYP n=500 0.94, GH n=1000 (lambda -1.5) 0.96 (82 of 100 fits had SEs); both Bowley versions unbiased to ~0.01, the data version 2.5-3.5x noisier than the model version. All these SEs assume iid observations -- for successive differences of a bout sequence they are TOO SMALL. Also fixed a third-revision bug: the Wald information was inverted in data units, so at data scale 1e-6 it failed the rcond test and every SE came back NaN; it is now inverted in z-units and mapped back. |
### Signal-level tools — fit no distribution at all

| File | Role |
| --- | --- |
| `chi2p.m` | Sokolove-Bushell chi-square periodogram against a block-permutation null. Makes no assumption about waveform shape, so a sharply peaked circadian profile scores on equal footing with a sinusoid. Answers *what period*, not *what distribution*. |
| `jsd_kde.m` | Jensen-Shannon distance between two samples by KDE, with bootstrap CIs and a noise-floor correction. Compares two empirical distributions; fits neither. Calls `ksdensity` (Statistics Toolbox). |
| `BoutIrreversibility.m` | Model-free test for a TIME ARROW in one bout type's sequence (added 2026-11-05). Forms `D_k = x(t+k) - x(t)` within blocks at lags 1, 2, 4, 8; for ANY reversible process `D_k` is exactly symmetric, whatever the bout distribution, so asymmetry can only come from order. Reports Bowley skewness (mid-quantiles), `UpFraction` (share of increases), moment skewness (DESCRIPTIVE ONLY) and Spearman memory `Rho`, per fly and as a group mean, against two within-block nulls: IAAFT surrogates (keep bout values and memory, reversible -- the null for the skew tests) and shuffles (the null for memory). **Measured** (200 simulated flies of 150 bouts per row, nominal 0.05): IAAFT rejects 0.025-0.065 for Bowley and UpFraction, with and without memory; the shuffle null over-rejects Bowley once memory is present (0.090 at lag 1); moment skewness is miscalibrated with BOTH nulls (IAAFT 0.06-0.115 even on iid bouts, shuffle 0.22 with memory), so it gets no p-value. Power against an AR(1) with exponential innovations, lag 1: UpFraction 0.875, Bowley 0.285. Group test on reversible data: 0.020 / 0.040. These calibrations are for the simulated processes only. **Strata option (profile-matched null, added 2026-11-05 for LD data):** each bout carries a label such as its 1-h onset ZT bin, and surrogates permute bouts only among positions with the same label, so the day/night profile is in the null. Needed because the profile alone fools IAAFT and shuffle (per-fly UpFraction rejection 0.070-0.135 on reversible LD-like simulated flies; strata 0.035-0.050). Its GROUP MEAN is biased (0.15-0.18 on 10-fly groups: permuting within ZT bins cuts the link between a bout's length and the next onset time), so no group strata p is reported; the group verdict is the COUNT of flies with per-fly UpFraction p < 0.05 against Binomial(n, 0.05), measured 0.030 on 40-fly groups. Bowley has per-fly strata p-values only (count would reject 0.130). Suite `test_BoutIrreversibility` (18 tests, 535 s; `'fast'` skips the 7 calibration tests). No toolbox. Claude-written, same self-consistency caveat. |
| `BoutAR1Noise.m` | Fits `x(t+1) = A + B*x(t) + eps` by least squares IN BOUT SPACE (seconds, no transform), per fly and one bout type, pairs within blocks (added 2026-11-05 to answer "can a linear AR(1) generate the bouts, and what is eps?"). Checks against shuffled (or, with `Strata`, profile-matched) sequences: `B` (memory), `Hetero` (Spearman of \|eps\| with x(t): the linear model says 0), `ShapeKS` (eps versus the no-dynamics eps of the same bouts), `BelowHalfB` (EAR(1) puts no pair below half of B; descriptive), and REGENERATION from the fitted A, B with resampled residuals, counting generated bouts below `Xmin` (300 s sleep, 2 s wake; never clipped). **Measured** (200 simulated flies per row): Hetero rejects 0.085 for EAR(1) and 0.695 for a Gaussian-copula AR(1); ShapeKS fires whenever there is memory; in LD the profile alone gives B ~0.18 and the shuffle null rejects B 0.465 (strata null 0.030). Regeneration below the floor: EAR(1) 0.000, copula AR(1) 0.086, and a negative B alone breaks a 2-s wake floor. For LD, call once per phase with each 12-h period of each day as a block. Suite `test_BoutAR1Noise` (10 tests, 1094 s). No toolbox. Claude-written. |
| `FitResidualFamilies.m` | Fits eight families to one residual vector by ML and returns comparable LogLike/AIC/BIC: shifted exponential, Normal, Laplace (closed forms), asymmetric Laplace, shifted lognormal, shifted Weibull, one-component GenNorm, and the GH family via `GeneralizedHyperbolic_MLE` `AUTO`/`NoProfile` (NIG or hyperbolic K = 4, full GH K = 5, read from its `Model` name -- matching "hyperbolic (" would wrongly catch the full GH). **Measured:** on 100 residual sets that ARE shifted exponential (n = 180), the shifted exponential won by BIC in 100/100. About 10 s per call, mostly GH. Suite `test_FitResidualFamilies` (4 tests, 789 s). Needs the six fitters' toolboxes. Claude-written. |
| `stressTest_jsd_kde.m`, `stressTest_chi2p.m` | Stress suites for those two: reported 36 and 44 tests (measured by another session in MATLAB R2026a at 37cf612; not re-measured here). **Not independent validation** -- per the git history in the analysis repo, the suites AND the code they test were both Claude-authored (every commit touching `chi2p`, the suites, and the recent `jsd_kde` / `GeneralizedHyperbolic_MLE` edits carries a Claude co-author line). A passing suite here means self-consistency, not a second opinion. Named `stressTest_*`, so anything globbing `test_*` misses them. |
| `stress_test_GH_MLE.m` | The GH suite: **137 tests, measured** in MATLAB R2026a on 2026-11-05 after the fourth revision: 137/137 passed in 8857 s (under concurrent load). Section 14, 14A-14F, tests the skewness outputs: NIG closed form, HYP against GIG moments by numerical integration, model Bowley against a brute-force CDF, data Bowley against `quantile()` and a hand value with ties, SkewnessSE by the delta method through `Cov`, and equivariance plus the z-unit Wald fix at scale 1e-6. The third revision measured 131/131 in 5156 s on a machine also running 4 coverage workers (the 120-test version took 1506 s unloaded). Section 13, 13A-13K, tests the hyperbolic model and the SEs against a closed-form hyperbolic density sharing no code with the fitter (alias, closed-form LogLike, an independent 30-start Nelder-Mead optimum, GH >= both nested fits, AUTO's restricted choice and p-values, SEs against an independent physical-coordinate Hessian to 1e-3, Wald formulas, NoProfile, SE equivariance, all-NaN SEs at the delta floor, no guard-made finite bounds). Before that: 86 original plus Section 11, 11A-11F, one guard per first-audit defect -- open bounds, `delta_at_bound` at the floor, unit equivariance, `single` input, the tiny-sample Normal limit; Section 1 now expects n < 10 to error; Section 12, one guard per second-audit fix plus a stress pass over seeds, sizes, heavy tails, near-Normal data, ties and outliers asserting invariants only). The pre-fix 94 was re-measured in MATLAB R2026a on 2026-11-04: 94/94 passed in 5046 s. The suite's own t sampler is toolbox-free (`trnd_rs`). Fixes 5 (`gamma = a/cosh(xi)`) and 6 (the n-dependent cap's formula) live in local functions and are covered by inspection plus indirect behaviour (12F checks the cap's noise bound at n=20000), not by a direct test. Claude-written, like the code it tests, so passing is self-consistency. Named `stress_test_*`, outside both the `test_*` and `stressTest_*` patterns. |
| `test_shiftlognormal_MLE.m` | The shifted-lognormal suite: **25 tests** -- 13 original smoke tests, 14-20 from the Nov 2026 audit, and 21-25 for the SEs and CIs (offset and scale equivariance, heavy-tail shift recovery, global optimum against an independent dense profile, `SmallSample`, a bound reached exactly, Fit/LogLike consistency). Measured in MATLAB R2026a on 2026-11-04: 25/25 in 8 s; on the pre-audit code tests 14-19 fail (6), as they should. A separate invariant-only stress pass (12 data shapes incl. Cauchy, uniform, ties, outliers, 1e7 offsets, 1e-7 scales; n 15-3000; both fitters, 264 cases) found 0 violations. Claude-written, like the code it tests, so passing is self-consistency; the ground truth inside it (a closed-form profile likelihood) shares no code with the fitter. |
| `test_GenNorm_MLE.m` | The GenNorm suite: **19 tests** (13-19 for the SEs and CIs), ground truth from an independent closed-form K=1 profile, a direct mixture density and an independent Nelder-Mead multistart. Measured 2026-11-04: 19/19 in 148 s; the pre-audit code passed 2 of the first 12 (the 10 audit guards fail). Claude-written; same self-consistency caveat. |
| `test_multigaussian_MLE.m` | The Gaussian-mixture suite: **16 tests**, ground truth from the closed-form K=1 MLE, a direct mixture density, numerical Hessians and an independent non-spurious multistart up to K=4. Measured 2026-11-04: 16/16 in 140 s. Claude-written; same self-consistency caveat. |
| `test_AsymLaplace_MLE.m` | The asymmetric Laplace suite: **15 tests**, ground truth from the density summed directly, a brute-force profile, independent Nelder-Mead profile bounds and a Monte Carlo Fisher information (2e6 draws). Measured 2026-11-04: 15/15 in 6.5 s. Claude-written; same self-consistency caveat. |
| `test_shiftweibull_MLE.m` | The shifted-Weibull suite: **16 tests**, ground truth from the density summed directly, a 40-start Nelder-Mead, independent Nelder-Mead profile bounds and a physical-coordinate Hessian. Measured 2026-11-04: 16/16 in 464 s (mostly the independent multistarts). Claude-written; same self-consistency caveat. |

## Analysing a new dataset

Works for any left-truncated durations or sizes, not only sleep and wake
bouts. **The order matters**: the hazard decides which families are
ADMISSIBLE, the likelihood decides which admissible one is BEST, and the
guards decide whether the winner's parameters mean anything. Run it the other
way round and you get a confident fit from a family the data had already
excluded.

**0. Establish `xmin` and `dt` first. They are facts about the protocol, not
choices.** `xmin` is the threshold below which durations do not exist -- 300 s
for a 5-minute sleep criterion, or just the scoring resolution if nothing else
censors. `dt` (`SamplingInterval`) is the recording grid. Everything conditions
on `T >= xmin` and bins on `dt`; a wrong `dt` silently destroys resolution, and
the grid-mismatch guard will name the spacing it detects instead.

**Build `xmin` as (integer number of bins) x `dt`**, using the same
multiplication that produced the durations, so `row >= xmin` is exact rather
than exact-to-rounding. From the sleep/wake pipeline (change notice
2026-09-21) that is `xmin = (floor(bridgeWakeSeconds/dt) + 1) * dt` for wake
bouts, and `ceil(300/dt) * dt` for the 5-minute sleep criterion. `xmin` is the
caller's to construct: the 2026-09-21 rewrite of
`extractBoutsForModelCompetition` removed its `xmin` output. (That rewrite is
independent of the same day's sleep=0/wake=1 recoding, which has NO bearing
here: this suite receives durations, never a trace, so the coding is resolved
upstream. The formula is convention-free in any case -- it follows from what
`ironout` BRIDGES, so the shortest surviving wake bout is `bridgeBins+1` bins
whichever digit means wake.)

**Nothing in this suite handles CENSORING.** This is a long-standing gap, not
a consequence of any recent pipeline change -- it was listed as an open issue
on 2026-09-21 but has been true throughout. The only mention of it in the repo
is a comment in `CompareBoutModels` noting its absence, and every likelihood
here treats each duration as observed in full. That has one known consequence
worth knowing BEFORE a phase-conditioned fit: the pipeline assigns a bout to a
phase by its START bin and does not clip it at the phase boundary, so a bout
that begins inside a ZT window and ends outside it enters the fit at its full
length, carrying out-of-phase time. **The affected bouts are the long ones** --
a long bout is likelier to span a boundary -- so the contamination concentrates
in exactly the tail where a hazard rise would be read. The direction of the
resulting bias on the hazard is NOT established; do not guess it, measure it
(and note that the per0 DD results in `FINDINGS.md` are whole-record, not
phase-split, so they are not affected). The principled fix is to treat a
boundary-crossing bout as right-censored at the boundary -- likelihood
contribution `S(c)/S(xmin)` instead of a bin probability -- which needs
censoring support added to the engine. Until then, a phase-conditioned
duration fit is reporting on a sample it has partly mismeasured, and rule 5
obliges you to say so.

**1. Look at the hazard BEFORE fitting anything.**

```matlab
H = BoutHazard(x, xmin, SamplingInterval=dt, NumBins=16, Plot=true);
```

A hazard that falls and then RISES excludes -- by construction, not by
evidence -- every hyperexponential at every order (a sum of decreasing
exponentials is decreasing), every other monotone family here, and the
exponentiated Weibull (one turning point at most). It also rules out
between-individual heterogeneity as the explanation, since mixing shifts weight
toward longer-lived components and can only make a hazard fall FASTER. No
information criterion will tell you any of this: they rank the candidates you
thought of. For a test rather than a picture, pass `NullSurvivalHandle` (a
fitted hyperexponential's `SurvivalHandle`) and `NullReplicates`, fix the
`NumBins` set before looking, and check `H.NullDropped` is 0.

**2. Fit the library.**

```matlab
R = CompareBoutModels(x, xmin, SamplingInterval=dt, MaxComponents=5);
```

Read it in this order: **`R.Skipped`** (a model that never competed is not a
model that lost), then `R.Table` (ranking, `Degenerate`, `Reason`), then
`R.Selected` and `R.SelectionPath` -- the walk-down stops at the first row that
passes the G test, so rows below the winner are UNTESTED, not accepted.

**3. Check the winner before believing it.**

| field | what it tells you |
| --- | --- |
| `Diagnostics.GuardOK` | false means parameters are unidentified, whatever the likelihood says |
| `Diagnostics.MixtureMinCount` | observations in the smallest branch; a handful means that branch is decoration |
| `Diagnostics.TailFraction` | `S(xmin)`, how much of the untruncated law survived truncation |
| `q` vs `Theta`'s `w` | report `q`. `w` is extrapolation below `xmin` |
| branch MEANS, not stage counts | `m/lambda` is identified; the integer `m` usually is not -- read `ShapeSweep` and see how flat it is |
| `k` vs `numel(Params)` | criteria use `k`; mixtures print one number more than they charge |

**4. When two orders are close, test rather than rank.** AICc and BIC can
disagree by design (they penalise differently), and neither answers "is this
component real". `HyperexponentialLRT` does, with a simulated null -- the
mixture-order LRT is non-regular, so no chi-square applies. `CompareBoutModels`
triggers it automatically when either criterion separates neighbouring orders
by less than `OrderLRTThreshold`.

**5. Prefer a comparison at MATCHED k.** When two candidates have the same
free-parameter count, the information-criterion penalties cancel and the
difference is pure likelihood -- immune to every argument about how `k` is
counted. On the per0 wake bouts, `hyper_erlang` against `hyperexp K=3` (both
k=5) differed by 18.7 nats, and `dAICc = 37.41` is exactly twice that.

**6. State the caveats that the numbers cannot.** If a family was added
BECAUSE the earlier library failed on these data, its goodness-of-fit p-value
is optimistic -- the family was chosen having seen the data. Say so, or re-test
on held-out individuals. Structural arguments (what a family's hazard CAN do)
do not carry that problem and are the stronger claim.

## Conventions that are not optional

**`SamplingInterval` is required everywhere**, in the same units as `xmin`.
It sets the discrete bin width and `n_min = max(1, round(xmin/dt))`. A wrong
value silently destroys resolution; the grid-mismatch guard warns and names
the spacing it detects in the data.

**`DistributionType` defaults to `"discrete"`.** Discrete log-likelihoods are
unit-free; continuous ones shift by `n*log(dt)`. **Never compare AIC/AICc/BIC
across modes** — different dominating measure, not different evidence.

**Report `q`, not `w`.** `q_j` is the observed-population weight (the fraction
of retained bouts from component *j*); `w_j` is the untruncated weight, related
by a factor `exp(+lambda_j * xmin)` that can reach 240x. When a time constant
sits below `xmin`, `w` is extrapolation into a region the protocol excluded.
Check `WeightsUntruncatedReliable`.

**`k` is not `numel(Params)`.** `Params` reports every parameter; `k` counts
the free ones. Hyperexponential at K=3 reports 6 numbers but has k=5, because
`sum(q)=1`. Information criteria use `k`. `nReported` records the difference.

**`q` is the OBSERVED weight in EVERY mixture** — hyperexponential, `hyper_erlang`
and `weibull_mix` alike. It is component *j*'s share of the bouts that were
actually retained, `q_j ~ w_j S_j(xmin)`. The fitted weights inside
`hyper_erlang` and `weibull_mix` are untruncated `w`; `M.report` converts them
at output. Never print a `w` under the name `q`: one fit had `w = 0.849` on a
component holding `2.5e-18` of the data.

**`Theta` is not `Params`.** `H.Theta` is the model's natural parameter vector,
what `M.sf`/`M.cdf`/`M.guard` are written in and the only thing safe to feed
back into the model. `H.Params` is what gets *reported*, and for a mixture it is
a transform of `Theta` with its own delta-method SEs. `nReport` exceeds `nPar`
there, because `sum(q)=1` leaves one weight determined but still worth printing.

**`pearson3` is opt-in: name it in `Models=` or it is not fitted.** It walked a
ridge to the optimizer's evaluation budget on every start -- 516 s for THREE
starts on 1500 bouts, so ~70 min for the engine's 26 -- and its own guard then
rejected the result, so the default library was spending most of its compute on
a row that could not be selected. `Models={'pearson3'}` still fits it, measured
at 35 s for n=120 in Octave, and `Degenerate` came back 0 there: the ridge
needs enough data to become the attractor, so do not read the guard as
unconditional. Fit it deliberately when you want to see the ridge -- on 1500
bouts the REJECTED fit outscored every legitimate model (`FINDINGS.md`).

**A model absent from `R.Table` did not necessarily lose.** `R.Skipped` lists
candidates that never entered the comparison, with the reason — usually a
mixture whose guard rejected every configuration, which is a statement about
the data rather than a malfunction. Absent from both is a bug; absent from the
table but present in `R.Skipped` is a finding. One case is in BOTH on purpose:
a `hyper_erlang` refused at `HyperErlangComponents=J` steps down to J-1..2,
so the refusal stays in `R.Skipped` while the first identified order enters
the table as `hyper_erlang J=j` — named for the order it was fitted at, never
the bare family name. `HyperErlangStepDown=false` turns it off.

**`S(xmin)`** is the *untruncated* survival at the threshold — the likelihood's
normalizer, reported as `Diagnostics.TailFraction`. **`SurvivalHandle`** is the
*truncated* survival `S(t)/S(xmin)`, which is 1 at `xmin`. Different objects;
one is a scalar, the other a curve.

**`BoutHazard`: only `RiseNullP` is a test, and it depends on `NumBins`.**
`RiseRatio` is a maximum taken after a minimum, so it is >= 1 for any curve;
`RiseCI` bootstraps the data, not a null, and is DESCRIPTIVE only (its lower
limit exceeded 1 on 9 of 50 no-rise datasets). A test needs
`NullSurvivalHandle` (a fitted monotone model, e.g. a hyperexponential's
`SurvivalHandle`) and `NullReplicates`; without them `NonMonotone` falls back
to `RiseDisjoint`, which is pointwise. `RiseNullP` moves with the bin count --
0.010 at 16 bins and 0.068 at 20 on the same wake bouts -- so **fix the set of
`NumBins` before looking, and report every one** (FINDINGS judgement call 8).

## Traps that have already bitten, with tests guarding them

- **A mixture guard must test the OBSERVED share, not the mixing weight.**
  `min(w)*n` misses the failure it exists to catch: under left truncation a
  component whose mass sits below `xmin` carries a large `w` and holds none of
  the data. A 3-component fit at `xmin=100` had `w = [0.055, 0.809, 0.136]`
  with rates `[0.019, 0.566, 0.0074]` — the LARGEST weight had mean 1.8 s, so
  `S(100)=3e-25` and it held zero of 1500 observations. `min(w)*n` read 82.3.
  Both `heGuard` and `weibullMixGuard` now test `w_j S_j(xmin)`, normalized.
  The phantom component had bought a stage count the data never supported.
- **The same truncated law has several parametrizations, and only some look
  degenerate.** The fit above and one with `w = [0.288, 1.3e-10, 0.712]` have
  identical log-likelihoods to ten digits — they ARE the same distribution
  above `xmin`. Judge identifiability on what the data see, never on where the
  optimizer happened to stop.
- **Never pass a name-value pair twice.** Octave's generated twins parse
  `varargin` in a loop and take the last; MATLAB's `arguments` block rejects a
  duplicate name outright. A one-platform failure the twins cannot catch,
  because the twins are what replace the `arguments` block.
- **Integer-shape Erlang needs no `gammainc`.** `S = exp(-x) sum_{j<m} x^j/j!`
  is exact, ~11x faster for the survival and ~6x for the CDF, and more accurate
  in the deep tail than Octave's `gammainc` (which drifts to 1.5e-3 relative at
  `x=0.1, m=8`).
- **Benchmark the function the hot path CALLS, not its sibling.** That 11x was
  measured on `erlangSFint`. The discrete likelihood calls the *CDF* for every
  bin edge, and `erlangCDFint` went unbenchmarked: it summed the upper tail
  directly whenever `x < m`, which is most of the fitted range for a slow
  branch, and was **0.7x** — 40% SLOWER than the `gammainc` it replaced. A
  full-library call that had cost 379 s was still running after 17 minutes,
  while the commit message advertised a speedup.
- **Order numerical branches by what threatens the answer, not by the obvious
  split.** `F = 1 - S` cancels with relative error about `eps/F`: harmless at
  `F = 0.1`, fatal at `F = 1e-300`. Branching on `x < m` put the slow path in
  the common case; branching on `F < 1e-6` puts it only where cancellation
  actually bites — and those are the small-`x` elements where the tail
  converges in a few terms, so the accurate branch is also the rare and cheap
  one.
- **The same seed is NOT the same data across implementations.** `rng(25)`
  then `rand` gives different draws in MATLAB and Octave, and a generator that
  rejection-samples diverges on the first rejection. So a log-likelihood
  measured in one and compared against the other is comparing two datasets, not
  two fitters. This produced a "77-nat" gap that was really 13 once both sides
  were measured on the same sample. Cross-platform claims need the same numbers
  computed twice, or the data written to a file and read by both.
- **Octave is not a proxy for MATLAB.** `*_oct.m` twins are local scaffolding
  (gitignored, regenerable) used for development. Octave has accepted three
  things MATLAB rejects: growing a struct array from `struct([])`,
  `struct2table` on a 1x1 struct with a `''` field, and it hid a false-positive
  guard because output was filtered to `[FAIL]` lines. **Anything touching
  struct arrays, tables or graphics needs a MATLAB run.**
- **A new MATLAB figure's theme is not predictable — pin it.** In R2026a
  batch mode, `figure('Name',..)` came out light while a `'Visible','off'`
  figure a moment later came out dark (axes `Color` 0.07, `XColor` 0.85), and
  another session got the opposite. `CompareBoutModels`' figure for the per0
  wake bouts rendered dark and its black data curve was invisible.
  `figure('Color','w')` alone does NOT fix it: the margin goes white while the
  axes stay dark. Create the figure, then `try theme(fh,'light'); catch, end`
  BEFORE any axes, so children inherit (`theme` is absent in Octave and before
  R2025a). `plotFit` does this and test 23 asserts a white figure and axes.
- **A statistical criterion is not validated by tests on single datasets —
  calibrate it on many under a null.** `BoutHazard`'s first rise criterion
  passed all 10 of its tests, including known-hazard laws, because they asked
  "does it detect a rise that exists, and decline one that doesn't" on one
  dataset each. Its error rate was never measured. Run on 50 datasets from a
  strictly decreasing hazard, it claimed a rise in 16-20% of them. Before
  quoting any p-value or accept/reject rule, simulate from the null it claims
  and count how often it fires.
- **A null replicate must pass through every data-dependent step the observed
  statistic did.** `RiseNullP` first binned each null dataset on the OBSERVED
  data's edges and at-risk/saturation mask, while the observed value had been
  computed after choosing those from itself. Not exchangeable, and it pushed p
  from 0.019 to 0.048 (16 bins) and from 0.057 to 0.173 (20 bins). A sampler one
  grid step early, found in the same code, moved p by at most 0.012 there --
  within Monte Carlo error -- on the real bouts with a fitted null. That is not
  a general reprieve for discretization (next entry): measure which defect
  matters IN THE SETTING AT HAND before fixing the one that looks worse. Fixed
  in 2371913: each replicate runs the whole estimator on its own edges, and the
  sampler uses `pmf(k) = S(g_k - dt) - S(g_k)`.
- **A calibration harness must discretize its simulated data exactly as the
  estimator's null sampler does.** Calibrating `RiseNullP` (pre-2371913) on
  synthetic data generated with `round()` -- nearest bin -- while the null
  sampler followed the engine's convention (a bout of k bins lies in
  `((k-1)dt, k dt]`, i.e. `ceil`) gave min p 0.080, median 0.610, 0/40 below
  0.05: a clean, confident "the test is conservative". Changing only the
  harness to `ceil()` gave min 0.025, median 0.542, 3/40 below 0.05 --
  approximately uniform. The conservatism was the HARNESS. The two traps above
  are about the estimator; this one is about the thing it is tested with,
  which nothing else here guards. The same half-bin question barely moved p on
  the real bouts, so which discretization defect dominates depends on the
  setting (a plausible, unmeasured reason: log-spaced edges follow `max(t)`,
  which moves little across synthetic draws and a lot between a real sample
  and a null draw). Simulate through the engine's own survival convention --
  `sampleFromSurvival`-style inversion of `SurvivalHandle` -- not `round()`.
- **Never floor a probability without its normalizer.** Flooring a bin
  probability at `realmin` while `S(xmin)` sat in the subnormals made `p/S`
  reach 4.5e15 — a conditional probability above 1, worth +36 log-likelihood
  per observation. A Weibull scored `logL = +39972` and won a comparison.
  Engine test 23 asserts `logL < 0`.
- **A small `S(xmin)` is not a broken fit.** For a family anchored at 0 it only
  means the untruncated interpretation is vacuous. A gamma at `shape -> 0` has
  `S(xmin) ~ 3e-12` and is sound: `Gamma(a,z) -> E1(z)`, so the truncated law
  tends to the power-law-with-cutoff limit. Gating on it excluded good fits.
  Only a **free location** produces a real ridge — see `pearson3Guard`.
- **Don't write tests that assert where the optimizer lands.** These
  likelihoods are non-concave; MATLAB and Octave find different legitimate
  optima. Assert the decision *rule* instead (engine test 24).
- **Some rows are the same model twice.** `powerlaw_cutoff` IS `gamma`
  (`alpha = 1-shape`, bit-identical logL); `erlang` at shape 1 IS `hyperexp K=1`
  IS the exponential. `R.Nesting` lists them. Equal log-likelihoods there are an
  implementation check, not a coincidence — disagreement beyond ~1e-9 is a bug.
- **An absolute optimizer tolerance is a units bug waiting to happen -- but
  check where it actually bites before calling it one.** A fit that is not
  scale-equivariant gives different answers for the same bouts scored in
  seconds and in minutes. Two mechanisms: an absolute `TolX` on a parameter
  carrying the data's units, and `fzero`, whose `TolX` it scales BY `|x|`, so
  at `|mu| ~ 1e4` the effective tolerance is ~1e4 and it returns a bracket end
  rather than a root. Measured here, on our own code:
    - **The bout-duration suite uses no root-finder at all** -- `fzero`,
      `fminbnd` and `fsolve` appear only in the two noise fitters -- so the
      `fzero` mechanism cannot reach the pipeline.
    - **The engine is structurally immune**, not lucky: every positive
      parameter is optimized in log space (`M.unpack = @(z) exp(z)`), so
      rescaling time shifts `z` by a constant and an absolute `TolX = 1e-10`
      on `z` is a RELATIVE tolerance on the rate. In discrete mode, scaling
      `x`, `xmin` and `dt` together leaves `n_min` and every bin count
      identical, so the likelihood is invariant by construction.
    - **`shiftlognormal_MLE` was measured, not assumed**: fitting the same
      sample at `c` from 1e-6 to 1e6, the recovered `shift/c`, `exp(mu)/c` and
      `sigma` agree to <= 6e-8 relative for `c >= 1e-2`, and the worst case is
      1e-4 at `c = 1e-6`. Its `fminbnd` refinement does use an absolute
      `TolX = 1e-8`, which is the mechanism -- but it only bites once the data
      span shrinks to ~1e-5, which no duration in seconds, minutes, hours or
      days reaches. Its boundary diagnostic is already relative
      (`max(10*TolX, 1e-6*intervalWidth)`). **That sweep missed OFFSETS**
      (Nov 2026 audit): `fminbnd` also has a RELATIVE tolerance,
      `sqrt(eps)*|shift|`, so data carrying a large constant (`x + 1e8`)
      moved the shift by 2.3% and LogLike by 0.049. Fixed by searching over
      `log(min(x) - shift)`; scale AND offset are now equivariant to ~1e-8.
  So the trap is real, the exposure here is not, and the difference took one
  grep and one 12-decade sweep to establish rather than an argument.
- **A profile-likelihood bound from a finite search errs in ONE direction,
  and that fixes which side a disagreement indicts.** The profile is
  `p(theta) = min_eta nll(theta, eta)`. A finite nuisance search returns
  `phat >= p`, so the computed deviance `phat - nllMin` is too LARGE, crosses
  the threshold too EARLY, and the reported bound is too NARROW. To use that,
  evaluate an independent profile AT the estimator's bound, where the
  estimator's own deviance equals the threshold by construction:
    - independent value **below** threshold: the checker found a better
      nuisance fit, so the true bound lies further out and **the ESTIMATOR was
      too narrow**. This was the audit's false finite `beta` bounds, on ~10% of
      NIG n=400 fits.
    - independent value **above** threshold: both numbers are upper bounds on
      the same truth, so the higher one is the worse fit and **the CHECKER
      failed**. This was GH C's `alpha` lower bound, +1.26 above threshold.
  In terms of intervals: a checker reporting a WIDER bound indicts the
  ESTIMATOR, a NARROWER one indicts the CHECKER. An earlier version of this
  entry had both backwards, and the mistake is instructive -- "the estimator
  can only err narrow" invites the conclusion that a wider independent bound
  must be wrong, when a wider bound is precisely the signature of that error.
  **The asymmetry needs two conditions, and both have failed here:**
    - **`nllMin` must be the true optimum.** Too high an `nllMin` understates
      every deviance, so the bound crosses LATE and comes out too WIDE --
      the one-sidedness reverses. On the Normal ridge the point estimate was
      measured 0.002-0.005 nats short (nll 1041.9695 against 1041.9672 in the
      final version), a small WIDENING effect. **Do not confuse it with GH's
      surviving ridge limitation, which points the other way**: the measured
      artifact there is a bound too NARROW -- `mu` lower `-9244` on a Normal
      n=500 NIG fit, where the independent profile was still 1.07 BELOW
      threshold -- caused by nuisance fits stalling on the flat ridge
      (interior-point `exitflag 2`), which OVERstates the profile in the
      ordinary direction and is about 500x larger than the `nllMin` effect.
      Both conditions are live on that ridge at once, in opposite directions,
      and the stalls are what survive.
    - **The checker must carry the SAME constraints.** An unconstrained
      checker can reach below `delta`'s floor, and then a wider bound reflects
      a larger feasible set rather than an error in the estimator.
- **`chisquared` is not unit-invariant** (its scale is fixed at 2), so its fit
  depends on whether you scored in seconds or minutes. Exclude it rather than
  interpret its defeat.

## Testing

```matlab
test_FitTruncatedDiscreteMLE      % 28
test_FitHyperexponentialMLE       % 26
test_FitExponentiatedWeibullMLE   % 23
test_CompareBoutModels            % 29
test_HyperexponentialLRT          % 13
test_BoutHazard                   % 10
```

**In MATLAB, run them by name from this directory.** Check `which
test_CompareBoutModels -all` first: a stale copy elsewhere on the path wins
over the repo unless you are `cd`'d into it, and that has already caused a
"pass" that ran 24 of 27 tests.

**Anything that moves `pwd` away from this repo defeats that check**, and
`run('some\path\script.m')` does exactly that: MATLAB changes into the
script's folder for the duration. The mechanism is NOT that the script sits
beside stale copies — the scratch folder may contain no fitters at all. It is
that leaving the repo removes the current-folder precedence that was making the
repo win, so resolution falls through to the MATLAB **path**, and whatever
other copy sits on the path takes over. `which -all` run from the repo
therefore reports the repo copy winning while `run()` quietly uses another one.
This has already produced a discarded 35-minute timing measured against a
pre-fix engine.

Either `addpath` the scratch folder and call the script by NAME, or open with

    assert(isequal(fileparts(which('CompareBoutModels')), pwd))

so a run against the wrong copy fails loudly instead of returning numbers.

**Do not "clean up" the other copies.** The folder that shadows this one is the
user's main analysis repo, not a pile of stale files: its
`FitHyperexponentialMLE.m` is tracked there and called by scripts in it, and
its copies of the newer fitters are UNTRACKED, so deleting them cannot be
undone. Fix the resolution, not the files — take that folder off the path while
working here, or refresh its copy of whatever you need.

**That folder now carries `Icode_sync_stamp.txt`**, naming the commit its
copies were synced from (`1397554`, 2026-09-14) with a SHA-256 prefix per file
and the backup location of what was replaced. A session working from there
should compare that file's Commit line against `git log -1` in Icode BEFORE
trusting the copies. And note the asymmetry: `addpath('...Icode','-begin')`
does NOT protect you, because MATLAB's current folder beats the path — the
user's analysis scripts live in that folder, so its copies win whenever it is
`pwd`. The only durable options are to re-sync after any `.m` change, or to
keep one copy of each function.

**In Octave they do NOT run by name.** Three of the suites call the bare
`Fit*MLE` names, which resolve to the MATLAB originals and die on the
`arguments` block. They need a generated test twin — see `verify_all.sh` in the
scratchpad. `test_CompareBoutModels` dispatches to `_oct` names itself and does
run directly.

**The BOUT SUITE needs no toolbox** — local Marsaglia-Tsang gamma sampler, and
chi-square tails via `gammainc(...,'upper')` rather than `chi2cdf`. Engine test
20 cross-checks against `gamcdf`/`wblcdf` and self-skips when absent, so it only
runs in MATLAB. (Verified: `chi2cdf`, `chi2inv`, `fmincon` and `ksdensity` appear
in the bout files only inside comments explaining what is avoided, never as
calls.) **That claim does NOT extend to the rest of the repo**:
five of the non-bout tools need toolboxes, each together with its own stress
suite (verified on `origin/main` by grepping for calls outside comments):

| tool | needs |
| --- | --- |
| `GeneralizedHyperbolic_MLE` + `stress_test_GH_MLE` | **Optimization** (`fmincon`) and **Statistics** (`chi2cdf`, `chi2inv`) |
| `GenNorm_MLE` + `test_GenNorm_MLE` | **Optimization** (`fmincon`) and **Statistics** (`kmeans`, `iqr`, `quantile`, `chi2cdf`; the test calls `iqr` itself). `shiftlognormal_MLE` and its test need no toolbox |
| `multigaussian_MLE` + `test_multigaussian_MLE` | **Optimization** (`fmincon`) only; no Statistics calls |
| `AsymLaplace_MLE`, `shiftweibull_MLE` + their tests | None, except **Statistics** (`prctile`) in `BOOTSTRAP` mode |
| `jsd_kde` + `stressTest_jsd_kde` | **Statistics** (`ksdensity`, `skewness`, `iqr`, `quantile`; the suite calls `skewness` itself) |
| `chi2p` + `stressTest_chi2p` | **Statistics** (`chi2cdf` line 271, `chi2inv` line 347). Also `parfor` at line 309, which runs SERIALLY without the Parallel Computing Toolbox -- correct results, just slower |

The two stress suites for `chi2p` and GH make no toolbox calls of their own;
they need them through the functions they test. So a machine that runs all 129
bout tests can still fail on ten of these files.
