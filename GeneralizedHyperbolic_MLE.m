function results = GeneralizedHyperbolic_MLE(data, modelFlag, mode)
% GENERALIZEDHYPERBOLIC_MLE  Maximum-likelihood estimation for the
% Generalized Hyperbolic (GH) family of distributions.
%
% SYNTAX
%   results = GeneralizedHyperbolic_MLE(data)
%   results = GeneralizedHyperbolic_MLE(data, modelFlag)
%   results = GeneralizedHyperbolic_MLE(data, modelFlag, mode)
%
% INPUTS
%   data      : real finite numeric vector of at least 10 observations
%               (n >= 10; error TooFewObservations below that, and a
%               SmallSample warning below 50), with at least two distinct
%               values (converted to double).
%   modelFlag : 'GH'   – fit the full 5-parameter Generalized Hyperbolic.
%               'NIG'  – fit the 4-parameter Normal-Inverse Gaussian
%                        (GH with lambda fixed at -0.5).
%               'HYP'  – fit the 4-parameter hyperbolic distribution (GH with
%                        lambda fixed at 1); 'HYPERBOLIC' is accepted too.
%               'AUTO' – (default) fit all three, take the better of NIG and
%                        hyperbolic (same number of parameters, so by
%                        likelihood, a VALID fit preferred), then choose
%                        between it and GH by a likelihood-ratio test AND
%                        BIC, again preferring a VALID fit over an invalid
%                        one whatever the criteria say. Profile CIs are
%                        computed only for the selected model.
%   mode      : text containing 'PLOT' displays a histogram with the fitted
%               density; text containing 'NOPROFILE' skips the profile CIs
%               (the dominant cost; 10-30x faster) and reports SEs and Wald
%               CIs only. Both may be combined, e.g. 'NoProfile,Plot'.
%
% PHYSICAL PARAMETERS
%   The GH distribution is parameterized by five physical quantities:
%     mu     : location (real)
%     lambda : tail index (real; controls the shape of the tails)
%     alpha  : tail heaviness (> |beta|)
%     beta   : skewness PARAMETER (real). It fixes the SIGN of the skewness
%              but not its size, and it has units of 1/x, so it cannot be
%              compared across samples of different scale. For a comparable
%              number use results.Skewness or results.SkewnessBowleyData (see
%              SKEWNESS below).
%     delta  : scale (> 0)
%   The NIG model fixes lambda = -0.5; the hyperbolic model fixes lambda = 1.
%
%   Special cases of the GH family:
%     lambda = -0.5              : Normal-Inverse Gaussian (NIG)
%     lambda =  1                : Hyperbolic, log f = -alpha*sqrt(delta^2+(x-mu)^2)
%                                  + beta*(x-mu) + const
%     delta -> 0 with lambda > 0 : Variance-Gamma limit. This model is NOT
%                                  inside the fitted family (delta has a floor):
%                                  VG-like data can drive delta to the floor,
%                                  and delta's lower confidence bound is then
%                                  open (reported as 0).
%     alpha -> infinity          : approaches the Normal distribution
%     gamma -> 0 with lambda < 0 : skewed Student-t limit, nu = -2*lambda
%                                  (alpha -> |beta|). The likelihood goes FLAT
%                                  in psi there (logL identical from gamma =
%                                  1e-8 to 1e-40), so alpha and gamma are not
%                                  identified; flagged by student_t_limit.
%                                  NOT every heavy-tailed sample lands there:
%                                  t(3) samples of n = 1000 fitted an interior
%                                  gamma (Nov 2026).
%   (Earlier versions of this header listed lambda = 0.5 as the hyperbolic and
%   lambda -> 0 as the variance-gamma limit; both were wrong.)
%
% INTERNAL OPTIMIZATION COORDINATES
%   To enforce positivity and ordering constraints without explicit bounds,
%   estimation is performed in a transformed coordinate system:
%
%     GH  theta = [lambda; beta; psi; log_delta; mu]   (5-vector)
%     NIG and hyperbolic theta = [beta; psi; log_delta; mu]   (4-vector)
%
%   where:
%     psi       = log(gamma),  gamma = sqrt(alpha^2 - beta^2) > 0
%     log_delta = log(delta)
%
%   Physical parameters are recovered as:
%     gamma = exp(psi)                          (always > 0)
%     alpha = hypot(beta, gamma)                (always > |beta|)
%     delta = exp(log_delta)                    (always > 0)
%
%   Advantage over the older phi = log(alpha - |beta|) parameterization:
%   the mapping beta -> alpha = hypot(beta, gamma) is smooth and
%   differentiable everywhere, including at beta = 0.  The older form
%   had a kink there that could stall quasi-Newton optimizers.
%
% ESTIMATION STRATEGY
%   0. The data are standardized, z = (data - median)/std, the model is fitted
%      to z, and every output is mapped back exactly (the GH family is closed
%      under location-scale maps, and logL changes by exactly -n*log(std)).
%      The fit, its intervals and every tolerance are therefore independent
%      of the units and location of the data.
%   1. Grid search over a coarse (lambda, beta, delta, mu) lattice to
%      identify promising starting regions.
%   2. fmincon (interior-point) is launched from the best 8 grid points.
%      A lower bound is applied to log_delta, because delta -> 0 makes the
%      density spike at mu and, on tied or discretized data, that spike can
%      drive the likelihood up without limit.  The same bound is applied to
%      every nuisance optimization in the profile CIs, so the fit and its
%      intervals describe one feasible region.  An active bound is reported
%      via results.delta_at_bound and a DeltaAtBound warning; such a point is
%      a constrained optimum, not an interior MLE.
%   3. For GH: two more warm starts are always attempted, from the polished
%      NIG and hyperbolic point estimates (lambda freed from -0.5 and from
%      1). Every model standardizes identically, so their thetas are reused
%      in z-units. A warm-start solution replaces the current one if it
%      converged when that did not, or if both converged and it is better.
%      The warm-start POINTS are candidates too, so the GH likelihood can
%      never come out below the NIG or hyperbolic one (GH nests both).
%   4. The best converged solution is selected; if none converged, the
%      best finite solution is used with a warning.
%   5. The selected solution is polished by a derivative-free (Nelder-Mead)
%      search within the same bounds, and replaced if that lowers the NLL.
%      Near the Normal limit the likelihood is a long flat ridge on which the
%      interior-point solver stops early (exitflag 2, step below tolerance);
%      measured there (Sep 2026), it ended 0.005 nats short of the optimum.
%      The polish is accepted only if its gain exceeds the NLL's own
%      floating-point noise (about n*max Bessel argument*eps).
%   6. Steps 0-5 are the POINT stage: they give LogLike, BIC and valid, which
%      is all AUTO needs to choose. Profile CIs (the CI stage, the dominant
%      cost) then run only for the model that is returned.
%   7. The CI stage tracks the lowest NLL found by ANY profile nuisance fit.
%      If it beats the point estimate by more than 1e-3 nats (or the NLL's
%      noise), the point estimate was not the optimum: the fit is
%      re-optimized from that theta (fmincon, then the polish of step 5) and
%      the CIs are recomputed, at most twice (results.profile_refits); after
%      that a ProfileImprovedOptimum warning is raised. The one-sidedness of
%      profile bounds (too narrow, never too wide) holds only when the point
%      estimate is the true optimum.
%
% LOG-LIKELIHOOD
%   The GH log-density is evaluated using the scaled Bessel function
%   besselk(nu, x, 1) = besselk(nu, x) * exp(x).  Scaling cancels the
%   dominant exponential growth and keeps intermediate values finite for
%   large arguments.  The log-likelihood is -sum(logpdf) over observations.
%
%   NUMERICAL-VALIDITY GUARD
%   The log-density forms (log K - arg) for both Bessel terms, so it adds and
%   subtracts quantities of order ARG.  In double precision that cancellation
%   carries an absolute error of about n*arg*eps in the summed NLL.  Parameter
%   vectors whose Bessel arguments exceed min(1e10, 1e-2/(n*eps)) are
%   therefore declared infeasible (NLL = Inf), which caps that error at
%   1e-2 nats: beyond it the error swamps genuine differences in likelihood
%   and an optimizer will "improve" the fit by chasing floating-point noise.
%   The cap is 1e10 for n <= ~4500. (It used to be a fixed 1e10, justified as
%   "below ~1e-2 even for n = 1e6"; in fact n*1e10*eps is 0.22 nats at
%   n = 1e5 and 2.2 at n = 1e6, above the 1.92 CI threshold.) Bessel
%   arguments are scale-invariant, so the cap means the same in any units.
%
%   This is not hypothetical.  Without the guard, NIG fits to n >= 1000 ran
%   away to alpha ~ 4e52 and delta ~ 1e36, reporting a log-likelihood of
%   +4.5e19 -- impossible for any density -- while results.Fit evaluated to
%   all NaN and results.converged was true.  Legitimate fits sit many orders
%   of magnitude below the cap.
%
% PROFILE LIKELIHOOD CONFIDENCE INTERVALS
%   95% CIs are computed by profile likelihood rather than the Wald
%   (Hessian-inversion) approach, because:
%     - GH parameters are strongly coupled, making Wald intervals
%       unreliable and potentially non-symmetric.
%     - Profile CIs remain valid even when the Hessian is ill-conditioned.
%
%   For each parameter p, the profile CI solves:
%     2 * (nll_profile(p) - nll_min) = chi2inv(0.95, 1)
%   where nll_profile(p) is the minimum NLL with p fixed and all other
%   parameters re-optimized (profiled out under the same box constraints).
%
%   Parameter-specific profiling strategy:
%     lambda, beta, mu : profiled directly in their (untransformed) coordinates.
%     delta            : profiled in log_delta space, then back-transformed.
%     alpha            : profiled in log(alpha) space with alpha fixed
%                        physically.  For each trial alpha = a, the
%                        constraint |beta| < a is enforced by substituting
%                        beta = a * tanh(xi) (xi free), which also gives
%                        gamma = a / cosh(xi) automatically.  This is an
%                        exact profile, not an approximation.
%
%   Root finding uses an exponentially expanding bracket search sized from the
%   Wald standard errors, followed by fzero on the located bracket, with
%   sign-only bisection as the fallback. A non-finite profile value (besselk
%   overflow, a failed nuisance fit) is never taken as a crossing; the search
%   backs off toward the last finite point instead. A crossing whose refinement
%   meets non-finite, mutually contradictory or discontinuous (cliff) profile
%   values is rejected as an artifact of numerical breakdown. At delta's floor,
%   where one value decides between an open and a finite bound, a value not
%   certified below the threshold is re-fitted from every nearby cached
%   solution with a larger Nelder-Mead budget.
%   log(delta) is searched no lower than the estimator's own delta floor.
%   Each profile point warm-starts its nuisance fit from the nearest point
%   already evaluated. Any value at or past the threshold is re-checked from
%   the MLE's nuisance values and, if still at or past it, by a derivative-free
%   Nelder-Mead search, before it is allowed to set a bound.
%
% STANDARD ERRORS AND WALD INTERVALS (added Nov 2026)
%   results.SE (and Cov) are delta-method standard errors of the physical
%   parameters, in the order of Params, from results.information: the
%   OBSERVED information, i.e. the Hessian of the NLL at the optimum in the
%   transformed coordinates of theta_names, taken by central differences with
%   Richardson extrapolation on the exact likelihood (lambda, the Bessel
%   order, has no closed-form derivative). results.CIWald are the matching
%   95% Wald intervals: symmetric for mu, lambda and beta, on the log scale
%   for alpha and delta (so they stay positive). lambda's SE is 0 when it is
%   fixed (NIG, hyperbolic).
%   - At the delta floor (delta_at_bound) no Wald SE is reported (NaN,
%     SEValid = false): conditional-on-delta intervals were measured to cover
%     0 of 12 times there. Use the profile intervals.
%   - A singular or indefinite information matrix (flat ridges: the Normal
%     limit, the Student-t limit) gives NaN SEs and SEValid = false.
%   - The profile intervals are the better summary whenever they differ from
%     the Wald ones: the GH likelihood is far from quadratic in alpha, beta and
%     lambda at ordinary sample sizes. See MEASURED COVERAGE.
%
% SKEWNESS (added Nov 2026)
%   Three dimensionless skewness summaries. All three are unchanged by a shift
%   or a positive rescaling of the data and change sign when the data are
%   negated, so they can be compared across samples, units and flies; beta
%   cannot (see PHYSICAL PARAMETERS).
%
%   results.Skewness - moment skewness E[(X - EX)^3]/SD^3 of the FITTED
%     distribution, in closed form. Writing GH as the normal variance-mean
%     mixture X = mu + beta*V + sqrt(V)*Z, with V ~ GIG(lambda, delta, gamma)
%     and k1, k2, k3 the mean, variance and third central moment of V:
%         Var X    = k1 + beta^2*k2
%         Skewness = (3*beta*k2 + beta^3*k3) / (Var X)^(3/2)
%     with E[V^r] = (delta/gamma)^r * K_{lambda+r}(delta*gamma) /
%     K_lambda(delta*gamma). k2 and k3 are positive, so the sign is beta's.
%     For NIG this reduces to 3*beta/(alpha*sqrt(delta*gamma)).
%     NaN on the flagged Student-t limit (student_t_limit), where gamma is not
%     identified and the moment skewness grows without bound as gamma -> 0.
%     It is a property of the fitted model: heavy tails make it sensitive to
%     gamma, so a fit with a small gamma can report a large value with a
%     large SkewnessSE.
%   results.SkewnessSE, SkewnessCIWald - delta-method SE from the observed
%     information, and the matching 95% Wald interval. NaN whenever SEValid
%     is false (e.g. delta at its floor). See MEASURED COVERAGE.
%   results.SkewnessBowleyData - Bowley's quartile skewness of the DATA,
%     (Q3 + Q1 - 2*Q2)/(Q3 - Q1), in [-1, 1]. Model-free and robust to the
%     tails. Quartiles are MID-QUANTILES (Parzen 2004): the mid-distribution
%     function Fmid(u) = (#{x < u} + #{x = u}/2)/n at each distinct value,
%     interpolated linearly. Without ties this is exactly quantile()'s
%     default. With ties (integer-minute data) ordinary quartiles can only take
%     lattice values, and the statistic jumps in coarse steps: measured (Nov
%     2026) on NIG data scaled to a few units and rounded to integers, plain
%     quartiles gave 0.333 against 0.067 for the unrounded data; the
%     mid-quantiles gave 0.056.
%   results.SkewnessBowleyModel - Bowley's skewness of the FITTED
%     distribution, from its quartiles (CDF integrated numerically). Exists for
%     every fit, Student-t limit included. Far from SkewnessBowleyData, it
%     says the fit misplaces the body of the data, whatever the likelihood.
%
%   MEASURED (Nov 2026, 'NoProfile' fits to exact GIG-mixture samples, mu = 0,
%   alpha = 2, beta = 0.8, delta = 1; true values from the closed form and a
%   4e6-point CDF):
%                      true   SkewnessCIWald   sd of est.  median    Bowley sd
%                      skew   coverage (valid) Skewness    SE        model  data
%   NIG n=150  (200)   0.886  0.923 (195)      0.40        0.41      0.034  0.10
%   NIG n=500  (200)   0.886  0.965 (200)      0.20        0.21      0.017  0.059
%   HYP n=500  (200)   0.892  0.940 (200)      0.16        0.16      0.023  0.059
%   GH  n=1000 (100)   0.749  0.963 (82)       0.26        0.21      0.016  0.038
%     (GH lambda = -1.5; the 18 without an SE are delta-floor fits or the
%     Student-t limit.) Both Bowley estimates were within about 0.01 of the
%     true value on average (true 0.081, 0.081, 0.104, 0.058); the data
%     version is 2.5-3.5 times as variable as the model version, the price of
%     being model-free. Wald is adequate from n ~ 500; at n = 150 it
%     undercovers slightly.
%
%   INDEPENDENCE. Every SE and interval in this function assumes independent
%   observations. Successive differences of a sequence are not independent
%   (adjacent differences share a term; for independent durations their
%   correlation is -0.5), so for such data these SEs are too small, and
%   uncertainty should come from resampling that respects the order (e.g. a
%   shuffle or block null computed by the caller).
%
% MEASURED COVERAGE (Nov 2026, nominal 95%, 100 simulated samples each from
%   mu = 0, alpha = 2, beta = 0.5, delta = 1; exact GH variates from the GIG
%   mixture; every fit's LogLike was at or above the log-likelihood at the
%   true parameters, so no fit stopped short of the optimum)
%                       mu      lambda  alpha   beta    delta
%   NIG n=500  profile  0.96    -       0.96    0.95    0.96
%              Wald     0.94    -       0.98    0.95    0.97
%   HYP n=500  profile  0.94    -       0.90    0.97    0.93
%              Wald     0.96    -       0.96    0.97    0.97
%   GH  n=1000 profile  0.95    0.96    0.98    0.96    0.95   (lambda -1.5)
%              Wald     0.96    0.80    0.79    0.94    0.84
%   GH, split: 88 interior fits, Wald 0.97/0.91/0.90/0.96/0.96 (profile
%   0.96-1.00); 12 fits on the delta floor -- each confirmed to be the global
%   optimum, the data's own variance-gamma corner (lambda 1.9 to 5.4) -- where
%   the profile covered 0.75-0.92 and Wald, then still conditional on delta,
%   covered lambda, alpha and delta 0 of 12 (now reported NaN). Binomial error
%   on each entry is about +/-0.02. The hyperbolic alpha profile at 0.90 was
%   checked replicate by replicate: each missed bound sits on the threshold
%   to 1e-3 under an independent profile, so the shortfall is the chi-square
%   approximation at n = 500, not the computation.
%
% results.hessian (kept for compatibility) is NOT the information: it is the
%   matrix fmincon returned, a quasi-Newton (BFGS) approximation of the
%   Hessian of the LAGRANGIAN, accumulated along the optimizer's path, mapped
%   to data units. Use it, if at all, only to diagnose ill-conditioning. SEs
%   used to be left to the caller to derive from it; they would have been
%   wrong.
%
% OUTPUTS
%   results : structure with the following fields:
%     Model       – descriptive string identifying the fitted model.
%     Params      – [mu; lambda; alpha; beta; delta] physical parameters.
%     LogLike     – maximized log-likelihood.
%     Fit         – fitted density evaluated at each DATA point (original order).
%     FitX        – sorted copy of DATA (for smooth density plots).
%     FitPDF      – fitted density evaluated at FitX.
%     AIC         – Akaike Information Criterion  (2k - 2*LogLike).
%     BIC         – Bayesian Information Criterion (k*log(n) - 2*LogLike).
%     theta_hat   – MLE in transformed coordinates.
%     theta_names – cell array naming each element of theta_hat.
%     nll         – minimized negative log-likelihood.
%     hessian     – fmincon's quasi-Newton Hessian (NOT the information; see
%                   STANDARD ERRORS AND WALD INTERVALS).
%     information – observed information (NLL Hessian) in the coordinates of
%                   theta_names, data units; the basis of SE.
%     SE, Cov     – delta-method standard errors and covariance of Params
%                   (order [mu; lambda; alpha; beta; delta]).
%     SEValid     – false when the information could not be inverted (SE NaN).
%     CIWald      – 95% Wald intervals, same layout as profile_ci.
%     CILevel     – 0.95.
%     Skewness    – moment skewness of the fitted distribution (see SKEWNESS).
%     SkewnessSE, SkewnessCIWald – its delta-method SE and 95% Wald interval.
%     SkewnessBowleyData  – Bowley quartile skewness of the data (mid-quantiles).
%     SkewnessBowleyModel – Bowley quartile skewness of the fitted distribution.
%     CIMethod    – 'profile', or 'wald' when 'NoProfile' skipped the
%                   profiles (profile_ci is then NaN, fixed lambda excepted).
%     hessian_coordinates – same as theta_names, documents what hessian rows/cols mean.
%     profile_ci  – struct with 95% profile CI [lo hi] for each physical
%                   parameter: lambda, beta, alpha, delta, mu.
%
%                   A bound is one of three things, and they mean different
%                   things:
%                     finite  – the profile crossed the threshold there.
%                     +/-Inf  – the profile was evaluated all the way out (for
%                               delta: down to its floor) and never crossed:
%                               the likelihood is flat in that
%                               direction, so the parameter is not identified
%                               and the interval is genuinely OPEN.  This is a
%                               result, not a failure.  For delta and alpha,
%                               which are back-transformed through exp(), an
%                               open bound appears as 0 or Inf.
%                     NaN     – the profile could not be evaluated, so the
%                               bound is unknown.  This IS a failure.
%                   Open and unknown were previously both reported as NaN and
%                   so could not be told apart.
%     exitflag    – fmincon exit flag for the selected solution.
%     output      – fmincon output struct for the selected solution (of the
%                   standardized problem).
%     converged   – logical; true when exitflag > 0.  NOTE this reports only
%                   what the optimizer said; see 'valid' for whether the
%                   answer is usable.
%     LogFit      – fitted LOG-density at each data point. Fit = exp(LogFit)
%                   can underflow to 0 for a gross outlier (log-density
%                   below about -745) without anything being wrong.
%     valid       – logical; true when the fitted parameters yield a finite
%                   LOG-density at every data point AND that log-density
%                   sums to the reported LogLike (to within its
%                   floating-point noise). Judged on LogFit, not Fit, so an
%                   underflowed density does not count against it.  A run can be
%                   converged = true and valid = false: the optimizer stopped
%                   cleanly at parameters where the density cannot be
%                   evaluated.  Check this before using any result.
%     gaussian_check – struct diagnosing the Normal limit, where GH shape
%                   parameters stop being identifiable:
%                     zeta           = delta*gamma; large values mean the
%                                      density is close to Gaussian.
%                     LogLike_normal, AIC_normal for a fitted Normal.
%                     beats_normal   = AIC < AIC_normal.
%                   When beats_normal is false the fitted density may be
%                   perfectly good while alpha, beta and delta are essentially
%                   arbitrary; a NoGainOverNormal warning is raised.  This is
%                   the usual explanation for a fit that looks wrong even
%                   though every numerical check passes. In this regime even
%                   FINITE profile bounds can be artifacts of the flat ridge:
%                   measured (Sep 2026) on normal data, a mu lower bound of
%                   -9244 sat where an independent profile was still 1.07
%                   below the threshold.
%     delta_at_bound – logical; true when delta sits on (within 1% of) its
%                   lower bound, i.e. the solution is a CONSTRAINED
%                   optimum.  Profile CIs and
%                   Hessian-based standard errors assume an interior solution
%                   and are unreliable in that case.
%     student_t_limit – logical; true when lambda < 0 and pushing gamma down
%                   by a factor exp(-20) at fixed (lambda, beta, delta, mu)
%                   costs less than 1e-3 nats: the fit sits on the flat
%                   skewed Student-t limit (nu = -2*lambda), where alpha and
%                   gamma are not identified. Report the t limit, not alpha.
%     profile_refits – number of times (0-2) the CI stage found a better
%                   optimum than the point estimate and re-fitted from it.
%                   Nonzero means the first point estimate was short of the
%                   optimum; LogLike, Params and CIs are the refitted ones.
%     ModelSelection – (AUTO only) string describing which model was chosen
%                   and why, including when validity overrode the criteria
%                   or the selected model was re-optimized during profiling.
%     LRT         – (AUTO only) struct: LRstat, df, pValue for GH against
%                   the better fixed-lambda model (named in Restricted);
%                   LogLike_GH, LogLike_NIG, LogLike_HYP; and pValue_NIG,
%                   pValue_HYP, GH against each. All on the final
%                   (post-refit) log-likelihoods.
%                   The model NOT selected has no profile CIs computed.
%
% WARNINGS RAISED
%   Each is raised for the RETURNED result only: in AUTO, the model that is
%   discarded is fitted silently.
%   GeneralizedHyperbolic_MLE:SmallSample
%       Fewer than 50 observations; shape parameters weakly identified.
%   GeneralizedHyperbolic_MLE:DensityNotEvaluable
%       Fitted parameters give a non-finite log-density at some data point.
%   GeneralizedHyperbolic_MLE:LikelihoodInconsistent
%       Reported LogLike disagrees with the density it implies.
%   GeneralizedHyperbolic_MLE:DeltaAtBound
%       delta reached its lower bound; constrained optimum.
%   GeneralizedHyperbolic_MLE:NoGainOverNormal
%       Fit does not beat a Normal by AIC; shape parameters unidentifiable.
%   GeneralizedHyperbolic_MLE:StudentTLimit
%       The fit sits on the flat Student-t limit; alpha and gamma not identified.
%   GeneralizedHyperbolic_MLE:ProfileImprovedOptimum
%       After 2 refits the profile search still beat the point estimate; the
%       fit is unconverged and its CIs may be too WIDE.
%   GeneralizedHyperbolic_MLE:BothModelsInvalid
%       (AUTO) neither fit is valid; the usual rule chose, untrustworthily.
%
% KNOWN LIMITATIONS (measured, Sep 2026)
%   - Near the Normal limit (a NoGainOverNormal warning), finite profile
%     bounds can still be artifacts of the flat likelihood ridge; see
%     gaussian_check below.
%   - Upper bounds for alpha and beta many orders of magnitude above the
%     estimate (e.g. 1.2e6 for an estimate of 8) can mark where the nuisance
%     fits stop converging rather than a true crossing: an independent
%     profile was still below the threshold at two such bounds. Read them as
%     effectively open.
%     Since Nov 2026 the commonest case, a crossing made by the Bessel-argument
%     guard (bounds of a few e8 at n = 500), is reported as Inf itself; a far
%     bound well short of the guard still needs this reading.
%   - The profile CIs cost more than before the Sep 2026 fixes. On GH samples
%     of n = 1000 the median fit took about 2.3 times as long (171 s versus
%     75 s, both measured under concurrent load), mostly from searching delta
%     down to its floor and re-checking values near the threshold.
%   - AUTO is cheaper than before (CIs for the selected model only, and the
%     GH fit warm-starts from the polished NIG fit): on NIG data of n = 600
%     it took 17 s, against 81 s for stand-alone GH (64 s) plus NIG (17 s)
%     (measured Nov 2026 in two suite runs: 17.3 and 17.1 s against 81.1 and 75.6 s; not under controlled load).
%   - student_t_limit rarely fires, even on Student-t data (measured Nov
%     2026). On t(3), t(1.5) and Cauchy samples of n = 1000, 11 of 12 GH and
%     NIG fits sat at an interior gamma of roughly 1/max|x| (1e-4 to 0.19)
%     and were not flagged; one t(3) GH fit (lambda = -1.55) was, at a cost
%     of 0.00085 nats against the 1e-3 tolerance. Damping the tail beyond the
%     largest observation costs almost nothing in-sample, so the optimum
%     usually keeps a small positive gamma, and whether the flag fires turns
%     on a fraction of a millinat. Read it as a rule check, not a detector of
%     t-distributed data; a strongly negative lambda with a small gamma is
%     the practical sign of a t-like sample.
%   - A profile bound is exact for the likelihood, but its 95% LEVEL rests on
%     the chi-square approximation, which is loose for the hyperbolic alpha
%     at n = 500 (coverage 0.90; see MEASURED COVERAGE). Checked independently
%     on every missed replicate: the bounds sit on the threshold to 1e-3.
%   - Along a flat ridge the profile can be unknown (NaN) rather than open:
%     one hyperbolic sample of n = 500 had mu = -1.95 (SE 3.5) and an
%     independent profile stayed below the threshold to mu = -6, crossing only
%     between -10 and -20 with beta ~ 1e5; the lower mu bound came back NaN.
%     A Wald SE many times the usual one is the warning sign.
%
% REQUIREMENTS
%   Optimization Toolbox  : fmincon
%   Statistics Toolbox    : chi2cdf, chi2inv
%
% REFERENCES
%   Barndorff-Nielsen, O.E. (1977). Exponentially decreasing distributions
%     for the logarithm of particle size. Proc. R. Soc. Lond. A, 353, 401-419.
%   Barndorff-Nielsen, O.E. (1997). Normal Inverse Gaussian distributions
%     and stochastic volatility modelling. Scand. J. Statist., 24, 1-13.
%   Bowley, A.L. (1901). Elements of Statistics. P.S. King, London.
%   Parzen, E. (2004). Quantile probability and statistical data modeling.
%     Statistical Science, 19, 652-662.
%
% Originally written by Sheyum but beefed up/stabilized/optimized 
% by Claude, Gemini and ChatGPT, Dec 2025.
% Reparameterization, exact alpha CI, and robustness improvements: Jul 2026.
% Updated to fmincon with bounds to prevent component collapse: Aug 2026.
% Audited against independent ground truth and fixed (Claude Code, Sep 2026):
% standardization; open bounds no longer read off numerical overflow; delta
% intervals confined to the delta floor; checks that reject crossings built
% from failed nuisance fits (non-finite, contradictory or cliff-like values);
% warm-started and derivative-free re-checked profile fits; an absolute root
% tolerance; interval bounds no longer swapped by sort(NaN); the
% delta_at_bound flag; double conversion of the input; a noise-aware validity
% check; and corrected special cases (hyperbolic, variance-gamma) in this
% header.
% Second audit, 8 findings fixed (Claude Code, 2026-11-03): AUTO prefers a
% valid fit; profile CIs check that the point estimate was the optimum and
% re-fit if not (profile_refits); AUTO computes CIs for the selected model
% only and GH reuses the polished NIG fit as its warm start; the Student-t
% limit is reported (student_t_limit); exact gamma = a/cosh(xi) in the alpha
% profile; an n-dependent Bessel-argument cap; n >= 10 with a SmallSample
% warning below 50; validity judged on the log-density (LogFit), with one
% shared log-density function.
% Third revision (Claude Code, 2026-11-04): the hyperbolic model ('HYP', lambda
% = 1), which AUTO now considers alongside NIG; GH warm-starts from both nested
% fits and can no longer return a likelihood below theirs; standard errors
% (SE, Cov) from a numerical observed information and Wald CIs (CIWald),
% replacing the unusable quasi-Newton Hessian; the 'NoProfile' mode; profile
% steps for every parameter, alpha included, sized from that information
% and capped at 1 z-unit;
% crossings made by the Bessel-argument guard reported open instead of as
% ~1e9 bounds; and coverage of every interval measured by simulation.
% Fourth revision (Claude Code, 2026-11-05): skewness outputs (Skewness with
% SE and Wald CI, SkewnessBowleyData, SkewnessBowleyModel; see SKEWNESS);
% beta described as the skewness parameter, not the skewness. Fix: Wald SEs
% are now inverted in z-units; inverting the data-unit information made
% SEValid depend on the units (data rescaled by 1e-6 returned all-NaN SEs).

% ---- Validate and normalize inputs -----------------------------------------
if nargin < 2 || isempty(modelFlag), modelFlag = 'AUTO'; end
if nargin < 3 || isempty(mode),      mode      = '';     end

flagMsg = ['MODELFLAG must be ''GH'', ''NIG'', ''HYP'' (or ''HYPERBOLIC''), ' ...
    'or ''AUTO''.'];
if ~(ischar(modelFlag) || (isstring(modelFlag) && isscalar(modelFlag)))
    error('GeneralizedHyperbolic_MLE:InvalidModelFlag', flagMsg);
end
modelFlag = upper(char(modelFlag));
if strcmp(modelFlag, 'HYPERBOLIC'), modelFlag = 'HYP'; end
if ~ismember(modelFlag, {'GH','NIG','HYP','AUTO'})
    error('GeneralizedHyperbolic_MLE:InvalidModelFlag', flagMsg);
end
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('GeneralizedHyperbolic_MLE:InvalidMode', ...
        'MODE must be a character vector or scalar string.');
end
mode = char(mode);
doProfile = ~contains(upper(mode), 'NOPROFILE');

validateattributes(data, {'numeric'}, ...
    {'vector','real','finite','nonempty'}, mfilename, 'data', 1);
% Work in double precision whatever the input class: single data would run the
% likelihood, the tolerances and besselk in single precision, and integer
% classes saturate on subtraction.
data = double(data(:));
% A 5-parameter model cannot be estimated from a handful of points: n = 2 used
% to be accepted. 10 is a hard floor (test 11F fits n = 10 near the Normal
% limit and must stay legal); below 50 the fit runs but warns.
if numel(data) < MIN_N_ERROR()
    error('GeneralizedHyperbolic_MLE:TooFewObservations', ...
        'DATA must contain at least %d observations (got %d).', ...
        MIN_N_ERROR(), numel(data));
end
if ~(std(data) > 0)
    error('GeneralizedHyperbolic_MLE:ConstantData', ...
        'DATA must contain at least two distinct values.');
end
if numel(data) < MIN_N_WARN()
    warning('GeneralizedHyperbolic_MLE:SmallSample', ...
        ['Only %d observations for a model with up to 5 parameters. The fit ' ...
        'runs, but shape parameters (lambda, alpha, beta) are weakly ' ...
        'identified below about %d observations; expect wide or open ' ...
        'intervals.'], numel(data), MIN_N_WARN());
end

% Check required toolbox functions are available.
requiredFunctions = {'fmincon','chi2cdf','chi2inv'};
for k = 1:numel(requiredFunctions)
    if isempty(which(requiredFunctions{k}))
        error('GeneralizedHyperbolic_MLE:MissingFunction', ...
            'Required function %s is unavailable.', requiredFunctions{k});
    end
end

% ---- Fit the requested model(s) --------------------------------------------
% Every fit has two stages: a POINT stage (multistart, polish, validity and
% boundary diagnostics, standard errors) and a CI stage (profile likelihood).
% The CI stage is the dominant cost, so AUTO runs it only for the model it
% selects, and 'NoProfile' skips it. The GH point stage warm-starts from the
% polished NIG and hyperbolic point estimates; every fit standardizes
% identically, so those z-unit thetas are reused as is. GH alone runs both
% nested point stages too, so a GH fit is the same whether it was asked for
% directly or reached through AUTO.
switch modelFlag
    case 'GH'
        fprintf('Fitting Generalized Hyperbolic (GH) model...\n');
        stNIG   = fit_point(data, -0.5);
        stHYP   = fit_point(data, 1);
        stGH    = fit_point(data, NaN, {warm_theta(stNIG), warm_theta(stHYP)});
        results = fit_ci(stGH, doProfile);

    case 'NIG'
        fprintf('Fitting Normal-Inverse Gaussian (NIG) model...\n');
        results = fit_ci(fit_point(data, -0.5), doProfile);

    case 'HYP'
        fprintf('Fitting hyperbolic (HYP, lambda = 1) model...\n');
        results = fit_ci(fit_point(data, 1), doProfile);

    case 'AUTO'
        fprintf('Running automatic model selection (GH versus NIG and hyperbolic)...\n');
        stNIG  = fit_point(data, -0.5);
        stHYP  = fit_point(data, 1);
        stGH   = fit_point(data, NaN, {warm_theta(stNIG), warm_theta(stHYP)});
        resGH  = stGH.results;
        resNIG = stNIG.results;
        resHYP = stHYP.results;

        % Step 1: the better of the two fixed-lambda models. Both have four
        % free parameters, so AIC, BIC and the likelihood rank them alike; a
        % VALID fit is preferred over an invalid one whatever the likelihood.
        if resHYP.valid ~= resNIG.valid
            useHYP = resHYP.valid;
        else
            useHYP = resHYP.LogLike > resNIG.LogLike;
        end
        if useHYP, stR = stHYP; nameR = 'HYP'; else, stR = stNIG; nameR = 'NIG'; end
        resR = stR.results;

        % Step 2: GH against that model. Likelihood-ratio test, H0 = lambda
        % fixed at the restricted model's value; GH has one extra free
        % parameter, so df = 1. Testing against the BETTER restricted model
        % gives the smaller of the two LR statistics, so the test is
        % conservative for "lambda is -0.5 or 1"; each separate test is in
        % results.LRT as well.
        LR   = max(2*(resGH.LogLike - resR.LogLike), 0);
        pval = chi2cdf(LR, 1, 'upper');  % upper tail; accurate for small p

        % Prefer GH only when the improvement is both statistically
        % significant (LRT) AND reflected in a lower BIC (penalizes the extra
        % parameter). A VALID fit is preferred over an invalid one whatever
        % the criteria say: an invalid fit's LogLike and BIC are not
        % trustworthy (see 'valid'), so comparing them decides nothing.
        ruleGH = pval < 0.05 && resGH.BIC < resR.BIC;
        if resGH.valid && ~resR.valid
            useGH = true;
            why = sprintf('the %s fit is invalid (valid = false), so the criteria were not used', nameR);
        elseif ~resGH.valid && resR.valid
            useGH = false;
            why = 'the GH fit is invalid (valid = false), so the criteria were not used';
        else
            useGH = ruleGH;
            why = '';
            if ~resGH.valid   % all invalid
                warning('GeneralizedHyperbolic_MLE:BothModelsInvalid', ...
                    ['Neither the GH nor the NIG or hyperbolic fit is valid; ' ...
                    'the usual LRT/BIC rule was applied, but the selected ' ...
                    'result is not trustworthy.']);
                why = 'no fit is valid; usual rule applied';
            end
        end

        if useGH
            results = fit_ci(stGH, doProfile);
            sel = 'GH';  bicSel = resGH.BIC;  bicOth = resR.BIC;  oth = nameR;
        else
            results = fit_ci(stR, doProfile);
            sel = nameR; bicSel = resR.BIC;   bicOth = resGH.BIC; oth = 'GH';
        end
        if isempty(why)
            results.ModelSelection = sprintf( ...
                '%s selected (LRT GH vs %s p=%.4g; BIC %.3f versus %s %.3f).', ...
                sel, nameR, pval, bicSel, oth, bicOth);
        else
            results.ModelSelection = sprintf( ...
                '%s selected: %s (LRT GH vs %s p=%.4g; BIC %.3f versus %s %.3f).', ...
                sel, why, nameR, pval, bicSel, oth, bicOth);
        end
        % The CI stage can re-optimize the selected model (profile_refits),
        % which raises its LogLike. Report the LRT on the final values; the
        % selection stands, since a higher LogLike for the selected model only
        % strengthens the rule that chose it.
        llGH = resGH.LogLike;  llNIG = resNIG.LogLike;  llHYP = resHYP.LogLike;
        switch sel
            case 'GH',  llGH  = results.LogLike;
            case 'NIG', llNIG = results.LogLike;
            case 'HYP', llHYP = results.LogLike;
        end
        llR = llNIG;  if useHYP, llR = llHYP; end
        if results.profile_refits > 0
            LR   = max(2*(llGH - llR), 0);
            pval = chi2cdf(LR, 1, 'upper');
            results.ModelSelection = sprintf( ...
                '%s The selected model was re-optimized %d time(s) during profiling; LRT p is now %.4g.', ...
                results.ModelSelection, results.profile_refits, pval);
        end
        results.LRT = struct('LRstat', LR, 'df', 1, 'pValue', pval, ...
            'Restricted', nameR, ...
            'LogLike_GH', llGH, 'LogLike_NIG', llNIG, 'LogLike_HYP', llHYP, ...
            'pValue_NIG', chi2cdf(max(2*(llGH - llNIG), 0), 1, 'upper'), ...
            'pValue_HYP', chi2cdf(max(2*(llGH - llHYP), 0), 1, 'upper'));
end

if contains(upper(mode), 'PLOT'), plot_results(data, results); end
end

function w = warm_theta(st)
% A fixed-lambda point estimate as a GH theta (lambda first), in z-units.
w = [st.fixLam; st.thetaZ(:)];
end

function tf = nfix(fixLam)
% 1 when lambda is fixed (NIG, hyperbolic), 0 for GH; doubles as the offset
% of every theta index, since a fixed-lambda theta has no lambda entry.
tf = double(~isnan(fixLam));
end

function s = model_name(fixLam)
% Descriptive model string for a fixed lambda (NaN = GH, lambda free).
if isnan(fixLam)
    s = 'Generalized Hyperbolic (mu, lambda, alpha, beta, delta)';
elseif fixLam == -0.5
    s = 'Normal-Inverse Gaussian (mu, lambda, alpha, beta, delta)';
elseif fixLam == 1
    s = 'Hyperbolic (mu, lambda = 1, alpha, beta, delta)';
else
    s = sprintf('GH with lambda fixed at %g (mu, lambda, alpha, beta, delta)', fixLam);
end
end

% =========================================================================
%                          TOP-LEVEL MODEL FITTERS
% =========================================================================

function st = fit_point(data, fixLam, warmZ)
% POINT stage: fit on z = (data - loc)/scl and assemble the data-unit results
% WITHOUT profile CIs (profile_ci is all NaN until fit_ci runs). Returns the
% z-unit state that fit_ci needs. fixLam is the fixed lambda (-0.5 NIG, 1
% hyperbolic) or NaN for GH, where lambda is free. warmZ (GH only) is a cell of
% polished z-unit point estimates of the fixed-lambda models, used as GH warm
% starts.
%
% WHY. The multistart grid, the optimizer and root-finder tolerances and the
% profile step sizes are absolute numbers, so results used to depend on the
% units of the data. Measured (Sep 2026) on one NIG sample: multiplying it by
% 1e4 moved the reported mu interval by 17%; on a GH sample, rescaling left the
% fit up to 3e-3 nats worse and delta 13 times different. The GH family is
% closed under y = loc + scl*z, with logL_y = logL_z - n*log(scl), so fitting z
% and transforming back is exact, and every output is equivariant under a
% shift or rescaling of the data. The delta floor is a fixed fraction of
% std(data), so it maps onto itself.
n   = numel(data);
loc = median(data);
scl = std(data);
z   = (data - loc) / scl;

if nfix(fixLam)
    [theta, nll, hess, exitflag, output] = gh_multistart(z, fixLam);
    model = model_name(fixLam);
else
    if nargin < 3, warmZ = {}; end
    [theta, nll, hess, exitflag, output] = fit_GH_core(z, warmZ);
    model = model_name(fixLam);
end
[theta, nll] = polish_point(theta, nll, z, fixLam);
infoZ = gh_observed_information(theta, z, fixLam);

st = struct('data', data, 'z', z, 'loc', loc, 'scl', scl, 'n', n, ...
    'fixLam', fixLam, 'model', model, 'thetaZ', theta(:), 'nllZ', nll, ...
    'hessZ', hess, 'infoZ', infoZ, 'exitflag', exitflag, 'output', output);
ciNaN = struct('lambda', [NaN NaN], 'beta', [NaN NaN], 'delta', [NaN NaN], ...
    'mu', [NaN NaN], 'alpha', [NaN NaN]);
% Silent here: fit_ci re-packages the model that is returned and raises its
% diagnostics then, so each returned result warns once, and a model AUTO
% discards does not warn about a result the caller never sees.
prevWarn = warning;
warning('off', 'all');
try
    st.results = package_results(st, ciNaN);
catch ME
    warning(prevWarn);
    rethrow(ME);
end
warning(prevWarn);
end

function [theta, nll] = polish_point(theta, nll, z, fixLam)
% Step 5 of ESTIMATION STRATEGY: derivative-free polish of the point estimate.
% The profile deviance is measured from nll, so an nll left short of the
% optimum also mislocates every confidence bound.
[lbP, ubP] = theta_bounds(z, fixLam);
[thetaP, nllP] = polish_nelder_mead(@(t) gh_nll_wrapper(t, z, fixLam), theta, lbP, ubP);
% Accept only an improvement larger than the likelihood's own floating-point
% noise, about n*max(Bessel argument)*eps (see NUMERICAL-VALIDITY GUARD).
% Without this the polish slid along the flat Normal-limit ridge toward the
% guard and kept "improvements" made of noise: measured (Sep 2026) on n = 10,
% a 2e-5 gain at Bessel argument 6.2e9, where the noise is ~1e-5, after which
% the fit failed its own LikelihoodInconsistent check (valid = false).
if nllP < nll - max(1e-9, nll_noise(thetaP, z, fixLam))
    theta = thetaP(:);
    nll   = nllP;
end
theta = theta(:);
end

function results = fit_ci(st, doProfile)
% CI stage: profile-likelihood CIs for a point-stage state ST, followed by a
% check that the point estimate really was the optimum. With doProfile false
% ('NoProfile' mode) the profiles are skipped: profile_ci is all NaN (fixed
% lambda excepted) and only the Wald intervals are reported.
if nargin < 2, doProfile = true; end
if ~doProfile
    ciZ = struct('lambda', [NaN NaN], 'beta', [NaN NaN], 'delta', [NaN NaN], ...
        'mu', [NaN NaN], 'alpha', [NaN NaN]);
    if nfix(st.fixLam), ciZ.lambda = [st.fixLam st.fixLam]; end
    results = package_results(st, ciZ);
    results.CIMethod = 'wald';
    return
end
%
% WHY THE CHECK. Every profile evaluation is a nuisance re-fit, and the
% profiles explore the likelihood far more thoroughly than the point
% multistart does. If any of them finds an NLL below nllMin, the point
% estimate was NOT the optimum, and CLAUDE.md's one-sidedness argument fails:
% too high an nllMin understates every deviance, so bounds cross late and come
% out too WIDE. So the lowest NLL seen anywhere in the profiles is tracked with
% the full theta it came from; if it beats nllMin by more than
% PROFILE_IMPROVE_TOL (or the NLL's own noise), the fit is re-optimized from
% that theta (fmincon plus the Nelder-Mead polish) and the CIs are recomputed,
% at most MAX_PROFILE_REFITS times.
refits = 0;
while true
    [ciZ, best] = gh_profile_ci(st.thetaZ, st.z, st.nllZ, st.fixLam, st.infoZ);
    tolImp = max(PROFILE_IMPROVE_TOL(), nll_noise(best.theta, st.z, st.fixLam));
    if ~(best.nll < st.nllZ - tolImp)
        break
    end
    if refits >= MAX_PROFILE_REFITS()
        warning('GeneralizedHyperbolic_MLE:ProfileImprovedOptimum', ...
            ['After %d re-optimization(s), the profile search still found an ' ...
            'NLL %.4g nats below the point estimate. The reported optimum is ' ...
            'not the true one, so the profile CIs may be too WIDE. Treat the ' ...
            'fit as unconverged.'], refits, st.nllZ - best.nll);
        break
    end
    refits = refits + 1;
    st = reoptimize_point(st, best.theta, best.nll);
end
results = package_results(st, ciZ);
results.profile_refits = refits;
results.CIMethod = 'profile';
end

function st = reoptimize_point(st, theta0, nll0)
% Re-fit the point estimate from theta0 (a profile evaluation that beat the
% current optimum, with value nll0), by fmincon then the Nelder-Mead polish.
% The better of the fmincon result and theta0 itself is kept, so the NLL can
% only go down.
[lb, ub] = theta_bounds(st.z, st.fixLam);
opts = fit_options(2500, 1e-7);
theta = theta0(:);  nll = nll0;
try
    [tr, nr, er, orr, ~, ~, hr] = fmincon( ...
        @(t) gh_nll_wrapper(t, st.z, st.fixLam), theta0(:), ...
        [], [], [], [], lb, ub, [], opts);
    st.hessZ = hr;  st.exitflag = er;  st.output = orr;
    if isfinite(nr) && nr < nll
        theta = tr(:);  nll = nr;
    end
catch
    % Keep theta0, which is already better than the old optimum.
end
[theta, nll] = polish_point(theta, nll, st.z, st.fixLam);
st.thetaZ = theta;
st.nllZ   = nll;
st.infoZ  = gh_observed_information(theta, st.z, st.fixLam);
end

function results = package_results(st, ciZ)
% Map a z-unit state and z-unit CIs back to the data's units and assemble the
% results structure with every point-estimate diagnostic.
[theta, hess] = unstandardize_theta(st.thetaZ, st.hessZ, st.loc, st.scl, st.fixLam);
nll = st.nllZ + st.n*log(st.scl);
ci  = ciZ;
ci.mu    = st.loc + st.scl*ci.mu;   % scl > 0 preserves order; no sort (NaN)
ci.alpha = ci.alpha / st.scl;
ci.beta  = ci.beta / st.scl;
ci.delta = ci.delta * st.scl;

p = theta_to_params(theta, st.fixLam);
results = assemble_results(st.model, p, theta, nll, hess, ci, st.exitflag, ...
    st.output, st.data, st.fixLam);
results = flag_delta_at_bound(results, st.data, st.fixLam);
results = flag_student_t_limit(results, st);
results.profile_refits = 0;
[~, info] = unstandardize_theta(st.thetaZ, st.infoZ, st.loc, st.scl, st.fixLam);
results = add_wald(results, theta, info, st.fixLam, st.infoZ, st.scl);
results = add_skewness(results, st);
results.CIMethod = 'none';
end

% =========================================================================
%                                SKEWNESS
% =========================================================================

function results = add_skewness(results, st)
% Dimensionless skewness summaries; see SKEWNESS in the header. Everything is
% computed in z-units: each measure is invariant under y = loc + scl*z with
% scl > 0, so the z-unit value is the data-unit value.
%
%   Skewness            moment skewness of the FITTED distribution, closed form
%   SkewnessSE          delta method from the observed information (infoZ)
%   SkewnessCIWald      Skewness +/- 1.96*SkewnessSE
%   SkewnessBowleyData  Bowley's quartile skewness of the DATA (mid-quantiles)
%   SkewnessBowleyModel Bowley's quartile skewness of the FITTED distribution
pZ = theta_to_params(st.thetaZ, st.fixLam);
g  = gh_moment_skewness(pZ);
% On the flagged Student-t limit gamma is not identified (any value down to
% 1e-40 fits equally well), and the moment skewness depends on gamma without
% bound there, so a number would be wherever the optimizer stopped.
if results.student_t_limit, g = NaN; end

se = NaN;
if isfinite(g) && results.SEValid
    % Five-point central differences of the closed form in theta (z-units);
    % the closed form is smooth, so h = 1e-4 relative is ample.
    th   = st.thetaZ(:);
    k    = numel(th);
    f    = @(t) gh_moment_skewness(theta_to_params(t, st.fixLam));
    grad = NaN(k, 1);
    for i = 1:k
        e = zeros(k, 1);  e(i) = 1e-4 * max(1, abs(th(i)));
        grad(i) = (-f(th + 2*e) + 8*f(th + e) - 8*f(th - e) + f(th - 2*e)) / (12*e(i));
    end
    Hz = (st.infoZ + st.infoZ.') / 2;
    [R, flag] = chol(Hz);
    if all(isfinite(grad)) && flag == 0
        v  = R.' \ grad;          % grad' * inv(Hz) * grad = ||R' \ grad||^2
        se = sqrt(v.' * v);
    end
end
zc = sqrt(2) * erfinv(0.95);

results.Skewness            = g;
results.SkewnessSE          = se;
results.SkewnessCIWald      = g + [-1 1]*zc*se;
results.SkewnessBowleyData  = bowley_mid(st.z);
results.SkewnessBowleyModel = gh_bowley_model(pZ);
end

function g = gh_moment_skewness(p)
% Moment skewness of GH(lambda, alpha, beta, delta, mu), from the normal
% variance-mean mixture X = mu + beta*V + sqrt(V)*Z, V ~ GIG(lambda, delta,
% gamma). The cumulant generating function is K_X(t) = mu*t + K_V(beta*t +
% t^2/2), which gives
%   Var X  = k1 + beta^2*k2
%   kappa3 = 3*beta*k2 + beta^3*k3
% with k1, k2, k3 the mean, variance and third central moment of V, and
%   E[V^r] = (delta/gamma)^r * K_{lambda+r}(w) / K_lambda(w),  w = delta*gamma.
% Scaled Bessel functions are used; the exp(w) factors cancel in each ratio.
% Returns NaN where the Bessel ratios are not finite.
w = p.delta * p.gamma;
g = NaN;
if ~(isfinite(w) && w > 0), return, end
K0 = besselk(p.lambda, w, 1);
r  = [besselk(p.lambda + 1, w, 1), besselk(p.lambda + 2, w, 1), ...
      besselk(p.lambda + 3, w, 1)] / K0;
s  = p.delta / p.gamma;
m1 = s   * r(1);
m2 = s^2 * r(2);
m3 = s^3 * r(3);
k2 = m2 - m1^2;
k3 = m3 - 3*m1*m2 + 2*m1^3;
v  = m1 + p.beta^2 * k2;
g  = (3*p.beta*k2 + p.beta^3*k3) / v^1.5;
if ~(isfinite(g) && isreal(g)), g = NaN; end
end

function b = bowley_mid(x)
% Bowley's quartile skewness (Q3 + Q1 - 2*Q2)/(Q3 - Q1), with quartiles read
% off the MID-DISTRIBUTION function: at each distinct value u_j with count
% c_j, Fmid(u_j) = (#{x < u_j} + c_j/2)/n, linearly interpolated (Parzen
% 2004). Without ties this is exactly quantile()'s default ((i - 0.5)/n,
% linear). With ties -- integer-minute data -- ordinary sample quartiles can
% only take a few lattice values, so the statistic jumps in coarse steps
% whose size depends on the spread of the data; the mid-quantiles move
% continuously with the counts. NaN when Q3 = Q1.
x = sort(x(:));
n = numel(x);
[u, ~, j] = unique(x);
b = NaN;
if numel(u) < 2, return, end
c  = accumarray(j, 1);
Fm = (cumsum(c) - c/2) / n;
P  = [0.25 0.5 0.75];
q  = interp1(Fm, u, P, 'linear');
q(P < Fm(1))   = u(1);
q(P > Fm(end)) = u(end);
d = q(3) - q(1);
if d > 0, b = (q(3) + q(1) - 2*q(2)) / d; end
end

function b = gh_bowley_model(p)
% Bowley's quartile skewness of the fitted GH (z-units), from quartiles found
% by root-finding on the CDF, itself integrated numerically from the density.
% The integral is split at mu, where a small delta concentrates the density
% into a near-cusp that an adaptive rule can step over. Unlike the moment
% skewness it exists for every fit, the Student-t limit included. NaN if the
% density cannot be integrated.
b = NaN;
pdf = @(x) reshape(exp(gh_logpdf(x(:), p)), size(x));
opt = {'AbsTol', 1e-12, 'RelTol', 1e-10};
prevW = warning('off', 'MATLAB:integral:NonFiniteValue');
cleanW = onCleanup(@() warning(prevW));
try
    lowMass = integral(pdf, -Inf, p.mu, opt{:});
    tot     = lowMass + integral(pdf, p.mu, Inf, opt{:});
    % A density that does not integrate to 1 means the integration failed
    % (a NaN from an unevaluable point, or a missed spike): report nothing.
    if ~(abs(tot - 1) < 1e-6), return, end
    F       = @(x) gh_mass_below(pdf, p.mu, lowMass, x, opt) / tot;
    q   = zeros(1, 3);
    P   = [0.25 0.5 0.75];
    for i = 1:3
        lo = p.mu - 1;  hi = p.mu + 1;
        while F(lo) > P(i) && lo > p.mu - 1e6, lo = p.mu - 2*(p.mu - lo); end
        while F(hi) < P(i) && hi < p.mu + 1e6, hi = p.mu + 2*(hi - p.mu); end
        q(i) = fzero(@(x) F(x) - P(i), [lo hi], optimset('TolX', 1e-10));
    end
catch
    return
end
d = q(3) - q(1);
if all(isfinite(q)) && d > 0, b = (q(3) + q(1) - 2*q(2)) / d; end
end

function v = gh_mass_below(pdf, mu, lowMass, x, opt)
% Unnormalized mass of PDF below x, split at mu (lowMass = mass below mu).
if x <= mu
    v = integral(pdf, -Inf, x, opt{:});
else
    v = lowMass + integral(pdf, mu, x, opt{:});
end
end

% =========================================================================
%                 OBSERVED INFORMATION, STANDARD ERRORS, WALD CIs
% =========================================================================

function H = gh_observed_information(theta, z, fixLam)
% Observed information (Hessian of the NLL) at theta, in the transformed
% coordinates and z-units, by central differences with Richardson
% extrapolation: H = (4*H(h/2) - H(h))/3, error O(h^4).
%
% WHY NOT fmincon's HESSIAN. The interior-point algorithm returns a BFGS
% approximation of the Lagrangian's Hessian, accumulated along the path the
% optimizer happened to take. It is not the observed information, and
% standard errors from it are not standard errors. The GH log-density has no
% closed-form derivative in lambda (the order of the Bessel function), so the
% information is taken numerically from the exact NLL.
%
% STEP. Every coordinate is O(1) in z-units (log scales for gamma and delta,
% data scaled to unit SD), so h = 2e-3*max(1, |theta_i|) is well inside the
% quadratic region of a regular optimum and far above rounding: the NLL's
% noise is about n*arg*eps (see nll_noise), divided by h^2 ~ 1e-6.
%
% Returns NaN(k) if any evaluation is not finite (e.g. a step crosses the
% Bessel-argument guard), so a caller never reads SEs off a partial matrix.
f  = @(t) gh_nll_wrapper(t, z, fixLam);
k  = numel(theta);
h  = 2e-3 * max(1, abs(theta(:)));
H1 = fd_hessian(f, theta(:), h);
H2 = fd_hessian(f, theta(:), h/2);
H  = (4*H2 - H1) / 3;
H  = (H + H.') / 2;
if ~all(isfinite(H(:))), H = NaN(k); end
end

function H = fd_hessian(f, x, h)
% Central-difference Hessian of f at x with per-coordinate steps h.
k  = numel(x);
H  = NaN(k);
f0 = f(x);
fp = zeros(k, 1);  fm = zeros(k, 1);
for i = 1:k
    e = zeros(k, 1);  e(i) = h(i);
    fp(i) = f(x + e);
    fm(i) = f(x - e);
    H(i,i) = (fp(i) - 2*f0 + fm(i)) / h(i)^2;
end
for i = 1:k
    for j = i+1:k
        ei = zeros(k, 1);  ei(i) = h(i);
        ej = zeros(k, 1);  ej(j) = h(j);
        v = (f(x+ei+ej) - f(x+ei-ej) - f(x-ei+ej) + f(x-ei-ej)) / (4*h(i)*h(j));
        H(i,j) = v;  H(j,i) = v;
    end
end
end

function results = add_wald(results, theta, info, fixLam, infoZ, scl)
% Standard errors and Wald CIs for the physical parameters, by the delta
% method from the observed information INFOZ (z-units; INFO, the same matrix
% in data units, is only stored), coordinates
% theta = [lambda; beta; psi; log_delta; mu], lambda absent if fixed).
%
%   SE, Cov  : in the order of Params, [mu; lambda; alpha; beta; delta].
%              lambda's SE is 0 when it is fixed (NIG, hyperbolic).
%   CIWald   : same layout as profile_ci. mu, lambda and beta are symmetric;
%              delta and alpha are formed on the log scale and
%              back-transformed, so both stay positive.
%
% At the delta floor (delta_at_bound) the optimum is constrained, the
% information is not the curvature of an interior maximum, and NO Wald SE is
% reported (all NaN, SEValid = false). Holding log_delta fixed and inverting
% the rest was tried and measured (Nov 2026): on GH samples (lambda = -1.5,
% n = 1000) whose optimum sat on the floor -- 12 of 100, each confirmed as the
% global optimum by refitting from the true parameters -- those conditional
% Wald intervals covered lambda, alpha and delta 0 times in 12, because fixing
% delta hides the lambda-delta trade-off that the data leave open. The
% profile intervals covered 0.75-0.92 there.
%
% A singular or indefinite information matrix (the flat ridges near the
% Normal and Student-t limits) gives NaN SEs and SEValid = false; profile
% intervals remain the honest summary there.
o    = nfix(fixLam);
k    = numel(theta);
p    = theta_to_params(theta, fixLam);
iLd  = 4 - o;                         % index of log_delta in theta
free = true(k, 1);

% Invert in z-units, where every coordinate is O(1), and map the covariance to
% data units with the Jacobian of unstandardize_theta, jt = d theta_y/d theta_z
% = diag(1, 1/scl, 1, 1, scl). Inverting the data-unit information instead
% made the rcond test depend on the units: rescaling the data by 1e-6 puts
% scl^-2 ~ 1e12 on the mu entry, rcond fell below 1e-14, and every SE came
% back NaN on a sample whose unscaled copy had valid SEs (measured Nov 2026).
C = NaN(k);
Hf = infoZ(free, free);
ok = all(isfinite(Hf(:)));
if ok
    [R, flag] = chol((Hf + Hf.')/2);
    ok = flag == 0 && rcond(Hf) > 1e-14;
    if ok
        Ri = R \ eye(size(R));
        jt = [ones(1, 1 - o), 1/scl, 1, 1, scl].';
        C(free, free) = (jt * jt.') .* (Ri * Ri.');
        C(~free, :) = 0;  C(:, ~free) = 0;
    end
end

% Jacobian of [mu; lambda; alpha; beta; delta] with respect to theta.
J = zeros(5, k);
iB = 2 - o;  iP = 3 - o;  iMu = 5 - o;
J(1, iMu) = 1;
if ~o, J(2, 1) = 1; end
J(3, iB)  = p.beta / p.alpha;          % alpha = hypot(beta, exp(psi))
J(3, iP)  = p.gamma^2 / p.alpha;
J(4, iB)  = 1;
J(5, iLd) = p.delta;

if ok
    Cov = J * C * J.';
    Cov = (Cov + Cov.') / 2;
    se  = sqrt(max(diag(Cov), 0));
    if results.delta_at_bound
        ok = false;  Cov = NaN(5);  se = NaN(5, 1);
    end
else
    Cov = NaN(5);  se = NaN(5, 1);
end
if o, se(2) = 0; if ok, Cov(2, :) = 0; Cov(:, 2) = 0; end, end

zc = sqrt(2) * erfinv(0.95);          % 1.959964, no toolbox needed
ci.lambda = p.lambda + [-1 1]*zc*se(2);
ci.beta   = p.beta   + [-1 1]*zc*se(4);
ci.delta  = p.delta * exp([-1 1]*zc*se(5)/p.delta);
ci.mu     = p.mu     + [-1 1]*zc*se(1);
ci.alpha  = p.alpha * exp([-1 1]*zc*se(3)/p.alpha);

results.SE          = se;
results.Cov         = Cov;
results.CIWald      = ci;
results.CILevel     = 0.95;
results.SEValid     = ok;
results.information = info;
end

function results = flag_student_t_limit(results, st)
% Report the skewed Student-t boundary. With lambda < 0, letting gamma -> 0
% (alpha -> |beta|) gives the skewed Student-t limit with nu = -2*lambda, and
% the likelihood goes FLAT there: measured (Nov 2026), logL was identical from
% gamma = 1e-8 down to 1e-40. Heavy-tailed data then drift along psi =
% log(gamma) until the optimizer stalls, so alpha and gamma are not
% identified: alpha's lower bound is effectively |beta|, and any reported
% gamma is where the search stopped, not an estimate.
%
% The test re-evaluates the NLL with psi pushed STUDENT_T_PSI_SHIFT lower
% (gamma times ~2e-9) at fixed (lambda, beta, delta, mu). If that costs less
% than STUDENT_T_TOL nats, the fit sits on the flat Student-t limit.
results.student_t_limit = false;
p = theta_to_params(st.thetaZ, st.fixLam);
if ~(p.lambda < 0)
    return
end
t   = st.thetaZ;
ip  = 3 - nfix(st.fixLam);               % position of psi in theta
t(ip) = t(ip) - STUDENT_T_PSI_SHIFT();
cost  = gh_nll_wrapper(t, st.z, st.fixLam) - st.nllZ;
if isfinite(cost) && cost < STUDENT_T_TOL()
    results.student_t_limit = true;
    warning('GeneralizedHyperbolic_MLE:StudentTLimit', ...
        ['The fit sits on the skewed Student-t limit (lambda = %.3g < 0, ' ...
        'gamma -> 0, nu = -2*lambda = %.3g): moving gamma down by a factor ' ...
        'of %.2g changes the NLL by only %.2g nats. alpha and gamma are not ' ...
        'identified here (alpha -> |beta|); report the Student-t limit ' ...
        'rather than the fitted alpha.'], ...
        p.lambda, -2*p.lambda, exp(STUDENT_T_PSI_SHIFT()), cost);
end
end

function v = STUDENT_T_PSI_SHIFT()
% How far psi = log(gamma) is pushed down to probe the Student-t limit.
v = 20;
end

function v = STUDENT_T_TOL()
% NLL change (nats) below which the Student-t limit is declared flat.
v = 1e-3;
end

function v = PROFILE_STEP_CAP()
% Largest first step of a profile bound search, in z-units (log units for
% delta and alpha). See gh_profile_ci.
v = 1;
end

function v = PROFILE_IMPROVE_TOL()
% A profile NLL must beat nllMin by more than this (nats) to trigger a refit.
v = 1e-3;
end

function v = MAX_PROFILE_REFITS()
% Cap on point-estimate re-optimizations triggered by the profile search.
v = 2;
end

function v = MIN_N_ERROR()
% Fewest observations accepted at all.
v = 10;
end

function v = MIN_N_WARN()
% Below this many observations, a SmallSample warning is raised.
v = 50;
end

function s = nll_noise(theta, data, fixLam)
% Rough floating-point noise in the summed NLL at theta: each log-density
% subtracts Bessel arguments of order ARG, so the sum carries an absolute error
% of about n*ARG*eps. Factor 10 for margin.
p = theta_to_params(theta, fixLam);
arg = max([p.delta*p.gamma; p.alpha*sqrt(p.delta^2 + (data - p.mu).^2)]);
s = 10 * numel(data) * max(arg, 1) * eps;
end

function [theta, hess] = unstandardize_theta(theta, hess, loc, scl, fixLam)
% Map theta and its Hessian from z-units to data units, y = loc + scl*z:
%   lambda unchanged, beta/scl, psi - log(scl), log_delta + log(scl), loc + scl*mu.
% The map is affine with Jacobian J = diag(1, 1/scl, 1, 1, scl) in the order
% (lambda, beta, psi, log_delta, mu), so H_y = inv(J)'*H_z*inv(J) = D*H_z*D
% with D = diag(1, scl, 1, 1, 1/scl). NIG has no lambda entry.
o = (1 - nfix(fixLam));
theta(1+o) = theta(1+o) / scl;
theta(2+o) = theta(2+o) - log(scl);
theta(3+o) = theta(3+o) + log(scl);
theta(4+o) = loc + scl*theta(4+o);
d = [ones(1,o), scl, 1, 1, 1/scl].';
if isequal(size(hess), [numel(d) numel(d)])
    hess = hess .* (d * d.');
end
end

function [theta, nll, hess, exitflag, output] = fit_GH_core(data, warmZ)
% GH point estimate: grid multistart, then a warm start from each fixed-lambda
% fit. warmZ is a cell of GH theta vectors (lambda first) built from the
% polished NIG and hyperbolic point estimates in the SAME standardized units
% (fit_point standardizes every model identically). Pass {} to fall back to a
% NIG multistart here.
% Step 1: multi-start search over a parameter grid.
[theta, nll, hess, exitflag, output] = gh_multistart(data, NaN);

% Step 2: always attempt a warm start from each nested fit. A fixed-lambda
% solution is a valid point in GH space, so it provides a useful starting
% point in a region where the GH likelihood is known to be reasonable. A warm
% start is accepted over the current solution if it converged when that did
% not, or if both have equal convergence status and it achieves a lower NLL.
%
% The warm-start POINT itself is also a candidate. GH nests both models, so
% its maximized likelihood can never be below theirs; without this, an
% fmincon run that wandered off a warm start could return a GH fit worse
% than the NIG or hyperbolic fit it started from, and the LRT statistic would
% be clipped at 0 for a reason that has nothing to do with the data.
if nargin < 2 || isempty(warmZ)
    [tN, ~, ~, ~, ~] = gh_multistart(data, -0.5);
    warmZ = {[-0.5; tN(:)]};
end
opts = fit_options(2500, 1e-7);

% Same delta floor as the grid search, so the warm start explores the
% identical feasible region.
[lb, ub] = theta_bounds(data, NaN);
objGH = @(t) gh_nll_wrapper(t, data, NaN);

for w = 1:numel(warmZ)
    theta0 = warmZ{w}(:);
    try
        [tr, nr, er, orr, ~, ~, hr] = fmincon(objGH, theta0, ...
            [], [], [], [], lb, ub, [], opts);
        oldConv  = exitflag > 0;
        newConv  = er > 0;
        useRescue = isfinite(nr) && ((newConv && ~oldConv) || ...
            (newConv == oldConv && nr < nll));
        if useRescue
            theta = tr;  nll = nr;  hess = hr;
            exitflag = er;  output = orr;
        end
    catch
        % Retain the best solution so far if the warm start crashes.
    end
    n0 = objGH(theta0);
    if isfinite(n0) && n0 < nll
        % The nested optimum beats every GH run: keep it (the polish in
        % fit_point can still move it off lambda's fixed value).
        theta = theta0;  nll = n0;  exitflag = max(exitflag, 1);
    end
end
end

% =========================================================================
%                         MULTI-START OPTIMIZATION
% =========================================================================

function [bestTheta, bestNLL, bestHess, bestExit, bestOut] = gh_multistart(data, fixLam)
muMed  = median(data);
muMean = mean(data);
[counts, edges] = histcounts(data, 'Normalization', 'pdf');
[~, im] = max(counts);
muMode = mean(edges(im:im+1));
muGrid = unique([muMed muMean muMode]);

s0 = std(data);
deltaGrid  = s0 * [0.5 1 1.5];
betaGrid   = [0  -0.5/s0  0.5/s0];
if nfix(fixLam)
    lambdaGrid = fixLam;
else
    lambdaGrid = [-1.5  -0.5  0.5  1.5];
end

[L, B, D, M] = ndgrid(lambdaGrid, betaGrid, deltaGrid, muGrid);
grid = [L(:) B(:) D(:) M(:)];

% psi = log(gamma); initialize gamma = 1/s0 so the distribution spread
% roughly matches the data spread.
psi0 = log(1/s0);
n = size(grid, 1);
starts = cell(n, 1);
vals   = inf(n, 1);

for i = 1:n
    if nfix(fixLam)
        t = [grid(i,2); psi0; log(grid(i,3)); grid(i,4)];
    else
        t = [grid(i,1); grid(i,2); psi0; log(grid(i,3)); grid(i,4)];
    end
    starts{i} = t;
    vals(i)   = gh_nll_wrapper(t, data, fixLam);
end

% Sort by NLL and keep the best finite starting points.
[~, ord] = sort(vals);
ord = ord(isfinite(vals(ord)));
if isempty(ord)
    error('GeneralizedHyperbolic_MLE:NoFiniteStart', ...
        'No finite initial likelihood value was found.');
end
ord  = ord(1:min(8, numel(ord)));
opts = fit_options(2000, 1e-7);

% Protective lower bound on log_delta. delta -> 0 makes the density spike at
% mu, and with tied or discretized data (integer bout durations, for example)
% that spike can drive the likelihood up without limit. The floor is a fixed
% fraction of the data spread so it scales with the problem. See
% delta_lower_bound for why the value is what it is.
[lb, ub] = theta_bounds(data, fixLam);

template = struct('theta',[], 'nll',Inf, 'hess',[], 'exitflag',-Inf, 'output',[]);
sol = repmat(template, numel(ord), 1);

% Suppress the singular-matrix chatter interior-point emits on flat regions of
% the likelihood. Capture and restore the caller's ENTIRE warning state rather
% than forcing these two identifiers back on: the caller may deliberately have
% had them off, and onCleanup guarantees restoration even if a start throws.
prevWarn    = warning;
restoreWarn = onCleanup(@() warning(prevWarn));
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:illConditionedMatrix');

for j = 1:numel(ord)
    try
        [t, v, ef, out, ~, ~, h] = fmincon( ...
            @(x) gh_nll_wrapper(x, data, fixLam), starts{ord(j)}, ...
            [], [], [], [], lb, ub, [], opts);
        sol(j) = struct('theta',t,'nll',v,'hess',h,'exitflag',ef,'output',out);
    catch ex
        sol(j).output = struct('message', ex.message);
    end
end

v  = [sol.nll];
ef = [sol.exitflag];
finite    = isfinite(v);
converged = finite & ef > 0;

if any(converged)
    candidates = find(converged);
elseif any(finite)
    candidates = find(finite);
    warning('GeneralizedHyperbolic_MLE:NoConvergedStart', ...
        'No run converged; returning the best finite solution.');
else
    error('GeneralizedHyperbolic_MLE:OptimizationFailed', ...
        'All multistart optimization attempts failed.');
end

[~, jj] = min(v(candidates));
best = sol(candidates(jj));
bestTheta = best.theta;  bestNLL = best.nll;
bestHess  = best.hess;   bestExit = best.exitflag;
bestOut   = best.output;
end

% =========================================================================
%                     OPTIMIZER OPTIONS (SHARED HELPER)
% =========================================================================

function opts = fit_options(maxIter, tol)
% Return fmincon options for interior-point minimization. Shared by the
% point estimate and every profile-CI nuisance optimization so both explore
% the same feasible region under the same tolerances.
opts = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', ...
    'MaxIterations', maxIter, 'MaxFunctionEvaluations', 20000, ...
    'OptimalityTolerance', tol, 'StepTolerance', 1e-10);
end

% =========================================================================
%                     LOG-LIKELIHOOD (CORE MATH)
% =========================================================================

function nll = gh_nll_wrapper(theta, data, fixLam)
% Convert transformed coordinates to physical params, then evaluate NLL.
p = theta_to_params(theta, fixLam);
nll = gh_nll_from_params(p, data);
end

function nll = gh_nll_from_params(p, data)
% Evaluate the GH negative log-likelihood given a physical parameter struct.
% The formula itself lives in gh_logpdf, shared with get_pdf_vals, so it
% exists once. This function adds the feasibility rules an OPTIMIZER needs.
[logpdf, argMax] = gh_logpdf(data, p);
if ~all(isfinite(logpdf))
    nll = Inf;  return
end

% ---- Numerical-validity guard on the Bessel arguments --------------------
% The log-density forms (log K - arg) for both Bessel terms, so it adds and
% subtracts quantities of order ARG. In double precision the absolute error of
% that cancellation is about arg*eps per observation, and roughly n*arg*eps in
% the summed NLL. Once that error is comparable to real differences in
% likelihood, an optimizer will happily "improve" the fit by chasing pure
% floating-point noise. MAX_BESSEL_ARG(n) caps it at ~1e-2 nats.
%
% This is not hypothetical. Without this guard, NIG fits to n >= 1000 ran away
% to alpha ~ 4e52, delta ~ 1e36, arg ~ 4e88, reporting a log-likelihood of
% +4.5e19 (impossible) while results.Fit evaluated to all NaN. Declaring such
% parameters infeasible confines the search to the region where the likelihood
% is actually computable. Legitimate fits stay far below the cap: even
% delta = 50 with alpha = 3 gives arg ~ 150.
if argMax > MAX_BESSEL_ARG(numel(data))
    nll = Inf;  return
end

nll = -sum(logpdf);
if ~isfinite(nll) || ~isreal(nll), nll = Inf; end
end

function [logpdf, argMax] = gh_logpdf(x, p)
% GH log-density at x for physical parameters p. NaN marks a point where it
% cannot be evaluated (invalid parameters, Bessel overflow or underflow), so a
% caller can tell "not evaluable" from "legitimately tiny": a log-density of
% -800 is a valid value even though exp() of it underflows to 0. argMax is the
% largest Bessel argument used, for the noise estimate and the guard.
%
% The GH log-density is (Barndorff-Nielsen 1977):
%
%   log f(x) = lambda*log(gamma) - lambda*log(delta) - 0.5*log(2*pi)
%              - (lambda - 0.5)*log(alpha)
%              - log K_lambda(delta*gamma)          [normalizing Bessel]
%              + (lambda/2 - 0.25)*log(delta^2 + (x-mu)^2)
%              + log K_{lambda-0.5}(alpha*sqrt(delta^2+(x-mu)^2))  [data Bessel]
%              + beta*(x - mu)
%
% Bessel functions are evaluated in scaled form, besselk(nu, x, 1) =
% besselk(nu, x)*exp(x), so log K(nu, x) = log(besselk(nu, x, 1)) - x.
logpdf = NaN(size(x));
argMax = NaN;
if ~all(isfinite([p.lambda p.beta p.gamma p.alpha p.delta p.mu])) || ...
        p.gamma <= 0 || p.alpha <= 0 || p.delta <= 0
    return
end
xm   = x - p.mu;
d2   = p.delta^2 + xm.^2;
arg0 = p.delta * p.gamma;     % argument of the normalizing Bessel K_lambda
argz = p.alpha * sqrt(d2);    % arguments of the data Bessel K_{lambda-0.5}
argMax = max([arg0; argz(:)]);
if ~isfinite(arg0) || arg0 <= 0
    return
end
K0 = besselk(p.lambda,       arg0, 1);
Kz = besselk(p.lambda - 0.5, argz, 1);
if ~isfinite(K0) || K0 <= 0
    return
end
logpdf = p.lambda*log(p.gamma) - p.lambda*log(p.delta) - 0.5*log(2*pi) ...
    - (p.lambda - 0.5)*log(p.alpha) - (log(K0) - arg0) ...
    + (p.lambda/2 - 0.25).*log(d2) + (log(Kz) - argz) + p.beta*xm;
bad = ~isfinite(argz) | argz <= 0 | ~isfinite(Kz) | Kz <= 0 | ...
    ~isfinite(logpdf) | imag(logpdf) ~= 0;
logpdf = real(logpdf);
logpdf(bad) = NaN;
end

function d = delta_lower_bound(data)
% Lower bound on delta, as a fraction of the data spread.
%
% Chosen at 1e-3 of std(data): far enough below any plausible fit that it never
% interferes (fitted deltas are typically within a factor of a few of the
% spread, i.e. ~1e3 times this bound), while still preventing delta from
% collapsing to zero on tied or discretized data. The value is arbitrary in the
% sense that no principled choice exists; what matters is that an active bound
% is REPORTED rather than silently returned as an interior optimum.
d = 1e-3 * std(data);
end

function [lb, ub] = theta_bounds(data, fixLam)
% Box constraints in transformed coordinates. Only log_delta is bounded; every
% other coordinate is free. Index of log_delta is 3 for NIG, 4 for GH.
mld = log(delta_lower_bound(data));
if nfix(fixLam)
    lb = [-Inf; -Inf; mld;  -Inf];
    ub = [ Inf;  Inf; Inf;   Inf];
else
    lb = [-Inf; -Inf; -Inf; mld;  -Inf];
    ub = [ Inf;  Inf;  Inf; Inf;   Inf];
end
end

function results = flag_delta_at_bound(results, data, fixLam)
% Report when the fit is sitting on the delta floor. Such a point is a
% CONSTRAINED optimum, not an MLE: the profile CIs are computed without the
% bound, and standard errors from the Hessian assume an interior solution.
% Returning it unannounced would misrepresent a boundary solution as a fit.
idx = 4 - nfix(fixLam);          % position of log_delta in theta
mld = log(delta_lower_bound(data));
% Tolerance in log(delta) units: 1e-2 means within 1% of the floor. fmincon's
% interior-point algorithm never lands ON a bound, it stops strictly inside
% one. Measured (Sep 2026) on data whose MLE is at the floor: the fit stopped
% 1.5e-3 above it on variance-gamma data and 7.5e-4 above it on Laplace data,
% so the old tolerance of 1e-8 reported delta_at_bound = false every time.
results.delta_at_bound = isfinite(results.theta_hat(idx)) && ...
    results.theta_hat(idx) <= mld + 1e-2;
if results.delta_at_bound
    warning('GeneralizedHyperbolic_MLE:DeltaAtBound', ...
        ['delta reached its lower bound (%.6g). This is a constrained ' ...
        'optimum, not an interior MLE: profile CIs and Hessian-based ' ...
        'standard errors assume an interior solution and are unreliable ' ...
        'here. Check for tied or heavily discretized data.'], ...
        delta_lower_bound(data));
end
end

function r = local_bisect(fun, a, b, tolX)
%LOCAL_BISECT Sign-only bisection, tolerant of non-finite values at the ends.
%   Assumes fun(a) <= 0 and fun(b) >= 0 (the caller established the bracket).
%   Only the sign of fun is used, so +Inf at an endpoint is fine. Returns NaN
%   if the bracket cannot be confirmed or a midpoint cannot be evaluated,
%   keeping NaN's meaning as "unknown" rather than "open".
r = NaN;
fa = fun(a);
if ~(isfinite(fa) && fa <= 0)
    return
end
for it = 1:200
    if abs(b - a) <= tolX
        break
    end
    m  = 0.5 * (a + b);
    fm = fun(m);
    if ~isfinite(fm)
        % Cannot evaluate here, so the bound is genuinely unknown. +Inf used to
        % count as "past the boundary", which is how a numerical overflow edge
        % was reported as a confidence bound; see find_root_robust.
        return
    elseif fm <= 0
        a = m;          % still inside the interval
    else
        b = m;          % past the boundary (finite or +Inf alike)
    end
end
r = 0.5 * (a + b);
end

function se = local_wald_se(H, k)
%LOCAL_WALD_SE Approximate standard errors from a Hessian, for step sizing only.
%   Returns NaN(k,1) whenever H is missing, the wrong shape, non-finite, not
%   symmetric positive definite, or badly conditioned. Callers must treat NaN
%   as "no information" and fall back to their own heuristic. Nothing computed
%   here reaches the reported confidence intervals.
se = nan(k, 1);
if nargin < 2 || isempty(H) || ~isequal(size(H), [k k]) || ~all(isfinite(H(:)))
    return
end
Hs = (H + H.') / 2;                       % symmetrize away round-off
[R, flag] = chol(Hs);                     % also tests positive definiteness
if flag ~= 0 || rcond(Hs) < 1e-12
    return
end
d = sum(inv(R).^2, 2);                    % diagonal of inv(Hs), without inv(Hs)
if all(isfinite(d)) && all(d > 0)
    se = sqrt(d(:));
end
end

function v = MAX_BESSEL_ARG(n)
% Largest Bessel argument for which the summed NLL over n points is
% numerically meaningful. The NLL carries an absolute error of roughly
% n*arg*eps (see gh_nll_from_params), so the cap keeps that at or below
% 1e-2 nats, two orders under the 1.92 profile threshold:
%   v = min(1e10, 1e-2/(n*eps)).
% A fixed 1e10 used to be justified as "below ~1e-2 even for n = 1e6"; in fact
% n*1e10*eps is 0.22 nats at n = 1e5 and 2.2 at n = 1e6, above the threshold.
% For n <= ~4500 the cap is still 1e10. Bessel arguments are products like
% alpha*delta and are scale-invariant, so the cap means the same in z-units
% and in data units.
v = min(1e10, 1e-2 / (n * eps));
end

% =========================================================================
%                     PROFILE LIKELIHOOD CONFIDENCE INTERVALS
% =========================================================================

function [ci, best] = gh_profile_ci(theta, data, nllMin, fixLam, hess)
% Compute 95% profile likelihood CIs for all physical parameters.
%
% For each parameter, the CI boundary satisfies:
%   nll_profile(param) - nllMin = chi2inv(0.95, 1) / 2  (~1.92)
%
% The search for each boundary is performed in a coordinate chosen to
% make the profile approximately parabolic and the step size well-scaled:
%   lambda, beta, mu : searched in their natural (untransformed) units.
%   delta            : searched in log(delta) space (scale-free).
%   alpha            : searched in log(alpha) space (scale-free, keeps alpha > 0).
crit = chi2inv(0.95, 1) / 2;   % ~1.9208; profile deviance threshold
% Nuisance optimizations run under the SAME box constraints as the point
% estimate, so fit and CI describe one feasible region (see regular_profile_dev).
opts = fit_options(1500, 1e-6);
p    = theta_to_params(theta, fixLam);
s    = std(data);               % used to set scale-appropriate step sizes

% Wald standard errors in the TRANSFORMED coordinates, used only to size the
% initial bracket step. The profile boundary sits near 1.96 SE for a locally
% quadratic profile, so starting there brackets the root in one or two probes
% instead of walking out from an arbitrary guess; every wasted probe is a full
% nuisance re-optimization. These SEs never enter the reported interval, which
% remains a pure profile result, so a poor Hessian costs speed and nothing
% else. theta's coordinates line up with the profiling coordinates exactly
% (log_delta is profiled in log space, beta/mu/lambda in their own units),
% which is what makes this usable directly.
seT = local_wald_se(hess, numel(theta));

% Lowest nuisance-refit NLL seen in ANY profile evaluation, with the full
% theta it came from (a handle, so every evaluation updates it). fit_ci uses
% it to check that nllMin really was the optimum.
tracker = containers.Map('KeyType', 'char', 'ValueType', 'any');
tracker('nll')   = nllMin;
tracker('theta') = theta(:);

if nfix(fixLam)
    ci.lambda = [fixLam fixLam];   % lambda is fixed (NIG, hyperbolic)
    names = {'beta','delta','mu'};
else
    ci = struct();
    names = {'lambda','beta','delta','mu'};
end

for i = 1:numel(names)
    name = names{i};
    switch name
        case 'lambda'
            idx = 1;
            q0   = p.lambda;
            step = max(0.1, 0.1*max(1, abs(q0)));
        case 'beta'
            idx  = 2 - nfix(fixLam);
            q0   = p.beta;
            step = max(0.05/s, 0.1*max(1/s, abs(q0)));
        case 'delta'
            idx  = 4 - nfix(fixLam);
            q0   = log(p.delta);
            step = 0.1;
        case 'mu'
            idx  = 5 - nfix(fixLam);
            q0   = p.mu;
            step = max(0.05*s, 0.1*max(s, abs(q0)));
    end

    % Prefer a Hessian-derived step when one is available; fall back to the
    % heuristic above when the Hessian is unusable (not positive definite,
    % ill-conditioned, or absent).
    if idx <= numel(seT) && isfinite(seT(idx)) && seT(idx) > 0
        step = 2 * seT(idx);
    end
    % Cap the first step (z-units). On a flat profile the Wald SE is huge,
    % and a first probe sized from it lands at the far end of the range,
    % where a cold-started nuisance fit fails; measured (Nov 2026) on a
    % near-Normal NIG sample of n = 60: a 7.1-unit first step in log(delta)
    % to 0.3 above the floor, a failed fit (+49.5), and a bound reported
    % UNKNOWN that walking out finds OPEN. The expanding bracket (1.5^k)
    % still reaches 100 units in about 12 probes.
    step = min(step, PROFILE_STEP_CAP());

    % One warm-start cache per parameter: each profile point starts its
    % nuisance optimization from the solution at the nearest point evaluated.
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    fun = @(q) regular_profile_dev(q, idx, theta, data, nllMin, crit, fixLam, opts, cache, tracker);
    % log(delta) is searched no lower than the estimator's own floor: the fit
    % may not return delta below it, so the interval must not claim support
    % there either. Without the limit the search walked down to log(delta)
    % ~ -700 and reported where besselk overflowed (delta ~ 1e-305) as a bound.
    if strcmp(name, 'delta')
        loLimit = log(delta_lower_bound(data));
    else
        loLimit = -Inf;
    end
    % The value AT the limit decides between an open bound and a finite one,
    % so it gets a harder nuisance fit (see find_root_robust).
    funHard = @(qq) regular_profile_dev(qq, idx, theta, data, nllMin, crit, fixLam, opts, cache, tracker, true);
    lo  = find_root_robust(fun, q0, -1, step, loLimit, funHard);
    hi  = find_root_robust(fun, q0,  1, step, Inf);
    % A crossing produced by the Bessel-argument guard, not the likelihood,
    % is an open bound (see open_if_capped).
    toP = @(qq, xx) theta_to_params(insert_param(setq(theta, idx, qq), xx, idx_mask(numel(theta), idx)), fixLam);
    lo  = open_if_capped(lo, -1, cache, toP, data);
    hi  = open_if_capped(hi,  1, cache, toP, data);

    % lo and hi come from searches below and above the MLE, so [lo hi] is
    % already ordered. Do NOT sort: sort moves NaN to the end, so an unknown
    % lower bound used to push a genuine upper bound into the lower slot
    % (delta reported as [0.2062, NaN] instead of [NaN, 0.2062]).
    if strcmp(name, 'delta')
        ci.(name) = exp([lo hi]);   % back-transform from log space
    else
        ci.(name) = [lo hi];
    end
end

% Alpha is profiled exactly by fixing it physically and re-optimizing
% all other parameters (see profile_alpha_ci for details).
% Its first step is 2 Wald SEs of log(alpha) when the information allows,
% as for the other parameters (it used to be a fixed 0.1).
stepA = 0.1;
if all(isfinite(seT))
    k  = numel(theta);  iB = 2 - nfix(fixLam);  iP = 3 - nfix(fixLam);
    g  = zeros(k, 1);
    g(iB) = p.beta / p.alpha^2;            % d log(alpha) / d beta
    g(iP) = p.gamma^2 / p.alpha^2;         % d log(alpha) / d psi
    Hs = (hess + hess.')/2;
    vA = g.' * (Hs \ g);
    if isfinite(vA) && vA > 0
        stepA = min(max(2*sqrt(vA), 1e-3), PROFILE_STEP_CAP());
    end
end
ci.alpha = profile_alpha_ci(theta, p.alpha, data, nllMin, crit, fixLam, opts, tracker, stepA);
best = struct('nll', tracker('nll'), 'theta', tracker('theta'));
end

function d = regular_profile_dev(q, idx, theta, data, nllMin, crit, fixLam, opts, cache, tracker, hard)
% Profile deviance for parameters other than alpha.
% Fix theta(idx) = q, re-optimize the remaining parameters,
% and return  nll_profile - nllMin - crit.
%
% The nuisance optimization is CONSTRAINED to the same feasible region as the
% point estimate. Without this the profile would be free to push delta below
% the floor that the fit itself was forbidden to cross, so the CI could report
% support for parameter values the estimator was never allowed to return.
%
% STARTS. The nuisance optimum drifts as q moves away from the MLE, and a
% single start at the MLE's nuisance values can stall short of it. A stalled
% optimization OVERSTATES the profile, an overstated profile crosses the
% threshold early, and the interval comes out too narrow. So each point starts
% from the solution at the nearest point already evaluated (CACHE), and a
% value at or past the threshold -- the only kind that can set a bound -- is
% re-checked from the MLE's nuisance values, keeping the lower of the two.
mask       = true(numel(theta), 1);
mask(idx)  = false;
base       = theta;
base(idx)  = q;
obj = @(x) gh_nll_wrapper(insert_param(base, x, mask), data, fixLam);

[lbFull, ubFull] = theta_bounds(data, fixLam);
lb = lbFull(mask);
ub = ubFull(mask);

if nargin < 11, hard = false; end
starts = {min(max(theta(mask), lb), ub)};   % fmincon needs a feasible start
if idx == 4 - nfix(fixLam)
    % Profiling log(delta): also start from the MLE's nuisance values with
    % psi = log(gamma) shifted by the same amount as log(delta). That keeps
    % delta/gamma, hence the mean shift beta*delta/gamma and the scale of the
    % variance, roughly where the data need them, instead of asking the
    % optimizer to rediscover them from values fitted at another delta.
    % Measured (Sep 2026), log-normal data at the delta floor: every other
    % start stalled at +252 above the threshold where an independent profile
    % was 1.92 below it. psi sits at position 3 - fixLam of the masked vector.
    xm = theta(mask);
    xm(3 - nfix(fixLam)) = xm(3 - nfix(fixLam)) + (q - theta(idx));
    starts = [starts, {min(max(xm, lb), ub)}];
end
if nargin >= 9 && isKey(cache, 'Q')
    % Warm start from the nearest evaluated point; a HARD evaluation (a value
    % that decides between an open and a finite bound) uses the 8 nearest.
    Q = cache('Q');  X = cache('X');
    [~, order] = sort(abs(Q - q));
    nWarm = 1;
    if hard, nWarm = min(8, numel(order)); end
    warm = arrayfun(@(j) min(max(X(:,j), lb), ub), order(1:nWarm), 'UniformOutput', false);
    starts = [warm, starts];
end
[v, xbest] = multistart_nuisance(obj, starts, lb, ub, opts, nllMin + crit, hard);
if nargin >= 9 && isfinite(v)
    cache_store(cache, q, xbest);
end
if isfinite(v)
    track_best(tracker, v, insert_param(base, xbest, mask));
end
d = v - nllMin - crit;
end

function track_best(tracker, v, thetaFull)
% Record (v, thetaFull) in TRACKER (a containers.Map handle) if v is the
% lowest NLL seen so far in any profile evaluation.
if v < tracker('nll')
    tracker('nll')   = v;
    tracker('theta') = thetaFull(:);   %#ok<NASGU> containers.Map is a handle
end
end

function [v, xbest] = multistart_nuisance(obj, starts, lb, ub, opts, vStop, hard)
% Minimize OBJ from each start in turn, stopping as soon as one lands below
% vStop: a value below the threshold cannot set a bound, so it needs no second
% opinion. Returns the lowest value found (Inf if every attempt failed).
%
% If every fmincon start ends at or past vStop, the answer decides a bound,
% and fmincon is not trusted with it alone. On the flat ridge near the Normal
% limit the interior-point solver stops early (exitflag 2, step below
% tolerance): measured (Sep 2026) at profile values 1,490-1,598 ABOVE the
% threshold where Nelder-Mead reached 1.13 BELOW it, and those stalls laid
% out a smooth, self-consistent false crossing that nothing else could catch.
% So the best point found and the first (warm) start are polished by
% derivative-free search before the value is returned.
v = Inf;  xbest = starts{end};
for s = 1:numel(starts)
    try
        [xs, vs] = fmincon(obj, starts{s}, [], [], [], [], lb, ub, [], opts);
    catch
        continue
    end
    if isfinite(vs) && vs < v
        v = vs;  xbest = xs;
    end
    if v < vStop
        return
    end
end
if nargin >= 7 && hard
    cand = [{xbest}, starts];   budget = 5000;   % every start, 5x the budget
else
    cand = {xbest, starts{1}};  budget = 1000;
end
for s = 1:numel(cand)
    [xs, vs] = polish_nelder_mead(obj, cand{s}, lb, ub, budget);
    if vs < v
        v = vs;  xbest = xs;
    end
    if v < vStop
        return
    end
end
end

function [x, v] = polish_nelder_mead(obj, x0, lb, ub, budget)
% Derivative-free minimization of OBJ within [lb, ub], bounds imposed by
% clamping, restarted from its own result because Nelder-Mead stalls. Returns
% x0 unchanged (with its value) if it cannot improve on it.
clampd = @(y) min(max(y, lb), ub);
f = @(y) obj(clampd(y));
x = clampd(x0);
v = f(x);
if ~isfinite(v), return; end
if nargin < 5, budget = 1000; end
o = optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8, ...
             'MaxFunEvals', budget*numel(x), 'MaxIter', budget*numel(x));
for r = 1:2
    try
        [xr, vr] = fminsearch(f, x, o);
    catch
        break
    end
    if isfinite(vr) && vr < v
        x = clampd(xr);  v = vr;
    else
        break
    end
end
end
function cache_store(cache, q, x)
% Append one evaluated profile point to a warm-start cache (a handle object).
if isKey(cache, 'Q')
    cache('Q') = [cache('Q'), q];
    cache('X') = [cache('X'), x(:)];   %#ok<NASGU> containers.Map is a handle
else
    cache('Q') = q;
    cache('X') = x(:);                 %#ok<NASGU> containers.Map is a handle
end
end
% ---- Exact profile CI for physical alpha -----------------------------------

function ci = profile_alpha_ci(theta, aHat, data, nllMin, crit, fixLam, opts, tracker, step)
% Profile alpha by searching in log(alpha) space.
% Using q = log(alpha) keeps alpha positive throughout the search and
% makes the step size scale-invariant regardless of alpha's magnitude.
if nargin < 9, step = 0.1; end
q0    = log(aHat);
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
fun   = @(q) alpha_profile_dev_log(q, theta, data, nllMin, crit, fixLam, opts, cache, tracker);
lo  = find_root_robust(fun, q0, -1, step, -Inf);
hi  = find_root_robust(fun, q0,  1, step,  Inf);
toP = @(qq, xx) alpha_params(exp(qq), xx, fixLam);
lo  = open_if_capped(lo, -1, cache, toP, data);
hi  = open_if_capped(hi,  1, cache, toP, data);
ci  = exp([lo hi]);   % ordered by construction; see gh_profile_ci on sort
end

function d = alpha_profile_dev_log(logA, theta, data, nllMin, crit, fixLam, opts, cache, tracker)
% Profile deviance as a function of log(alpha).
a = exp(logA);
if ~isfinite(a) || a <= 0, d = Inf;  return;  end
v = alpha_fixed_nll(a, theta, data, fixLam, opts, cache, logA, nllMin + crit, tracker);
if isfinite(v), d = v - nllMin - crit;  else, d = Inf;  end
end

function nll = alpha_fixed_nll(a, theta, data, fixLam, opts, cache, logA, vStop, tracker)
% Minimize NLL over all parameters except alpha, which is fixed at a.
%
% Reparameterize beta = a*tanh(xi), xi free in (-inf, inf).
% This enforces |beta| < a automatically without any inequality constraint.
% Consequently gamma = sqrt(a^2 - beta^2) = a*sqrt(1 - tanh(xi)^2)
%                    = a / cosh(xi),  always > 0.
%
% Free optimization variables:
%   NIG : [xi; log_delta; mu]
%   GH  : [lambda; xi; log_delta; mu]
%
% log_delta carries the same floor as the point estimate, so this profile
% explores the identical feasible region (see regular_profile_dev, which also
% explains the warm start and the re-check against vStop).
p = theta_to_params(theta, fixLam);

% Convert the current beta estimate to the xi parameterization.
% Clamp to avoid atanh(+/-1) = +/-inf at the boundary.
r   = max(min(p.beta/a, 1-1e-8), -1+1e-8);
xi0 = atanh(r);

mld = log(delta_lower_bound(data));
if nfix(fixLam)
    x0  = [xi0; log(p.delta); p.mu];
    lb  = [-Inf; mld; -Inf];
    ub  = [ Inf; Inf;  Inf];
    obj = @(x) alpha_core(a, fixLam, x(1), x(2), x(3), data);
else
    x0  = [p.lambda; xi0; log(p.delta); p.mu];
    lb  = [-Inf; -Inf; mld; -Inf];
    ub  = [ Inf;  Inf; Inf;  Inf];
    obj = @(x) alpha_core(a, x(1), x(2), x(3), x(4), data);
end
starts = {min(max(x0, lb), ub)};   % fmincon needs a feasible start
if isKey(cache, 'Q')
    Q = cache('Q');  X = cache('X');
    [~, j] = min(abs(Q - logA));
    starts = [{min(max(X(:,j), lb), ub)}, starts];
end
[nll, xbest] = multistart_nuisance(obj, starts, lb, ub, opts, vStop);
if isfinite(nll)
    cache_store(cache, logA, xbest);
    % Back to theta: beta = a*tanh(xi), psi = log(a/cosh(xi)), with
    % log(cosh(xi)) = |xi| + log1p(exp(-2|xi|)) - log(2), which cannot
    % overflow. Re-evaluated in theta coordinates, so the tracked value is
    % exactly what a refit from that theta would start at.
    xi   = xbest(2 - nfix(fixLam));
    beta = a * tanh(xi);
    psi  = log(a) - (abs(xi) + log1p(exp(-2*abs(xi))) - log(2));
    if nfix(fixLam)
        thetaFull = [beta; psi; xbest(2); xbest(3)];
    else
        thetaFull = [xbest(1); beta; psi; xbest(3); xbest(4)];
    end
    track_best(tracker, gh_nll_wrapper(thetaFull, data, fixLam), thetaFull);
end
end
function root = open_if_capped(root, direction, cache, toP, data)
% Report a profile crossing as OPEN (+/-Inf) when the nuisance fits that
% located it sit against the Bessel-argument guard (MAX_BESSEL_ARG), i.e.
% when the feasible region, not the likelihood, ended the search.
%
% WHY. Near a flat ridge the nuisance optimum runs out to ever larger alpha
% and beta until the guard declares it infeasible; the best FEASIBLE value
% then rises, and the root finder records that rise as a crossing. Measured
% (Nov 2026, hyperbolic samples of n = 500): alpha upper bounds of 7.7e8 to
% 8.5e8 against estimates of 3 to 13, in 7 of the 9 replicates examined --
% always the same few e8, the place where alpha*sqrt(delta^2 + x^2) reaches
% the cap of 1e10, whatever the data. CLAUDE.md already listed such bounds as
% "effectively open"; this makes them so. A bound is only ever WIDENED here,
% and only where the likelihood was not what stopped it.
if ~isfinite(root) || ~isKey(cache, 'Q')
    return
end
Q = cache('Q');  X = cache('X');
[~, ord] = sort(abs(Q - root));
cap = MAX_BESSEL_ARG(numel(data));
for j = ord(1:min(2, numel(ord)))
    [~, argMax] = gh_logpdf(data, toP(Q(j), X(:, j)));
    if isfinite(argMax) && argMax > CAPPED_FRACTION() * cap
        root = direction * Inf;
        return
    end
end
end

function v = CAPPED_FRACTION()
% A nuisance fit whose largest Bessel argument is within this fraction of the
% guard counts as stopped by the guard. Legitimate fits sit near 1e1-1e3;
% only fits on the Normal ridge approach the guard, where a bound cannot be
% confirmed either way and open is the honest report.
v = 0.1;
end

function t = setq(theta, idx, q)
% theta with coordinate idx set to q.
t = theta;  t(idx) = q;
end

function m = idx_mask(k, idx)
% Logical mask of the k theta coordinates other than idx.
m = true(k, 1);  m(idx) = false;
end

function p = alpha_params(a, x, fixLam)
% Physical parameters from the alpha-profile coordinates (alpha = a fixed):
% x = [xi; log_delta; mu] (fixed lambda) or [lambda; xi; log_delta; mu].
o = nfix(fixLam);
if o, p.lambda = fixLam; else, p.lambda = x(1); end
xi = x(2 - o);
p.beta  = a * tanh(xi);
p.gamma = a / cosh(xi);
p.alpha = a;
p.delta = exp(x(3 - o));
p.mu    = x(4 - o);
end

function nll = alpha_core(a, lambda, xi, logDelta, mu, data)
% NLL with alpha = a fixed and beta = a*tanh(xi).
% gamma = a/cosh(xi) exactly. The old a*sqrt(max(1 - tanh(xi)^2, realmin))
% cancels: relative error 8e-5 at xi = 15 and 100% from xi = 19, where tanh
% rounds to 1. cosh overflow gives gamma = 0, which the density treats as
% invalid (NLL = Inf).
p.lambda   = lambda;
p.beta     = a * tanh(xi);
p.gamma    = a / cosh(xi);
p.alpha    = a;
p.delta    = exp(logDelta);
p.mu       = mu;
nll        = gh_nll_from_params(p, data);
end

% =========================================================================
%                         ROOT FINDER (PROFILE BOUNDARY)
% =========================================================================

function root = find_root_robust(fun, q0, direction, scale, limit, funHard)
% Find where the profile deviance FUN crosses zero, moving from q0 in
% DIRECTION, by an exponentially expanding bracket search that never passes
% LIMIT (default: unbounded).
%
% Algorithm:
%   1. Evaluate f0 = fun(q0). If it is non-finite or > 1e-6, return NaN: the
%      deviance should be about -1.92 at the MLE, so anything else is a
%      structural problem (inconsistent nllMin/theta, a failed nuisance
%      optimization), not noise.
%   2. Probe q0 + direction*scale*1.5^k, k = 0..30, clipped to LIMIT.
%   3. A NON-FINITE value is never read as a crossing. The search backs off
%      toward the last finite point (back_off_to_finite), because a long step
%      can jump over a real crossing into a region that cannot be evaluated.
%   4. At the first sign change between two FINITE values, fzero on that
%      bracket, with sign-only bisection as the fallback. Every evaluation
%      used to refine it is recorded, and the crossing is REJECTED if any was
%      non-finite or they contradict a single crossing (see
%      profile_evaluations_suspect); a rejected crossing is handled as in 5.
%   5. No crossing: +/-Inf (OPEN) when there is positive evidence the profile
%      stays below the threshold -- a finite value at LIMIT, finite values all
%      the way out, or finite values at least 10*min(step, 1) standardized
%      units from the MLE before evaluation broke down. Otherwise NaN
%      (UNKNOWN).
%
% WHY A NON-FINITE VALUE IS NOT A CROSSING. For this family the density is
% positive and finite wherever the parameters are valid, so an infinite
% profile NLL is never the likelihood genuinely rising: it is besselk
% overflowing or a nuisance fit failing. This routine used to treat +Inf as
% "past the boundary" and bisect onto the point where it first appeared.
% Measured (Sep 2026): delta intervals [2.4e-305, 0.91] on variance-gamma data
% and [3.7e-297, 1.76] on Laplace data, where the profile is flat all the way
% down and the honest lower bound is 0; lambda [-292.75, 292.76] on a weakly
% identified GH sample, where 292.76 was simply where besselk gave up.
if nargin < 5, limit = direction*Inf; end
if nargin < 6, funHard = []; end
root = NaN;
f0 = fun(q0);
if ~isfinite(f0) || f0 > 1e-6, return; end

previous     = q0;
fp           = f0;
farthest     = 0;       % distance of the farthest finite below-threshold probe
reachedLimit = false;   % finite below-threshold value obtained AT the limit
brokeDown    = false;   % any probe failed to evaluate
for k = 0:30
    current = q0 + direction * scale * (1.5^k);
    atLimit = direction * (current - limit) >= 0;
    if atLimit, current = limit; end
    fc = fun(current);
    if atLimit && ~isempty(funHard) && ~(isfinite(fc) && fc < 0)
        % The limit decides between OPEN and a finite bound, so one failed
        % nuisance fit must not decide it. Measured (Sep 2026) at the delta
        % floor on log-normal data: the ordinary fit returned +252 where an
        % independent profile was 1.92 below the threshold.
        fcH = funHard(current);
        if isfinite(fcH) && (~isfinite(fc) || fcH < fc)
            fc = fcH;
        end
    end
    if ~isfinite(fc)
        brokeDown = true;
        % The step may have jumped straight past a real crossing into a region
        % where the profile cannot be evaluated. Measured (Sep 2026): a Wald
        % step of 27.6 in log(delta) went from the MLE to log(delta) = 26,
        % where every value is Inf (Bessel guard), although the profile crosses
        % the threshold far inside that. So back off toward the last finite
        % point until a finite value is found.
        [previous, fp, farthest, current, fc, found] = ...
            back_off_to_finite(fun, q0, previous, fp, farthest, current);
        if ~found
            break   % the evaluable region ends here without a crossing in it
        end
    end
    if fc >= 0 && fp <= 0
        % Sign change found: bracket is [previous, current].
        %
        % Give fzero a tolerance matched to what a confidence bound is worth.
        % By default fzero refines to machine precision, and here every one of
        % its iterations costs a full nuisance re-optimization, so the default
        % spends most of the fit computing digits that are then reported to
        % four significant figures. Profiling an n=2000 NIG fit put 6.87 s of
        % 14.96 s inside fzero alone. A tolerance of 1e-4 of the local scale
        % leaves the bound far more accurate than the 95% interval is
        % meaningful to.
        %
        % The tolerance is ABSOLUTE, 1e-4 of the local step. fzero does not
        % treat TolX as absolute: it stops once the bracket is narrower than
        % 2*TolX*max(|b|,1), so TolX is scaled down by the bracket's
        % magnitude here. This used to pass 1e-4*max(|q0|,scale) straight in,
        % which at |q| ~ 1e4 is an effective tolerance of ~1e4: fzero stopped
        % at once and returned a bracket END as the bound (a mu interval off
        % by 17%, Sep 2026). On standardized data lambda can still reach ~50.
        tolX  = max(1e-12, 1e-4 * scale);
        tolFz = tolX / max([abs(previous), abs(current), 1]);
        rec = containers.Map('KeyType', 'char', 'ValueType', 'any');
        rec('q') = [previous, current];
        rec('f') = [fp, fc];
        frec = @(q) recorded_eval(fun, q, rec);
        try
            root = fzero(frec, sort([previous current]), optimset('TolX', tolFz));
        catch
            root = NaN;
        end
        if ~isfinite(root)
            % fzero needs finite values at both bracket ends and throws
            % otherwise; bisection needs only the sign at each midpoint.
            root = local_bisect(frec, previous, current, tolX);
        end
        % A crossing is only as good as the evaluations that located it. Near a
        % numerical breakdown the nuisance fit can return a FINITE but grossly
        % wrong value, which fakes a crossing as readily as Inf did. Measured
        % (Sep 2026) on a weakly identified GH sample: the lambda profile read
        % -0.76 at 290.30 on one evaluation and +10.4 at the same point on the
        % next, an earlier probe returned +245, and an independent profile was
        % -0.76 throughout. So a crossing whose refinement met a non-finite or
        % contradictory value is rejected and treated as a breakdown.
        if profile_evaluations_suspect(rec('q'), rec('f'), direction)
            % Credit the farthest point the refinement itself evaluated below
            % the threshold as distance covered.
            Qr = rec('q');  Fr = rec('f');
            okBelow = isfinite(Fr) & Fr <= -0.1;
            if any(okBelow)
                farthest = max(farthest, max(abs(Qr(okBelow) - q0)));
            end
            root = open_or_unknown(direction, farthest, scale);
        end
        return
    end
    previous = current;
    fp       = fc;
    farthest = abs(current - q0);
    if atLimit
        reachedLimit = true;
        break
    end
end

% No crossing. Open and unknown mean opposite things, so report them apart.
% A finite below-threshold value AT the limit is open even when farthest is
% 0, which happens when the MLE itself sits on the limit (delta on its floor).
if reachedLimit || (farthest > 0 && ~brokeDown)
    root = direction * Inf;
elseif farthest > 0
    root = open_or_unknown(direction, farthest, scale);
end
end

function root = open_or_unknown(direction, farthest, scale)
% No usable crossing. OPEN (+/-Inf) if the profile had been evaluated below
% the threshold at least 10*min(scale, 1) from the MLE; UNKNOWN (NaN) otherwise.
% The profiles run on standardized data, where every coordinate is O(1), so 10
% units is far beyond where an identified parameter crosses. The min() matters
% because the Wald step can be absurd: measured (Sep 2026) steps of 27.6 and
% 50.4 in log(delta) from a poor quasi-Newton Hessian, which would have
% demanded evidence out to log(delta) = 500. Measured the other way: a lambda
% profile still at -1.85 at 243 units from its MLE before besselk overflowed is
% flat, not unknown.
if farthest >= 10*min(scale, 1)
    root = direction * Inf;
else
    root = NaN;
end
end

function [previous, fp, farthest, current, fc, found] = ...
    back_off_to_finite(fun, q0, previous, fp, farthest, current)
% Bisect, on finiteness, between the last finite point PREVIOUS (below the
% threshold) and a non-finite point CURRENT. Every finite value below the
% threshold moves PREVIOUS (and FARTHEST) outward. Returns found = true, with
% CURRENT and fc set, at the first finite value at or past the threshold: a
% crossing bracket [previous, current]. Returns found = false once the edge of
% the evaluable region is located to 0.1% of its distance from the MLE.
found = false;  fc = NaN;
lo = previous;  hi = current;
tol = 1e-3 * max(1, abs(hi - q0));
for it = 1:60
    if abs(hi - lo) <= tol
        break
    end
    mid = 0.5 * (lo + hi);
    fm = fun(mid);
    if ~isfinite(fm)
        hi = mid;
    elseif fm >= 0
        current = mid;  fc = fm;  found = true;
        return
    else
        lo = mid;  previous = mid;  fp = fm;  farthest = abs(mid - q0);
    end
end
current = hi;
end

function v = recorded_eval(fun, q, rec)
% Evaluate FUN at q and append (q, value) to REC, a containers.Map (a handle,
% so the caller sees every evaluation fzero and local_bisect make).
v = fun(q);
rec('q') = [rec('q'), q];
rec('f') = [rec('f'), v];   %#ok<NASGU> containers.Map is a handle
end

function bad = profile_evaluations_suspect(Q, F, direction)
% True when the evaluations that located a crossing cannot all belong to one
% profile crossing the threshold once: any non-finite value; a point at or
% INSIDE another that is clearly above the threshold while the farther one is
% clearly below it; or a CLIFF, where even the nearest positive value sits
% more than 1 deviance unit above the threshold. Once fzero or bisection has
% converged, the nearest positive value lies within the root tolerance and is
% a small fraction of a unit; a larger one means the crossing never refined,
% because the values jump -- a nuisance fit falling over, not the likelihood
% rising. Measured (Sep 2026): -0.70 then +252 within 0.015 in log(delta)
% beside the delta floor on log-normal data; and alpha/beta upper bounds of
% 1.2e6, where an independent profile was still 1.30 below the threshold.
% Rejecting is safe in one direction only, which is the right one: profile
% evaluations over-estimate the profile, so they can make a bound too narrow
% but never too wide. The 0.1 of slack absorbs optimizer noise at the root.
margin = 0.1;
bad = any(~isfinite(F));
if bad, return; end
s = direction * Q(:);            % larger s = farther from the MLE
above = F(:) >= margin;
below = F(:) <= -margin;
pos = F(:) > 0;
if any(pos)
    Fp = F(pos);  Sp = s(pos);
    [~, j] = min(Sp);
    if Fp(j) > 1
        bad = true;
        return
    end
end
for i = find(above).'
    if any(below & s >= s(i))
        bad = true;
        return
    end
end
end

% =========================================================================
%                              UTILITIES
% =========================================================================

function p = theta_to_params(theta, fixLam)
% Convert transformed optimization coordinates to physical GH parameters.
%
% Transformation:
%   gamma = exp(psi)              ensures gamma > 0
%   alpha = hypot(beta, gamma)    ensures alpha > |beta| (smooth at beta=0)
%   delta = exp(log_delta)        ensures delta > 0
%
% hypot is used instead of sqrt(beta^2+gamma^2) for numerical stability
% when either beta or gamma is very large or very small.
if nfix(fixLam)
    p.lambda = fixLam;
    p.beta   = theta(1);
    p.gamma  = exp(theta(2));
    p.delta  = exp(theta(3));
    p.mu     = theta(4);
else
    p.lambda = theta(1);
    p.beta   = theta(2);
    p.gamma  = exp(theta(3));
    p.delta  = exp(theta(4));
    p.mu     = theta(5);
end
p.alpha = hypot(p.beta, p.gamma);
end

function theta = insert_param(base, sub, mask)
% Reconstruct a full theta vector by inserting free-parameter values
% (sub) back into the positions indicated by mask.
theta        = base;
theta(mask)  = sub;
end

function results = assemble_results(model, p, theta, nll, hess, ci, ef, out, data, fixLam)
% Package all estimation outputs into a single results structure.
n    = numel(data);
k    = numel(theta);
fitX = sort(data);

results.Model   = model;
results.Params  = [p.mu; p.lambda; p.alpha; p.beta; p.delta];
results.LogLike = -nll;
[logFit, argMax] = gh_logpdf(data, p);
results.Fit     = exp(logFit);              % density at original data points
results.LogFit  = logFit;                   % log-density there (NaN = not evaluable)
results.FitX    = fitX;                     % sorted x for smooth density plots
results.FitPDF  = get_pdf_vals(fitX, p);    % density at sorted points
results.AIC     = 2*k + 2*nll;
results.BIC     = k*log(n) + 2*nll;
results.theta_hat = theta;

% ---- Validate that the reported solution is a real fit -------------------
% A returned parameter vector is only meaningful if the density can actually
% be evaluated there AND reproduces the reported likelihood. Nothing used to
% look, so a degenerate run could report converged=true alongside an all-NaN
% density and a nonsensical log-likelihood. Check both and say so loudly.
%
% The check works on the LOG-density. It used to test exp(logpdf) > 0, and a
% legitimate log-density below about -745 (a gross outlier) underflows to 0,
% which set valid = false with a misleading "optimizer ran away" warning.
results.valid = true;
badFit = ~all(isfinite(logFit));
if badFit
    results.valid = false;
    warning('GeneralizedHyperbolic_MLE:DensityNotEvaluable', ...
        ['The fitted parameters do not yield an evaluable density (%d of %d ' ...
        'points have a non-finite log-density). results.LogLike is not ' ...
        'trustworthy. This usually means the optimizer ran away to extreme ' ...
        'parameters.'], nnz(~isfinite(logFit)), n);
else
    llCheck = sum(logFit);
    % The tolerance includes the log-density's own floating-point noise (see
    % NUMERICAL-VALIDITY GUARD). Near the Normal limit, Bessel arguments of
    % 1e9-1e10 are legitimate and the summed log-likelihood carries about
    % n*arg*eps of noise; a fixed 1e-6 failed a sound fit (n = 10, argument
    % 6.2e9, discrepancy 2e-5; Sep 2026).
    tol = max([1e-6, 1e-8 * abs(nll), 10 * n * max(argMax, 1) * eps]);
    if ~isfinite(llCheck) || abs(llCheck + nll) > tol
        results.valid = false;
        warning('GeneralizedHyperbolic_MLE:LikelihoodInconsistent', ...
            ['Reported log-likelihood (%.6g) disagrees with the density ' ...
            'evaluated at the fitted parameters (%.6g). The optimum is not ' ...
            'numerically reliable.'], -nll, llCheck);
    end
end

if nfix(fixLam)
    results.theta_names = {'beta','psi_log_gamma','log_delta','mu'};
else
    results.theta_names = {'lambda','beta','psi_log_gamma','log_delta','mu'};
end
% ---- Gaussian-limit / identifiability diagnostic -------------------------
% GH approaches the Normal distribution as zeta = delta*gamma grows. Deep in
% that limit the shape parameters stop being identifiable: wildly different
% (alpha, beta, delta) give near-identical densities, the likelihood develops
% a flat ridge, and the fit can return parameters orders of magnitude from the
% truth while remaining a perfectly valid density. Observed on data generated
% with alpha=2, delta=50: the fit returned alpha=157, delta=4204 with a beta
% profile flat to 1e-4 and a CI spanning [1.3, 7373] for alpha.
%
% Nothing is wrong with the arithmetic in that situation, so none of the
% validity checks above will fire. The honest signal is a comparison against
% the 2-parameter Normal the family is collapsing onto: if the extra
% parameters buy nothing by AIC, they should not be interpreted.
mu_n  = mean(data);
s2_n  = mean((data - mu_n).^2);              % ML variance (1/n)
if s2_n > 0
    ll_n  = -0.5*n*log(2*pi*s2_n) - 0.5*n;
else
    ll_n  = NaN;
end
aic_n = 2*2 - 2*ll_n;
results.gaussian_check = struct( ...
    'zeta',            p.delta * p.gamma, ...
    'LogLike_normal',  ll_n, ...
    'AIC_normal',      aic_n, ...
    'beats_normal',    results.AIC < aic_n);
if isfinite(aic_n) && ~results.gaussian_check.beats_normal
    warning('GeneralizedHyperbolic_MLE:NoGainOverNormal', ...
        ['The %d-parameter fit does not improve on a 2-parameter Normal by ' ...
        'AIC (%.2f versus %.2f, zeta = delta*gamma = %.3g). The data are ' ...
        'close to Gaussian, where GH shape parameters are not identifiable: ' ...
        'the density may be fine while alpha, beta and delta are essentially ' ...
        'arbitrary. Treat the parameter values and their CIs as unreliable.'], ...
        k, results.AIC, aic_n, p.delta * p.gamma);
end

results.nll                  = nll;
results.hessian              = hess;
results.hessian_coordinates  = results.theta_names;
results.profile_ci           = ci;
results.exitflag             = ef;
results.output               = out;
results.converged            = ef > 0;
end

function pdf = get_pdf_vals(x, p)
% Evaluate the GH density at points x given physical parameters p. NaN flags
% a point where the log-density cannot be evaluated; a legitimately tiny
% density underflows to 0, which is not a failure.
pdf = exp(gh_logpdf(x, p));
end

function plot_results(data, results)
% Plot histogram of data with superimposed fitted density.
% Pin the light theme BEFORE creating axes: in R2025a and later a new figure
% follows the desktop theme, and 'Color','w' alone can leave dark axes.
fh = figure('Color', 'w');
try, theme(fh, 'light'); catch, end   %#ok<NOCOMMA> theme() does not exist before R2025a
histogram(data, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'k');
hold on
plot(results.FitX, results.FitPDF, 'r-', 'LineWidth', 2);
xlabel('Data Value');
ylabel('Density');
title(['Fit: ' results.Model]);
grid on;  box on;
legend('Data', 'Fitted Model', 'Location', 'best');
end
