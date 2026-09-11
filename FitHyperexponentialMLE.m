function H = FitHyperexponentialMLE(eventseries, xmin, options)
%FITHYPEREXPONENTIALMLE Maximum-likelihood fit of a left-truncated mixture
%of exponentials (continuous mode) or geometric distributions (discrete
%mode), with sample-size-corrected AIC (AICc) model-order selection and
%identifiability gating.
%
%   H = FITHYPEREXPONENTIALMLE(eventseries, xmin)
%   H = FITHYPEREXPONENTIALMLE(eventseries, xmin, Name=Value, ...)
%
%   MODEL
%   (T denotes bout duration as a random variable; t denotes a specific
%   value/observation it takes -- so "f(t)" is a density evaluated at t,
%   and "T >= xmin" is an event/condition on the random variable itself.)
%
%   Continuous mode: each component j is Exponential(lambda_j). The
%   UNTRUNCATED mixture density and survival function, for K components
%   with weights w_j (w_j>=0, sum_j w_j = 1) and rates lambda_j>0, are
%
%       f(t) = sum_j  w_j * lambda_j * exp(-lambda_j * t),    t > 0
%       S(t) = sum_j  w_j * exp(-lambda_j * t)
%
%   and the density actually fit -- conditioned on T >= xmin, since that
%   is all you can ever observe -- is
%
%       f(t | T >= xmin) = f(t) / S(xmin),    t >= xmin.
%
%   Discrete mode: each component j is Geometric on n = 1,2,3,... (number
%   of SamplingInterval-sized steps until the event), with per-step
%   "success" probability p_j = 1 - exp(-lambda_j * SamplingInterval) --
%   the exact discretization of a rate-lambda_j Poisson clock observed
%   only once every SamplingInterval. lambda_j HERE IS THE SAME
%   continuous-time rate as in continuous mode (tau_j=1/lambda_j the same
%   time constant, same units as xmin/SamplingInterval) -- it is what's
%   actually fitted and reported in Rates/Tau in both modes. p_j is a
%   derived, dt-dependent quantity used only inside the likelihood below,
%   not itself the reported parameter, which is exactly why Rates/Tau
%   from continuous- and discrete-mode fits stay directly comparable to
%   each other. Writing dt = SamplingInterval and n_min = max(1,
%   round(xmin/dt)), the untruncated mixture pmf and survival are
%
%       P(N=n)  = sum_j  w_j * (1-p_j)^(n-1) * p_j,    n = 1,2,3,...
%       P(N>=n) = sum_j  w_j * (1-p_j)^(n-1)
%
%   and the truncated pmf actually fit is
%
%       P(N=n | N>=n_min) = P(N=n) / P(N>=n_min),    n >= n_min,
%
%   with observed duration recovered as t = n*dt (each eventseries value
%   is likewise mapped to n = round(t/dt) before evaluating this).
%
%   This fits the RAW, UNBINNED durations directly by maximum likelihood
%   -- using the mixture density itself, not a histogram-count
%   least-squares surrogate -- and CONDITIONS the likelihood on T>=xmin
%   via the mixture's own survival function above, rather than shifting
%   the data by xmin before fitting.
%
%   The shift-by-xmin trick (subtract xmin, fit an unshifted exponential
%   to what's left) is exact for a SINGLE exponential, by the memoryless
%   property. It is NOT exact for a mixture: a hyperexponential mixture's
%   hazard is decreasing, not constant, so conditioning on T>=xmin
%   reweights the components (the faster component is disproportionately
%   depleted by xmin) rather than simply re-centering. This function
%   normalizes by S(xmin) directly instead, which is exact for any number
%   of components.
%
%   ======================================================================
%   TWO KINDS OF WEIGHT -- READ THIS BEFORE REPORTING ANY NUMBER
%   ======================================================================
%   The truncated density above can be rewritten EXACTLY as a mixture over
%   the *observable* population:
%
%       f(t | T>=xmin) = sum_j  q_j * lambda_j * exp(-lambda_j*(t-xmin))
%       q_j = w_j * exp(-lambda_j*xmin) / S(xmin)      (continuous)
%       q_j = w_j * (1-p_j)^(n_min-1) / P(N>=n_min)    (discrete)
%
%   so the data constrain (q, lambda) directly, and w is only recovered
%   afterwards by the back-transform w_j ~ q_j * exp(+lambda_j*xmin).
%
%   H.Selected.WeightsObserved     = q  (fraction of the OBSERVED, i.e.
%                                        >= xmin, bouts from component j)
%   H.Selected.WeightsUntruncated  = w  (mixing weights over ALL bouts the
%                                        process produced, those below
%                                        xmin included)
%
%   There is deliberately NO field called "Weights". The unqualified name
%   is what a reader reaches for expecting "the fraction of bouts in this
%   component" -- which is q -- whereas the mixture's own weight
%   PARAMETER is w, and on real data the two are different numbers by a
%   factor of five or more. Rather than let ".Weights" quietly return
%   whichever one, it does not exist: H.AllFits(K).Weights raises
%   "Reference to non-existent field", and you have to say which you
%   meant. Their standard errors follow the same names,
%   WeightsObservedSE and WeightsUntruncatedSE.
%
%   REPORT q, NOT w, when describing your data. The two differ by orders
%   of magnitude whenever lambda_j*xmin is large: a component with
%   tau_j << xmin is almost entirely truncated away, so it can carry a
%   large untruncated weight while contributing almost nothing observable.
%   The back-transform gain exp(+lambda_j*xmin) also destroys the
%   floating-point precision of w: at tau_j/xmin = 1/160 the gain is
%   ~1e69, which pushes the remaining weights to ~1e-67 and eventually to
%   literal zeros. All-but-one weights printing as 0.0000 is that
%   underflow, not a fitting error -- check WeightsObserved, which stays
%   well scaled.
%
%   Two further traps when writing q up:
%
%     - tau_j is NOT the typical duration of that component's bouts. Every
%       observation is conditioned on T >= xmin and each component is
%       individually memoryless, so component j's OBSERVED durations have
%       mean xmin + tau_j. A tau_j = 99 s component under xmin = 300 s
%       gives bouts averaging 399 s, not 99 s. tau_j is the decay constant
%       of the excess over the threshold -- which is precisely what makes
%       it comparable across different choices of xmin.
%
%     - q_j is a proportion over a LATENT component label, not a rule for
%       classifying individual bouts. The components overlap heavily; one
%       observation's posterior is
%           P(j | t) ~ q_j * lambda_j * exp(-lambda_j*(t-xmin))
%       and on a real 3-component sleep-bout fit only about a third of
%       bouts could be assigned to a single component with >90 percent
%       confidence. n*q_j is an expected count, not a subset you can point
%       at. (As an internal check, those posteriors summed over all
%       observations equal n*q_j exactly at the MLE -- the EM stationarity
%       identity, and a useful confirmation that the optimizer really is
%       at a stationary point rather than merely somewhere with a good
%       function value.)
%
%   So what is w for? In principle it is the xmin-INVARIANT quantity: q is
%   defined relative to your threshold, so it is not comparable between
%   studies using different cutoffs, nor between sleep bouts (xmin=300)
%   and wake bouts (xmin=2), whereas w is a property of the process alone.
%   In practice that advantage does not cash out. Refitting one real
%   dataset at xmin = 300, 350, ..., 500 held tau_2 and tau_3 to ~1.5
%   percent while w for the fast component wandered over 0.61..1.00 and
%   collapsed entirely at the top of that range -- recovering w means
%   dividing q by an ever-smaller survival probability, which amplifies
%   noise faster than the added truncation removes it. Use w for
%   simulating the untruncated process, and for estimating what fraction
%   of the population fell below xmin (1 - S(xmin)), an extrapolation
%   resting entirely on the exponential form holding below xmin where
%   there is no data. To compare against a study using a different cutoff
%   x, do NOT pass w around; predict their weights from your own fit as
%       q_j(x) ~ w_j * exp(-lambda_j*x),
%   keeping the compared quantity one the data actually constrain.
%
%   ======================================================================
%   THE CONTINUOUS LIKELIHOOD IS UNBOUNDED IF ANY OBSERVATION EQUALS xmin
%   ======================================================================
%   In continuous mode the truncated density at the truncation point is
%   f(xmin | T>=xmin) = sum_j q_j*lambda_j, which is UNBOUNDED in
%   lambda_j. So if m = #{i : t_i == xmin} > 0, then letting one
%   component's rate diverge at fixed q gives
%
%       logL  ->  const + m*log(lambda_j)  ->  +Inf,
%
%   and NO maximum likelihood estimate exists for K >= 2. The optimizer
%   climbs this ridge until the weight logits saturate in floating point,
%   and then stops at an arbitrary point -- reporting a spurious
%   sub-second component, all-but-one weights at ~0, and a
%   non-positive-definite Hessian. This is the left-truncation analogue of
%   the classical unbounded-likelihood problem for normal mixtures (Day
%   1969, Biometrika 56:463-474; McLachlan & Peel 2000, Finite Mixture
%   Models, ch. 2) and of the unbounded three-parameter shifted-lognormal
%   likelihood handled in SHIFTLOGNORMAL_MLE.
%
%   Only ties AT xmin do this. Repeated values in the interior are
%   harmless: the exponential components are anchored at xmin and decrease
%   from there, so no component can concentrate on an interior point.
%
%   CALLING THIS WITH xmin = min(eventseries) GUARANTEES THE PATHOLOGY,
%   because the minimum is attained by at least one observation by
%   definition. Do not do it. xmin is a pre-specified truncation threshold
%   set by your recording/scoring protocol (e.g. 300 s for the standard
%   5-minute sleep rule), not a quantity to estimate from the data.
%
%   Three ways out, in order of preference:
%     1. Use DistributionType="discrete" (with SamplingInterval matching
%        your acquisition rate). Integer-valued durations ARE discrete,
%        the geometric pmf is bounded by 1, and no divergence is possible.
%        This is the statistically correct model for gridded data, not
%        just a numerical workaround.
%     2. Pass the protocol's true xmin rather than min(eventseries).
%     3. Keep continuous mode but rely on the MaxRate ceiling and the
%        degeneracy gating below, which stop the runaway and refuse to
%        SELECT a model order that exploits it. The fits are still
%        reported, flagged, in H.AllFits.
%
%   IDENTIFIABILITY GATING
%   A model order K is excluded from AICc selection (but still returned in
%   H.AllFits, with Degenerate=true and a reason string) when any of:
%     - (n - k - 1) <= 0, so AICc is undefined;
%     - a component's expected observed count n*q_j < MinExpectedCount:
%       the component is not supported by the observable data;
%     - two components' time constants agree to within DuplicateRateTol
%       in log space: the mixture has collapsed to fewer components;
%     - a component sits at the MaxRate ceiling: the optimizer was pushing
%       tau below the sampling resolution;
%     - lambda_j*xmin > MaxTruncationExponent: the component is so far below
%       the truncation point that exp(-lambda_j*xmin) underflows relative to
%       the other components, so its untruncated weight carries no
%       significant digits at all (this is the deterministic guard against
%       the unbounded-likelihood runaway described above, and does not
%       depend on the Hessian being detected as singular);
%     - logL(K) - logL(K-1) < DeadComponentTol: the extra component is
%       dead and contributes nothing;
%     - the Hessian is not positive definite (CovValid=false), i.e. the
%       parameters are not locally identified.
%   This is why an unguarded fit can "select" K=4 whose AICc is lower only
%   because the extra parameters are unidentified: AICc's penalty assumes
%   k estimable parameters (Hurvich & Tsai 1989, Biometrika 76:297-307),
%   which fails exactly here.
%
%   REVISION NOTES: likelihood evaluation is done entirely in log space
%   via logsumexp (avoids the flat, zero-gradient region a linear-space
%   floor-at-realmin creates for fminsearch), and the discrete branch's
%   log(p) is computed via log(-expm1(-x)) rather than log(1-exp(-x))
%   (avoids catastrophic cancellation for slow rates relative to the
%   sampling interval -- confirmed numerically: the naive form is already
%   off by ~0.1 in log-probability at rate*dt=1e-16). logsumexp_mat
%   carries the same isinf guard as logsumexp_vec -- without it, a single
%   observation whose density underflows to exactly zero under every
%   component simultaneously poisons the entire log-likelihood to NaN for
%   that evaluation instead of correctly contributing -Inf for just that
%   one observation.
%
%   INPUTS
%   eventseries : numeric vector of positive event durations, same units
%                 as xmin (e.g. seconds)
%   xmin        : positive scalar, left-truncation cutoff (e.g. 300 for
%                 sleep bouts under the standard 5-minute rule; your
%                 chosen minimum for wake bouts). A PROTOCOL constant --
%                 never min(eventseries), see above.
%
%   NAME-VALUE OPTIONS
%   MaxComponents      positive integer, default 4. Fits K = 1..MaxComponents
%                       and selects among them by AICc.
%   DistributionType    "continuous" (default) or "discrete". Use
%                       "discrete" whenever the durations live on a grid
%                       (integer multiples of the sampling interval),
%                       which is the usual case for bout data, and ALWAYS
%                       when any observation equals xmin. The relevant
%                       criterion is grid-valued data / ties at xmin, NOT
%                       the size of xmin/SamplingInterval: a fit with
%                       xmin/dt = 300 is still broken by a single tie at
%                       xmin. Same concern Clauset/Shalizi/Newman (2009,
%                       SIAM Rev. 51:661-703) raise for discrete vs.
%                       continuous power laws.
%   SamplingInterval    REQUIRED positive scalar, same time units as
%                       xmin. There is no default -- omitting it raises
%                       SamplingIntervalRequired. Sets the discrete-mode
%                       step size, the default MaxRate ceiling, and the
%                       rate range used
%                       to seed the multistart -- so it is consulted in
%                       BOTH modes, not only in discrete mode. Set it to
%                       your true acquisition interval, IN THE SAME UNITS
%                       as eventseries and xmin. If you rescale the
%                       durations you must rescale this too: converting
%                       seconds to minutes means xmin=5 AND
%                       SamplingInterval=1/60, not the default 1. Getting
%                       this wrong is silent and costly -- see the
%                       grid-spacing guard, reported in
%                       H.Diagnostics.GridMismatch.
%   MaxRate             positive scalar, default 1/SamplingInterval. Hard
%                       ceiling on fitted rates, i.e. tau_j >=
%                       SamplingInterval: a time constant shorter than the
%                       sampling interval is not resolvable in principle.
%                       Enforced by a smooth quadratic barrier in log-rate
%                       so fminsearch keeps a restoring gradient rather
%                       than hitting a flat plateau. Pass Inf to disable
%                       (not recommended in continuous mode).
%   MinExpectedCount    nonnegative scalar, default 5. A component whose
%                       expected observed count n*q_j falls below this is
%                       treated as unsupported (degenerate).
%   MaxTruncationExponent  positive scalar, default log(1/eps) ~ 36.04. A
%                       component with lambda_j*xmin above this is treated
%                       as degenerate: the back-transform gain
%                       exp(+lambda_j*xmin) then exceeds 1/eps, so the
%                       untruncated weights lose every significant digit
%                       (they print as 0.0000 while one weight prints as
%                       1.0000). At the default the component is also
%                       contributing less than eps of the observable
%                       probability relative to a component at tau ~ xmin.
%   DuplicateRateTol    nonnegative scalar, default 0.05. Two components
%                       whose |log(tau_i)-log(tau_j)| is below this are
%                       treated as collapsed (degenerate).
%   DeadComponentTol    nonnegative scalar, default 1e-4. Minimum logL
%                       improvement over K-1 required for K to count as a
%                       genuine extra component.
%   ErrorOnNoValidFit   logical, default true. When true (the default), a
%                       sample with fewer than 5 usable observations raises
%                       TooFewData and a sample where no model order passes
%                       the identifiability gating raises NoValidFit. Set
%                       it FALSE for batch work -- fitting one animal at a
%                       time -- where a single pathological sample would
%                       otherwise abort the whole loop. The call then
%                       returns normally with H.Failed = true,
%                       H.FailureReason naming the cause, H.SelectedK =
%                       NaN, H.Selected an empty (unfitted) template, and
%                       whatever per-K fits were obtained still in
%                       H.AllFits for inspection. The output struct has the
%                       same fields in the same order either way, so
%                       results(i) = H works across a mix of successes and
%                       failures; filter afterwards on [results.Failed].
%   SortComponents      logical, default true. Return components sorted by
%                       ascending Tau, so component indices mean the same
%                       thing across K, across animals, and across
%                       bootstrap replicates. The unconstrained mixture
%                       likelihood is invariant to relabelling, so an
%                       unsorted fit returns components in an arbitrary
%                       order (the label-switching problem: Stephens 2000,
%                       JRSS-B 62:795-809; Celeux, Hurn & Robert 2000,
%                       JASA 95:957-970).
%   nStartsBase, nStartsPerComponent, maxStarts, RandomSeed, MaxIter,
%   MaxFunEvals, TolX, TolFun, Verbose  -- multistart optimization
%                       settings.
%
%   OUTPUT
%   H.SelectedK        chosen number of components (minimum AICc among
%                       model orders that pass the identifiability gating)
%   H.AllFits          1 x MaxComponents struct array, one per K, with
%                       fields K, k, n, WeightsObserved, WeightsUntruncated,
%                       Rates, Tau, RateSE, TauSE, WeightsObservedSE,
%                       WeightsUntruncatedSE, PointwiseLogLik,
%                       CovValid, LogLik, AIC, AICc, Success, Converged,
%                       Degenerate, DegenerateReason, ExitFlag,
%                       BestParamVector
%   H.Selected         convenience copy of H.AllFits(H.SelectedK)
%   PointwiseLogLik    n x 1 per-observation log-likelihood contributions,
%                       summing to LogLik. Use these to compare this fit
%                       against a different model family (a power law, an
%                       exponentiated Weibull) by a Vuong (1989,
%                       Econometrica 57:307-333) test, which is the right
%                       non-nested test since the mixture likelihood-ratio
%                       has no chi-square null:
%                         D = H1.Selected.PointwiseLogLik - L2;
%                         V = sqrt(n)*mean(D)/std(D,1);   % ~N(0,1)
%                       Both sides must be fitted to the SAME observations
%                       under the SAME measure -- so the same
%                       DistributionType, SamplingInterval and xmin.
%   H.Failed           true if the fit could not be completed (only
%                       reachable with ErrorOnNoValidFit=false; otherwise
%                       those cases throw). When true, SelectedK is NaN.
%   H.FailureReason    '' on success, else 'TooFewData: ...' or
%                       'NoValidFit: ...'
%   H.Diagnostics      struct: nAtXmin, XminEqualsDataMin, LooksGridded,
%                       UnboundedContinuousLikelihood, GridSpacing,
%                       GridMismatch, nDistinctData, nDistinctOnGrid,
%                       nDropped fields, ExcludedK and Recommendation
%   H.DistributionType, H.xmin, H.n
%
%   Tau = 1./Rates (time units, matching eventseries/xmin). RateSE, TauSE,
%   WeightsUntruncatedSE, WeightsObservedSE are asymptotic standard errors, obtained
%   post-hoc (after fminsearch converges, not part of the optimization
%   itself) from a numerical Hessian of the negative log-likelihood at the
%   MLE -- the observed information matrix, which is the preferable
%   variance estimator here (Efron & Hinkley 1978, Biometrika 65:457-487)
%   -- inverted to get the covariance of the underlying log-rate/logit
%   parameters, then propagated to Rates, Tau, WeightsUntruncated and
%   WeightsObserved via the delta method. CovValid is false (and the SE
%   fields are NaN) when the Hessian is not usable as a covariance matrix
%   (non-positive-definite, e.g. near label-switching or a poorly
%   identified component) -- point estimates are still reported in that
%   case, just without uncertainty. These are asymptotic (large-n)
%   standard errors and inherit the same "don't fully trust on sparse
%   per-fly fits" caveat as everything else here; a cluster/parametric
%   bootstrap would be more robust but far more expensive, and isn't what
%   this computes.
%
%   IMPORTANT -- VALIDATE BEFORE TRUSTING: run this on simulated data with
%   known K, weights, rates, and xmin (matching your real sample sizes)
%   before trusting it on sparse per-fly fits. Multi-exponential
%   likelihoods are multimodal and only weakly identified (Redner & Walker
%   1984, SIAM Rev. 26:195-239), and AICc itself becomes unreliable once n
%   is only a few times k=2K-1.
%
%   A second check, cheaper than simulation and available on your own
%   data: refit at a higher threshold x > xmin and compare the refitted q
%   against what the ORIGINAL fit predicts there via
%   q_j(x) ~ w_j*exp(-lambda_j*x). That is a genuine out-of-sample test of
%   the mixture form -- a model fitted to bouts >= 300 s predicting the
%   composition of bouts >= 500 s -- and is stronger evidence for K than
%   the AICc ranking, which only ever compares in-sample fit. Agreement to
%   <= 0.01 in q is what a well-specified 3-component fit looks like. Do
%   not expect the FAST component's tau to survive this: it is exactly
%   what the raised threshold removes, so its estimate legitimately
%   degrades (99 s -> 36 s over that range on the dataset above) while the
%   slow components stay put to ~1.5 percent. AICc over-selection is real
%   too: in simulation at n=400 with a true K=2, AICc chose K=3 in 4 of 40
%   replicates, so prefer a parametric bootstrap likelihood-ratio test
%   (McLachlan 1987, Appl. Stat. 36:318-324) when the model ORDER is
%   itself the scientific claim -- the mixture LRT has no chi-square null
%   distribution (Hartigan 1985, Proc. Berkeley Conf. II:807-810; Chen,
%   Chen & Kalbfleisch 2001, JRSS-B 63:19-29).
%
%   See also SHIFTLOGNORMAL_MLE, SIGWORTHSINEV3,
%   CALCBOUTSIZESONSETSOFFSETS.

arguments
    eventseries double {mustBeReal}
    xmin (1,1) double {mustBePositive}
    options.MaxComponents (1,1) double {mustBeInteger,mustBePositive} = 4
    options.DistributionType (1,1) string {mustBeMember(options.DistributionType,["continuous","discrete"])} = "continuous"
    options.SamplingInterval (1,1) double = NaN
    options.MaxRate (1,1) double = NaN
    options.MinExpectedCount (1,1) double {mustBeNonnegative} = 5
    options.MaxTruncationExponent (1,1) double {mustBePositive} = log(1/eps)
    options.DuplicateRateTol (1,1) double {mustBeNonnegative} = 0.05
    options.DeadComponentTol (1,1) double {mustBeNonnegative} = 1e-4
    options.SortComponents (1,1) logical = true
    options.ErrorOnNoValidFit (1,1) logical = true
    options.nStartsBase (1,1) double {mustBeInteger,mustBePositive} = 10
    options.nStartsPerComponent (1,1) double {mustBeInteger,mustBePositive} = 10
    options.maxStarts (1,1) double {mustBeInteger,mustBePositive} = 80
    options.RandomSeed = []
    options.MaxIter (1,1) double {mustBeInteger,mustBePositive} = 4000
    options.MaxFunEvals (1,1) double {mustBeInteger,mustBePositive} = 20000
    options.TolX (1,1) double {mustBePositive} = 1e-9
    options.TolFun (1,1) double {mustBePositive} = 1e-9
    options.Verbose (1,1) logical = true
end

if isnan(options.SamplingInterval)
    error('FitHyperexponentialMLE:SamplingIntervalRequired', ...
        ['SamplingInterval is required. It is load-bearing in BOTH modes: ' ...
         'in discrete mode it is the bin width and sets n_min, and in ' ...
         'continuous mode it sets the default MaxRate, i.e. a hard floor ' ...
         'of tau >= SamplingInterval that can silently clamp a fast ' ...
         'component. Pass your acquisition interval in the SAME UNITS as ' ...
         'xmin (seconds data at 1 Hz: SamplingInterval=1; the same data ' ...
         'rescaled to minutes: SamplingInterval=1/60).']);
elseif ~(options.SamplingInterval > 0)
    error('FitHyperexponentialMLE:InvalidSamplingInterval', ...
        'SamplingInterval must be a positive scalar.');
end

if isnan(options.MaxRate)
    % tau_j >= SamplingInterval: you cannot resolve a time constant
    % shorter than the interval at which you looked. NaN is the "auto"
    % sentinel, so MaxRate carries no mustBePositive validator in the
    % arguments block (MATLAB validates defaults, and NaN would fail it)
    % and is range-checked here instead.
    options.MaxRate = 1 / options.SamplingInterval;
elseif ~(options.MaxRate > 0)
    error('FitHyperexponentialMLE:InvalidMaxRate', ...
        'MaxRate must be a positive scalar (or NaN to use 1/SamplingInterval).');
end

if ~isempty(options.RandomSeed)
    rng(options.RandomSeed);
end

isDiscrete = strcmp(options.DistributionType, "discrete");
dt = options.SamplingInterval;

% ---------------------------------------------------------------- data prep
raw = eventseries(:);
nRaw = numel(raw);
finiteMask = isfinite(raw);
nNonFinite = nnz(~finiteMask);

% Truncation must be applied on the same scale the likelihood is evaluated
% on. In discrete mode the support is round(t/dt) >= n_min, so filtering on
% the raw t >= xmin would wrongly drop observations in
% [xmin - dt/2, xmin) that round up into the support.
n_min = max(1, round(xmin / dt));
if isDiscrete
    keep = finiteMask & (round(raw / dt) >= n_min);
else
    keep = finiteMask & (raw >= xmin);
end
nBelowXmin = nnz(finiteMask & ~keep);
data = raw(keep);
n = numel(data);

if options.Verbose && (nNonFinite > 0 || nBelowXmin > 0)
    fprintf(['FitHyperexponentialMLE: dropped %d of %d input value(s) ' ...
        '(%d non-finite/NaN, %d finite but below the truncation point ' ...
        'xmin=%.4g).\n'], nNonFinite + nBelowXmin, nRaw, nNonFinite, ...
        nBelowXmin, xmin);
end

% Diagnostics are seeded with every field present from the outset, so that
% a soft return below has exactly the same shape as a successful one.
diag_ = initDiagnostics();
diag_.nNonFinite = nNonFinite;
diag_.nBelowXmin = nBelowXmin;

if n < 5
    msg = sprintf(['Fewer than 5 usable observations (n=%d) after applying ' ...
        'xmin.'], n);
    if options.ErrorOnNoValidFit
        error('FitHyperexponentialMLE:TooFewData', '%s', msg);
    end
    if options.Verbose
        fprintf('FitHyperexponentialMLE: %s Returning Failed=true.\n', msg);
    end
    H = assembleOutput(NaN, repmat(initFitStruct(), 1, options.MaxComponents), ...
        options, xmin, n, diag_, true, ['TooFewData: ' msg]);
    return
end

% ------------------------------------------------------- likelihood hazards
diag_.nAtXmin = nnz(data == xmin);
diag_.XminEqualsDataMin = (min(data) == xmin);
diag_.LooksGridded = all(abs(data/dt - round(data/dt)) < 1e-9);
diag_.UnboundedContinuousLikelihood = ~isDiscrete && diag_.nAtXmin > 0;
diag_.Recommendation = '';

% ---- grid-spacing guard -------------------------------------------------
% eventseries, xmin and SamplingInterval must share units. Rescaling the
% durations (seconds -> minutes, say) while leaving SamplingInterval at its
% default silently destroys resolution: round(t/dt) then collapses distinct
% durations onto a coarse grid, n_min lands on the wrong step, and the fit
% still returns plausible-looking numbers. Measured on one real dataset,
% minutes-valued durations left at SamplingInterval=1 discarded 92% of the
% resolution and inflated the fastest tau by 116% (99 s -> 215 s) with no
% error raised. Nothing else here can notice it: at the wrong dt the data
% genuinely are not on the dt grid, so LooksGridded is false and the "use
% discrete mode" hint above goes quiet too.
uData = unique(data);
diag_.nDistinctData = numel(uData);
if numel(uData) >= 2
    diag_.GridSpacing = min(diff(uData));
else
    diag_.GridSpacing = NaN;
end
diag_.nDistinctOnGrid = numel(unique(round(data / dt)));

g = diag_.GridSpacing;
% Only claim a grid if every observation really is a multiple of g -- for
% genuinely continuous durations the smallest gap is arbitrary and means
% nothing. The tolerance is loose enough for the float error in t/g at
% ratios of ~1e4 (~1e-12) and tight enough to reject non-multiples.
griddedAtG = isfinite(g) && g > 0 && all(abs(data/g - round(data/g)) < 1e-6);
diag_.GridMismatch = griddedAtG && abs(dt - g) > 1e-6 * g;

if isDiscrete && diag_.nDistinctOnGrid < diag_.nDistinctData
    lost = 100 * (1 - diag_.nDistinctOnGrid / diag_.nDistinctData);
    if griddedAtG
        hint = sprintf([' The durations lie on a grid of spacing %.6g, so ' ...
            'SamplingInterval=%.6g is very likely what you meant (n_min ' ...
            'would then be %d rather than %d).'], ...
            g, g, max(1, round(xmin/g)), n_min);
    else
        hint = [' The durations are not multiples of any single spacing, ' ...
            'so check what your acquisition interval actually is.'];
    end
    gridMsg = sprintf(['SamplingInterval=%.6g discards resolution: ' ...
        'round(t/SamplingInterval) collapses %d distinct durations into ' ...
        '%d (%.1f%% lost).%s Remember that eventseries, xmin and ' ...
        'SamplingInterval must share units -- rescaling the durations ' ...
        'requires rescaling SamplingInterval by the same factor.'], ...
        dt, diag_.nDistinctData, diag_.nDistinctOnGrid, lost, hint);
    diag_.Recommendation = gridMsg;
    warning('FitHyperexponentialMLE:GridMismatch', '%s', gridMsg);
elseif diag_.GridMismatch && options.Verbose
    fprintf(['FitHyperexponentialMLE: durations lie on a grid of spacing ' ...
        '%.6g but SamplingInterval=%.6g. SamplingInterval also sets the ' ...
        'MaxRate ceiling and the multistart seeding, so set it to your ' ...
        'true acquisition interval (same units as xmin).\n'], g, dt);
end

if diag_.UnboundedContinuousLikelihood
    if diag_.XminEqualsDataMin
        xminHint = ', or pass the protocol xmin instead of min(eventseries)';
    else
        xminHint = '';
    end
    rec = sprintf(['%d observation(s) sit exactly at xmin=%.6g, so the ' ...
        'CONTINUOUS truncated log-likelihood is unbounded (logL -> const ' ...
        '+ %d*log(lambda) as any one rate diverges) and no MLE exists for ' ...
        'K>=2. Fits are still returned but model orders that exploit the ' ...
        'divergence are excluded from selection. Use ' ...
        'DistributionType="discrete" (SamplingInterval=%.6g)%s.'], ...
        diag_.nAtXmin, xmin, diag_.nAtXmin, dt, ...
        xminHint);
    diag_.Recommendation = rec;
    warning('FitHyperexponentialMLE:UnboundedLikelihood', '%s', rec);
elseif ~isDiscrete && diag_.LooksGridded && options.Verbose
    fprintf(['FitHyperexponentialMLE: durations are all multiples of ' ...
        'SamplingInterval=%.6g, i.e. genuinely discrete. ' ...
        'DistributionType="discrete" is the correct model here.\n'], dt);
end

if diag_.XminEqualsDataMin && options.Verbose
    fprintf(['FitHyperexponentialMLE: xmin equals min(eventseries). xmin is ' ...
        'a pre-specified protocol threshold, not a statistic of the data; ' ...
        'estimating it from the sample guarantees ties at xmin.\n']);
end

% ------------------------------------------------------------------- fitting
allFits = repmat(initFitStruct(), 1, options.MaxComponents);
for K = 1:options.MaxComponents
    allFits(K) = fitOneK(data, xmin, n_min, K, isDiscrete, options);
end

% ------------------------------------- identifiability gating and selection
valid = false(1, options.MaxComponents);
aicc = inf(1, options.MaxComponents);
for K = 1:options.MaxComponents
    f = allFits(K);
    if ~f.Success
        continue
    end
    reasons = {};
    if (f.n - f.k - 1) <= 0
        reasons{end+1} = 'AICc undefined (n-k-1<=0)'; %#ok<AGROW>
    end
    if any(f.n * f.WeightsObserved < options.MinExpectedCount)
        reasons{end+1} = sprintf('component with expected observed count %.2f < %g', ...
            min(f.n * f.WeightsObserved), options.MinExpectedCount); %#ok<AGROW>
    end
    if K > 1
        lt = sort(log(f.Tau));
        if any(diff(lt) < options.DuplicateRateTol)
            reasons{end+1} = 'duplicate/collapsed time constants'; %#ok<AGROW>
        end
    end
    if any(f.Rates * xmin > options.MaxTruncationExponent)
        reasons{end+1} = sprintf(['component almost entirely truncated away ' ...
            '(lambda*xmin = %.1f > %.1f): its untruncated weight is ' ...
            'numerically meaningless'], max(f.Rates)*xmin, ...
            options.MaxTruncationExponent); %#ok<AGROW>
    end
    if f.AtRateBound
        reasons{end+1} = sprintf('component at the MaxRate ceiling (tau = %.6g)', ...
            1/options.MaxRate); %#ok<AGROW>
    end
    if K > 1 && allFits(K-1).Success && ...
            (f.LogLik - allFits(K-1).LogLik) < options.DeadComponentTol
        reasons{end+1} = 'dead component (no logL gain over K-1)'; %#ok<AGROW>
    end
    if ~f.CovValid
        reasons{end+1} = 'Hessian not positive definite (not locally identified)'; %#ok<AGROW>
    end

    if isempty(reasons)
        valid(K) = true;
        aicc(K) = f.AICc;
        allFits(K).Degenerate = false;
        allFits(K).DegenerateReason = '';
    else
        allFits(K).Degenerate = true;
        allFits(K).DegenerateReason = strjoin(reasons, '; ');
    end
end

if ~any(valid)
    msg = sprintf(['No model order passed the identifiability gating. %s ' ...
        'Inspect H.AllFits(K).DegenerateReason via a Verbose call, relax ' ...
        'MinExpectedCount/DuplicateRateTol, reduce MaxComponents, or ' ...
        'collect more data.'], diag_.Recommendation);
    if options.ErrorOnNoValidFit
        error('FitHyperexponentialMLE:NoValidFit', '%s', msg);
    end
    if options.Verbose
        fprintf('FitHyperexponentialMLE: %s Returning Failed=true.\n', msg);
    end
    diag_.ExcludedK = find([allFits.Success]);
    H = assembleOutput(NaN, allFits, options, xmin, n, diag_, true, ...
        ['NoValidFit: ' msg]);
    return
end

[~, selectedK] = min(aicc);
diag_.ExcludedK = find(~valid & [allFits.Success]);

H = assembleOutput(selectedK, allFits, options, xmin, n, diag_, false, '');

if options.Verbose
    fprintf('FitHyperexponentialMLE: selected K=%d by AICc (n=%d, %s).\n', ...
        selectedK, n, options.DistributionType);
    for K = 1:options.MaxComponents
        f = allFits(K);
        if f.Success
            flag = '';
            if f.Degenerate
                flag = sprintf('  [EXCLUDED: %s]', f.DegenerateReason);
            end
            fprintf('  K=%d: logL=%.4f  k=%d  AIC=%.4f  AICc=%.4f%s\n', ...
                K, f.LogLik, f.k, f.AIC, f.AICc, flag);
        else
            fprintf('  K=%d: fit failed\n', K);
        end
    end

    sel = H.Selected;
    fprintf('  Selected fit (K=%d) parameter estimates:\n', selectedK);
    if ~sel.CovValid
        fprintf('    [Hessian-based standard errors unavailable/unreliable for this fit -- point estimates only]\n');
    end
    fprintf('    (q = observed-population weight: the fraction of the n=%d bouts >= xmin\n', n);
    fprintf('     from each component. Report q, not the untruncated w.)\n');
    for j = 1:selectedK
        if sel.CovValid
            fprintf(['    Component %d: tau=%.4g (SE %.4g)   q=%.4f (SE %.4f)' ...
                '   w_untrunc=%.4g   rate=%.4g (SE %.4g)\n'], ...
                j, sel.Tau(j), sel.TauSE(j), sel.WeightsObserved(j), ...
                sel.WeightsObservedSE(j), sel.WeightsUntruncated(j), sel.Rates(j), sel.RateSE(j));
        else
            fprintf('    Component %d: tau=%.4g   q=%.4f   w_untrunc=%.4g   rate=%.4g\n', ...
                j, sel.Tau(j), sel.WeightsObserved(j), sel.WeightsUntruncated(j), sel.Rates(j));
        end
    end
end

end

% ========================================================================
function H = assembleOutput(selectedK, allFits, options, xmin, n, diag_, ...
    failed, failureReason)
%ASSEMBLEOUTPUT Build the output struct. Every return path goes through
%here so the field names AND their order are identical whether the fit
%succeeded or soft-failed -- MATLAB refuses "results(i) = H" between
%structs whose fields differ in name or order, which would otherwise break
%the batch loops that ErrorOnNoValidFit=false exists to support.
H = struct();
H.SelectedK = selectedK;
H.AllFits = allFits;
if isnan(selectedK)
    % An unfitted template rather than [], so that downstream code reading
    % H.Selected.Tau gets an empty value (which propagates as empty)
    % instead of erroring on a field reference into [].
    H.Selected = initFitStruct();
else
    H.Selected = allFits(selectedK);
end
H.DistributionType = options.DistributionType;
H.xmin = xmin;
H.n = n;
H.Failed = failed;
H.FailureReason = failureReason;
H.Diagnostics = diag_;
end

function d = initDiagnostics()
%INITDIAGNOSTICS All diagnostic fields, in a fixed order, at defaults.
d = struct('nNonFinite', 0, 'nBelowXmin', 0, 'nAtXmin', 0, ...
    'XminEqualsDataMin', false, 'LooksGridded', false, ...
    'UnboundedContinuousLikelihood', false, 'GridSpacing', NaN, ...
    'GridMismatch', false, 'nDistinctData', 0, 'nDistinctOnGrid', 0, ...
    'ExcludedK', [], 'Recommendation', '');
end

function s = initFitStruct()
s = struct('K', NaN, 'k', NaN, 'n', NaN, 'WeightsObserved', [], ...
    'WeightsUntruncated', [], 'Rates', [], 'Tau', [], 'RateSE', [], ...
    'TauSE', [], 'WeightsObservedSE', [], 'WeightsUntruncatedSE', [], ...
    'CovValid', false, 'LogLik', NaN, 'AIC', NaN, 'AICc', NaN, ...
    'PointwiseLogLik', [], ...
    'Success', false, 'Converged', false, 'Degenerate', true, ...
    'DegenerateReason', 'not fitted', 'AtRateBound', false, ...
    'ExitFlag', NaN, 'BestParamVector', []);
end

function out = fitOneK(data, xmin, n_min, K, isDiscrete, options)
PENALTY = 1e12;

out = initFitStruct();
out.K = K;
out.k = 2*K - 1;
out.n = numel(data);

nStarts = min(options.maxStarts, options.nStartsBase + options.nStartsPerComponent*K);
starts = makeInitialPoints(K, nStarts, data, xmin, options);

optOptions = optimset('Display', 'off', 'MaxIter', options.MaxIter, ...
    'MaxFunEvals', options.MaxFunEvals, 'TolX', options.TolX, 'TolFun', options.TolFun);

objfun = @(z) negLogLik(z, K, data, xmin, n_min, isDiscrete, ...
    options.SamplingInterval, options.MaxRate);

bestNegLL = Inf;
bestZ = [];
bestExitFlag = NaN;

for si = 1:size(starts,1)
    z0 = starts(si,:);
    try
        [zFit, fval, exitflag] = fminsearch(objfun, z0, optOptions);
        if isfinite(fval) && fval < bestNegLL
            bestNegLL = fval;
            bestZ = zFit;
            bestExitFlag = exitflag;
        end
    catch
        % Ignore failed starts.
    end
end

% A best value still sitting at the numerical penalty means every start
% landed outside the representable region: that is not a fit.
if isempty(bestZ) || bestNegLL >= PENALTY
    return
end

[~, weights, rates, extra] = negLogLik(bestZ, K, data, xmin, n_min, ...
    isDiscrete, options.SamplingInterval, options.MaxRate);

out.WeightsUntruncated = weights;
out.WeightsObserved = extra.WeightsObserved;
out.PointwiseLogLik = extra.PointwiseLogLik;
out.Rates = rates;
out.Tau = 1 ./ rates;
out.LogLik = -bestNegLL;
out.AIC = 2*out.k - 2*out.LogLik;
if (out.n - out.k - 1) > 0
    out.AICc = out.AIC + (2*out.k*(out.k+1)) / (out.n - out.k - 1);
else
    out.AICc = Inf;
end
out.Success = true;
out.Converged = (bestExitFlag == 1);
out.ExitFlag = bestExitFlag;
out.BestParamVector = bestZ;
out.AtRateBound = any(rates >= options.MaxRate * (1 - 1e-6));

% --- Post-hoc asymptotic standard errors, via a numerical Hessian of the
% negative log-likelihood AT the already-found MLE (the optimization
% above is untouched by any of this -- this only runs once, after). ---
rateSE = nan(1, K);
tauSE = nan(1, K);
if K == 1
    weightUntruncSE = 0;            % w == 1 identically, no uncertainty
    weightObsSE = 0;
else
    weightUntruncSE = nan(1, K);
    weightObsSE = nan(1, K);
end
covValid = false;

try
    Hess = computeNumericalHessian(objfun, bestZ);

    % Test positive-definiteness via Cholesky rather than inverting and
    % checking the diagonal: a matrix can have an all-positive diagonal
    % without being positive definite (e.g. [[1,2],[2,1]], eigenvalues
    % -1 and 3), so a diagonal-only check can pass on a Hessian that
    % isn't actually a valid covariance matrix -- exactly in the
    % borderline cases (near label-switching, poorly identified
    % components) this check exists to catch. chol tests this directly,
    % and inverting through the Cholesky factor is also the numerically
    % preferred route once PD-ness is confirmed.
    [R, cholFlag] = chol(Hess);

    if cholFlag == 0
        % Restore the state of exactly the ids being touched, via
        % onCleanup so an error in between cannot leak them either.
        % Capturing the global state with a bare "warning" and restoring
        % that does NOT work: an id already sitting at its default state
        % is not recorded in the returned struct, so the restore silently
        % leaves it switched off for the rest of the session -- the same
        % leak as turning two ids off and restoring only the first.
        wsNear = warning('query', 'MATLAB:nearlySingularMatrix');
        wsSing = warning('query', 'MATLAB:singularMatrix');
        restoreWarnings = onCleanup(@() restoreWarningStates(wsNear, wsSing));
        warning('off', 'MATLAB:nearlySingularMatrix');
        warning('off', 'MATLAB:singularMatrix');
        % R is triangular with a strictly positive diagonal because chol
        % just succeeded, so backslash is an exact triangular solve here
        % and never needs the general-purpose inv().
        Rinv = R \ eye(out.k);
        clear restoreWarnings

        covZ = Rinv * Rinv';
        dvec = diag(covZ);

        if all(isfinite(dvec)) && all(dvec > 0)
            covValid = true;
            seZ = sqrt(dvec)'; % 1 x k, one SE per z-parameter (log-rates then weight logits)

            % log_rates = z(1:K), so rate_j = exp(z_j): SE(rate_j) = rate_j*SE(z_j).
            rateSE = rates .* seZ(1:K);
            % tau_j = 1/rate_j = exp(-z_j): SE(tau_j) = tau_j*SE(z_j) (delta method,
            % same |derivative| magnitude as for rate since d(exp(-z))/dz = -exp(-z)).
            tauSE = out.Tau .* seZ(1:K);

            if K > 1
                Jw = softmaxJacobianFreeParams(weights);
                covAlpha = covZ(K+1:end, K+1:end);
                covW = Jw * covAlpha * Jw';
                weightUntruncSE = sqrt(max(diag(covW), 0))';
            end

            % q depends on BOTH the logits and the rates, so propagate it
            % through a numerical Jacobian of the full map z -> q rather
            % than the softmax Jacobian alone.
            Jq = observedWeightJacobian(@(zz) getObservedWeights(zz, K, data, ...
                xmin, n_min, isDiscrete, options.SamplingInterval, options.MaxRate), bestZ, K);
            covQ = Jq * covZ * Jq';
            weightObsSE = sqrt(max(diag(covQ), 0))';
        end
    end
catch
    % Leave SEs as NaN / covValid = false -- point estimates above are
    % unaffected, they just won't have uncertainty attached.
end

out.RateSE = rateSE;
out.TauSE = tauSE;
out.WeightsUntruncatedSE = weightUntruncSE;
out.WeightsObservedSE = weightObsSE;
out.CovValid = covValid;

% Canonical component ordering. The mixture likelihood is invariant to
% relabelling, so without this the returned order is whatever the
% optimizer happened to end on -- which makes "component 1" mean different
% things across K, across animals, and across bootstrap replicates. All
% per-component vectors are permuted together, so the pairing between
% WeightsUntruncated/WeightsObserved/Rates/Tau and their SEs is preserved.
if options.SortComponents && K > 1
    [~, ord] = sort(out.Tau, 'ascend');
    out.WeightsUntruncated = out.WeightsUntruncated(ord);
    out.WeightsObserved = out.WeightsObserved(ord);
    out.Rates = out.Rates(ord);
    out.Tau = out.Tau(ord);
    out.RateSE = out.RateSE(ord);
    out.TauSE = out.TauSE(ord);
    out.WeightsUntruncatedSE = out.WeightsUntruncatedSE(ord);
    out.WeightsObservedSE = out.WeightsObservedSE(ord);
end

end

function starts = makeInitialPoints(K, nStarts, data, xmin, options)
%MAKEINITIALPOINTS Multistart seeds in (log-rate, weight-logit) space.
%
% A uniform draw over the whole admissible log-rate box wastes most starts
% on implausibly fast components: for sleep bouts (xmin=300 s, max ~1e4 s)
% only ~46% of a uniform draw lands at tau in [50, 5000] s, so the chance
% that all K rates are simultaneously plausible falls to ~4% at K=4 --
% about two usable seeds out of fifty. Half the seeds are therefore drawn
% around the data's own mean excess duration, which is a consistent
% estimate of the mixture mean and lands in the right decade by
% construction, and one deterministic geometric ladder is always included.
k = 2*K - 1;
starts = zeros(nStarts, k);

posData = data(data > 0);
lo = log(1 / max(posData));
hi = log(min(options.MaxRate, 2 / options.SamplingInterval));
if hi <= lo
    hi = lo + 2; % Safety floor
end

% Mean excess over the truncation point: a scale for the slow bulk.
meanExcess = max(mean(data) - xmin, eps);
centre = min(max(log(1/meanExcess), lo), hi);

for si = 1:nStarts
    if si == 1
        % Deterministic ladder spanning the admissible range.
        rateGuess = linspace(lo, hi, K);
        wGuess = ones(1, K) / K;
    elseif mod(si, 2) == 0
        % Data-driven: log-rates scattered around the mean-excess scale.
        rateGuess = sort(min(max(centre + 1.5*randn(1, K), lo), hi));
        wGuess = -log(max(rand(1, K), realmin));
        wGuess = wGuess / sum(wGuess);
    else
        % Broad uniform sweep over the admissible box (original scheme).
        rateGuess = sort(lo + (hi - lo) * rand(1, K));
        wGuess = -log(max(rand(1, K), realmin));
        wGuess = wGuess / sum(wGuess);
    end

    if K == 1
        starts(si,:) = rateGuess;
    else
        logits = log(wGuess(1:K-1) ./ wGuess(K));
        starts(si,:) = [rateGuess, logits];
    end
end
end

function [negLL, weights, rates, extra] = negLogLik(z, K, data, xmin, ...
    n_min, isDiscrete, samplingInterval, maxRate)
z = z(:).';
log_rates_raw = z(1:K);

% Smooth barrier at the rate ceiling. Clamping alone would create a flat
% plateau (the same zero-gradient trap the linear-space realmin floor used
% to create for fminsearch), so the clamped density is evaluated AND a
% quadratic penalty in log-rate is added on the excess, giving the
% simplex a restoring gradient that points back into the feasible region.
% Inside the feasible region the penalty is exactly zero, so the reported
% MLE is unaffected.
log_rate_max = log(maxRate);
excess = max(0, log_rates_raw - log_rate_max);
log_rates = min(log_rates_raw, log_rate_max);
rates = exp(log_rates);

if K == 1
    log_w = 0;
    weights = 1;
else
    logits = [z(K+1:end), 0];
    log_w = logits - logsumexp_vec(logits);
    weights = exp(log_w);
end

T = data(:); % N x 1 vector

if ~isDiscrete
    % log f(t) = logsumexp_j (log_w_j + log_rate_j - rate_j * t)
    log_comp = log_w + log_rates - (T * rates); % N x K
    log_fObs = logsumexp_mat(log_comp, 2);      % N x 1

    % log S(xmin) = logsumexp_j (log_w_j - rate_j * xmin)
    log_surv_j = log_w - rates * xmin;          % 1 x K
    log_Sxmin = logsumexp_vec(log_surv_j);

    negLL = -( sum(log_fObs) - numel(T) * log_Sxmin );
else
    dt = samplingInterval;
    n_obs = round(T / dt);

    log_q = -rates * dt;                 % log(1 - p) = -lambda * dt
    log_p = log(-expm1(-rates * dt));    % Numerically stable log(p) = log(1 - exp(-lambda*dt))

    % log P(N = n) = logsumexp_j (log_w_j + (n - 1)*log_q_j + log_p_j)
    log_comp = log_w + (n_obs - 1) * log_q + log_p; % N x K
    log_fObs = logsumexp_mat(log_comp, 2);

    % log S(n_min) = log P(N >= n_min) = logsumexp_j (log_w_j + (n_min - 1)*log_q_j)
    % n_min is clamped to >= 1 by the caller: at n_min = 0 this survival
    % would exceed 1 and the "normalizer" would inflate the likelihood.
    log_surv_j = log_w + (n_min - 1) * log_q;
    log_Sxmin = logsumexp_vec(log_surv_j);

    negLL = -( sum(log_fObs) - numel(n_obs) * log_Sxmin );
end

negLL = negLL + numel(T) * sum(excess.^2);

if ~isfinite(negLL)
    negLL = 1e12; % Penalty value for numerical non-convergence
end

if nargout > 3
    % Observed-population (truncated) weights: the fraction of the bouts
    % that actually clear xmin contributed by each component. This is the
    % quantity the data identify; the untruncated weights above are its
    % exp(+lambda_j*xmin)-amplified back-transform.
    % Per-observation log-likelihood contributions, summing to LogLik.
    % Needed to compare this fit against another model family by a Vuong
    % (1989) test, which works on the pointwise log-ratios rather than the
    % totals.
    extra = struct('WeightsObserved', exp(log_surv_j - log_Sxmin), ...
        'PointwiseLogLik', log_fObs - log_Sxmin);
end
end

function q = getObservedWeights(z, K, data, xmin, n_min, isDiscrete, dt, maxRate)
%GETOBSERVEDWEIGHTS Thin wrapper returning only q, for the Jacobian below.
% Evaluated on a single observation: q depends on the parameters alone, not
% on the sample, so this avoids re-sweeping all n durations per finite
% difference.
[~, ~, ~, extra] = negLogLik(z, K, data(1), xmin, n_min, isDiscrete, dt, maxRate);
q = extra.WeightsObserved(:);
end

function J = observedWeightJacobian(qfun, z0, K)
%OBSERVEDWEIGHTJACOBIAN Central-difference dq_i/dz_m, K x k.
z0 = z0(:).';
p = numel(z0);
h = max(1e-6, 1e-6*abs(z0));
J = zeros(K, p);
for m = 1:p
    zp = z0; zp(m) = zp(m) + h(m);
    zm = z0; zm(m) = zm(m) - h(m);
    J(:,m) = (qfun(zp) - qfun(zm)) / (2*h(m));
end
end

% Helpers for numerically stable log-sum-exp calculations
function s = logsumexp_vec(x)
max_x = max(x);
if isinf(max_x)
    s = max_x;
else
    s = max_x + log(sum(exp(x - max_x)));
end
end

function s = logsumexp_mat(X, dim)
% Row-wise (or dim-wise) logsumexp. Carries the SAME isinf guard as
% logsumexp_vec: if every entry along DIM is -Inf (a single observation
% with zero density under every mixture component), max_X is -Inf and
% X - max_X would otherwise be -Inf - (-Inf) = NaN, poisoning that row's
% result -- and, one layer up, the entire summed log-likelihood -- to NaN
% instead of the correct -Inf.
max_X = max(X, [], dim);
out_size = size(max_X);
s = zeros(out_size);
isInfMask = isinf(max_X);
if any(isInfMask(:))
    s(isInfMask) = max_X(isInfMask);
end
if any(~isInfMask(:))
    idx = ~isInfMask;
    % Broadcast-safe: only recompute the finite entries.
    diffX = X - max_X;
    expX = exp(diffX);
    sumExp = sum(expX, dim);
    finiteResult = max_X + log(sumExp);
    s(idx) = finiteResult(idx);
end
end

function H = computeNumericalHessian(fun, z0)
%COMPUTENUMERICALHESSIAN Central-difference Hessian of scalar function FUN
%at point Z0. Post-hoc only: used once after fminsearch has already
%converged, to estimate the observed information matrix for asymptotic
%standard errors. Step size h ~ eps^(1/4), the standard heuristic
%balancing truncation error against floating-point round-off for
%central-difference second derivatives.
z0 = z0(:).';
p = numel(z0);
h = max(1e-4, 1e-4*abs(z0));
H = zeros(p, p);
f0 = fun(z0);

fPlus = zeros(1, p);
fMinus = zeros(1, p);
for i = 1:p
    zp = z0; zp(i) = zp(i) + h(i);
    zm = z0; zm(i) = zm(i) - h(i);
    fPlus(i) = fun(zp);
    fMinus(i) = fun(zm);
    H(i,i) = (fPlus(i) - 2*f0 + fMinus(i)) / h(i)^2;
end

for i = 1:p-1
    for j = i+1:p
        zpp = z0; zpp(i) = zpp(i)+h(i); zpp(j) = zpp(j)+h(j);
        zpm = z0; zpm(i) = zpm(i)+h(i); zpm(j) = zpm(j)-h(j);
        zmp = z0; zmp(i) = zmp(i)-h(i); zmp(j) = zmp(j)+h(j);
        zmm = z0; zmm(i) = zmm(i)-h(i); zmm(j) = zmm(j)-h(j);
        val = (fun(zpp) - fun(zpm) - fun(zmp) + fun(zmm)) / (4*h(i)*h(j));
        H(i,j) = val;
        H(j,i) = val;
    end
end
end

function restoreWarningStates(varargin)
%RESTOREWARNINGSTATES Put each captured warning id back to its own state.
for ii = 1:numel(varargin)
    w = varargin{ii};
    warning(w.state, w.identifier);
end
end

function J = softmaxJacobianFreeParams(weights)
%SOFTMAXJACOBIANFREEPARAMS d(weight_i)/d(alpha_m) for the K-1 free logits
%alpha (weights = softmax([alpha, 0]), the same fixed-last-logit
%parametrization used throughout this file). Standard softmax Jacobian
%w_i*(delta_im - w_m), restricted to the free coordinates, used to
%delta-method-propagate log-rate/logit covariance into weight standard
%errors.
K = numel(weights);
J = zeros(K, K-1);
for m = 1:K-1
    for i = 1:K
        J(i,m) = weights(i) * ((i==m) - weights(m));
    end
end
end
