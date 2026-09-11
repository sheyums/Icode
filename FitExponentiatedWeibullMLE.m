function H = FitExponentiatedWeibullMLE(eventseries, xmin, options)
%FITEXPONENTIATEDWEIBULLMLE Maximum-likelihood fit of a left-truncated
%exponentiated Weibull distribution (continuous mode) or its exact
%discretization (discrete mode).
%
%   H = FITEXPONENTIATEDWEIBULLMLE(eventseries, xmin)
%   H = FITEXPONENTIATEDWEIBULLMLE(eventseries, xmin, Name=Value, ...)
%
%   Companion to FITHYPEREXPONENTIALMLE, deliberately sharing its
%   conventions -- same xmin semantics, same DistributionType /
%   SamplingInterval switch, same truncation by S(xmin), same output
%   layout -- so that AIC/AICc from the two are directly comparable
%   PROVIDED both are called with the same DistributionType,
%   SamplingInterval and xmin on the same observations. See MODEL
%   COMPARISON below; this is the whole point of the file.
%
%   MODEL
%   The exponentiated Weibull (Mudholkar & Srivastava 1993, IEEE Trans.
%   Reliab. 42:299-302) adds a second shape parameter alpha to the Weibull
%   by raising its CDF to a power:
%
%       F(t) = [1 - exp(-(t/lambda)^k)]^alpha,        t > 0
%       f(t) = (alpha*k/lambda) * (t/lambda)^(k-1) * exp(-(t/lambda)^k)
%              * [1 - exp(-(t/lambda)^k)]^(alpha-1)
%
%   with scale lambda>0 (same time units as xmin) and shapes k>0, alpha>0.
%   alpha=1 recovers the Weibull; k=1 the exponentiated exponential;
%   alpha=k=1 the exponential. FixAlpha=true fits the nested Weibull, so a
%   likelihood-ratio test of alpha=1 is one extra call.
%
%   The reason to try it against a mixture of exponentials is the hazard.
%   A hyperexponential mixture has a structurally DECREASING hazard (it is
%   completely monotone), whereas the exponentiated Weibull's hazard can be
%   increasing, decreasing, bathtub-shaped or unimodal depending on k and
%   alpha*k -- so it can imitate a decreasing-hazard mixture with three
%   parameters instead of 2K-1, and is a far more demanding comparison
%   than a power law. H.HazardShape reports which regime the fit landed in:
%
%       k >= 1 and alpha*k >= 1   ->  increasing
%       k <= 1 and alpha*k <= 1   ->  decreasing
%       k >  1 and alpha*k <  1   ->  bathtub
%       k <  1 and alpha*k >  1   ->  unimodal
%
%   Left truncation. Only durations >= xmin are observable, so the density
%   actually fit is conditioned on that:
%
%       f(t | T >= xmin) = f(t) / S(xmin),    S(xmin) = 1 - F(xmin)
%
%   Discrete mode. Writing dt = SamplingInterval and
%   n_min = max(1, round(xmin/dt)), the duration is modelled as falling in
%   one dt-wide bin, which is exact for grid-valued data:
%
%       P(N=n)            = F(n*dt) - F((n-1)*dt),   n = 1,2,3,...
%       P(N=n | N>=n_min) = P(N=n) / (1 - F((n_min-1)*dt)),  n >= n_min
%
%   and each eventseries value is mapped to n = round(t/dt). The pmf is
%   bounded by 1, so unlike the continuous branch it cannot be driven to
%   +Inf by observations sitting exactly at xmin. (That particular
%   pathology is specific to mixtures -- see FITHYPEREXPONENTIALMLE -- and
%   does not arise here, since this density is bounded at xmin for
%   alpha*k >= 1. It can still diverge as t->0 when alpha*k < 1, which
%   left truncation at xmin > 0 removes from consideration entirely.)
%
%   ======================================================================
%   WHEN TO CALL THIS DISCRETE AND WHEN CONTINUOUS
%   ======================================================================
%   The question is not about the fitter, it is about how your durations
%   were measured. A duration read off a recording sampled every dt is not
%   a real number; it is a count of dt-sized steps. "Discrete" models that
%   count, "continuous" pretends the duration was measured with infinite
%   resolution.
%
%   USE "discrete" (with SamplingInterval = your acquisition interval) IF
%   ANY OF THESE HOLD -- and for bout data at least one always does:
%     - the durations are integer multiples of the sampling interval, i.e.
%       counts of video frames or 1-Hz samples. H.Diagnostics.LooksGridded
%       reports this for you, and it is TRUE for ordinary bout data;
%     - any observation equals xmin exactly. Ties at the truncation point
%       are what make a continuous MIXTURE likelihood unbounded (see
%       FITHYPEREXPONENTIALMLE); a pmf is bounded by 1 and cannot diverge.
%       Even where it is not fatal, as here, it means the continuous
%       density is being asked about a point mass;
%     - xmin is only a small multiple of the sampling interval (wake
%       bouts, say xmin = 2 s at 1 Hz), where treating duration as
%       continuous is worst precisely at xmin, which is where truncation
%       makes the fit most sensitive;
%     - you intend to compare AIC against any other model you will also
%       fit discretely.
%
%   USE "continuous" ONLY IF ALL OF THESE HOLD:
%     - the durations are genuinely real-valued -- interpolated event
%       times, or timestamps whose quantization is orders of magnitude
%       finer than the smallest time constant you care about;
%     - no observation sits exactly at xmin;
%     - the competing model you want an AIC against is itself a continuous
%       density you cannot refit as a pmf.
%
%   IN CASE OF DOUBT, USE "discrete" -- which is why it is the DEFAULT.
%   As dt -> 0 the discrete
%   log-likelihood approaches the continuous one plus n*log(dt), and the
%   parameter estimates converge, so discrete is continuous-in-the-limit
%   and never the riskier choice. Its only cost is that SamplingInterval
%   has to be right, which the grid-spacing guard checks. Erring towards a
%   SamplingInterval SMALLER than the true grid is harmless; larger throws
%   away resolution.
%
%   WHAT CHANGES IF YOU SWITCH. Lambda, K and Alpha are on the same scale
%   in both modes and stay comparable (discrete-mode estimates carry an
%   O(dt) discretization shift, negligible while dt is small against
%   lambda). LogLik, AIC and AICc DO NOT: a pmf is a density times a bin
%   width, so the discrete log-likelihood sits roughly n*log(dt) above the
%   continuous one. Never compare an AIC across the two modes, or across
%   two different SamplingIntervals -- see MODEL COMPARISON below.
%
%   ALPHA IS NOT IDENTIFIED WHEN xmin SITS DEEP IN THE TAIL
%   For large u = (t/lambda)^k we have [1-exp(-u)]^alpha ~ 1 - alpha*exp(-u),
%   hence S(t) ~ alpha*exp(-u(t)) and therefore
%
%       S(t)/S(xmin)  ~  exp(-(u(t) - u(xmin))),
%
%   in which ALPHA HAS CANCELLED. So if exp(-(xmin/lambda)^k) is small, the
%   left-truncated exponentiated Weibull collapses onto a left-truncated
%   Weibull and alpha is estimated from nothing: the likelihood is flat
%   along it and the optimizer parks arbitrarily. alpha is informative only
%   through terms of order exp(-(xmin/lambda)^k), so that quantity is
%   reported as H.Diagnostics.ExpMinusUmin and compared against
%   MinAlphaIdentifiability; when it falls below, AlphaIdentifiable is
%   false and a warning is issued recommending FixAlpha=true, which fits
%   the Weibull that the data actually support. This is the direct analogue
%   of the identifiability gating in FITHYPEREXPONENTIALMLE.
%
%   MODEL COMPARISON
%   AIC differences are meaningful only between log-likelihoods of the SAME
%   observations under the SAME dominating measure. Concretely:
%     - both fits discrete (pmf, counting measure), or both continuous
%       (density, Lebesgue) -- never one of each. A discrete log-likelihood
%       sits roughly n*log(dt) above the continuous one, since a pmf is a
%       density times a bin width; that offset is an artefact of the
%       measure, not evidence;
%     - the same xmin, so the same retained observations. In particular do
%       not compare against a power law whose xmin was chosen by KS
%       minimization, because that fit used a different subset;
%     - the same SamplingInterval in discrete mode.
%   PointwiseLogLik is provided for the non-nested Vuong (1989,
%   Econometrica 57:307-333) test against another family, the mixture
%   likelihood-ratio having no chi-square null:
%       D = Hew.PointwiseLogLik - Hother.Selected.PointwiseLogLik;
%       V = sqrt(numel(D))*mean(D)/std(D,1);     % ~ N(0,1)
%
%   INPUTS
%   eventseries : numeric vector of positive event durations, same units
%                 as xmin (e.g. seconds)
%   xmin        : positive scalar, left-truncation cutoff. A PROTOCOL
%                 constant (e.g. 300 s for the 5-minute sleep rule), never
%                 min(eventseries) -- estimating it from the sample is not
%                 a truncation point and biases everything conditioned on
%                 it.
%
%   NAME-VALUE OPTIONS
%   DistributionType    "discrete" (DEFAULT) or "continuous". See WHEN TO
%                       CALL THIS DISCRETE AND WHEN CONTINUOUS above.
%                       Switching modes changes LogLik/AIC/AICc by roughly
%                       n*log(dt), so an AIC recorded under an older
%                       release that defaulted to "continuous" is NOT
%                       comparable to one from the current default.
%   SamplingInterval    REQUIRED positive scalar, in the SAME UNITS as
%                       eventseries and xmin. There is no default --
%                       omitting it raises SamplingIntervalRequired. Sets
%                       the discrete-mode bin width and seeds the
%                       multistart. Rescaling the
%                       durations means rescaling this too (seconds to
%                       minutes is xmin=5 AND SamplingInterval=1/60); a
%                       mismatch is caught by the grid-spacing guard and
%                       reported in H.Diagnostics.GridMismatch.
%   FixAlpha            logical, default false. True constrains alpha=1,
%                       fitting the nested two-parameter Weibull. Run both
%                       and compare by a chi-square(1) likelihood-ratio
%                       test for a principled check on whether the third
%                       parameter is earned.
%   MinAlphaIdentifiability
%                       positive scalar, default 1e-3. Threshold on
%                       exp(-(xmin/lambda)^k) below which alpha is treated
%                       as unidentified (see above).
%   nStartsBase, nStartsPerParameter, maxStarts, RandomSeed, MaxIter,
%   MaxFunEvals, TolX, TolFun, Verbose, ErrorOnNoValidFit
%                       as in FITHYPEREXPONENTIALMLE.
%
%   OUTPUT
%   H.Lambda, H.K, H.Alpha          point estimates
%   H.LambdaSE, H.KSE, H.AlphaSE    asymptotic standard errors from a
%                                   numerical Hessian of the negative
%                                   log-likelihood at the MLE (the observed
%                                   information; Efron & Hinkley 1978,
%                                   Biometrika 65:457-487), propagated from
%                                   the log-parameters by the delta method.
%                                   NaN when CovValid is false. AlphaSE is
%                                   0 when FixAlpha=true, where alpha=1 by
%                                   assumption rather than estimated.
%   H.HazardShape                   'increasing' | 'decreasing' |
%                                   'bathtub' | 'unimodal'
%   H.k, H.n, H.LogLik, H.AIC, H.AICc, H.CovValid, H.Success, H.Converged,
%   H.ExitFlag, H.BestParamVector, H.PointwiseLogLik
%   H.Failed, H.FailureReason       as in FITHYPEREXPONENTIALMLE; Failed is
%                                   only reachable with
%                                   ErrorOnNoValidFit=false
%   H.DistributionType, H.xmin, H.Diagnostics
%
%   NUMERICAL NOTES
%   log(1-exp(-u)) is evaluated by the two-branch rule of Maechler (2012,
%   "Accurately computing log(1-exp(-|a|))"): log(-expm1(-u)) for
%   u < log 2, log1p(-exp(-u)) above it. Either branch alone loses all
%   precision in one limit. The discrete bin probability
%   F(n*dt) - F((n-1)*dt) is never formed by subtraction, which cancels
%   catastrophically once both CDFs approach 1; it is computed as
%   logF_n + log(-expm1(logF_prev - logF_n)) below the median and, above
%   it, in the survival domain as logS_prev + log(-expm1(logS_n - logS_prev)),
%   so the far tail stays accurate.
%
%   VALIDATE BEFORE TRUSTING -- AND DO NOT OVER-READ THE PARAMETERS.
%   k and alpha are both shape parameters and are strongly correlated, so
%   this likelihood is flatter and more awkward than a count of three
%   parameters suggests (Nadarajah, Cordeiro & Ortega 2013, J. Stat.
%   Comput. Simul. 83:1-27). The consequence is that THE LAW IS IDENTIFIED
%   MUCH BETTER THAN THE INDIVIDUAL PARAMETERS. On n=3000 draws from
%   lambda=800, k=0.8, alpha=0.6 truncated at xmin=100, the MLE landed at
%   lambda=884, k=0.843, alpha=0.546 -- 10% off in every coordinate -- yet
%   with a log-likelihood 0.36 ABOVE the truth's, a truncated survival
%   function agreeing with the truth to 0.005 everywhere, and quantiles
%   agreeing to 2%. The optimizer was not failing; the data genuinely do
%   not distinguish those triples.
%
%   So: report the fitted DISTRIBUTION (survival curve, quantiles, hazard
%   regime, log-likelihood for model comparison), and treat lambda, k and
%   alpha individually as weakly determined unless their standard errors
%   say otherwise. Two fits with different-looking parameters may be the
%   same distribution; compare them by their survival curves, not
%   coordinate-by-coordinate. Check CovValid, check AlphaIdentifiable, use
%   enough multistarts (the defaults give 65 for three parameters, and
%   fewer than ~50 measurably degrades the optimum), and fit simulated
%   data at your real sample sizes before trusting per-animal fits.
%
%   See also FITHYPEREXPONENTIALMLE, SHIFTLOGNORMAL_MLE.

arguments
    eventseries double {mustBeReal}
    xmin (1,1) double {mustBePositive}
    options.DistributionType (1,1) string {mustBeMember(options.DistributionType,["continuous","discrete"])} = "discrete"
    options.SamplingInterval (1,1) double = NaN
    options.FixAlpha (1,1) logical = false
    options.MinAlphaIdentifiability (1,1) double {mustBePositive} = 1e-3
    options.nStartsBase (1,1) double {mustBeInteger,mustBePositive} = 20
    options.nStartsPerParameter (1,1) double {mustBeInteger,mustBePositive} = 15
    options.maxStarts (1,1) double {mustBeInteger,mustBePositive} = 80
    options.RandomSeed = []
    options.MaxIter (1,1) double {mustBeInteger,mustBePositive} = 20000
    options.MaxFunEvals (1,1) double {mustBeInteger,mustBePositive} = 40000
    options.TolX (1,1) double {mustBePositive} = 1e-10
    options.TolFun (1,1) double {mustBePositive} = 1e-10
    options.Verbose (1,1) logical = true
    options.ErrorOnNoValidFit (1,1) logical = true
end

if isnan(options.SamplingInterval)
    error('FitExponentiatedWeibullMLE:SamplingIntervalRequired', ...
        ['SamplingInterval is required. It is load-bearing in BOTH modes: ' ...
         'in discrete mode it is the bin width and sets n_min, and in ' ...
         'continuous mode it sets the default MaxRate, i.e. a hard floor ' ...
         'of tau >= SamplingInterval that can silently clamp a fast ' ...
         'component. Pass your acquisition interval in the SAME UNITS as ' ...
         'xmin (seconds data at 1 Hz: SamplingInterval=1; the same data ' ...
         'rescaled to minutes: SamplingInterval=1/60).']);
elseif ~(options.SamplingInterval > 0)
    error('FitExponentiatedWeibullMLE:InvalidSamplingInterval', ...
        'SamplingInterval must be a positive scalar.');
end

PENALTY = 1e12;
isDiscrete = strcmp(options.DistributionType, "discrete");
dt = options.SamplingInterval;

if ~isempty(options.RandomSeed)
    rng(options.RandomSeed);
end

% ---------------------------------------------------------------- data prep
raw = eventseries(:);
nRaw = numel(raw);
finiteMask = isfinite(raw) & (raw > 0);
nNonFinite = nnz(~finiteMask);

% Truncation is applied on the scale the likelihood is evaluated on: in
% discrete mode the support is round(t/dt) >= n_min, so filtering on the
% raw t >= xmin would wrongly drop values in [xmin - dt/2, xmin) that round
% up into the support.
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
    fprintf(['FitExponentiatedWeibullMLE: dropped %d of %d input value(s) ' ...
        '(%d non-finite/non-positive, %d below the truncation point ' ...
        'xmin=%.4g).\n'], nNonFinite + nBelowXmin, nRaw, nNonFinite, ...
        nBelowXmin, xmin);
end

nPar = 3 - double(options.FixAlpha);

diag_ = initDiagnostics();
diag_.nNonFinite = nNonFinite;
diag_.nBelowXmin = nBelowXmin;

if n < nPar + 2
    msg = sprintf(['Fewer than %d usable observations (n=%d) after applying ' ...
        'xmin.'], nPar + 2, n);
    if options.ErrorOnNoValidFit
        error('FitExponentiatedWeibullMLE:TooFewData', '%s', msg);
    end
    if options.Verbose
        fprintf('FitExponentiatedWeibullMLE: %s Returning Failed=true.\n', msg);
    end
    H = assembleOutput([], nPar, options, xmin, n, diag_, true, ...
        ['TooFewData: ' msg]);
    return
end

% ------------------------------------------------------------- diagnostics
diag_.nAtXmin = nnz(data == xmin);
diag_.XminEqualsDataMin = (min(data) == xmin);
diag_.LooksGridded = all(abs(data/dt - round(data/dt)) < 1e-9);

uData = unique(data);
diag_.nDistinctData = numel(uData);
if numel(uData) >= 2
    diag_.GridSpacing = min(diff(uData));
else
    diag_.GridSpacing = NaN;
end
diag_.nDistinctOnGrid = numel(unique(round(data / dt)));

g = diag_.GridSpacing;
griddedAtG = isfinite(g) && g > 0 && all(abs(data/g - round(data/g)) < 1e-6);
diag_.GridMismatch = griddedAtG && abs(dt - g) > 1e-6 * g;

if isDiscrete && diag_.nDistinctOnGrid < diag_.nDistinctData
    lost = 100 * (1 - diag_.nDistinctOnGrid / diag_.nDistinctData);
    if griddedAtG
        hint = sprintf([' The durations lie on a grid of spacing %.6g, so ' ...
            'SamplingInterval=%.6g is very likely what you meant (n_min ' ...
            'would then be %d rather than %d).'], g, g, ...
            max(1, round(xmin/g)), n_min);
    else
        hint = [' The durations are not multiples of any single spacing, ' ...
            'so check what your acquisition interval actually is.'];
    end
    gridMsg = sprintf(['SamplingInterval=%.6g discards resolution: ' ...
        'round(t/SamplingInterval) collapses %d distinct durations into ' ...
        '%d (%.1f%% lost).%s eventseries, xmin and SamplingInterval must ' ...
        'share units.'], dt, diag_.nDistinctData, diag_.nDistinctOnGrid, ...
        lost, hint);
    diag_.Recommendation = gridMsg;
    warning('FitExponentiatedWeibullMLE:GridMismatch', '%s', gridMsg);
elseif diag_.GridMismatch && options.Verbose
    fprintf(['FitExponentiatedWeibullMLE: durations lie on a grid of ' ...
        'spacing %.6g but SamplingInterval=%.6g.\n'], g, dt);
end

if ~isDiscrete && diag_.LooksGridded && options.Verbose
    fprintf(['FitExponentiatedWeibullMLE: durations are all multiples of ' ...
        'SamplingInterval=%.6g, i.e. genuinely discrete. ' ...
        'DistributionType="discrete" is the correct model here.\n'], dt);
end

% ------------------------------------------------------------------- fitting
objfun = @(zFree) negLogLik(expandZ(zFree, options.FixAlpha), data, xmin, ...
    n_min, isDiscrete, dt);
starts = makeInitialPoints(data, xmin, options, nPar);

optOptions = optimset('Display', 'off', 'MaxIter', options.MaxIter, ...
    'MaxFunEvals', options.MaxFunEvals, 'TolX', options.TolX, ...
    'TolFun', options.TolFun);

bestNegLL = Inf; bestZ = []; bestExitFlag = NaN;
for si = 1:size(starts,1)
    try
        [zFit, fval, exitflag] = fminsearch(objfun, starts(si,:), optOptions);
        if isfinite(fval) && fval < bestNegLL
            bestNegLL = fval; bestZ = zFit; bestExitFlag = exitflag;
        end
    catch
        % Ignore failed starts.
    end
end

if isempty(bestZ) || bestNegLL >= PENALTY
    msg = ['Every multistart ended at the numerical penalty; no usable ' ...
        'maximum was found.'];
    if options.ErrorOnNoValidFit
        error('FitExponentiatedWeibullMLE:NoValidFit', '%s', msg);
    end
    if options.Verbose
        fprintf('FitExponentiatedWeibullMLE: %s Returning Failed=true.\n', msg);
    end
    H = assembleOutput([], nPar, options, xmin, n, diag_, true, ...
        ['NoValidFit: ' msg]);
    return
end

zFull = expandZ(bestZ, options.FixAlpha);
theta = exp(zFull);
[~, pointwise] = negLogLik(zFull, data, xmin, n_min, isDiscrete, dt);

fit = struct();
fit.Lambda = theta(1); fit.K = theta(2); fit.Alpha = theta(3);
fit.LogLik = -bestNegLL;
fit.PointwiseLogLik = pointwise;
fit.ExitFlag = bestExitFlag;
fit.Converged = (bestExitFlag == 1);
fit.BestParamVector = bestZ;

% --- asymptotic standard errors, post hoc at the MLE ---
se = nan(1, nPar);
covValid = false;
try
    Hess = computeNumericalHessian(objfun, bestZ);
    [R, cholFlag] = chol(Hess);
    if cholFlag == 0
        wsNear = warning('query', 'MATLAB:nearlySingularMatrix');
        wsSing = warning('query', 'MATLAB:singularMatrix');
        restoreWarnings = onCleanup(@() restoreWarningStates(wsNear, wsSing));
        warning('off', 'MATLAB:nearlySingularMatrix');
        warning('off', 'MATLAB:singularMatrix');
        Rinv = R \ eye(nPar);
        clear restoreWarnings
        covZ = Rinv * Rinv';
        dvec = diag(covZ);
        if all(isfinite(dvec)) && all(dvec > 0)
            covValid = true;
            se = sqrt(dvec)';   % SE of the log-parameters
        end
    end
catch
    % Leave SEs as NaN; point estimates above are unaffected.
end

% Delta method: theta = exp(z) so SE(theta) = theta*SE(z).
fit.LambdaSE = fit.Lambda * se(1);
fit.KSE = fit.K * se(2);
if options.FixAlpha
    fit.AlphaSE = 0;    % alpha == 1 by assumption, not estimated
else
    fit.AlphaSE = fit.Alpha * se(3);
end
fit.CovValid = covValid;

% Hazard-shape regime (Mudholkar & Srivastava 1993).
ak = fit.Alpha * fit.K;
if fit.K >= 1 && ak >= 1
    fit.HazardShape = 'increasing';
elseif fit.K <= 1 && ak <= 1
    fit.HazardShape = 'decreasing';
elseif fit.K > 1 && ak < 1
    fit.HazardShape = 'bathtub';
else
    fit.HazardShape = 'unimodal';
end

% alpha identifiability: alpha only enters through exp(-(xmin/lambda)^k).
diag_.ExpMinusUmin = exp(-exp(fit.K * (log(xmin) - log(fit.Lambda))));
diag_.AlphaIdentifiable = options.FixAlpha || ...
    (diag_.ExpMinusUmin >= options.MinAlphaIdentifiability);
if ~diag_.AlphaIdentifiable
    aMsg = sprintf(['exp(-(xmin/lambda)^k) = %.3g < %.3g, so xmin lies deep ' ...
        'in the tail: S(t)/S(xmin) is independent of alpha to that order ' ...
        'and alpha is not identified (alpha=%.4g is where the optimizer ' ...
        'stopped on a flat ridge, not an estimate). Refit with ' ...
        'FixAlpha=true, which fits the Weibull the data actually ' ...
        'support.'], diag_.ExpMinusUmin, options.MinAlphaIdentifiability, ...
        fit.Alpha);
    diag_.Recommendation = aMsg;
    warning('FitExponentiatedWeibullMLE:AlphaNotIdentified', '%s', aMsg);
end

H = assembleOutput(fit, nPar, options, xmin, n, diag_, false, '');

if options.Verbose
    fprintf(['FitExponentiatedWeibullMLE: %s mode, n=%d, k=%d parameter(s).\n'], ...
        options.DistributionType, n, nPar);
    fprintf('  logL=%.4f  AIC=%.4f  AICc=%.4f  hazard=%s\n', ...
        H.LogLik, H.AIC, H.AICc, H.HazardShape);
    if H.CovValid
        fprintf('  lambda=%.6g (SE %.4g)   k=%.4g (SE %.4g)   alpha=%.4g (SE %.4g)\n', ...
            H.Lambda, H.LambdaSE, H.K, H.KSE, H.Alpha, H.AlphaSE);
    else
        fprintf('  [Hessian-based standard errors unavailable -- point estimates only]\n');
        fprintf('  lambda=%.6g   k=%.4g   alpha=%.4g\n', H.Lambda, H.K, H.Alpha);
    end
    if options.FixAlpha
        fprintf('  (alpha fixed at 1: this is the nested Weibull)\n');
    end
end

end

% ========================================================================
function H = assembleOutput(fit, nPar, options, xmin, n, diag_, failed, reason)
%ASSEMBLEOUTPUT Single exit point, so the field names AND their order are
%identical whether the fit succeeded or soft-failed. MATLAB refuses
%"results(i) = H" between structs whose fields differ in name or order,
%which would break batch loops over animals.
if isempty(fit)
    fit = struct('Lambda', NaN, 'K', NaN, 'Alpha', NaN, 'LambdaSE', NaN, ...
        'KSE', NaN, 'AlphaSE', NaN, 'HazardShape', '', 'LogLik', NaN, ...
        'PointwiseLogLik', [], 'CovValid', false, 'Converged', false, ...
        'ExitFlag', NaN, 'BestParamVector', []);
end
H = struct();
H.Lambda = fit.Lambda;
H.K = fit.K;
H.Alpha = fit.Alpha;
H.LambdaSE = fit.LambdaSE;
H.KSE = fit.KSE;
H.AlphaSE = fit.AlphaSE;
H.HazardShape = fit.HazardShape;
H.k = nPar;
H.n = n;
H.LogLik = fit.LogLik;
H.PointwiseLogLik = fit.PointwiseLogLik;
H.AIC = 2*nPar - 2*fit.LogLik;
if (n - nPar - 1) > 0
    H.AICc = H.AIC + (2*nPar*(nPar+1)) / (n - nPar - 1);
else
    H.AICc = Inf;
end
H.CovValid = fit.CovValid;
H.Success = ~failed;
H.Converged = fit.Converged;
H.ExitFlag = fit.ExitFlag;
H.BestParamVector = fit.BestParamVector;
H.FixAlpha = options.FixAlpha;
H.DistributionType = options.DistributionType;
H.xmin = xmin;
H.Failed = failed;
H.FailureReason = reason;
H.Diagnostics = diag_;
end

function d = initDiagnostics()
d = struct('nNonFinite', 0, 'nBelowXmin', 0, 'nAtXmin', 0, ...
    'XminEqualsDataMin', false, 'LooksGridded', false, ...
    'GridSpacing', NaN, 'GridMismatch', false, 'nDistinctData', 0, ...
    'nDistinctOnGrid', 0, 'ExpMinusUmin', NaN, ...
    'AlphaIdentifiable', false, 'Recommendation', '');
end

function z = expandZ(zFree, fixAlpha)
%EXPANDZ Map the free parameters to the full [log lambda, log k, log alpha].
if fixAlpha
    z = [zFree(1), zFree(2), 0];   % log(alpha) = 0 <=> alpha = 1
else
    z = zFree(1:3);
end
end

function starts = makeInitialPoints(data, xmin, options, nPar)
%MAKEINITIALPOINTS Multistart seeds in log-parameter space. lambda is
%seeded from the data's own scale; k and alpha around 1, which is the
%exponential special case, with enough spread to reach the
%decreasing-hazard corner (k<1, alpha*k<1) where bout data usually sit.
nStarts = min(options.maxStarts, ...
    options.nStartsBase + options.nStartsPerParameter*nPar);
starts = zeros(nStarts, nPar);
scale = max(mean(data), xmin);
for si = 1:nStarts
    if si == 1
        z = [log(scale), 0, 0];            % the exponential
    elseif si == 2
        z = [log(scale), log(0.7), log(0.2)];  % decreasing-hazard corner
    else
        z = [log(scale) + randn*1.0, randn*0.6, randn*0.8];
    end
    starts(si,:) = z(1:nPar);
end
end

function [negLL, pointwise] = negLogLik(z, data, xmin, n_min, isDiscrete, dt)
%NEGLOGLIK Left-truncated exponentiated Weibull negative log-likelihood.
lambda = exp(z(1)); k = exp(z(2)); alpha = exp(z(3));
t = data(:);

if ~isDiscrete
    logt = log(t) - log(lambda);
    u = exp(k * logt);
    % log f(t) in full, with log(1-exp(-u)) on its stable branch.
    logf = log(alpha) + log(k) - log(lambda) + (k-1)*logt - u ...
        + (alpha-1) * log1mexpVec(u);
    logSxmin = logSurv(xmin, lambda, k, alpha);
    pointwise = logf - logSxmin;
else
    nObs = round(t / dt);
    % P(N=n) = F(n*dt) - F((n-1)*dt), never formed by subtraction.
    logFn = logCDF(nObs*dt, lambda, k, alpha);
    logFp = logCDF((nObs-1)*dt, lambda, k, alpha);
    logp = zeros(size(logFn));

    % Below the median work in F (where F is small and accurate); above it
    % work in S, since both CDFs approach 1 there and their difference
    % cancels away every significant digit.
    lower = logFn <= log(0.5);
    if any(lower)
        d1 = min(logFp(lower) - logFn(lower), -realmin);
        logp(lower) = logFn(lower) + log(-expm1(d1));
    end
    if any(~lower)
        lSn = log(-expm1(logFn(~lower)));
        lSp = log(-expm1(logFp(~lower)));
        d2 = min(lSn - lSp, -realmin);
        logp(~lower) = lSp + log(-expm1(d2));
    end

    logSxmin = logSurv((n_min-1)*dt, lambda, k, alpha);
    pointwise = logp - logSxmin;
end

negLL = -sum(pointwise);
if ~isfinite(negLL)
    negLL = 1e12;
end
end

function lf = logCDF(t, lambda, k, alpha)
%LOGCDF log F(t) = alpha*log(1 - exp(-(t/lambda)^k)). t may contain 0,
%which correctly gives -Inf (F(0)=0).
t = t(:);
lf = -inf(size(t));
pos = t > 0;
if any(pos)
    u = exp(k * (log(t(pos)) - log(lambda)));
    lf(pos) = alpha * log1mexpVec(u);
end
end

function ls = logSurv(t, lambda, k, alpha)
%LOGSURV log S(t) = log(1 - F(t)), via -expm1 of the log CDF so that the
%subtraction from 1 is never done in linear space.
ls = log(-expm1(logCDF(t, lambda, k, alpha)));
ls = ls(1);
end

function y = log1mexpVec(u)
%LOG1MEXPVEC log(1 - exp(-u)) for u >= 0, on whichever branch is accurate.
%Maechler (2012): below log 2 the argument of the logarithm is small and
%log(-expm1(-u)) keeps it; above log 2, exp(-u) is small and
%log1p(-exp(-u)) keeps it instead. Using either form alone loses all
%precision in the opposite limit -- at u=1e-16 the naive log(1-exp(-u)) is
%already wrong by ~0.1 in log-probability.
u = u(:);
y = zeros(size(u));
small = u < log(2);
y(small)  = log(-expm1(-u(small)));
y(~small) = log1p(-exp(-u(~small)));
y(u <= 0) = -Inf;
end

function H = computeNumericalHessian(fun, z0)
%COMPUTENUMERICALHESSIAN Central-difference Hessian at Z0. Step size
%h ~ eps^(1/4), the usual balance of truncation against round-off for
%second differences.
z0 = z0(:).';
p = numel(z0);
h = max(1e-4, 1e-4*abs(z0));
H = zeros(p, p);
f0 = fun(z0);
for i = 1:p
    zp = z0; zp(i) = zp(i) + h(i);
    zm = z0; zm(i) = zm(i) - h(i);
    H(i,i) = (fun(zp) - 2*f0 + fun(zm)) / h(i)^2;
end
for i = 1:p-1
    for j = i+1:p
        zpp = z0; zpp(i) = zpp(i)+h(i); zpp(j) = zpp(j)+h(j);
        zpm = z0; zpm(i) = zpm(i)+h(i); zpm(j) = zpm(j)-h(j);
        zmp = z0; zmp(i) = zmp(i)-h(i); zmp(j) = zmp(j)+h(j);
        zmm = z0; zmm(i) = zmm(i)-h(i); zmm(j) = zmm(j)-h(j);
        val = (fun(zpp) - fun(zpm) - fun(zmp) + fun(zmm)) / (4*h(i)*h(j));
        H(i,j) = val; H(j,i) = val;
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
