function H = FitTruncatedDiscreteMLE(eventseries, xmin, modelName, options)
%FITTRUNCATEDDISCRETEMLE Maximum-likelihood fit of a left-truncated
%distribution, discretized onto the sampling grid (default) or treated as
%continuous, for any family in the model registry below.
%
%   H = FITTRUNCATEDDISCRETEMLE(eventseries, xmin, modelName, Name=Value, ...)
%
%   This is the shared engine behind FITGAMMAMLE, FITERLANGMLE,
%   FITCHISQUAREDMLE, FITPEARSON3MLE, FITWEIBULLMLE and FITBETAMLE. It
%   exists so that the truncation, the discretization, the cancellation-safe
%   tail arithmetic, the multistart, the Hessian standard errors, the
%   grid-spacing guard and the soft-return are written ONCE. Adding a family
%   means adding a CDF and a start rule to the registry, and it inherits
%   every guard and every test.
%
%   FITHYPEREXPONENTIALMLE and FITEXPONENTIATEDWEIBULLMLE predate this
%   engine and remain standalone (a mixture needs order selection and
%   identifiability gating that do not generalise), but they use the same
%   conventions, so AIC/AICc/BIC from all of them are comparable -- see
%   MODEL COMPARISON.
%
%   MODEL
%   Every family is specified by its CDF F(t; theta). Writing
%   dt = SamplingInterval and n_min = max(1, round(xmin/dt)), the discrete
%   model says a duration falls in one dt-wide bin, and conditions on
%   clearing the truncation point:
%
%       P(N=n)            = F(n*dt) - F((n-1)*dt),          n = 1,2,3,...
%       P(N=n | N>=n_min) = P(N=n) / (1 - F((n_min-1)*dt)),  n >= n_min
%
%   with each observation mapped to n = round(t/dt). This is exact for
%   grid-valued durations, and the pmf is bounded by 1 so the likelihood
%   cannot diverge. In continuous mode the density is conditioned instead,
%   f(t)/S(xmin).
%
%   WHEN TO CALL THIS DISCRETE AND WHEN CONTINUOUS
%   Use "discrete" -- the DEFAULT -- whenever the durations are integer
%   multiples of the acquisition interval, which is the normal case for
%   bout data scored from a sampled recording, and always when any
%   observation equals xmin. Use "continuous" only for genuinely
%   real-valued durations whose quantization is far finer than any fitted
%   time scale. As dt -> 0 the discrete log-likelihood approaches the
%   continuous one plus n*log(dt) and the estimates converge, so discrete
%   is continuous-in-the-limit and never the riskier choice. See
%   FITEXPONENTIATEDWEIBULLMLE for the fuller discussion.
%
%   REGISTERED MODELS (modelName, parameters in order)
%     "gamma"           [shape k, scale]        k<1 gives a decreasing
%                                               hazard, which bout data
%                                               generally want
%     "gamma_fixedshape" [scale], with Shape supplied -- the building block
%                                               FITERLANGMLE sweeps over
%                                               integer shapes
%     "chisquared"      [nu]                    Gamma(nu/2, scale 2). The
%                                               scale is FIXED at 2, which
%                                               for durations in seconds is
%                                               an arbitrary constant, so
%                                               this is the gamma above with
%                                               one parameter deleted for no
%                                               physical reason. Provided
%                                               for completeness of a model
%                                               sweep; do not expect it to
%                                               compete.
%     "pearson3"        [shape k, scale, location]  Gamma shifted by a
%                                               location. The location is
%                                               parametrized as
%                                               loc = xmin - exp(z), so it
%                                               stays below xmin by
%                                               construction, and the gap
%                                               xmin-loc is guarded: as it
%                                               closes the likelihood
%                                               diverges, the same pathology
%                                               SHIFTLOGNORMAL_MLE documents
%                                               for its shift.
%     "weibull"         [scale lambda, shape k] Identical to
%                                               FITEXPONENTIATEDWEIBULLMLE
%                                               with FixAlpha=true; provided
%                                               as its own row for sweeps.
%     "powerlaw"        [alpha]                 Pareto with its scale PINNED
%                                               at the first bin's lower
%                                               edge. A free scale cancels
%                                               out of S(t)/S(xmin) and is
%                                               unidentified. NOTE this is
%                                               the BINNED CONTINUOUS power
%                                               law, not the zeta/Zipf
%                                               distribution Clauset et al.
%                                               use for discrete data; the
%                                               two are different models
%                                               with different
%                                               log-likelihoods. Binned
%                                               Pareto is used here because
%                                               every other family in the
%                                               registry is discretized the
%                                               same way, and mixing
%                                               conventions is exactly the
%                                               AIC error this engine is
%                                               built to avoid.
%     "beta"            [p, q], on [0, UpperBound]
%                                               Beta's support is bounded,
%                                               so UpperBound is REQUIRED
%                                               and must be fixed, not
%                                               fitted: estimating it makes
%                                               the likelihood unbounded as
%                                               it approaches max(data). Use
%                                               the recording length, which
%                                               a bout genuinely cannot
%                                               exceed.
%
%   INPUTS
%   eventseries : numeric vector of positive durations, same units as xmin
%   xmin        : positive scalar, the left-truncation threshold. A PROTOCOL
%                 constant (e.g. 300 s for the 5-minute sleep rule), not
%                 min(eventseries). Ties at xmin are fine here -- the pmf is
%                 bounded -- but see FITHYPEREXPONENTIALMLE for why they are
%                 fatal to a continuous mixture.
%   modelName   : one of the registry names above
%
%   NAME-VALUE OPTIONS
%   SamplingInterval    REQUIRED positive scalar, SAME UNITS as xmin. Sets
%                       the discrete bin width and n_min. Rescaling the
%                       durations means rescaling this too.
%   DistributionType    "discrete" (DEFAULT) or "continuous".
%   Shape               required for "gamma_fixedshape"; the fixed shape.
%   UpperBound          required for "beta"; the fixed upper support limit.
%   MinLocationGap      positive scalar, default 1e-3 (as a fraction of
%                       xmin). For "pearson3", warns when (xmin-loc)/xmin
%                       falls below it, i.e. the location has run up against
%                       the truncation point and its estimate is a stopping
%                       place on a divergence, not a maximum.
%   nStarts, MaxIter, MaxFunEvals, TolX, TolFun, RandomSeed, Verbose,
%   ErrorOnNoValidFit   as in FITHYPEREXPONENTIALMLE, but the optimizer
%                       caps are far lower here (24 starts, 4000 function
%                       evaluations) because these families have 1-3
%                       smooth parameters, where fminsearch converges in a
%                       few hundred evaluations. The mixture fitter needs
%                       its much larger budget because it has up to 2K-1
%                       parameters and a multimodal likelihood; copying
%                       those caps here only wastes time.
%
%   OUTPUT
%   H.Model, H.ParamNames, H.Params, H.ParamSE, H.CovValid
%   H.k, H.n, H.LogLik, H.PointwiseLogLik, H.AIC, H.AICc, H.BIC
%   H.SurvivalHandle   @(t) P(T > t | T >= xmin) for the fitted model. A
%                       bin's probability is Sh((n-1)*dt) - Sh(n*dt), so this
%                       one handle serves both goodness-of-fit expected counts
%                       and plotting. Empty on a soft-failed fit.
%   H.Success, H.Converged, H.ExitFlag, H.BestParamVector
%   H.Failed, H.FailureReason, H.DistributionType, H.xmin, H.Diagnostics
%
%   BIC = k*log(n) - 2*LogLik is reported alongside AICc because they
%   disagree by design: BIC's log(n) penalty is harsher than AIC's 2 for
%   n > 7, so BIC prefers simpler models. Report both and say which you
%   used; do not pick the one that favours your preferred model.
%
%   MODEL COMPARISON
%   AIC/AICc/BIC differences mean something only between fits of the SAME
%   observations under the SAME dominating measure. So compare only fits
%   made with the same DistributionType, the same SamplingInterval and the
%   same xmin. A discrete log-likelihood sits roughly n*log(dt) above a
%   continuous one, since a pmf is a density times a bin width; that gap is
%   an artefact of the measure, not evidence. PointwiseLogLik is provided
%   for the non-nested Vuong (1989, Econometrica 57:307-333) test:
%       D = Ha.PointwiseLogLik - Hb.PointwiseLogLik;
%       V = sqrt(numel(D))*mean(D)/std(D,1);      % ~ N(0,1)
%
%   NUMERICAL NOTES
%   The bin probability F(n*dt) - F((n-1)*dt) is never formed in whichever
%   tail it would cancel in. Below the median it is differenced from the
%   CDF; above it, from the SURVIVAL function, which MATLAB and Octave both
%   supply directly (gammainc(...,'upper'), betainc(...,'upper')). The
%   difference is not cosmetic: at gammainc(60,2), 1-lower underflows to
%   exactly 0 while upper returns 5.34e-25.
%
%   See also FITGAMMAMLE, FITERLANGMLE, FITCHISQUAREDMLE, FITPEARSON3MLE,
%   FITWEIBULLMLE, FITBETAMLE, FITHYPEREXPONENTIALMLE,
%   FITEXPONENTIATEDWEIBULLMLE.

arguments
    eventseries double {mustBeReal}
    xmin (1,1) double {mustBePositive}
    modelName (1,1) string
    options.SamplingInterval (1,1) double = NaN
    options.DistributionType (1,1) string {mustBeMember(options.DistributionType,["continuous","discrete"])} = "discrete"
    options.Shape (1,1) double = NaN
    options.UpperBound (1,1) double = NaN
    options.MinLocationGap (1,1) double {mustBePositive} = 1e-3
    options.MinTailFraction (1,1) double {mustBeNonnegative} = 1e-10
    options.nStarts (1,1) double {mustBeInteger,mustBePositive} = 24
    options.RandomSeed = []
    options.MaxIter (1,1) double {mustBeInteger,mustBePositive} = 2000
    options.MaxFunEvals (1,1) double {mustBeInteger,mustBePositive} = 4000
    options.TolX (1,1) double {mustBePositive} = 1e-10
    options.TolFun (1,1) double {mustBePositive} = 1e-10
    options.Verbose (1,1) logical = true
    options.ErrorOnNoValidFit (1,1) logical = true
end

PENALTY = 1e12;

if isnan(options.SamplingInterval)
    error('FitTruncatedDiscreteMLE:SamplingIntervalRequired', ...
        ['SamplingInterval is required, in the SAME UNITS as xmin. It is ' ...
         'the discrete bin width and sets n_min.']);
elseif ~(options.SamplingInterval > 0)
    error('FitTruncatedDiscreteMLE:InvalidSamplingInterval', ...
        'SamplingInterval must be a positive scalar.');
end

if ~isempty(options.RandomSeed)
    rng(options.RandomSeed);
end

isDiscrete = strcmp(options.DistributionType, "discrete");
dt = options.SamplingInterval;
n_min_pre = max(1, round(xmin / dt));
M  = getModel(modelName, xmin, options, (n_min_pre-1)*dt);

% ---------------------------------------------------------------- data prep
raw = eventseries(:);
nRaw = numel(raw);
finiteMask = isfinite(raw) & (raw > 0);
nNonFinite = nnz(~finiteMask);

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
    fprintf(['FitTruncatedDiscreteMLE(%s): dropped %d of %d input value(s) ' ...
        '(%d non-finite/non-positive, %d below the truncation point ' ...
        'xmin=%.4g).\n'], M.Name, nNonFinite + nBelowXmin, nRaw, ...
        nNonFinite, nBelowXmin, xmin);
end

diag_ = initDiagnostics();
diag_.nNonFinite = nNonFinite;
diag_.nBelowXmin = nBelowXmin;

if n < M.nPar + 2
    msg = sprintf(['Fewer than %d usable observations (n=%d) after applying ' ...
        'xmin.'], M.nPar + 2, n);
    if options.ErrorOnNoValidFit
        error('FitTruncatedDiscreteMLE:TooFewData', '%s', msg);
    end
    if options.Verbose
        fprintf('FitTruncatedDiscreteMLE(%s): %s Returning Failed=true.\n', M.Name, msg);
    end
    H = assembleOutput([], M, options, xmin, n, diag_, true, ['TooFewData: ' msg]);
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

% For bounded-support families, how much of the support the data actually
% occupy. A parameter governing the density near the upper bound is
% extrapolation wherever no observations reach it.
if isfinite(options.UpperBound) && options.UpperBound > 0
    diag_.SupportCoverage = max(data) / options.UpperBound;
end

g = diag_.GridSpacing;
griddedAtG = isfinite(g) && g > 0 && all(abs(data/g - round(data/g)) < 1e-6);
diag_.GridMismatch = griddedAtG && abs(dt - g) > 1e-6 * g;

if isDiscrete && diag_.nDistinctOnGrid < diag_.nDistinctData
    lost = 100 * (1 - diag_.nDistinctOnGrid / diag_.nDistinctData);
    if griddedAtG
        hint = sprintf([' The durations lie on a grid of spacing %.6g, so ' ...
            'SamplingInterval=%.6g is very likely what you meant.'], g, g);
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
    warning('FitTruncatedDiscreteMLE:GridMismatch', '%s', gridMsg);
end

% ------------------------------------------------------------------- fitting
objfun = @(z) negLogLik(z, M, data, xmin, n_min, isDiscrete, dt, PENALTY);
starts = M.starts(data, xmin, options.nStarts);

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
        error('FitTruncatedDiscreteMLE:NoValidFit', '%s', msg);
    end
    if options.Verbose
        fprintf('FitTruncatedDiscreteMLE(%s): %s Returning Failed=true.\n', M.Name, msg);
    end
    H = assembleOutput([], M, options, xmin, n, diag_, true, ['NoValidFit: ' msg]);
    return
end

[~, pointwise] = negLogLik(bestZ, M, data, xmin, n_min, isDiscrete, dt, PENALTY);

fit = struct();
fit.Params = M.unpack(bestZ);
fit.LogLik = -bestNegLL;
fit.PointwiseLogLik = pointwise;
fit.ExitFlag = bestExitFlag;
fit.Converged = (bestExitFlag == 1);
fit.BestParamVector = bestZ;

% --- asymptotic standard errors at the MLE, delta-methoded to natural scale
se = nan(1, M.nPar);
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
        Rinv = R \ eye(M.nPar);
        clear restoreWarnings
        covZ = Rinv * Rinv';
        dvec = diag(covZ);
        if all(isfinite(dvec)) && all(dvec > 0)
            covValid = true;
            J = numericalJacobian(@(zz) M.unpack(zz), bestZ, M.nPar);
            covTheta = J * covZ * J';
            se = sqrt(max(diag(covTheta), 0))';
        end
    end
catch
    % Leave SEs as NaN; point estimates above are unaffected.
end
fit.ParamSE = se;
fit.CovValid = covValid;

% Truncated survival function, P(T > t | T >= xmin), as a handle. One
% handle serves both the goodness-of-fit expected bin counts (a bin's
% probability is Sh((n-1)*dt) - Sh(n*dt)) and the plotted curve, so a
% comparison wrapper never has to re-derive each family's CDF.
if isDiscrete
    Sref = M.sf((n_min-1)*dt, fit.Params);
    fit.SurvivalHandle = @(t) min(M.sf(round(t/dt)*dt, fit.Params) / Sref, 1);
else
    Sref = M.sf(xmin, fit.Params);
    fit.SurvivalHandle = @(t) min(M.sf(t, fit.Params) / Sref, 1);
end

% --- model-specific guard (currently: pearson3's location running to xmin)
diag_ = M.guard(fit.Params, xmin, options, diag_);

% --- general guard: is the fit describing a TRUNCATED POPULATION at all?
% S(xmin) is the fraction of the untruncated model lying in the observed
% range, so n/S(xmin) is the population size the fit implies. A free
% location or bound lets the optimizer slide the whole distribution far
% below xmin and use the family purely as a tail SHAPE; shape, scale and
% location are then jointly unidentified along that ridge, and the
% likelihood can beat an honest 1- or 2-parameter fit.
%
% Observed on simulated exponential data: Pearson III converged to
% location = -5.34e5 (six days before zero) with S(xmin) = 8.8e-320,
% implying 1e322 underlying bouts from 400 observed. It won the comparison
% outright and passed the goodness-of-fit test.
%
% The threshold is deliberately far below anything a real truncation
% produces -- a sleep criterion retaining a few percent of bouts gives
% S(xmin) of 0.01 to 0.1 -- so this fires only on the degenerate ridge,
% never on legitimate heavy truncation. Set MinTailFraction=0 to disable.
xminEff = xmin;
if isDiscrete, xminEff = (n_min-1)*dt; end
diag_.TailFraction = M.sf(xminEff, fit.Params);
if options.MinTailFraction > 0 && ~(diag_.TailFraction >= options.MinTailFraction)
    diag_.GuardOK = false;
    % log10 of the implied population, because n/S overflows to Inf for
    % the S values this actually fires on (6.5e-320 gives 1e322).
    logPop = log10(n) - log10(max(diag_.TailFraction, realmin*1e-16));
    msg = sprintf(['the fit has slid almost entirely below the truncation ' ...
        'point: S(xmin) = %.3g, so only that fraction of the fitted ' ...
        'distribution lies in the observed range and the %d observations ' ...
        'imply an underlying population of about 1e%.0f. The family is ' ...
        'being used as a bare tail shape, along which its parameters are ' ...
        'jointly unidentified, so this is not a maximum in any useful ' ...
        'sense. Treat this fit as unusable and prefer a model with fewer ' ...
        'free parameters.'], diag_.TailFraction, n, logPop);
    if isempty(diag_.Recommendation)
        diag_.Recommendation = msg;
    else
        diag_.Recommendation = [diag_.Recommendation ' ALSO: ' msg];
    end
    warning('FitTruncatedDiscreteMLE:FitBelowTruncation', '%s', msg);
end

H = assembleOutput(fit, M, options, xmin, n, diag_, false, '');

if options.Verbose
    fprintf('FitTruncatedDiscreteMLE(%s): %s mode, n=%d, k=%d parameter(s).\n', ...
        M.Name, options.DistributionType, n, M.nPar);
    fprintf('  logL=%.4f  AICc=%.4f  BIC=%.4f\n', H.LogLik, H.AICc, H.BIC);
    for j = 1:M.nPar
        if H.CovValid
            fprintf('  %-10s = %.6g (SE %.4g)\n', M.ParamNames{j}, H.Params(j), H.ParamSE(j));
        else
            fprintf('  %-10s = %.6g\n', M.ParamNames{j}, H.Params(j));
        end
    end
    if ~H.CovValid
        fprintf('  [Hessian-based standard errors unavailable -- point estimates only]\n');
    end
end

end

% ========================================================================
function H = assembleOutput(fit, M, options, xmin, n, diag_, failed, reason)
%ASSEMBLEOUTPUT Single exit point, so field names and their order are
%identical whether the fit succeeded or soft-failed. MATLAB refuses
%"results(i) = H" between structs whose fields differ in name or order,
%which would break sweeps over models and over animals.
if isempty(fit)
    fit = struct('Params', nan(1, M.nPar), 'ParamSE', nan(1, M.nPar), ...
        'CovValid', false, 'LogLik', NaN, 'PointwiseLogLik', [], ...
        'SurvivalHandle', [], 'Converged', false, 'ExitFlag', NaN, ...
        'BestParamVector', []);
end
H = struct();
H.Model = M.Name;
H.ParamNames = M.ParamNames;
H.Params = fit.Params;
H.ParamSE = fit.ParamSE;
H.CovValid = fit.CovValid;
H.k = M.nPar;
H.n = n;
H.LogLik = fit.LogLik;
H.PointwiseLogLik = fit.PointwiseLogLik;
H.SurvivalHandle = fit.SurvivalHandle;
H.AIC = 2*M.nPar - 2*fit.LogLik;
if (n - M.nPar - 1) > 0
    H.AICc = H.AIC + (2*M.nPar*(M.nPar+1)) / (n - M.nPar - 1);
else
    H.AICc = Inf;
end
H.BIC = M.nPar*log(n) - 2*fit.LogLik;
H.Success = ~failed;
H.Converged = fit.Converged;
H.ExitFlag = fit.ExitFlag;
H.BestParamVector = fit.BestParamVector;
H.DistributionType = options.DistributionType;
H.SamplingInterval = options.SamplingInterval;
H.xmin = xmin;
H.Failed = failed;
H.FailureReason = reason;
H.Diagnostics = diag_;
end

function d = initDiagnostics()
d = struct('nNonFinite', 0, 'nBelowXmin', 0, 'nAtXmin', 0, ...
    'XminEqualsDataMin', false, 'LooksGridded', false, ...
    'GridSpacing', NaN, 'GridMismatch', false, 'nDistinctData', 0, ...
    'nDistinctOnGrid', 0, 'LocationGap', NaN, 'SupportCoverage', NaN, ...
    'TailFraction', NaN, 'GuardOK', true, ...
    'Recommendation', '');
end

function [negLL, pointwise] = negLogLik(z, M, data, xmin, n_min, isDiscrete, dt, PENALTY)
%NEGLOGLIK Left-truncated negative log-likelihood for any registry model.
th = M.unpack(z);
t = data(:);

if ~M.valid(th)
    negLL = PENALTY; pointwise = -inf(size(t)); return
end

if ~isDiscrete
    logf = M.logpdf(t, th);
    Sx = M.sf(xmin, th);
    if ~(Sx > 0), negLL = PENALTY; pointwise = -inf(size(t)); return; end
    pointwise = logf - log(Sx);
else
    hi = round(t / dt) * dt;
    lo = hi - dt;
    % Bin probability, differenced in whichever tail does not cancel.
    Fhi = M.cdf(hi, th);
    p = zeros(size(hi));
    lower = Fhi <= 0.5;
    if any(lower)
        p(lower) = Fhi(lower) - M.cdf(lo(lower), th);
    end
    if any(~lower)
        p(~lower) = M.sf(lo(~lower), th) - M.sf(hi(~lower), th);
    end
    Sx = M.sf((n_min-1)*dt, th);
    if ~(Sx > 0) || any(~isfinite(p))
        negLL = PENALTY; pointwise = -inf(size(t)); return
    end
    % Floor the CONDITIONAL bin probability p/Sx, and cap it at 1. Flooring
    % p alone -- what an earlier version did -- MANUFACTURES likelihood.
    % Deep in a bad region p and Sx both underflow, but Sx can reach the
    % smallest SUBNORMAL (4.9e-324, which still passes Sx>0) while the
    % floor pins p at realmin (2.2e-308). The ratio is then 4.5e15: a
    % conditional probability fifteen orders of magnitude above 1, worth
    % +36 of log-likelihood per observation. A Weibull with scale 10.7
    % against data starting at 300 scored logL = +39972 and won a model
    % comparison outright.
    %
    % The floor itself has to stay. Returning PENALTY on an underflowing
    % bin instead removes the ramp fminsearch needs to walk out of that
    % region: with no gradient anywhere on the plateau, every multistart
    % from a poor scale dies there and the fit fails outright. Flooring the
    % ratio keeps a finite objective while putting the plateau at
    % log(realmin) = -708 per observation -- catastrophically BAD rather
    % than spuriously good, so the optimizer leaves it instead of climbing
    % into it. Capping at 1 makes the sign of the log-likelihood an
    % invariant, not an accident of which quantity underflowed first.
    r = min(p / Sx, 1);
    r(r <= 0) = realmin;
    pointwise = log(r);
end

negLL = -sum(pointwise);
if ~isfinite(negLL)
    negLL = PENALTY;
end
end

% ---------------------------------------------------------------- registry
function M = getModel(name, xmin, options, binFloor)
%GETMODEL The model registry. Each entry supplies a CDF, a survival
%function, a log density, an unconstrained->natural parameter map, start
%values, a validity test and an optional guard. Everything else -- the
%truncation, the discretization, the optimizer, the standard errors -- is
%the engine's.
name = char(name);
M = struct('Name', name, 'guard', @(th,xm,o,d) d);

switch name
    case 'gamma'
        M.ParamNames = {'shape','scale'};
        M.nPar = 2;
        M.unpack = @(z) exp(z(1:2));
        M.valid  = @(th) all(isfinite(th)) && all(th > 0);
        M.cdf    = @(t,th) gammainc(max(t,0)/th(2), th(1), 'lower');
        M.sf     = @(t,th) gammainc(max(t,0)/th(2), th(1), 'upper');
        M.logpdf = @(t,th) (th(1)-1)*log(t) - t/th(2) - th(1)*log(th(2)) - gammaln(th(1));
        M.starts = @(d,xm,ns) gammaStarts(d, xm, ns);

    case 'gamma_fixedshape'
        if isnan(options.Shape) || ~(options.Shape > 0)
            error('FitTruncatedDiscreteMLE:ShapeRequired', ...
                'Model "gamma_fixedshape" requires a positive Shape.');
        end
        k0 = options.Shape;
        M.ParamNames = {'scale'};
        M.nPar = 1;
        M.unpack = @(z) exp(z(1));
        M.valid  = @(th) isfinite(th) && th > 0;
        M.cdf    = @(t,th) gammainc(max(t,0)/th(1), k0, 'lower');
        M.sf     = @(t,th) gammainc(max(t,0)/th(1), k0, 'upper');
        M.logpdf = @(t,th) (k0-1)*log(t) - t/th(1) - k0*log(th(1)) - gammaln(k0);
        M.starts = @(d,xm,ns) log(max((mean(d))/max(k0,eps), eps)) + ...
                              [0; randn(max(ns-1,1),1)*0.8];

    case 'chisquared'
        % chi2(nu) == Gamma(nu/2, scale 2). The scale is fixed at 2, which
        % for second-valued durations is an arbitrary constant; see the
        % header. One free parameter.
        M.ParamNames = {'nu'};
        M.nPar = 1;
        M.unpack = @(z) exp(z(1));
        M.valid  = @(th) isfinite(th) && th > 0;
        M.cdf    = @(t,th) gammainc(max(t,0)/2, th(1)/2, 'lower');
        M.sf     = @(t,th) gammainc(max(t,0)/2, th(1)/2, 'upper');
        M.logpdf = @(t,th) (th(1)/2-1)*log(t) - t/2 - (th(1)/2)*log(2) - gammaln(th(1)/2);
        M.starts = @(d,xm,ns) log(max(mean(d), eps)) + [0; randn(max(ns-1,1),1)*0.5];

    case 'pearson3'
        % Gamma shifted by a location. loc = xmin - exp(z3) keeps loc < xmin
        % by construction; the gap is guarded because the likelihood
        % diverges as it closes.
        M.ParamNames = {'shape','scale','location'};
        M.nPar = 3;
        M.unpack = @(z) [exp(z(1)), exp(z(2)), xmin - exp(z(3))];
        M.valid  = @(th) all(isfinite(th)) && th(1) > 0 && th(2) > 0 && th(3) < xmin;
        M.cdf    = @(t,th) gammainc(max(t-th(3),0)/th(2), th(1), 'lower');
        M.sf     = @(t,th) gammainc(max(t-th(3),0)/th(2), th(1), 'upper');
        M.logpdf = @(t,th) (th(1)-1)*log(max(t-th(3),realmin)) - (t-th(3))/th(2) ...
                           - th(1)*log(th(2)) - gammaln(th(1));
        M.starts = @(d,xm,ns) pearson3Starts(d, xm, ns);
        M.guard  = @pearson3Guard;

    case 'weibull'
        M.ParamNames = {'scale','shape'};
        M.nPar = 2;
        M.unpack = @(z) exp(z(1:2));
        M.valid  = @(th) all(isfinite(th)) && all(th > 0);
        M.cdf    = @(t,th) -expm1(-(max(t,0)/th(1)).^th(2));
        M.sf     = @(t,th) exp(-(max(t,0)/th(1)).^th(2));
        M.logpdf = @(t,th) log(th(2)) - log(th(1)) + (th(2)-1)*(log(t)-log(th(1))) ...
                           - (t/th(1)).^th(2);
        M.starts = @(d,xm,ns) [log(mean(d)) 0; ...
                               log(mean(d))+randn(max(ns-1,1),1)*0.8, randn(max(ns-1,1),1)*0.6];

    case 'beta'
        if isnan(options.UpperBound) || ~(options.UpperBound > xmin)
            error('FitTruncatedDiscreteMLE:UpperBoundRequired', ...
                ['Model "beta" requires UpperBound > xmin, fixed rather ' ...
                 'than fitted (the recording length is the natural ' ...
                 'choice). Estimating it makes the likelihood unbounded ' ...
                 'as it approaches max(eventseries).']);
        end
        B = options.UpperBound;
        M.ParamNames = {'p','q'};
        M.nPar = 2;
        M.unpack = @(z) exp(z(1:2));
        M.valid  = @(th) all(isfinite(th)) && all(th > 0);
        M.cdf    = @(t,th) betainc(min(max(t,0)/B,1), th(1), th(2), 'lower');
        M.sf     = @(t,th) betainc(min(max(t,0)/B,1), th(1), th(2), 'upper');
        M.logpdf = @(t,th) (th(1)-1)*log(t/B) + (th(2)-1)*log(max(1-t/B,realmin)) ...
                           - betaln(th(1), th(2)) - log(B);
        M.starts = @(d,xm,ns) [0 0; randn(max(ns-1,1),2)*0.8];

    case 'powerlaw'
        % Pareto. Its scale is PINNED, not fitted: under left truncation
        % S(t)/S(xmin) = (xmin/t)^alpha, in which a free scale cancels
        % exactly, so it is unidentified -- the same cancellation that
        % un-identifies the exponentiated Weibull's alpha deep in the tail.
        % Pinning it at the lower edge of the first retained bin also keeps
        % that bin's probability positive; pinning at xmin itself would give
        % the first bin zero mass and assign probability 0 to every
        % observation sitting at xmin.
        if binFloor > 0
            s0 = binFloor;
        else
            s0 = xmin;
        end
        M.ParamNames = {'alpha'};
        M.nPar = 1;
        M.unpack = @(z) exp(z(1));
        M.valid  = @(th) isfinite(th) && th > 0;
        M.cdf    = @(t,th) 1 - min((s0./max(t,s0)).^th(1), 1);
        M.sf     = @(t,th) min((s0./max(t,s0)).^th(1), 1);
        M.logpdf = @(t,th) log(th(1)) + th(1)*log(s0) - (th(1)+1)*log(max(t,s0));
        M.starts = @(d,xm,ns) log(max(1/max(mean(log(d/s0)),eps), 1e-3)) + ...
                              [0; randn(max(ns-1,1),1)*0.5];

    otherwise
        error('FitTruncatedDiscreteMLE:UnknownModel', ...
            ['Unknown model "%s". Registered: gamma, gamma_fixedshape, ' ...
             'chisquared, pearson3, weibull, beta, powerlaw.'], name);
end
end

function s = gammaStarts(d, xmin, ns)
% Method of moments on the excess over xmin, plus scatter. The excess is
% used because the fit conditions on T >= xmin, so its scale is what the
% likelihood actually sees.
m = max(mean(d) - xmin, eps);
v = max(var(d), eps);
k0 = max(m^2/v, 1e-3);
s = [log(k0), log(max(v/m, eps)); ...
     log(0.5)+zeros(1,1), log(max(m,eps)); ...
     bsxfun(@plus, [log(k0) log(max(v/m,eps))], randn(max(ns-2,1),2)*0.9)];
end

function s = pearson3Starts(d, xmin, ns)
m = max(mean(d) - xmin, eps);
v = max(var(d), eps);
k0 = max(m^2/v, 1e-3);
% z3 parametrizes the gap below xmin: loc = xmin - exp(z3).
s = [log(k0), log(max(v/m,eps)), log(max(xmin,eps)); ...
     log(0.5), log(max(m,eps)),  log(max(xmin/2,eps)); ...
     bsxfun(@plus, [log(k0) log(max(v/m,eps)) log(max(xmin,eps))], ...
            randn(max(ns-2,1),3).*[0.9 0.9 1.2])];
end

function diag_ = pearson3Guard(th, xmin, options, diag_)
gap = (xmin - th(3)) / xmin;
diag_.LocationGap = gap;
diag_.GuardOK = gap >= options.MinLocationGap;
if ~diag_.GuardOK
    msg = sprintf(['the fitted location has run up against the truncation ' ...
        'point: (xmin-location)/xmin = %.3g < %.3g. The Pearson III ' ...
        'likelihood diverges as the location approaches the smallest ' ...
        'retained duration, so location=%.6g is where the optimizer ' ...
        'stopped on that divergence, not a maximum. Treat this fit as ' ...
        'unusable and prefer the 2-parameter gamma.'], gap, ...
        options.MinLocationGap, th(3));
    diag_.Recommendation = msg;
    warning('FitTruncatedDiscreteMLE:LocationAtBoundary', '%s', msg);
end
end

% ---------------------------------------------------------------- numerics
function J = numericalJacobian(fun, z0, m)
%NUMERICALJACOBIAN Central-difference d(natural params)/d(z), m x p.
z0 = z0(:).';
p = numel(z0);
h = max(1e-6, 1e-6*abs(z0));
J = zeros(m, p);
for i = 1:p
    zp = z0; zp(i) = zp(i) + h(i);
    zm = z0; zm(i) = zm(i) - h(i);
    % Assigned to temporaries first: fun(zp)(:) is legal Octave but a
    % syntax error in MATLAB, which forbids indexing a call result.
    a = fun(zp); b = fun(zm);
    J(:,i) = (a(:) - b(:)) / (2*h(i));
end
end

function H = computeNumericalHessian(fun, z0)
%COMPUTENUMERICALHESSIAN Central-difference Hessian at Z0, step ~ eps^(1/4).
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
        zpp = z0; zpp(i)=zpp(i)+h(i); zpp(j)=zpp(j)+h(j);
        zpm = z0; zpm(i)=zpm(i)+h(i); zpm(j)=zpm(j)-h(j);
        zmp = z0; zmp(i)=zmp(i)-h(i); zmp(j)=zmp(j)+h(j);
        zmm = z0; zmm(i)=zmm(i)-h(i); zmm(j)=zmm(j)-h(j);
        val = (fun(zpp) - fun(zpm) - fun(zmp) + fun(zmm)) / (4*h(i)*h(j));
        H(i,j) = val; H(j,i) = val;
    end
end
end

function restoreWarningStates(varargin)
for ii = 1:numel(varargin)
    w = varargin{ii};
    warning(w.state, w.identifier);
end
end
