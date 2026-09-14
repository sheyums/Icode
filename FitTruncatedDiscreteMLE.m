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
%   StartZ              extra starting points for the optimizer, one per
%                       row, in its INTERNAL coordinates rather than the
%                       reported parameters. Advanced: it exists so a
%                       caller that has already fitted a neighbouring
%                       model can say where the optimum is. They are tried
%                       first and the model's own nStarts still run, so a
%                       supplied start can only add a candidate optimum,
%                       never suppress one. FITHYPERERLANGMLE uses it to
%                       walk its shape sweep. A row of the wrong width is
%                       an error; non-finite rows are dropped.
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
    options.Shapes double = []
    options.UpperBound (1,1) double = NaN
    options.MinLocationGap (1,1) double {mustBePositive} = 1e-3
    options.MinLocationSpanRatio (1,1) double {mustBeNonnegative} = 0.05
    options.MinMixtureCount (1,1) double {mustBeNonnegative} = 5
    options.nStarts (1,1) double {mustBeInteger,mustBePositive} = 24
    options.StartZ double = []
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
    % Pinned to 'twister' rather than bare rng(seed): rng keeps whatever
    % generator is current, and a parallel worker's default generator is
    % not the client's, so a bare call makes RandomSeed reproduce only
    % within one execution mode. In the client this is identical to the
    % default, so nothing changes serially.
    rng(options.RandomSeed, 'twister');
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
diag_.DataSpan = max(data) - xmin;   % used by the pearson3 location guard
diag_.n = n;                         % used by the mixture weight guard
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

% StartZ: caller-supplied starting points, tried FIRST and in addition to
% the model's own. The coordinates are the optimizer's internal ones, not
% the reported parameters -- this is a hook for a caller that has already
% fitted a neighbouring model and knows where the optimum is, not a user
% dial. FitHyperErlangMLE uses it to walk its shape sweep: the fit at m
% stages is a short step from the fit at m-1 once the branch mean is held
% fixed, so starting there converges in a fraction of the iterations a
% cold start needs. Nothing is removed by supplying it -- the cold starts
% still run -- so a warm start can only add a candidate optimum, never
% hide one.
if ~isempty(options.StartZ)
    sz = options.StartZ;
    if size(sz, 2) ~= size(starts, 2)
        error('FitTruncatedDiscreteMLE:StartZWidth', ...
            ['StartZ has %d columns but model "%s" optimizes over %d. ' ...
             'StartZ is in the internal coordinates, one row per ' ...
             'starting point.'], size(sz,2), M.Name, size(starts,2));
    end
    sz = sz(all(isfinite(sz), 2), :);
    starts = [sz; starts];
end

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
% THETA is the model's natural parameter vector -- what M.sf, M.cdf and
% M.guard are written in. PARAMS is what gets REPORTED, which for a
% mixture is not the same thing: the fitted weights are untruncated, and
% the weight worth reporting is the OBSERVED one. M.report maps between
% them and is the identity for every non-mixture family.
fit.Theta  = M.unpack(bestZ);
fit.Params = M.report(fit.Theta, xmin);
fit.LogLik = -bestNegLL;
fit.PointwiseLogLik = pointwise;
fit.ExitFlag = bestExitFlag;
fit.Converged = (bestExitFlag == 1);
fit.BestParamVector = bestZ;

% --- asymptotic standard errors at the MLE, delta-methoded to natural scale
se = nan(1, M.nReport);
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
            % Delta method onto the REPORTED parameters, so a mixture's
            % observed weights carry their own standard errors rather
            % than the untruncated weights' -- the transform is part of
            % the estimator, not a relabelling after it.
            J = numericalJacobian(@(zz) M.report(M.unpack(zz), xmin), ...
                                  bestZ, M.nReport);
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
    Sref = M.sf((n_min-1)*dt, fit.Theta);
    fit.SurvivalHandle = @(t) min(M.sf(round(t/dt)*dt, fit.Theta) / Sref, 1);
else
    Sref = M.sf(xmin, fit.Theta);
    fit.SurvivalHandle = @(t) min(M.sf(t, fit.Theta) / Sref, 1);
end

% --- model-specific guard (currently: pearson3's location running to xmin)
diag_ = M.guard(fit.Theta, xmin, options, diag_);

% --- S(xmin) as a DIAGNOSTIC, not a gate.
% The fraction of the untruncated model lying in the observed range. Worth
% reporting; must NOT be gated on, and an earlier version of this file made
% exactly that mistake.
%
% A tiny S(xmin) is normal and harmless for a family whose support is
% ANCHORED at zero. A gamma driven to shape -> 0 has S(xmin) -> 0 because
% Gamma(a) -> infinity, yet the TRUNCATED law it defines is perfectly well
% behaved: Gamma(a,z) -> E1(z), a finite positive number, so the truncated
% density tends to t^-1 exp(-t/theta) / E1(xmin/theta), a proper
% distribution identified in theta. That is the power-law-with-cutoff
% limit, a legitimate boundary member of the family. Gating on
% S(xmin) = 3.1e-12 excluded a sound gamma fit. Nothing in the analysis
% ever evaluates the model below xmin, so the untruncated interpretation
% being vacuous costs nothing.
%
% The pathology that DOES matter needs a free LOCATION to slide, and is
% caught where it belongs, in pearson3Guard.
xminEff = xmin;
if isDiscrete, xminEff = (n_min-1)*dt; end
diag_.TailFraction = M.sf(xminEff, fit.Theta);

H = assembleOutput(fit, M, options, xmin, n, diag_, false, '');

if options.Verbose
    fprintf('FitTruncatedDiscreteMLE(%s): %s mode, n=%d, k=%d parameter(s).\n', ...
        M.Name, options.DistributionType, n, M.nPar);
    fprintf('  logL=%.4f  AICc=%.4f  BIC=%.4f\n', H.LogLik, H.AICc, H.BIC);
    for j = 1:M.nReport
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
    fit = struct('Theta', nan(1, M.nPar), ...
        'Params', nan(1, M.nReport), 'ParamSE', nan(1, M.nReport), ...
        'CovValid', false, 'LogLik', NaN, 'PointwiseLogLik', [], ...
        'SurvivalHandle', [], 'Converged', false, 'ExitFlag', NaN, ...
        'BestParamVector', []);
end
H = struct();
H.Model = M.Name;
H.ParamNames = M.ParamNames;
H.Params = fit.Params;
% The natural vector, for a caller that needs to re-enter the model --
% FITHYPERERLANGMLE warm-starts its sweep from it. Reported Params may be
% a transform of this and must not be fed back in.
H.Theta = fit.Theta;
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
    'TailFraction', NaN, 'DataSpan', NaN, 'LocationSpanRatio', NaN, ...
    'n', 0, 'MixtureWeight', NaN, 'MixtureMinCount', NaN, ...
    'MixtureWeightUntruncated', NaN, ...
    'GuardOK', true, ...
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
% report: natural parameters -> the vector that gets REPORTED. The
% identity for every family whose parameters are already the ones worth
% printing; the mixtures override it so their weights are reported as
% OBSERVED shares, which is the convention across every mixture in this
% library (see FITHYPEREXPONENTIALMLE). nReport is that vector's length,
% which for a mixture EXCEEDS nPar: sum(q)=1 means one weight is not free,
% but all of them are worth printing.
M = struct('Name', name, 'guard', @(th,xm,o,d) d, ...
           'report', @(th,xm) th, 'nReport', NaN);

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

    case 'hyper_erlang'
        % HYPER-ERLANG: sum_j w_j * Erlang(m_j, lambda_j), integer m_j.
        %
        %   f_j(t) = lam^m t^(m-1) exp(-lam t) / (m-1)!
        %   S_j(t) = exp(-lam t) * sum_{i<m} (lam t)^i / i!     [erlangSFint]
        %   S(t)   = sum_j w_j S_j(t),  F = 1 - S               [heSF, heCDF]
        %
        % mean m/lam per branch; the branch hazard rises from 0 to lam for
        % m > 1 and is flat at lam for m = 1. The weights here are the
        % UNTRUNCATED w; heReport converts them to observed q for output.
        %
        % The exponential-native way to get a NON-MONOTONE hazard. A
        % hyperexponential arranges its phases in PARALLEL -- enter one of
        % K states, leave at that state's own rate -- and that topology
        % forces a completely monotone hazard at every K. An Erlang
        % arranges phases in SERIES, must pass through m stages, and its
        % hazard RISES. Allow both and the mixture hazard can fall, rise
        % and fall again while every parameter is still a rate.
        %
        % It NESTS the hyperexponential exactly: Shapes = ones(1,K) IS
        % FitHyperexponentialMLE at order K, to the digit. So "are
        % memoryless states enough?" becomes a constraint on this model
        % rather than a comparison between two of them.
        %
        % Phase-type of order sum(m_j), so it keeps a Markov reading: an
        % Erlang branch with m=2 says a bout in that branch passes through
        % two sequential sub-stages before it can end -- a refractory or
        % cumulative process, not a memoryless one. Hyper-Erlangs are dense
        % in the distributions on [0,inf) (Tijms 1994), so a hyper-Erlang
        % that still will not fit is telling you something structural.
        %
        % NOTE an order-2 phase-type cannot hump at all -- PH(2) hazards
        % are monotone -- so the smallest useful shape vector for a humped
        % hazard has sum(m_j) >= 3.
        if isempty(options.Shapes) || any(options.Shapes < 1) ...
                || any(options.Shapes ~= round(options.Shapes))
            error('FitTruncatedDiscreteMLE:ShapesRequired', ...
                ['Model "hyper_erlang" requires Shapes, a vector of ' ...
                 'positive INTEGER stage counts, one per component. ' ...
                 'Shapes=ones(1,K) reproduces the K-component ' ...
                 'hyperexponential; give one component m>=2 for a rising ' ...
                 'hazard contribution.']);
        end
        mm = sort(options.Shapes(:).');          % non-decreasing, see unpack
        J  = numel(mm);
        M.ParamNames = heParamNames(mm);
        M.nPar = 2*J - 1;                        % J rates + J weights, sum=1
        M.nReport = 2*J;                         % ...but all J are printed
        M.report = @(th,xm) heReport(th, mm, xm);
        M.unpack = @(z) heUnpack(z, mm);
        M.valid  = @(th) heValid(th, J);
        M.cdf    = @(t,th) heCDF(t, th, mm);
        M.sf     = @(t,th) heSF(t, th, mm);
        M.logpdf = @(t,th) heLogPdf(t, th, mm);
        M.starts = @(d,xm,ns) heStarts(d, mm, ns);
        M.guard  = @(th,xm,o,dg) heGuard(th, mm, xm, o, dg);

    case 'weibull_mix'
        % TWO-COMPONENT WEIBULL MIXTURE, for a NON-MONOTONE HAZARD.
        %
        %   S_i(t) = exp(-(t/a_i)^b_i),  a = scale, b = shape
        %   f_i(t) = (b/a)(t/a)^(b-1) S_i(t)
        %   S(t)   = w S_1(t) + (1-w) S_2(t),  F = 1 - S    [weibullMixSF]
        %
        % b = 1 is the exponential, b < 1 a decreasing component hazard,
        % b > 1 an increasing one. w here is the UNTRUNCATED weight;
        % wmReport converts to the observed q1/q2 for output.
        %
        % Every other family here has a monotone hazard, and a
        % hyperexponential has a strictly decreasing one at ANY order --
        % f(t) = sum q_j lam_j exp(-lam_j t) is a sum of decreasing terms,
        % so no number of exponential components can produce a hump. Data
        % whose hazard falls, rises, then falls again are therefore outside
        % all of it by construction rather than by evidence, and that is
        % what per0 DD WAKE bouts do: hazard 1.2e-2 at 2 s falling to
        % 1.1e-3 by 70 s, rising to 1.65e-3 near 500 s, then falling away
        % past 1000 s.
        %
        % This reaches all three regimes with two components: shape1 < 1
        % gives the steep early decline, shape2 > 1 contributes the rising
        % middle, and whichever component is heavier-tailed dominates again
        % at the top. An exponential first component cannot do it -- its
        % hazard is flat, so the mixture would START flat rather than
        % falling tenfold.
        %
        % It NESTS the useful special cases, which makes each a testable
        % hypothesis rather than a separate fit: shape1=shape2=1 is the
        % 2-component hyperexponential, shape1=1 is exponential + Weibull,
        % and identical components are a single Weibull. Because shape=1 is
        % an INTERIOR point of shape>0 and the model stays identified
        % there, testing shape2=1 is a REGULAR hypothesis -- ordinary
        % chi2(1) applies, unlike the mixture-order question, which needs
        % HYPEREXPONENTIALLRT.
        % q1/q2 are the OBSERVED weights, the same symbol every mixture
        % in this library uses; both are printed although only one is
        % free, so k stays 5.
        M.ParamNames = {'q1','scale1','shape1','q2','scale2','shape2'};
        M.nPar = 5;
        M.nReport = 6;
        M.report = @(th,xm) wmReport(th, xm);
        M.unpack = @(z) weibullMixUnpack(z);
        M.valid  = @(th) all(isfinite(th)) && th(1) > 0 && th(1) < 1 ...
                         && all(th(2:5) > 0);
        M.cdf    = @(t,th) 1 - weibullMixSF(t, th);
        M.sf     = @(t,th) weibullMixSF(t, th);
        M.logpdf = @(t,th) weibullMixLogPdf(t, th);
        M.starts = @(d,xm,ns) weibullMixStarts(d, ns);
        M.guard  = @weibullMixGuard;

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
             'chisquared, pearson3, weibull, hyper_erlang, weibull_mix, ' ...
             'beta, powerlaw.'], name);
end

% A family that did not override report/nReport reports its natural
% parameters unchanged, one name each.
if ~isfinite(M.nReport)
    M.nReport = M.nPar;
end
if numel(M.ParamNames) ~= M.nReport
    error('FitTruncatedDiscreteMLE:ParamNameCount', ...
        ['Model "%s" declares %d reported parameters but names %d. Every ' ...
         'reported number must carry its own name, or the table prints ' ...
         'one parameter''s value under another''s label.'], ...
        name, M.nReport, numel(M.ParamNames));
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

function nm = heParamNames(mm)
%HEPARAMNAMES One name per REPORTED number: J rates, then J weights.
% All J weights are named and printed even though only J-1 are free --
% sum(q)=1 makes the last one determined, not uninteresting. k stays at
% 2J-1 and the information criteria are unaffected.
%
% The weights are q, the OBSERVED weights, as in every mixture here.
J = numel(mm);
nm = cell(1, 2*J);
for j = 1:J, nm{j} = sprintf('rate%d_m%d', j, mm(j)); end
for j = 1:J, nm{J+j} = sprintf('q%d', j); end
end

function pr = heReport(th, mm, xmin)
%HEREPORT Natural [rates, w] -> reported [rates, q].
% The fitted weights are UNTRUNCATED: heSF forms sum_j w_j S_j(t) and the
% engine divides by S(xmin). What a reader wants, and what every other
% mixture in this library reports, is the OBSERVED weight -- component
% j's share of the bouts that were actually RETAINED:
%
%       q_j  =  w_j S_j(xmin) / sum_i w_i S_i(xmin)
%
% The two differ by a factor that reaches 240x in these data, and w is
% extrapolation into a region the protocol excluded whenever a component's
% timescale sits below xmin. Reporting w under the name q, as this did,
% put a number in the table that meant something else entirely.
J = numel(mm);
rates = th(1:J); w = th(J+1:2*J);
Sx = zeros(1, J);
for j = 1:J
    Sx(j) = erlangSFint(rates(j)*xmin, mm(j));
end
contrib = w(:).' .* Sx;
total = sum(contrib);
if total > 0 && all(isfinite(contrib))
    q = contrib / total;
else
    q = nan(1, J);
end
pr = [rates(:).', q];
end

function th = heUnpack(z, mm)
%HEUNPACK z -> [rate1..rateJ, q1..qJ]. Weights by softmax with the last
%logit pinned at 0, so sum(q)=1 holds by construction and only J-1 of them
%are free -- the same device the hyperexponential uses in q coordinates.
%
% Components are ordered by (shape, rate). Shapes arrive sorted, so this
% only permutes within equal-shape blocks -- which is exactly where the
% label switching lives, since two components sharing a shape are
% exchangeable and the optimizer would otherwise return whichever of the
% equivalent optima it reached.
J = numel(mm);
rates = exp(z(1:J));
if J == 1
    q = 1;
else
    logits = [z(J+1:2*J-1), 0];
    mx = max(logits);
    e = exp(logits - mx);
    q = e / sum(e);
end
[~, ord] = sortrows([mm(:), rates(:)], [1 2]);
th = [rates(ord), q(ord)];
end

function ok = heValid(th, J)
ok = all(isfinite(th)) && all(th(1:J) > 0) && all(th(J+1:end) >= 0) ...
     && abs(sum(th(J+1:end)) - 1) < 1e-8;
end

function S = heSF(t, th, mm)
J = numel(mm); t = max(t, 0);
S = zeros(size(t));
for j = 1:J
    S = S + th(J+j) * erlangSFint(th(j)*t, mm(j));
end
end

function F = heCDF(t, th, mm)
% Computed from the lower branch directly rather than as 1-S, so the
% engine's median split gets an accurate small-probability value at the
% bottom of the range instead of a cancellation.
J = numel(mm); t = max(t, 0);
F = zeros(size(t));
for j = 1:J
    F = F + th(J+j) * erlangCDFint(th(j)*t, mm(j));
end
end

function S = erlangSFint(x, m)
%ERLANGSFINT Erlang survival at INTEGER shape, with no special function.
%   S = P(Poisson(x) <= m-1) = exp(-x) * sum_{j=0}^{m-1} x^j / j!, where
%   x = lambda*t. This is the Erlang-Poisson identity, EXACT rather than
%   an approximation, and it costs m terms -- m <= MaxShape, so at most a
%   handful. The recursion p_j = p_{j-1}*x/j is the standard stable
%   Poisson-pmf recursion: every term is positive, so nothing cancels and
%   the relative error stays at rounding.
%
%   Why not gammainc: a hyper-Erlang's shapes are integers BY
%   CONSTRUCTION -- that is the premise of the shape sweep -- while
%   gammainc is a general incomplete gamma that cannot know it and pays
%   for the general case on every call. In Octave it is interpreted, and
%   measurement put one hyper-Erlang row at 379 s against 49 s for seven
%   other families combined, with the cost nearly independent of n: it
%   was per-CALL overhead, not per-observation work.
x = max(x, 0);
p = exp(-x);                    % the j = 0 term
S = p;
for j = 1:m-1
    p = p .* x / j;
    S = S + p;
end
S = min(max(S, 0), 1);
end

function F = erlangCDFint(x, m)
%ERLANGCDFINT Erlang CDF at INTEGER shape: F = P(Poisson(x) >= m).
%   1 - S almost always, and the direct upper-tail sum only where that
%   subtraction would actually lose the answer.
%
%   S is the m-term sum and costs nothing. Forming F = 1 - S cancels, and
%   the relative error of the result is about eps/F -- harmless at
%   F = 0.1, fatal at F = 1e-300. So the direct sum is reserved for
%   elements where F has fallen below TINY, and those are exactly the
%   elements where x is small, where the tail sum's terms decay by x/j
%   per step and converge in a handful of iterations. The expensive
%   branch is therefore also the rare and cheap one.
%
%   Doing it the other way round -- direct sum whenever x < m -- is
%   correct but SLOWER THAN GAMMAINC: measured 0.7x at m=4 with a third
%   of the range below m, against 5x faster when the sum is skipped. The
%   discrete likelihood calls the CDF for every bin edge, so that branch
%   sat in the hot path and cost more than the special function it
%   replaced.
TINY = 1e-6;
x = max(x, 0);
S = erlangSFint(x, m);
F = 1 - S;
small = F < TINY;
if any(small(:))
    xs = x(small);
    p = exp(-xs);
    for j = 1:m
        p = p .* xs / j;        % p is now the j = m term
    end
    s = p;
    j = m;
    % max(s(:), realmin) is the two-argument elementwise max, so each
    % term is judged against ITS OWN running sum. Measured against an
    % independent series expansion this branch holds ~1e-14 relative
    % error down to F ~ 1e-69 -- better than Octave's own gammainc,
    % which drifts to 1.5e-3 relative at x = 0.1, m = 8.
    while j <= m + 1000
        j = j + 1;
        p = p .* xs / j;
        s = s + p;
        if all(p(:) <= eps * max(s(:), realmin)), break; end
    end
    F(small) = s;
end
F = min(max(F, 0), 1);
end

function lp = heLogPdf(t, th, mm)
%HELOGPDF Log density by logsumexp. Components of a hyper-Erlang differ by
%orders of magnitude over most of the range -- a fast branch contributes
%nothing at long durations and vice versa -- so summing densities directly
%underflows one away.
J = numel(mm); t = max(t(:), realmin);
L = zeros(numel(t), J);
for j = 1:J
    qj = th(J+j);
    if qj <= 0
        L(:,j) = -Inf;
    else
        % gammaln(m) for integer m is log((m-1)!), a sum of at most
        % MaxShape-1 logs. Computed directly so the density, like the
        % survival, makes no special-function call at all.
        L(:,j) = log(qj) + mm(j)*log(th(j)) + (mm(j)-1)*log(t) ...
                 - th(j)*t - sum(log(1:mm(j)-1));
    end
end
mx = max(L, [], 2);
lp = mx + log(sum(exp(L - repmat(mx, 1, J)), 2));
lp(~isfinite(mx)) = -Inf;
lp = reshape(lp, size(t));
end

function st = heStarts(d, mm, ns)
%HESTARTS Rates spread across the data's own timescales.
% An Erlang branch with m stages has mean m/rate, so its rate start is
% scaled by m -- otherwise a 3-stage branch starts three times too slow and
% the multistart wastes its budget walking there.
ds = sort(d(:));
q = @(p) max(ds(max(1, min(numel(ds), round(p*numel(ds))))), realmin);
J = numel(mm);
anchors = zeros(1, J);
for j = 1:J
    p = (j - 0.5) / J;                    % spread over the quantiles
    anchors(j) = log(mm(j) / q(p));       % rate = stages / typical duration
end
base = [anchors, zeros(1, J-1)];
extra = max(ns-1, 1);
sc = [repmat(0.9, 1, J), repmat(0.7, 1, J-1)];
st = [base; repmat(base, extra, 1) + randn(extra, 2*J-1) .* repmat(sc, extra, 1)];
end

function diag_ = heGuard(th, mm, xmin, options, diag_)
%HEGUARD The same two failures as any mixture: a component nobody is in,
%and two components that have become one.
%
% "Nobody is in it" must be judged on the OBSERVED share, not on the
% mixing weight. These weights are UNTRUNCATED -- heSF forms
% sum_j w_j S_j(t) and the engine divides by S(xmin) -- so a component
% whose mass lies below the truncation point can carry a large w while
% contributing nothing to the data. Its share of what was actually
% RETAINED is
%
%       obs_j  proportional to  w_j * S_j(xmin)
%
% and that is the quantity its rate is estimated from. Testing min(w)
% instead misses the case entirely: a 3-component fit at xmin=100 was
% seen with w = [0.055, 0.809, 0.136] and rates [0.019, 0.566, 0.0074],
% where the middle component -- the largest weight of the three -- has
% mean 1.8 s, so S(100) = 3e-25 and it holds ZERO of the 1500
% observations. The guard passed it, and the phantom component bought a
% shape the data never supported. The same fit reached from a different
% starting point sends that weight to 1e-10 instead, where min(w) does
% catch it: the two are the same truncated law, one of them visibly
% degenerate and one of them not, which is precisely why the test has to
% be on the observed share.
%
% This is the q-versus-w distinction the hyperexponential fitter makes
% everywhere, applied to the guard rather than to the report.
J = numel(mm);
w = th(J+1:end);
Sx = zeros(1, J);
for j = 1:J
    Sx(j) = erlangSFint(th(j)*xmin, mm(j));
end
contrib = w(:).' .* Sx;
total = sum(contrib);
if ~(total > 0) || ~all(isfinite(contrib))
    diag_.GuardOK = false;
    diag_.MixtureWeight = 0;
    diag_.MixtureMinCount = 0;
    msg = ['no hyper-Erlang component retains any mass above xmin, so ' ...
           'the truncated likelihood is not identified at all.'];
    diag_.Recommendation = msg;
    warning('FitTruncatedDiscreteMLE:MixtureComponentEmpty', '%s', msg);
    return
end
share = contrib / total;
diag_.MixtureWeight = min(share);               % OBSERVED share
diag_.MixtureWeightUntruncated = min(w);        % the mixing weight
nEff = diag_.n * min(share);
diag_.MixtureMinCount = nEff;
if nEff < options.MinMixtureCount
    diag_.GuardOK = false;
    [~, jm] = min(share);
    msg = sprintf(['a hyper-Erlang component holds only %.2f of the %d ' ...
        'observations (observed share %.4g, mixing weight %.4g, rate ' ...
        '%.4g at %d stage(s)), so its rate is not estimated and the ' ...
        'parameter count the information criteria charge is wrong. Drop a ' ...
        'component, or use the shape vector that fits.'], nEff, diag_.n, ...
        min(share), w(jm), th(jm), mm(jm));
    diag_.Recommendation = msg;
    warning('FitTruncatedDiscreteMLE:MixtureComponentEmpty', '%s', msg);
    return
end
for a = 1:J-1
    if mm(a) == mm(a+1) && abs(log(th(a+1)/th(a))) < 0.05
        diag_.GuardOK = false;
        msg = sprintf(['two hyper-Erlang components share a shape (m=%d) ' ...
            'and have converged on the same rate (%.6g vs %.6g), so their ' ...
            'weights are unidentified -- mass slides between the twins ' ...
            'with no change in likelihood and the optimum is a ridge, not ' ...
            'a point. Use fewer components at that shape.'], ...
            mm(a), th(a), th(a+1));
        diag_.Recommendation = msg;
        warning('FitTruncatedDiscreteMLE:MixtureCollapsed', '%s', msg);
        return
    end
end
end

function th = weibullMixUnpack(z)
%WEIBULLMIXUNPACK Unconstrained z to [w1 scale1 shape1 scale2 shape2].
% Components are SORTED by scale. A mixture is invariant to relabelling
% its components, so without a rule the optimizer drifts between two
% equivalent optima, the multistart returns whichever it happened to land
% on, and the standard errors describe a parameter whose identity changes
% between runs.
w = 1 / (1 + exp(-z(1)));
sc = [exp(z(2)), exp(z(4))];
sh = [exp(z(3)), exp(z(5))];
if sc(1) > sc(2)
    sc = sc([2 1]); sh = sh([2 1]); w = 1 - w;
end
th = [w, sc(1), sh(1), sc(2), sh(2)];
end

function S = weibullMixSF(t, th)
t = max(t, 0);
S = th(1) * exp(-(t/th(2)).^th(3)) + (1-th(1)) * exp(-(t/th(4)).^th(5));
end

function lp = weibullMixLogPdf(t, th)
%WEIBULLMIXLOGPDF Log density by logsumexp, not log(sum(exp)).
% One component is routinely many orders of magnitude below the other over
% part of the range -- that is the point of a heterogeneous mixture -- so
% summing the densities directly lets the smaller one underflow the larger
% away and loses the very structure being fitted.
t = max(t, realmin);
l1 = log(th(1))     + log(th(3)) - log(th(2)) ...
     + (th(3)-1)*(log(t)-log(th(2))) - (t/th(2)).^th(3);
l2 = log(1-th(1))   + log(th(5)) - log(th(4)) ...
     + (th(5)-1)*(log(t)-log(th(4))) - (t/th(4)).^th(5);
m = max(l1, l2);
lp = m + log(exp(l1-m) + exp(l2-m));
lp(~isfinite(m)) = -Inf;      % both components dead here
end

function st = weibullMixStarts(d, ns)
%WEIBULLMIXSTARTS Starts must BRACKET shape 1 in both components.
% The whole reason for this family is a hazard that falls and then rises,
% which needs shape < 1 in one component and shape > 1 in the other. Starts
% clustered on one side of 1 leave the optimizer to cross a region where
% the likelihood barely moves, and it frequently does not.
ds = sort(d(:));
q = @(p) ds(max(1, min(numel(ds), round(p*numel(ds)))));
lo = log(max(q(0.25), realmin));
hi = log(max(q(0.90), realmin));
base = [0, lo, log(0.7), hi, log(1.6)];
extra = max(ns-1, 1);
jit = randn(extra, 5) .* repmat([0.8 0.9 0.45 0.9 0.45], extra, 1);
st = [base; repmat(base, extra, 1) + jit];
end

function pr = wmReport(th, xmin)
%WMREPORT Natural [w, scale1, shape1, scale2, shape2] -> reported
%[q1, scale1, shape1, q2, scale2, shape2], with q the OBSERVED weights.
% Same argument as HEREPORT: the fitted w is untruncated, and under left
% truncation a component whose mass lies below xmin can hold a large w
% while contributing nothing to the retained data.
w = [th(1), 1-th(1)];
sc = [th(2), th(4)]; sh = [th(3), th(5)];
Sx = exp(-(max(xmin,0)./sc).^sh);
contrib = w .* Sx;
total = sum(contrib);
if total > 0 && all(isfinite(contrib))
    q = contrib / total;
else
    q = [NaN NaN];
end
pr = [q(1), th(2), th(3), q(2), th(4), th(5)];
end

function diag_ = weibullMixGuard(th, xmin, options, diag_)
%WEIBULLMIXGUARD A mixture that has stopped being a mixture.
% Two ways for five parameters to describe a one-component model, and in
% both the extra parameters are unidentified and the information matrix is
% singular, so the standard errors are meaningless and the parameter count
% the information criteria charge is wrong:
%   the weight has gone to a boundary, leaving one component carrying
%   everything -- judged against sample size, since a component holding
%   3 of 4000 observations is not estimated either;
%   the two components have converged on the same scale AND shape, which
%   is a single Weibull along a flat ridge, mass sliding freely between
%   the twins.
% Judged on the OBSERVED share, not on the untruncated weight -- see
% HEGUARD for the case that forced this: a component sitting below xmin
% keeps a healthy w while holding none of the data, and min(w) passes it.
pr = wmReport(th, xmin);
q = [pr(1), pr(4)];
w = th(1);
diag_.MixtureWeight = min(q);                   % OBSERVED share
diag_.MixtureWeightUntruncated = min(w, 1-w);   % the mixing weight
nEff = diag_.n * min(q);
diag_.MixtureMinCount = nEff;
sameScale = abs(log(th(4)/th(2))) < 0.05;
sameShape = abs(log(th(5)/th(3))) < 0.05;
if nEff < options.MinMixtureCount
    diag_.GuardOK = false;
    msg = sprintf(['one mixture component holds only %.2f of the %d ' ...
        'observations (observed share %.4g, mixing weight %.4g), so it ' ...
        'is not estimated: this is a single Weibull carrying five ' ...
        'parameters, three of them unidentified. Prefer the 1-component ' ...
        'weibull, whose parameter count the information criteria will ' ...
        'charge correctly.'], nEff, diag_.n, min(q), min(w, 1-w));
    diag_.Recommendation = msg;
    warning('FitTruncatedDiscreteMLE:MixtureComponentEmpty', '%s', msg);
elseif sameScale && sameShape
    diag_.GuardOK = false;
    msg = sprintf(['the two mixture components have converged on the same ' ...
        'law (scale %.6g vs %.6g, shape %.4g vs %.4g), so the weight is ' ...
        'unidentified -- mass slides between the twins with no change in ' ...
        'likelihood, and the optimum is a flat ridge rather than a ' ...
        'point. Prefer the 1-component weibull.'], ...
        th(2), th(4), th(3), th(5));
    diag_.Recommendation = msg;
    warning('FitTruncatedDiscreteMLE:MixtureCollapsed', '%s', msg);
end
end

function diag_ = pearson3Guard(th, xmin, options, diag_)
%PEARSON3GUARD The location fails in BOTH directions, and both failures are
%the same one: the location stops being identified.
%
% UP against xmin: the likelihood diverges as the location approaches the
% smallest retained duration, so the optimizer slides into the divergence
% and stops there rather than at a maximum.
%
% DOWN towards minus infinity: with (xmin - location) far larger than the
% SPAN of the data, (t - location) is nearly constant across every
% observation, the density flattens to const * exp(-t/scale), and shape
% and location stop being separately identifiable. Observed at
% location = -5.34e5 on data spanning 100 to 7060: an offset of 534000
% against a span of 6960, 1.3% of it. That fit beat the true model and won
% a comparison outright.
%
% The second check is on the SPAN RATIO, not on S(xmin). A small S(xmin)
% is a symptom shared with families that cannot slide at all -- a gamma at
% shape -> 0 has it too and is perfectly sound -- so gating on it excluded
% good fits. Only a free location produces this ridge.
gap = (xmin - th(3)) / xmin;
diag_.LocationGap = gap;
diag_.GuardOK = gap >= options.MinLocationGap;

offset = xmin - th(3);
span = diag_.DataSpan;
if diag_.GuardOK && offset > 0 && span > 0
    ratio = span / offset;
    diag_.LocationSpanRatio = ratio;
    if ratio < options.MinLocationSpanRatio
        diag_.GuardOK = false;
        msg = sprintf(['the fitted location has run far BELOW the data: ' ...
            'xmin-location = %.6g while the data span only %.6g, a ratio ' ...
            'of %.4g < %.4g. Over that range (t-location) is nearly ' ...
            'constant, so the density is indistinguishable from an ' ...
            'exponential and shape and location are not separately ' ...
            'identified -- the optimizer has stopped on a ridge, not at a ' ...
            'maximum. Treat this fit as unusable and prefer the ' ...
            '2-parameter gamma.'], offset, span, ratio, ...
            options.MinLocationSpanRatio);
        diag_.Recommendation = msg;
        warning('FitTruncatedDiscreteMLE:LocationBelowData', '%s', msg);
        return
    end
end

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
