function ModelResults = GenNorm_MLE(X, models_to_test, mode)
% GENNORM_MLE Fit generalized Normal (generalized Gaussian) mixtures by MLE.
%
%   ModelResults = GenNorm_MLE(X)
%   ModelResults = GenNorm_MLE(X, MODELS_TO_TEST)
%   ModelResults = GenNorm_MLE(X, MODELS_TO_TEST, MODE)
%
%   One component has density
%
%       f(x) = alpha / (2*sigma*Gamma(1/alpha)) * exp(-(|x - mu|/sigma)^alpha)
%
%   alpha = 2 is a Normal (with sigma = sqrt(2)*SD), alpha = 1 a Laplace,
%   alpha -> Inf a uniform on [mu - sigma, mu + sigma]. A K-term model is
%   sum_j w_j * f_j(x) with weights summing to 1.
%
%   INPUTS
%       X              Real numeric vector. NaN and Inf are removed (with a
%                      warning).
%       MODELS_TO_TEST Vector of component counts K (positive integers),
%                      default 1. A K needs at least 4K + 1 observations
%                      (one more than its 4K - 1 parameters, so that AICc
%                      is finite); larger K are dropped with a warning, and
%                      an error is raised if none remain.
%       MODE           Character vector or string; keywords may be combined
%                      ('Fit Select'), case-insensitive:
%                        'Fit'    - fill Fit and FitX.
%                        'Plot'   - histogram with the fitted densities.
%                        'Select' - keep only the model chosen by a greedy
%                                   likelihood-ratio test (see LIMITATIONS).
%                        'Profile'- profile-likelihood CIs for mu, sigma and
%                                   alpha (slower: a K = 2 fit at n = 500 took
%                                   25 s instead of 4 s). Recommended when
%                                   any alpha < 1.5.
%
%   OUTPUT
%       ModelResults   Structure array, one element per K (one element
%                      after 'Select'), with fields:
%         Model        Description string.
%         Params       K x 4 matrix [mu, sigma, alpha, weight], rows sorted
%                      by mu.
%         LogLike      Maximized log-likelihood.
%         AIC          2p - 2*LogLike, p = 4K - 1.  (Before Nov 2026 this
%                      field held AICc; AICc now has its own field.)
%         AICc         AIC + 2p(p+1)/(n - p - 1).
%         BIC          p*log(n) - 2*LogLike.
%         Fit          Fitted density at FitX ('Fit' mode only, else []).
%         FitX         Sorted X ('Fit' mode only, else []).
%         N            Number of observations used.
%         NParameters  p = 4K - 1.
%         Converged    true when the start that produced the returned fit
%                      ended with a positive fmincon exit flag.
%         ExitFlag     That exit flag.
%         AlphaAtBound K x 1 logical, rows as in Params: alpha within 1e-4
%                      (relative) of its bound 0.5 or 10. The true shape is
%                      then outside the range searched (e.g. uniform data
%                      want alpha -> Inf, Cauchy-like data alpha < 0.5).
%         SigmaAtFloor K x 1 logical: sigma on its floor 1e-3*IQR(X), i.e.
%                      a component collapsed onto a few points.
%         SE           K x 4 standard errors, laid out like Params. From the
%                      observed information (central differences of the
%                      analytic score), mapped to [mu, sigma, alpha, weight]
%                      by the delta method -- EXCEPT mu of a component with
%                      alpha < 1.5, which uses the outer product of the
%                      scores (OPG): there the likelihood has a cusp at each
%                      data point and the Hessian is meaningless. NaN for a
%                      parameter on its bound or floor (held fixed), for the
%                      weight when K = 1, and everywhere if the information
%                      matrix is not positive definite.
%         SEHessian,   The two SE estimates separately (diagnostic).
%         SEOPG
%         Cov          4K x 4K covariance (observed information) in the
%                      order [mu_1..K, sigma_1..K, alpha_1..K, w_1..K], rows
%                      as in Params.
%         CILower,     K x 4 confidence limits, laid out like Params. Wald by
%         CIUpper      default; with 'Profile', profile likelihood for mu,
%                      sigma and alpha (weights stay Wald). See CIMethod.
%         CIOpen       K x 4 logical: a profile end never crossed the
%                      threshold within the parameter's range, so the
%                      interval is open there (mu -> -/+Inf, sigma -> its
%                      floor or Inf, alpha -> 0.5 or 10).
%         CIMethod     'Wald' or 'profile (mu, sigma, alpha); Wald (weight)'.
%         CIWaldLower, Wald limits, always: mu +/- z*SE; sigma and alpha on
%         CIWaldUpper  the log scale; weights on the logit scale.
%         CILevel      0.95.
%
%   WARNINGS
%       GenNorm_MLE:NonfiniteDataRemoved   NaN/Inf dropped from X.
%       GenNorm_MLE:TooFewObservations     some K dropped (too few points).
%       GenNorm_MLE:AlphaAtBound           a shape parameter on its bound.
%       GenNorm_MLE:SigmaAtFloor           a scale parameter on its floor.
%       GenNorm_MLE:NotConverged           no start converged for some K.
%       GenNorm_MLE:CuspSE                 alpha < 1.5 without 'Profile': the
%                                          Wald interval for mu under-covers.
%       GenNorm_MLE:ProfileImprovedOptimum a profile fit beat the estimate.
%
%   ESTIMATION
%       The data are standardized, z = (x - median(x))/std(x), and fitted in
%       z-units; mu, sigma and LogLike are transformed back exactly, so the
%       fit is equivariant under x -> a*x + b (LogLike shifts by -n*log|a|).
%       Each K is fitted by fmincon (interior-point) from several starts:
%       k-means centres crossed with alpha in {1, 2, 5}, quantile-spread
%       centres, one broad start and, for K > 1, random starts; the best is
%       polished by a second fmincon run with tighter tolerances.
%       Coordinates: mu free; log(sigma) >= log(1e-3*IQR(z)); log(alpha) in
%       [log(0.5), log(10)]; K - 1 weight logits in [-20, 20].
%       Starts use a private random stream (seed 0); the caller's global
%       stream is neither used nor advanced, so results are reproducible.
%
%   LIMITATIONS (measured Nov 2026; see the history below)
%       - 'Select' compares nested K by a likelihood-ratio test with a
%         chi-square reference on 4*dK degrees of freedom. The number of
%         mixture components is a NON-REGULAR test (a component weight can
%         go to 0 at the boundary), so that reference is not valid and the
%         p-value is only a heuristic. Prefer BIC across K.
%       - Mixtures (K > 1) have many local optima; the multistart does not
%         guarantee the global one. Measured against an independent
%         40-start Nelder-Mead search on three-population data: K = 2 and
%         K = 3, 4 samples each, the fit was never more than 2e-6 below
%         the independent best (and above it in 3 of 8); before the random
%         starts the old fitter was 11.8, 5.97 and 1.09 below in 3 of 4
%         K = 2 samples.
%       - With alpha <= 1 the log-likelihood has a cusp at every data point
%         in mu, where a gradient-based optimizer is not exact (handled by
%         a data-point polish, also inside every profile fit).
%       - Confidence intervals: measured coverage of the nominal 95%
%         (K = 1: n = 300, 300 samples for Wald, 100 for profile; K = 2:
%         n = 500, 100 and 50; Monte Carlo SE about 0.01-0.03). Columns
%         mu / sigma / alpha:
%           alpha0 = 0.7  default Wald 0.80 / 0.96 / 0.96   profile 0.94 / 0.97 / 0.96
%           alpha0 = 1.0  default Wald 0.91 / 0.95 / 0.96   profile 0.93 / 0.95 / 0.94
%           alpha0 = 1.3  default Wald 0.94 / 0.94 / 0.96   profile 0.95 / 0.95 / 0.93
%           alpha0 = 2    default Wald 0.94 / 0.94 / 0.97   profile 0.93 / 0.94 / 0.97
%           alpha0 = 4    default Wald 0.94 / 0.96 / 0.96   profile 0.96 / 0.95 / 0.94
%           K = 2 (alpha 2 and 1.5, well separated): Wald 0.93-0.97 for
%           every parameter including the weights; profile 0.92-0.98.
%         The Hessian Wald interval for mu alone would cover 0.04, 0.33
%         and 0.88 at alpha0 = 0.7, 1 and 1.3 -- why mu switches to OPG
%         below alpha = 1.5. Below alpha ~ 1 use 'Profile' for mu. Before
%         the profile fits were polished at the cusps, profile sigma and
%         alpha covered only 0.63 at alpha0 = 0.7.
%       - Intervals for overlapping mixture components were not measured;
%         label switching (components trading places along a profile) can
%         make them unreliable there.
%
%   REQUIREMENTS
%       Optimization Toolbox (fmincon); Statistics and Machine Learning
%       Toolbox (kmeans, iqr, quantile, chi2cdf).
%
%   Originally written by Sheyum. Audited against an independent profile
%   likelihood (Claude Code, Nov 2026), see the history at the end.

%% 1. Input validation
if nargin < 2 || isempty(models_to_test)
    models_to_test = 1;
    fprintf('Default: Fitting 1 term.\n');
end
if nargin < 3 || isempty(mode)
    mode = '';
end
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('GenNorm_MLE:InvalidMode', ...
        'MODE must be a character vector or string scalar.');
end
mode = upper(char(mode));

if ~isnumeric(X) || ~isreal(X) || isempty(X)
    error('GenNorm_MLE:InvalidData', ...
        'X must be a non-empty real numeric vector.');
end
X = double(X(:));
nRemoved = sum(~isfinite(X));
X = X(isfinite(X));
if nRemoved > 0
    warning('GenNorm_MLE:NonfiniteDataRemoved', ...
        '%d nonfinite observation(s) were removed.', nRemoved);
end
if isempty(X)
    error('GenNorm_MLE:InvalidData', 'X has no finite observations.');
end
if max(X) == min(X)
    error('GenNorm_MLE:ConstantData', ...
        'X is constant; a generalized Normal cannot be fitted.');
end

if ~isnumeric(models_to_test) || any(models_to_test(:) < 1) || ...
        any(models_to_test(:) ~= round(models_to_test(:)))
    error('GenNorm_MLE:InvalidModels', ...
        'MODELS_TO_TEST must contain positive integers.');
end
models_to_test = unique(models_to_test(:))';

nObs    = numel(X);
maxK    = floor((nObs - 1) / 4);       % needs n >= 4K + 1
dropped = models_to_test(models_to_test > maxK);
models_to_test = models_to_test(models_to_test <= maxK);
if isempty(models_to_test)
    error('GenNorm_MLE:TooFewObservations', ...
        ['%d observations are too few: a K-term model needs at least ' ...
         '4K + 1 (5 for K = 1).'], nObs);
end
if ~isempty(dropped)
    warning('GenNorm_MLE:TooFewObservations', ...
        ['K = %s dropped: %d observations support at most K = %d ' ...
         '(a K-term model needs 4K + 1).'], ...
        mat2str(dropped), nObs, maxK);
end

%% 2. Private random stream and warning state, both restored on exit
globalStream = RandStream.getGlobalStream;
savedState   = globalStream.State;
savedWarn    = warning;
cleanup      = onCleanup(@() restore_state(globalStream, savedState, savedWarn));
RandStream.setGlobalStream(RandStream('twister', 'Seed', 0));
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'stats:kmeans:FailedToConverge');

%% 3. Standardize
loc = median(X);
scl = std(X);
Z   = (X - loc) / scl;

ModelResults = repmat(struct('Model', [], 'Params', [], 'LogLike', [], ...
    'AIC', [], 'AICc', [], 'BIC', [], 'Fit', [], 'FitX', [], 'N', [], ...
    'NParameters', [], 'Converged', [], 'ExitFlag', [], ...
    'AlphaAtBound', [], 'SigmaAtFloor', [], 'SE', [], 'Cov', [], ...
    'CILower', [], 'CIUpper', [], 'CIOpen', [], 'CIMethod', [], ...
    'CIWaldLower', [], 'CIWaldUpper', [], 'CILevel', [], ...
    'SEHessian', [], 'SEOPG', []), numel(models_to_test), 1);
ciLevel   = 0.95;
doProfile = contains(mode, 'PROFILE');

%% 4. Fit each K
for m = 1:numel(models_to_test)
    kComp = models_to_test(m);
    fprintf('Fitting %d-term Generalized Normal model...\n', kComp);

    stats = FitGenNorm_Engine(Z, kComp);

    P = stats.Params;                         % z-units
    P(:, 1) = loc + scl * P(:, 1);            % mu
    P(:, 2) = scl * P(:, 2);                  % sigma
    LL = stats.LogLike - nObs * log(scl);

    nParams = 4*kComp - 1;
    ModelResults(m).Model  = sprintf( ...
        '%d-term Generalized Normal (mu, sigma, alpha, weight)', kComp);
    ModelResults(m).Params      = P;
    ModelResults(m).LogLike     = LL;
    ModelResults(m).AIC         = 2*nParams - 2*LL;
    ModelResults(m).AICc        = 2*nParams - 2*LL + ...
        2*nParams*(nParams + 1) / (nObs - nParams - 1);
    ModelResults(m).BIC         = nParams*log(nObs) - 2*LL;
    ModelResults(m).N           = nObs;
    ModelResults(m).NParameters = nParams;
    ModelResults(m).Converged   = stats.Converged;
    ModelResults(m).ExitFlag    = stats.ExitFlag;
    ModelResults(m).AlphaAtBound = stats.AlphaAtBound;
    ModelResults(m).SigmaAtFloor = stats.SigmaAtFloor;

    U = gn_uncertainty(Z, stats, kComp, loc, scl, ciLevel, doProfile);
    ModelResults(m).SE          = U.SE;
    ModelResults(m).Cov         = U.Cov;
    ModelResults(m).CILower     = U.CILower;
    ModelResults(m).CIUpper     = U.CIUpper;
    ModelResults(m).CIOpen      = U.CIOpen;
    ModelResults(m).CIMethod    = U.CIMethod;
    ModelResults(m).CIWaldLower = U.CIWaldLower;
    ModelResults(m).CIWaldUpper = U.CIWaldUpper;
    ModelResults(m).CILevel     = ciLevel;
    ModelResults(m).SEHessian   = U.SEHessian;
    ModelResults(m).SEOPG       = U.SEOPG;
    if any(U.CuspMu) && ~doProfile
        warning('GenNorm_MLE:CuspSE', ...
            ['K = %d: component(s) %s have alpha < 1.5, where the Wald ' ...
             'interval for mu under-covers (measured 0.80 at alpha = 0.7, ' ...
             '0.91 at alpha = 1 for a nominal 0.95). Use the ''Profile'' ' ...
             'mode for mu there.'], kComp, mat2str(find(U.CuspMu)'));
    end
    if U.ProfileBest > stats.LogLike + 1e-6
        warning('GenNorm_MLE:ProfileImprovedOptimum', ...
            ['K = %d: a profile fit found a log-likelihood %.3g above ' ...
             'the estimate; the estimate is not the optimum and the ' ...
             'CIs are approximate.'], kComp, U.ProfileBest - stats.LogLike);
    end

    if ~stats.Converged
        warning('GenNorm_MLE:NotConverged', ...
            ['K = %d: no start converged (best exit flag %d). The ' ...
             'returned fit is the best finite one found.'], ...
            kComp, stats.ExitFlag);
    end
    if any(stats.AlphaAtBound)
        warning('GenNorm_MLE:AlphaAtBound', ...
            ['K = %d: alpha of component(s) %s is on its bound (0.5 or ' ...
             '10); the best shape lies outside the range searched.'], ...
            kComp, mat2str(find(stats.AlphaAtBound)'));
    end
    if any(stats.SigmaAtFloor)
        warning('GenNorm_MLE:SigmaAtFloor', ...
            ['K = %d: sigma of component(s) %s is on its floor ' ...
             '(1e-3*IQR); a component has collapsed onto a few points.'], ...
            kComp, mat2str(find(stats.SigmaAtFloor)'));
    end
end

%% 5. Model selection (greedy likelihood-ratio test; see LIMITATIONS)
if contains(mode, 'SELECT') && numel(models_to_test) > 1
    best_idx   = 1;
    alpha_crit = 0.01;
    stopP      = NaN;
    for i = 1:numel(models_to_test) - 1
        df      = 4 * (models_to_test(i+1) - models_to_test(i));
        LR_stat = max(0, 2 * (ModelResults(i+1).LogLike - ...
                              ModelResults(i).LogLike));
        pval    = chi2cdf(LR_stat, df, 'upper');
        if pval < alpha_crit
            best_idx = i + 1;
        else
            stopP = pval;
            break
        end
    end
    ModelResults = ModelResults(best_idx);
    if isnan(stopP)
        fprintf(['>> Auto-Selected: %s (largest K tested; every step ' ...
                 'had p < %g)\n'], ModelResults.Model, alpha_crit);
    else
        fprintf(['>> Auto-Selected: %s (the next larger K gave p = %.4g ' ...
                 '>= %g; chi-square reference is approximate)\n'], ...
            ModelResults.Model, stopP, alpha_crit);
    end
end

%% 6. Plot
if contains(mode, 'PLOT')
    plot_results(X, ModelResults);
end

%% 7. Fitted density at the sorted data
if contains(mode, 'FIT')
    sorted_x = sort(X);
    for i = 1:numel(ModelResults)
        ModelResults(i).FitX = sorted_x;
        ModelResults(i).Fit  = exp(GenNorm_LogPDF(sorted_x, ...
            ModelResults(i).Params));
    end
end
end

%% ============================================================
%  CORE OPTIMIZATION ENGINE (z-units)
% ============================================================
function stats = FitGenNorm_Engine(X, kComp)

range_iqr = iqr(X);
if range_iqr == 0
    range_iqr = std(X);
end
sig_min   = 1e-3 * range_iqr;      % relative floor; X is standardized
max_logit = 20;
alphaLo   = 0.5;
alphaHi   = 10;

LB = [-inf(1, kComp), log(sig_min)*ones(1, kComp), ...
      log(alphaLo)*ones(1, kComp), -max_logit*ones(1, kComp-1)];
UB = [ inf(1, kComp),  inf(1, kComp), ...
      log(alphaHi)*ones(1, kComp),  max_logit*ones(1, kComp-1)];

seeds = GenerateSeeds(X, kComp, sig_min, max_logit);

opt = optimoptions('fmincon', 'Display', 'off', ...
    'Algorithm', 'interior-point', 'MaxIterations', 2000, ...
    'MaxFunctionEvaluations', 2e4, 'OptimalityTolerance', 1e-6);
optPolish = optimoptions(opt, 'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12);

obj = @(p) GenNorm_NLL(X, p, kComp);

bestNLL = Inf;  bestP = [];  bestFlag = -Inf;
bestConvNLL = Inf;  bestConvP = [];  bestConvFlag = -Inf;
for i = 1:size(seeds, 1)
    p0 = min(max(seeds(i, :), LB), UB);
    try
        [popt, nll, flag] = fmincon(obj, p0, [], [], [], [], LB, UB, [], opt);
    catch
        continue                       % a start that errors is skipped
    end
    if ~isfinite(nll)
        continue
    end
    if nll < bestNLL
        bestNLL = nll;  bestP = popt;  bestFlag = flag;
    end
    if flag > 0 && nll < bestConvNLL
        bestConvNLL = nll;  bestConvP = popt;  bestConvFlag = flag;
    end
end

% Prefer the best converged start; fall back to the best finite one.
if ~isempty(bestConvP)
    bestNLL = bestConvNLL;  bestP = bestConvP;  bestFlag = bestConvFlag;
end

if isempty(bestP)
    stats.Params       = NaN(kComp, 4);
    stats.LogLike      = -Inf;
    stats.Converged    = false;
    stats.ExitFlag     = NaN;
    stats.AlphaAtBound = false(kComp, 1);
    stats.SigmaAtFloor = false(kComp, 1);
    stats.p = [];
    return
end

% Polish the winner; keep the polish only if it does not get worse.
try
    [pp, nllp, flagp] = fmincon(obj, bestP, [], [], [], [], LB, UB, [], ...
        optPolish);
    if isfinite(nllp) && nllp <= bestNLL
        bestP = pp;  bestNLL = nllp;
        if flagp > 0 || bestFlag <= 0
            bestFlag = flagp;
        end
    end
catch
end

% Cusp polish: with alpha <= 1 the optimum mu of a component sits ON a
% data point, where the gradient does not exist (see cusp_polish).
if any(exp(bestP(2*kComp+1:3*kComp)) < cusp_alpha())
    [bestP, bestNLL] = cusp_polish(X, bestP, bestNLL, kComp, LB, UB, ...
        optPolish, false(size(bestP)));
end

[stats.Params, order] = ReconstructParams(bestP, kComp);
stats.LogLike   = -bestNLL;
stats.Converged = bestFlag > 0;
stats.ExitFlag  = bestFlag;

la = bestP(2*kComp+1:3*kComp);
ls = bestP(kComp+1:2*kComp);
tolA = 1e-4;
atA  = abs(la - log(alphaLo)) <= tolA | abs(la - log(alphaHi)) <= tolA;
atS  = abs(ls - log(sig_min)) <= tolA;
stats.AlphaAtBound = atA(order)';
stats.SigmaAtFloor = atS(order)';

% For the uncertainty stage: the optimizer coordinates, the sort order,
% the box, and which coordinates sit on an active bound (held fixed there).
stats.p     = bestP;
stats.order = order;
stats.LB    = LB;
stats.UB    = UB;
stats.opt   = optPolish;
stats.fixed = abs(bestP - LB) <= tolA | abs(bestP - UB) <= tolA;
end

function [p, nll] = cusp_polish(X, p, nll, kComp, LB, UB, opt, held)
% For alpha <= 1, sum_i w_i*|x_i - mu|^alpha is concave between data points,
% so the likelihood's optimum in mu_j lies AT a data point -- a cusp, where
% fmincon's finite-difference gradient is meaningless and it stops short
% (measured: up to 0.08 log-likelihood units at K = 1, alpha = 0.6).
% Slightly above 1 the same stall happens: once mu sits on a data point the
% mu gradient is garbage and alpha stops early (measured: alpha 1.04495
% instead of 1.04465, 3.2e-6 short), so components with alpha below
% cusp_alpha() = 1.5 are polished. Only improvements are kept.
% For each such component: screen the 41 data points nearest mu_j with the
% other parameters fixed, re-optimize the other parameters with mu_j held at
% the 3 best, keep any improvement; repeat until nothing improves (max 3).
% HELD marks coordinates that must not move (a profile's fixed parameter):
% a held mu is not screened, and held coordinates stay out of every re-fit.
nCand = 41;  nKeep = 3;
xs = sort(X);
for iter = 1:3
    improved = false;
    for j = 1:kComp
        if exp(p(2*kComp + j)) >= cusp_alpha() || held(j)
            continue
        end
        [~, near] = sort(abs(xs - p(j)));
        cand = xs(near(1:min(nCand, numel(xs))));
        screen = zeros(numel(cand), 1);
        for c = 1:numel(cand)
            q = p;  q(j) = cand(c);
            screen(c) = GenNorm_NLL(X, q, kComp);
        end
        [~, ord] = sort(screen);
        free = ~held;  free(j) = false;
        for c = ord(1:min(nKeep, end))'
            q = p;  q(j) = cand(c);
            f = @(r) GenNorm_NLL(X, fill_free(q, free, r), kComp);
            try
                [r, fr] = fmincon(f, q(free), [], [], [], [], ...
                    LB(free), UB(free), [], opt);
            catch
                continue
            end
            if isfinite(fr) && fr < nll - 1e-12
                p = fill_free(q, free, r);  nll = fr;  improved = true;
            end
        end
    end
    if ~improved
        break
    end
end
end

function q = fill_free(q, free, r)
q(free) = r;
end

function a = cusp_alpha()
a = 1.5;
end

% ============================================================
%  LIKELIHOOD
% ============================================================
function nll = GenNorm_NLL(x, p, kComp)
P   = unpack(p, kComp);
nll = -sum(GenNorm_LogPDF(x, P));
if ~isfinite(nll)
    nll = Inf;
end
end

function logf = GenNorm_LogPDF(x, Params)
% Log of the mixture density at x (column), by log-sum-exp.
% Params: K x 4 [mu, sigma, alpha, weight].
mu = Params(:, 1)';  sigma = Params(:, 2)';
alpha = Params(:, 3)';  w = Params(:, 4)';
log_const = log(w) + log(alpha) - log(2) - log(sigma) - gammaln(1 ./ alpha);
L = log_const - (abs(x(:) - mu) ./ sigma) .^ alpha;      % n x K
Lmax = max(L, [], 2);
logf = Lmax + log(sum(exp(L - Lmax), 2));
logf(isnan(logf)) = -Inf;
end

function P = unpack(p, kComp)
mu    = p(1:kComp);
sigma = exp(p(kComp+1:2*kComp));
alpha = exp(p(2*kComp+1:3*kComp));
if kComp > 1
    v = [p(3*kComp+1:end), 0];
    v = v - max(v);
    w = exp(v) / sum(exp(v));
else
    w = 1;
end
P = [mu(:), sigma(:), alpha(:), w(:)];
end

% ============================================================
%  STARTING POINTS
% ============================================================
function seeds = GenerateSeeds(x, kComp, sig_min, max_logit)
% K-means centres x alpha in {1, 2, 5}; quantile spread; one broad start;
% and, for K > 1, random starts. Uses the global stream, which the caller
% has replaced with a private seeded one.
seeds = [];
std_global = std(x);
v_equal = zeros(1, kComp - 1);

try
    [idx, C] = kmeans(x, kComp, 'Replicates', 3, 'MaxIter', 100, ...
        'Display', 'off');
    [mu_k, ord] = sort(C'); %#ok<TRSRT> C is K x 1; a row is wanted
    sig_k = zeros(1, kComp);
    w_k   = zeros(1, kComp);
    for j = 1:kComp
        cluster_data = x(idx == ord(j));
        if numel(cluster_data) > 1
            sig_k(j) = std(cluster_data);
        else
            sig_k(j) = std_global / kComp;
        end
        w_k(j) = numel(cluster_data) / numel(x);
    end
    sig_k = max(sig_k, sig_min);
    w_k   = max(w_k, 0.01);  w_k = w_k / sum(w_k);
    if kComp > 1
        v_k = log(w_k(1:end-1)) - log(w_k(end));
        v_k = max(min(v_k, max_logit), -max_logit);
    else
        v_k = [];
    end
    for a_val = [1.0, 2.0, 5.0]
        seeds = [seeds; mu_k, log(sig_k), log(a_val)*ones(1, kComp), v_k]; %#ok<AGROW>
    end
catch
    % k-means failed; the remaining starts still apply
end

if kComp > 1
    mu_q  = quantile(x, linspace(0.1, 0.9, kComp));
    sig_q = repmat(std_global / kComp, 1, kComp);
    seeds = [seeds; mu_q, log(sig_q), log(2)*ones(1, kComp), v_equal];
end

mu_g = median(x) + linspace(-1, 1, kComp) * std_global * 0.1;
seeds = [seeds; sort(mu_g), log(std_global)*ones(1, kComp), ...
    log(1.5)*ones(1, kComp), v_equal];

if kComp > 1
    nRandom = 4 * kComp;
    n = numel(x);
    for r = 1:nRandom
        mu_r  = sort(x(randi(n, 1, kComp)))';
        sig_r = std_global * (0.1 + 0.9*rand(1, kComp));
        a_r   = exp(log(0.7) + (log(5) - log(0.7)) * rand(1, kComp));
        seeds = [seeds; mu_r, log(max(sig_r, sig_min)), log(a_r), ...
            randn(1, kComp - 1)]; %#ok<AGROW>
    end
end
end

% ============================================================
%  UTILITIES
% ============================================================
function [Params, order] = ReconstructParams(p, kComp)
P = unpack(p, kComp);
[~, order] = sortrows(P, 1);
Params = P(order, :);
end

% ============================================================
%  STANDARD ERRORS AND CONFIDENCE INTERVALS
% ============================================================
function U = gn_uncertainty(Z, stats, k, loc, scl, level, doProfile)
% Standard errors, Wald CIs and (optionally) profile-likelihood CIs, all in
% the units of X and in the row order of Params.
zc = sqrt(2) * erfinv(level);
nanK = NaN(k, 4);
U = struct('SE', nanK, 'Cov', NaN(4*k), 'CILower', nanK, ...
    'CIUpper', nanK, 'CIOpen', false(k, 4), 'CIMethod', 'Wald', ...
    'CIWaldLower', nanK, 'CIWaldUpper', nanK, 'SEHessian', nanK, ...
    'SEOPG', nanK, 'ProfileBest', -Inf, 'CuspMu', false(k, 1));
if isempty(stats.p)
    return
end
p = stats.p;  np = numel(p);  order = stats.order;
free = ~stats.fixed;

% Per-observation scores (analytic) and the observed information (central
% differences of the analytic total score).
S = gn_scores(Z, p, k);
g = @(q) sum(gn_scores(Z, q, k), 1)';
H = zeros(np);
for c = 1:np
    h = 1e-5 * max(1, abs(p(c)));
    e = zeros(size(p));  e(c) = h;
    H(:, c) = (g(p + e) - g(p - e)) / (2*h);
end
H = (H + H') / 2;
covH = restricted_inverse(-H, free);
covO = restricted_inverse(S' * S, free);

J = phys_jacobian(p, k, scl);
[seH, CovH] = phys_se(J, covH, k, order);
seO = phys_se(J, covO, k, order);
U.SEHessian = seH;
U.SEOPG     = seO;
U.SE        = seH;
U.Cov       = CovH;
% mu of a component in the cusp regime: the Hessian is meaningless there
% (Wald coverage of mu measured 0.04 at alpha = 0.7 and 0.33 at alpha = 1),
% so SE takes the OPG value, which is less wrong (0.80 and 0.91). Profile
% intervals are the reliable choice for those mu (see the header).
P0 = ReconstructParams(p, k);
cusp = P0(:, 3) < cusp_alpha();
U.SE(cusp, 1) = seO(cusp, 1);
U.CuspMu = cusp;

P = ReconstructParams(p, k);
P(:, 1) = loc + scl * P(:, 1);
P(:, 2) = scl * P(:, 2);
[lo, hi] = wald_ci(P, U.SE, zc);
U.CIWaldLower = lo;
U.CIWaldUpper = hi;
U.CILower = lo;
U.CIUpper = hi;

if doProfile
    [plo, phi, popen, best] = gn_profile(Z, stats, k, loc, scl, ...
        zc^2, covH);
    U.CILower(:, 1:3) = plo;
    U.CIUpper(:, 1:3) = phi;
    U.CIOpen(:, 1:3)  = popen;
    U.CIMethod = 'profile (mu, sigma, alpha); Wald (weight)';
    U.ProfileBest = best;
end
end

function S = gn_scores(x, p, k)
% n x (4k-1) per-observation scores of log f in the optimizer coordinates
% [mu, log sigma, log alpha, weight logits]. With u = |x - mu|/sigma and
% responsibilities r:
%   d/dmu_j        = r_j * alpha_j u^(alpha_j - 1) sign(x - mu_j) / sigma_j
%   d/dlog sigma_j = r_j * (alpha_j u^alpha_j - 1)
%   d/dlog alpha_j = r_j * (1 + psi(1/alpha_j)/alpha_j - alpha_j u^alpha_j log u)
%   d/dv_m         = r_m - w_m
% At u = 0 (a data point exactly at mu) the mu score is set to 0, the
% subgradient at the cusp.
P = unpack(p, k);
mu = P(:, 1)';  sg = P(:, 2)';  al = P(:, 3)';  w = P(:, 4)';
x = x(:);
u = abs(x - mu) ./ sg;
L = log(w) + log(al) - log(2) - log(sg) - gammaln(1 ./ al) - u .^ al;
r = exp(L - max(L, [], 2));
r = r ./ sum(r, 2);
ua  = u .^ al;
um1 = u .^ (al - 1);
um1(u == 0) = 0;
ulog = ua .* log(u);
ulog(u == 0) = 0;
dmu = r .* al .* um1 .* sign(x - mu) ./ sg;
dls = r .* (al .* ua - 1);
dla = r .* (1 + psi(1 ./ al) ./ al - al .* ulog);
dv  = r(:, 1:k-1) - w(1:k-1);
S = [dmu, dls, dla, dv];
end

function C = restricted_inverse(I, free)
% Inverse of the information over the free coordinates; coordinates on an
% active bound get zero variance (they are reported as NaN later).
C = NaN(size(I));
If = I(free, free);
[R, flag] = chol(If);
if flag > 0 || any(~isfinite(If(:)))
    return
end
Ri = R \ eye(size(If));
C = zeros(size(I));
C(free, free) = Ri * Ri';
C(~free, :) = NaN;
C(:, ~free) = NaN;
end

function J = phys_jacobian(p, k, scl)
% d[mu; sigma; alpha; w] / dp, physical units, components in fitted order.
np = numel(p);
P = unpack(p, k);
J = zeros(4*k, np);
for j = 1:k
    J(j, j)             = scl;          % mu = loc + scl*p_j
    J(k + j, k + j)     = P(j, 2)*scl;  % sigma = scl*exp(p)
    J(2*k + j, 2*k + j) = P(j, 3);      % alpha = exp(p)
end
w = P(:, 4);
for j = 1:k
    for m = 1:k-1
        J(3*k + j, 3*k + m) = w(j) * ((j == m) - w(m));
    end
end
end

function [SE, Cov] = phys_se(J, covP, k, order)
% Physical covariance J*covP*J' with NaN where a coordinate is fixed,
% reordered to the rows of Params. Cov is 4k x 4k in the order
% [mu_1..mu_k, sigma_1..k, alpha_1..k, w_1..k] (rows of Params).
fixedCols = isnan(diag(covP))';     % a fixed coordinate's own variance
Cz = covP;
Cz(isnan(Cz)) = 0;
Cov = J * Cz * J';
touched = any(J(:, fixedCols) ~= 0, 2);
Cov(touched, :) = NaN;
Cov(:, touched) = NaN;
if all(isnan(covP(:)))
    Cov(:) = NaN;
end
idx = [order(:); k + order(:); 2*k + order(:); 3*k + order(:)];
Cov = Cov(idx, idx);
se  = sqrt(max(diag(Cov), 0));
se(isnan(diag(Cov))) = NaN;
SE  = reshape(se, k, 4);
if k == 1
    SE(4) = NaN;            % a single weight is fixed at 1
end
end

function [lo, hi] = wald_ci(P, SE, zc)
% mu directly; sigma and alpha on the log scale; weights on the logit scale.
lo = NaN(size(P));  hi = lo;
lo(:, 1) = P(:, 1) - zc*SE(:, 1);
hi(:, 1) = P(:, 1) + zc*SE(:, 1);
for c = 2:3
    lo(:, c) = P(:, c) .* exp(-zc*SE(:, c) ./ P(:, c));
    hi(:, c) = P(:, c) .* exp( zc*SE(:, c) ./ P(:, c));
end
w = P(:, 4);
sl = SE(:, 4) ./ (w .* (1 - w));            % SE of logit(w)
lw = log(w ./ (1 - w));
lo(:, 4) = 1 ./ (1 + exp(-(lw - zc*sl)));
hi(:, 4) = 1 ./ (1 + exp(-(lw + zc*sl)));
end

function [lo, hi, isOpen, bestSeen] = gn_profile(Z, stats, k, loc, scl, ...
    chiCrit, covH)
% Profile-likelihood CIs for mu, sigma and alpha of each component: fix one
% optimizer coordinate, re-fit the others (fmincon, warm-started), and find
% where the log-likelihood falls chiCrit/2 below the maximum. Rows follow
% Params. An end that reaches the coordinate's box without crossing is OPEN:
% mu -> -/+Inf, sigma -> its floor or Inf, alpha -> 0.5 or 10.
p0 = stats.p;  LB = stats.LB;  UB = stats.UB;  opt = stats.opt;
order = stats.order;
llMax = stats.LogLike;
thr = llMax - chiCrit/2;
nll = @(q) GenNorm_NLL(Z, q, k);
lo = NaN(k, 3);  hi = lo;  isOpen = false(k, 3);
bestSeen = -Inf;
warmP = p0;

for blk = 1:3
    for j = 1:k
        c = (blk - 1)*k + j;
        step = sqrt(covH(c, c));
        if ~(isfinite(step) && step > 0)
            step = 0.1;
        end
        ends = NaN(1, 2);  opens = false(1, 2);
        for sd = [-1, 1]
            warmP = p0;
            prev = p0(c);  fPrev = llMax - thr;     % > 0 at the estimate
            found = false;
            for kk = 0:30
                v = p0(c) + sd*step*2^kk;
                v = min(max(v, LB(c)), UB(c));
                pl = prof(v);
                if pl < thr
                    ends((sd + 3)/2) = crossing(prev, fPrev, v, pl - thr);
                    found = true;
                    break
                end
                if v == LB(c) || v == UB(c)
                    break
                end
                prev = v;  fPrev = pl - thr;
            end
            if ~found
                ends((sd + 3)/2) = v;
                opens((sd + 3)/2) = true;
            end
        end
        row = find(order == j);
        switch blk
            case 1   % mu: the box is unbounded, so open means infinite
                val = loc + scl*ends;
                if opens(1), val(1) = -Inf; end
                if opens(2), val(2) =  Inf; end
            case 2   % sigma: open below = its floor, open above = Inf
                val = scl*exp(ends);
                if opens(2), val(2) = Inf; end
            case 3   % alpha: open = the bound 0.5 or 10 it reached
                val = exp(ends);
        end
        lo(row, blk) = val(1);
        hi(row, blk) = val(2);
        isOpen(row, blk) = any(opens);
    end
end

    function ll = prof(v)
        % Maximum log-likelihood with coordinate c held at v: the better of
        % two nuisance fits, one started from the estimate and one from the
        % last profile point, so a poor warm start cannot leave a value
        % too low.
        fr = true(size(p0));  fr(c) = false;
        ll = -Inf;  bestQ = warmP;
        for start = {p0, warmP}
            q = start{1};  q(c) = v;
            f = @(r) nll(fill_free(q, fr, r));
            try
                [r, fv] = fmincon(f, q(fr), [], [], [], [], LB(fr), UB(fr), ...
                    [], opt);
            catch
                continue
            end
            if isfinite(fv) && -fv > ll
                ll = -fv;
                bestQ = fill_free(q, fr, r);
            end
        end
        % The nuisance fit stalls at the mu cusps exactly as the estimate
        % does (measured: without this, profile CIs for sigma and alpha at
        % alpha0 = 0.7 covered 0.63 instead of 0.95), so polish it too.
        if isfinite(ll) && any(exp(bestQ(2*k+1:3*k)) < cusp_alpha())
            [bestQ, nq] = cusp_polish(Z, bestQ, -ll, k, LB, UB, opt, ~fr);
            ll = -nq;
        end
        if isfinite(ll)
            warmP = bestQ;
        end
        bestSeen = max(bestSeen, ll);
    end

    function e = crossing(a, fa, b, fb)
        % Root of prof - thr between a (fa >= 0) and b (fb < 0) by the
        % Illinois false-position method, using only the values it has
        % evaluated, so a profile that is not exactly reproducible on a
        % second evaluation cannot break the bracket.
        side = 0;
        e = b;
        for it = 1:40
            e  = b - fb*(b - a)/(fb - fa);
            fe = prof(e) - thr;
            if abs(fe) < 1e-7 || abs(b - a) < 1e-9*max(1, abs(e))
                return
            end
            if fe >= 0
                a = e;  fa = fe;
                if side == -1, fb = fb/2; end
                side = -1;
            else
                b = e;  fb = fe;
                if side == 1, fa = fa/2; end
                side = 1;
            end
        end
    end
end

function restore_state(stream, state, warnState)
RandStream.setGlobalStream(stream);
stream.State = state;
warning(warnState);
end

function plot_results(x, res)
figure('Color', 'w');  hold on;
histogram(x, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
    'EdgeColor', [0.2 0.2 0.2]);
xx = linspace(min(x), max(x), 500)';
colors = lines(numel(res));
for i = 1:numel(res)
    yy = exp(GenNorm_LogPDF(xx, res(i).Params));
    plot(xx, yy, 'LineWidth', 2, 'Color', colors(i, :), ...
        'DisplayName', res(i).Model);
end
legend('show');  grid on;  xlabel('Data');  ylabel('PDF');
title('Generalized Normal Fit');
end

% HISTORY
%   Dec 2025 - Jan 2026  Written by Sheyum.
%   Nov 2026  Audit against an independent profile likelihood (Claude
%     Code). Fixed: AIC held AICc (now AIC is plain and AICc is separate);
%     tiny samples returned garbage (n = 2..7 at K = 2 gave finite AICc of
%     -25 to -97, n = 8 gave Inf; now K needs n >= 4K + 1); alpha sat on its
%     bounds silently (uniform data at 10, Cauchy at 0.5; now flagged);
%     sigma floor had an absolute 1e-6 part and the fit ran in raw units
%     (now z-units with exact back-transforms); kmeans and the starts used
%     and advanced the caller's global random stream (now a private seeded
%     one, restored on exit); the warning state was forced on rather than
%     restored; the 'Select' message printed the p-value of the rejected
%     step as if it belonged to the selected model; failures were silent
%     (now Converged/ExitFlag and a NotConverged warning); random starts
%     and a polish run added for K > 1. Measured results are in the
%     LIMITATIONS above and in test_GenNorm_MLE.m.
%   Nov 2026  Standard errors (observed information; OPG for mu in the cusp
%     regime), Wald CIs, and optional profile-likelihood CIs ('Profile')
%     added, with coverage measured by simulation (LIMITATIONS).
