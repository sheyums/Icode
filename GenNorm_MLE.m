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
%
%   WARNINGS
%       GenNorm_MLE:NonfiniteDataRemoved   NaN/Inf dropped from X.
%       GenNorm_MLE:TooFewObservations     some K dropped (too few points).
%       GenNorm_MLE:AlphaAtBound           a shape parameter on its bound.
%       GenNorm_MLE:SigmaAtFloor           a scale parameter on its floor.
%       GenNorm_MLE:NotConverged           no start converged for some K.
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
%         in mu, where a gradient-based optimizer is not exact.
%       - No standard errors or confidence intervals are returned.
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
    'AlphaAtBound', [], 'SigmaAtFloor', []), numel(models_to_test), 1);

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
        optPolish);
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
end

function [p, nll] = cusp_polish(X, p, nll, kComp, LB, UB, opt)
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
nCand = 41;  nKeep = 3;
xs = sort(X);
for iter = 1:3
    improved = false;
    for j = 1:kComp
        if exp(p(2*kComp + j)) >= cusp_alpha()
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
        free = true(size(p));  free(j) = false;
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
