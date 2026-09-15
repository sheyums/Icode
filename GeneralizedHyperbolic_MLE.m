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
%   data      : real finite numeric vector, at least two distinct values.
%   modelFlag : 'GH'   – fit the full 5-parameter Generalized Hyperbolic.
%               'NIG'  – fit the 4-parameter Normal-Inverse Gaussian
%                        (GH with lambda fixed at -0.5).
%               'AUTO' – (default) fit both models, then select the better
%                        one via a likelihood-ratio test AND BIC.
%   mode      : pass 'PLOT' to display a histogram with the fitted density.
%
% PHYSICAL PARAMETERS
%   The GH distribution is parameterized by five physical quantities:
%     mu     : location (real)
%     lambda : tail index (real; controls the shape of the tails)
%     alpha  : tail heaviness (> |beta|)
%     beta   : skewness (real)
%     delta  : scale (> 0)
%   The NIG model fixes lambda = -0.5.
%
%   Special cases of the GH family:
%     lambda = -0.5          : Normal-Inverse Gaussian (NIG)
%     lambda =  0.5          : Hyperbolic
%     lambda -> 0            : approaches the Variance-Gamma limit
%     alpha -> infinity      : approaches the Normal distribution
%
% INTERNAL OPTIMIZATION COORDINATES
%   To enforce positivity and ordering constraints without explicit bounds,
%   estimation is performed in a transformed coordinate system:
%
%     GH  theta = [lambda; beta; psi; log_delta; mu]   (5-vector)
%     NIG theta = [beta; psi; log_delta; mu]            (4-vector)
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
%   3. For GH: an additional warm start is always attempted using the
%      NIG MLE as the initial point (with lambda freed to -0.5).
%      The warm-start solution replaces the grid solution if it converged
%      when the grid did not, or if both converged and it is better.
%   4. The best converged solution is selected; if none converged, the
%      best finite solution is used with a warning.
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
%   vectors whose Bessel arguments exceed 1e10 are therefore declared
%   infeasible (NLL = Inf): beyond that point the error swamps any genuine
%   difference in likelihood and an optimizer will "improve" the fit by
%   chasing floating-point noise.
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
%   Wald standard errors, with NaN/Inf guards, followed by fzero on the located
%   bracket and sign-only bisection as a fallback when an endpoint is infinite.
%
% HESSIAN WARNING
%   results.hessian is the Hessian returned by fmincon in TRANSFORMED
%   coordinates [lambda, beta, psi, log_delta, mu] (or the NIG subset).
%   Standard errors derived directly from this matrix are WRONG for
%   physical alpha and delta because those involve nonlinear transforms.
%   To obtain physical-space standard errors, apply the delta method:
%     Var(g(theta)) ~ J * inv(H) * J'
%   where J is the Jacobian of (alpha, delta) with respect to theta.
%   The profile CIs do not have this limitation.
%
%   SECOND CAVEAT, specific to the interior-point algorithm: fmincon returns
%   a Hessian of the LAGRANGIAN (a quasi-Newton approximation), not the
%   observed information matrix of the objective.  When a bound is active it
%   includes the constraint contribution.  It is adequate for diagnosing
%   ill-conditioning but should not be treated as the information matrix.
%   Prefer the profile CIs, which make no such assumption.
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
%     hessian     – Hessian in transformed coordinates (see warning above).
%     hessian_coordinates – same as theta_names, documents what hessian rows/cols mean.
%     profile_ci  – struct with 95% profile CI [lo hi] for each physical
%                   parameter: lambda, beta, alpha, delta, mu.
%
%                   A bound is one of three things, and they mean different
%                   things:
%                     finite  – the profile crossed the threshold there.
%                     +/-Inf  – the profile was evaluated all the way out and
%                               never crossed: the likelihood is flat in that
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
%     output      – fmincon output struct for the selected solution.
%     converged   – logical; true when exitflag > 0.  NOTE this reports only
%                   what the optimizer said; see 'valid' for whether the
%                   answer is usable.
%     valid       – logical; true when the fitted parameters yield a finite,
%                   positive density at every data point AND that density
%                   reproduces the reported LogLike.  A run can be
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
%                   though every numerical check passes.
%     delta_at_bound – logical; true when delta sits on its lower bound, i.e.
%                   the solution is a CONSTRAINED optimum.  Profile CIs and
%                   Hessian-based standard errors assume an interior solution
%                   and are unreliable in that case.
%     ModelSelection – (AUTO only) string describing which model was chosen and why.
%     LRT         – (AUTO only) struct with LRstat, df, pValue, LogLike_GH, LogLike_NIG.
%
% WARNINGS RAISED
%   GeneralizedHyperbolic_MLE:DensityNotEvaluable
%       Fitted parameters give NaN or non-positive density values.
%   GeneralizedHyperbolic_MLE:LikelihoodInconsistent
%       Reported LogLike disagrees with the density it implies.
%   GeneralizedHyperbolic_MLE:DeltaAtBound
%       delta reached its lower bound; constrained optimum.
%   GeneralizedHyperbolic_MLE:NoGainOverNormal
%       Fit does not beat a Normal by AIC; shape parameters unidentifiable.
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
%
% Originally written by Sheyum but beefed up/stabilized/optimized 
% by Claude, Gemini and ChatGPT, Dec 2025.
% Reparameterization, exact alpha CI, and robustness improvements: Jul 2026.
% Updated to fmincon with bounds to prevent component collapse: Aug 2026.

% ---- Validate and normalize inputs -----------------------------------------
if nargin < 2 || isempty(modelFlag), modelFlag = 'AUTO'; end
if nargin < 3 || isempty(mode),      mode      = '';     end

if ~(ischar(modelFlag) || (isstring(modelFlag) && isscalar(modelFlag)))
    error('GeneralizedHyperbolic_MLE:InvalidModelFlag', ...
        'MODELFLAG must be ''GH'', ''NIG'', or ''AUTO''.');
end
modelFlag = upper(char(modelFlag));
if ~ismember(modelFlag, {'GH','NIG','AUTO'})
    error('GeneralizedHyperbolic_MLE:InvalidModelFlag', ...
        'MODELFLAG must be ''GH'', ''NIG'', or ''AUTO''.');
end
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('GeneralizedHyperbolic_MLE:InvalidMode', ...
        'MODE must be a character vector or scalar string.');
end
mode = char(mode);

validateattributes(data, {'numeric'}, ...
    {'vector','real','finite','nonempty'}, mfilename, 'data', 1);
data = data(:);
if numel(data) < 2
    error('GeneralizedHyperbolic_MLE:TooFewObservations', ...
        'DATA must contain at least two observations.');
end
if ~(std(data) > 0)
    error('GeneralizedHyperbolic_MLE:ConstantData', ...
        'DATA must contain at least two distinct values.');
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
switch modelFlag
    case 'GH'
        fprintf('Fitting Generalized Hyperbolic (GH) model...\n');
        results = fit_GH(data);

    case 'NIG'
        fprintf('Fitting Normal-Inverse Gaussian (NIG) model...\n');
        results = fit_NIG(data);

    case 'AUTO'
        fprintf('Running automatic model selection (GH versus NIG)...\n');
        resGH  = fit_GH(data);
        resNIG = fit_NIG(data);

        % Likelihood-ratio test: H0 = NIG (lambda = -0.5).
        % GH has one extra free parameter (lambda), so df = 1.
        LR   = max(2*(resGH.LogLike - resNIG.LogLike), 0);
        pval = chi2cdf(LR, 1, 'upper');  % upper tail; accurate for small p

        % Prefer GH only when the improvement is both statistically
        % significant (LRT) AND reflected in a lower BIC (penalizes
        % the extra parameter).

        if pval < 0.05 && resGH.BIC < resNIG.BIC
            results = resGH;
            results.ModelSelection = sprintf( ...
                'GH selected (LRT p=%.4g; BIC %.3f versus %.3f).', ...
                pval, resGH.BIC, resNIG.BIC);
        else
            results = resNIG;
            results.ModelSelection = sprintf( ...
                'NIG selected (LRT p=%.4g; BIC %.3f versus %.3f).', ...
                pval, resNIG.BIC, resGH.BIC);
        end
        results.LRT = struct('LRstat', LR, 'df', 1, 'pValue', pval, ...
            'LogLike_GH', resGH.LogLike, 'LogLike_NIG', resNIG.LogLike);
end

if contains(upper(mode), 'PLOT'), plot_results(data, results); end
end

% =========================================================================
%                          TOP-LEVEL MODEL FITTERS
% =========================================================================

function results = fit_GH(data)
% Fit the full 5-parameter GH model.

% Step 1: multi-start search over a parameter grid.
[theta, nll, hess, exitflag, output] = gh_multistart(data, false);

% Step 2: always attempt a warm start from the NIG MLE.
% The NIG solution is a valid point in GH space (lambda = -0.5), so it
% provides a useful starting point in a region where the GH likelihood
% is known to be reasonable.  The warm start is accepted over the grid
% solution if it converged when the grid did not, or if both have equal
% convergence status and the warm start achieves a lower NLL.
[thetaNIG, ~, ~, ~, ~] = gh_multistart(data, true);
theta0 = [-0.5; thetaNIG];   % prepend lambda = -0.5 to NIG theta
opts = fit_options(2500, 1e-7);

% Same delta floor as the grid search, so the warm start explores the
% identical feasible region.
[lb, ub] = theta_bounds(data, false);

try
    [tr, nr, er, orr, ~, ~, hr] = fmincon( ...
        @(t) gh_nll_wrapper(t, data, false), theta0, ...
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
    % Retain the best multistart solution if the warm start crashes.
end

p  = theta_to_params(theta, false);
ci = gh_profile_ci(theta, data, nll, false, hess);
results = assemble_results( ...
    'Generalized Hyperbolic (mu, lambda, alpha, beta, delta)', ...
    p, theta, nll, hess, ci, exitflag, output, data, false);
results = flag_delta_at_bound(results, data, false);
end

function results = fit_NIG(data)
% Fit the 4-parameter NIG model (lambda fixed at -0.5).
[theta, nll, hess, exitflag, output] = gh_multistart(data, true);
p  = theta_to_params(theta, true);
ci = gh_profile_ci(theta, data, nll, true, hess);
results = assemble_results( ...
    'Normal-Inverse Gaussian (mu, lambda, alpha, beta, delta)', ...
    p, theta, nll, hess, ci, exitflag, output, data, true);
results = flag_delta_at_bound(results, data, true);
end

% =========================================================================
%                         MULTI-START OPTIMIZATION
% =========================================================================

function [bestTheta, bestNLL, bestHess, bestExit, bestOut] = gh_multistart(data, isNIG)
muMed  = median(data);
muMean = mean(data);
[counts, edges] = histcounts(data, 'Normalization', 'pdf');
[~, im] = max(counts);
muMode = mean(edges(im:im+1));
muGrid = unique([muMed muMean muMode]);

s0 = std(data);
deltaGrid  = s0 * [0.5 1 1.5];
betaGrid   = [0  -0.5/s0  0.5/s0];
if isNIG
    lambdaGrid = -0.5;
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
    if isNIG
        t = [grid(i,2); psi0; log(grid(i,3)); grid(i,4)];
    else
        t = [grid(i,1); grid(i,2); psi0; log(grid(i,3)); grid(i,4)];
    end
    starts{i} = t;
    vals(i)   = gh_nll_wrapper(t, data, isNIG);
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
[lb, ub] = theta_bounds(data, isNIG);

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
            @(x) gh_nll_wrapper(x, data, isNIG), starts{ord(j)}, ...
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

function nll = gh_nll_wrapper(theta, data, isNIG)
% Convert transformed coordinates to physical params, then evaluate NLL.
p = theta_to_params(theta, isNIG);
nll = gh_nll_from_params(p, data);
end

function nll = gh_nll_from_params(p, data)
% Evaluate the GH negative log-likelihood given a physical parameter struct.
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
% Bessel functions are evaluated in scaled form:
%   besselk(nu, x, 1) = besselk(nu, x) * exp(x)
% so  log K(nu, x) = log(besselk(nu, x, 1)) - x.
% This avoids overflow/underflow for large arguments.

% Guard against any non-finite or physically invalid parameter values.
if ~all(isfinite([p.lambda p.beta p.gamma p.alpha p.delta p.mu])) || ...
        p.gamma <= 0 || p.alpha <= 0 || p.delta <= 0
    nll = Inf;  return
end

xm   = data - p.mu;
d2   = p.delta^2 + xm.^2;
z    = sqrt(d2);

% Arguments to the two Bessel evaluations.
arg0 = p.delta * p.gamma;   % argument of the normalizing Bessel K_lambda
argz = p.alpha * z;         % arguments of the data Bessel K_{lambda-0.5}

if ~isfinite(arg0) || arg0 <= 0 || any(~isfinite(argz) | argz <= 0)
    nll = Inf;  return
end

% ---- Numerical-validity guard on the Bessel arguments --------------------
% The log-density below forms (log K - arg) for both Bessel terms, so it adds
% and subtracts quantities of order ARG. In double precision the absolute
% error of that cancellation is about arg*eps per observation, and roughly
% n*arg*eps in the summed NLL. Once arg exceeds ~1e10 that error swamps any
% real difference in likelihood, and an optimizer will happily "improve" the
% fit by chasing pure floating-point noise.
%
% This is not hypothetical. Without this guard, NIG fits to n >= 1000 ran away
% to alpha ~ 4e52, delta ~ 1e36, arg ~ 4e88, reporting a log-likelihood of
% +4.5e19 (impossible) while results.Fit evaluated to all NaN. Declaring such
% parameters infeasible confines the search to the region where the likelihood
% is actually computable. Legitimate fits stay far below the cap: even
% delta = 50 with alpha = 3 gives arg ~ 150.
if arg0 > MAX_BESSEL_ARG() || any(argz > MAX_BESSEL_ARG())
    nll = Inf;  return
end

% Scaled Bessel evaluations; check for non-positive results (unphysical).
K0 = besselk(p.lambda,       arg0, 1);
Kz = besselk(p.lambda - 0.5, argz, 1);
if ~isfinite(K0) || K0 <= 0 || any(~isfinite(Kz) | Kz <= 0)
    nll = Inf;  return
end

% Log-density: scaled Bessel logs are (log(K_scaled) - argument).
logpdf = p.lambda*log(p.gamma) - p.lambda*log(p.delta) - 0.5*log(2*pi) ...
    - (p.lambda - 0.5)*log(p.alpha) - (log(K0) - arg0) ...
    + (p.lambda/2 - 0.25).*log(d2) + (log(Kz) - argz) + p.beta*xm;

nll = -sum(logpdf);
if ~isfinite(nll) || ~isreal(nll), nll = Inf; end
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

function [lb, ub] = theta_bounds(data, isNIG)
% Box constraints in transformed coordinates. Only log_delta is bounded; every
% other coordinate is free. Index of log_delta is 3 for NIG, 4 for GH.
mld = log(delta_lower_bound(data));
if isNIG
    lb = [-Inf; -Inf; mld;  -Inf];
    ub = [ Inf;  Inf; Inf;   Inf];
else
    lb = [-Inf; -Inf; -Inf; mld;  -Inf];
    ub = [ Inf;  Inf;  Inf; Inf;   Inf];
end
end

function results = flag_delta_at_bound(results, data, isNIG)
% Report when the fit is sitting on the delta floor. Such a point is a
% CONSTRAINED optimum, not an MLE: the profile CIs are computed without the
% bound, and standard errors from the Hessian assume an interior solution.
% Returning it unannounced would misrepresent a boundary solution as a fit.
idx = 4 - double(isNIG);          % position of log_delta in theta
mld = log(delta_lower_bound(data));
results.delta_at_bound = isfinite(results.theta_hat(idx)) && ...
    results.theta_hat(idx) <= mld + 1e-8;
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
    if isnan(fm)
        return          % cannot evaluate here; the bound is genuinely unknown
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

function v = MAX_BESSEL_ARG()
% Largest Bessel argument for which the log-density is numerically meaningful.
% See the derivation in gh_nll_from_params: the summed NLL carries an absolute
% error of roughly n*arg*eps, so 1e10 keeps that below ~1e-2 even for n = 1e6,
% while sitting many orders of magnitude above any legitimate fit.
v = 1e10;
end

% =========================================================================
%                     PROFILE LIKELIHOOD CONFIDENCE INTERVALS
% =========================================================================

function ci = gh_profile_ci(theta, data, nllMin, isNIG, hess)
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
p    = theta_to_params(theta, isNIG);
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

if isNIG
    ci.lambda = [-0.5 -0.5];   % lambda is fixed for NIG
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
            idx  = 2 - double(isNIG);
            q0   = p.beta;
            step = max(0.05/s, 0.1*max(1/s, abs(q0)));
        case 'delta'
            idx  = 4 - double(isNIG);
            q0   = log(p.delta);
            step = 0.1;
        case 'mu'
            idx  = 5 - double(isNIG);
            q0   = p.mu;
            step = max(0.05*s, 0.1*max(s, abs(q0)));
    end

    % Prefer a Hessian-derived step when one is available; fall back to the
    % heuristic above when the Hessian is unusable (not positive definite,
    % ill-conditioned, or absent).
    if idx <= numel(seT) && isfinite(seT(idx)) && seT(idx) > 0
        step = 2 * seT(idx);
    end

    fun = @(q) regular_profile_dev(q, idx, theta, data, nllMin, crit, isNIG, opts);
    lo  = find_root_robust(fun, q0, -1, step);
    hi  = find_root_robust(fun, q0,  1, step);

    if strcmp(name, 'delta')
        ci.(name) = exp(sort([lo hi]));   % back-transform from log space
    else
        ci.(name) = sort([lo hi]);
    end
end

% Alpha is profiled exactly by fixing it physically and re-optimizing
% all other parameters (see profile_alpha_ci for details).
ci.alpha = profile_alpha_ci(theta, p.alpha, data, nllMin, crit, isNIG, opts);
end

function d = regular_profile_dev(q, idx, theta, data, nllMin, crit, isNIG, opts)
% Profile deviance for parameters other than alpha.
% Fix theta(idx) = q, re-optimize the remaining parameters,
% and return  nll_profile - nllMin - crit.
%
% The nuisance optimization is CONSTRAINED to the same feasible region as the
% point estimate. Without this the profile would be free to push delta below
% the floor that the fit itself was forbidden to cross, so the CI could report
% support for parameter values the estimator was never allowed to return.
mask       = true(numel(theta), 1);
mask(idx)  = false;
base       = theta;
base(idx)  = q;
obj = @(x) gh_nll_wrapper(insert_param(base, x, mask), data, isNIG);

[lbFull, ubFull] = theta_bounds(data, isNIG);
lb = lbFull(mask);
ub = ubFull(mask);

x0 = theta(mask);
x0 = min(max(x0, lb), ub);   % fmincon needs a feasible start
try
    [~, v] = fmincon(obj, x0, [], [], [], [], lb, ub, [], opts);
catch
    v = Inf;
end
d = v - nllMin - crit;
end

% ---- Exact profile CI for physical alpha -----------------------------------

function ci = profile_alpha_ci(theta, aHat, data, nllMin, crit, isNIG, opts)
% Profile alpha by searching in log(alpha) space.
% Using q = log(alpha) keeps alpha positive throughout the search and
% makes the step size scale-invariant regardless of alpha's magnitude.
q0  = log(aHat);
fun = @(q) alpha_profile_dev_log(q, theta, data, nllMin, crit, isNIG, opts);
lo  = find_root_robust(fun, q0, -1, 0.1);
hi  = find_root_robust(fun, q0,  1, 0.1);
ci  = exp(sort([lo hi]));
end

function d = alpha_profile_dev_log(logA, theta, data, nllMin, crit, isNIG, opts)
% Profile deviance as a function of log(alpha).
a = exp(logA);
if ~isfinite(a) || a <= 0, d = Inf;  return;  end
v = alpha_fixed_nll(a, theta, data, isNIG, opts);
if isfinite(v), d = v - nllMin - crit;  else, d = Inf;  end
end

function nll = alpha_fixed_nll(a, theta, data, isNIG, opts)
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
% explores the identical feasible region (see regular_profile_dev).
p = theta_to_params(theta, isNIG);

% Convert the current beta estimate to the xi parameterization.
% Clamp to avoid atanh(±1) = ±inf at the boundary.
r   = max(min(p.beta/a, 1-1e-8), -1+1e-8);
xi0 = atanh(r);

mld = log(delta_lower_bound(data));
if isNIG
    x0  = [xi0; log(p.delta); p.mu];
    lb  = [-Inf; mld; -Inf];
    ub  = [ Inf; Inf;  Inf];
    obj = @(x) alpha_core(a, -0.5, x(1), x(2), x(3), data);
else
    x0  = [p.lambda; xi0; log(p.delta); p.mu];
    lb  = [-Inf; -Inf; mld; -Inf];
    ub  = [ Inf;  Inf; Inf;  Inf];
    obj = @(x) alpha_core(a, x(1), x(2), x(3), x(4), data);
end
x0 = min(max(x0, lb), ub);   % fmincon needs a feasible start
try
    [~, nll] = fmincon(obj, x0, [], [], [], [], lb, ub, [], opts);
catch
    nll = Inf;
end
end

function nll = alpha_core(a, lambda, xi, logDelta, mu, data)
% NLL with alpha = a fixed and beta = a*tanh(xi).
t          = tanh(xi);
p.lambda   = lambda;
p.beta     = a * t;
p.gamma    = a * sqrt(max(1 - t^2, realmin));
p.alpha    = a;
p.delta    = exp(logDelta);
p.mu       = mu;
nll        = gh_nll_from_params(p, data);
end

% =========================================================================
%                         ROOT FINDER (PROFILE BOUNDARY)
% =========================================================================

function root = find_root_robust(fun, q0, direction, scale)
% Find the root of fun in the given direction from q0 using an
% exponentially expanding bracket search.
%
% Algorithm:
%   1. Evaluate f0 = fun(q0).  If f0 > 1e-6, bail immediately:
%      the profile deviance should be approximately -1.92 at the MLE,
%      so a positive value at q0 indicates a structural problem
%      (inconsistent nllMin/theta, failed nuisance optimization, etc.)
%      rather than numerical noise.
%   2. Expand the search point exponentially: current = q0 + direction*scale*1.5^k.
%   3. Skip NaN/Inf returns (optimizer failure at extreme parameter values).
%   4. When a sign change is detected between the last finite point and
%      the current point, call fzero to pin down the root.
%   5. Return NaN if no bracket is found within 31 steps (~scale * 1.9e5).

root = NaN;
f0 = fun(q0);
% f0 should be around -1.92 at the MLE; anything positive signals trouble.
if ~isfinite(f0) || f0 > 1e-6, return; end

previous = q0;
fp = f0;
sawFiniteBelow = false;   % did the profile stay below the threshold throughout?
for k = 0:30
    current = q0 + direction * scale * (1.5^k);
    fc = fun(current);
    if isnan(fc), continue; end   % skip true optimizer crashes (NaN)
    % Inf is not skipped: it means the profile has gone to infinity, which
    % is a valid signal that we have passed the CI boundary. fzero itself
    % cannot cope with an infinite endpoint, so local_bisect takes over below.
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
        tolX = max(1e-10, 1e-4 * max(abs(q0), scale));
        try
            root = fzero(fun, sort([previous current]), ...
                optimset('TolX', tolX));
        catch
            root = NaN;
        end
        if ~isfinite(root)
            % fzero needs FINITE values at both bracket ends. The profile is
            % legitimately +Inf wherever the nuisance optimization cannot be
            % evaluated, and an earlier comment here wrongly asserted that
            % fzero could bracket against Inf: it throws, the catch above
            % swallows it, and the bound is reported as NaN even though the
            % root is known to lie inside this bracket. Bisection needs only
            % the SIGN at each end, so it succeeds exactly where fzero cannot.
            root = local_bisect(fun, previous, current, tolX);
        end
        return
    end
    if isfinite(fc) && fc < 0
        sawFiniteBelow = true;
    end
    previous = current;
    fp = fc;
end

% No bracket within the search range.  Distinguish the two reasons, because
% they mean opposite things and used to be reported identically as NaN:
%
%   * The profile was evaluated successfully the whole way out and never rose
%     to the threshold.  The likelihood is flat in this direction, so the
%     parameter is not identified and the interval is genuinely OPEN.  Report
%     -Inf or +Inf, which says exactly that.
%   * Evaluation kept failing (every fc was NaN).  Nothing was learned, so the
%     bound is unknown and NaN remains the honest answer.
if sawFiniteBelow
    root = direction * Inf;
end
end

% =========================================================================
%                              UTILITIES
% =========================================================================

function p = theta_to_params(theta, isNIG)
% Convert transformed optimization coordinates to physical GH parameters.
%
% Transformation:
%   gamma = exp(psi)              ensures gamma > 0
%   alpha = hypot(beta, gamma)    ensures alpha > |beta| (smooth at beta=0)
%   delta = exp(log_delta)        ensures delta > 0
%
% hypot is used instead of sqrt(beta^2+gamma^2) for numerical stability
% when either beta or gamma is very large or very small.
if isNIG
    p.lambda = -0.5;
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

function results = assemble_results(model, p, theta, nll, hess, ci, ef, out, data, isNIG)
% Package all estimation outputs into a single results structure.
n    = numel(data);
k    = numel(theta);
fitX = sort(data);

results.Model   = model;
results.Params  = [p.mu; p.lambda; p.alpha; p.beta; p.delta];
results.LogLike = -nll;
results.Fit     = get_pdf_vals(data, p);    % density at original data points
results.FitX    = fitX;                     % sorted x for smooth density plots
results.FitPDF  = get_pdf_vals(fitX, p);    % density at sorted points
results.AIC     = 2*k + 2*nll;
results.BIC     = k*log(n) + 2*nll;
results.theta_hat = theta;

% ---- Validate that the reported solution is a real fit -------------------
% A returned parameter vector is only meaningful if the density can actually
% be evaluated there AND reproduces the reported likelihood. get_pdf_vals
% already marks failed evaluations as NaN, but nothing used to look, so a
% degenerate run could report converged=true alongside an all-NaN density and
% a nonsensical log-likelihood. Check both and say so loudly.
results.valid = true;
badFit = ~all(isfinite(results.Fit)) || any(results.Fit <= 0);
if badFit
    results.valid = false;
    warning('GeneralizedHyperbolic_MLE:DensityNotEvaluable', ...
        ['The fitted parameters do not yield an evaluable density (%d of %d ' ...
        'points are NaN or non-positive). results.LogLike is not trustworthy. ' ...
        'This usually means the optimizer ran away to extreme parameters.'], ...
        nnz(~isfinite(results.Fit) | results.Fit <= 0), n);
else
    llCheck = sum(log(results.Fit));
    tol = max(1e-6, 1e-8 * abs(nll));
    if ~isfinite(llCheck) || abs(llCheck + nll) > tol
        results.valid = false;
        warning('GeneralizedHyperbolic_MLE:LikelihoodInconsistent', ...
            ['Reported log-likelihood (%.6g) disagrees with the density ' ...
            'evaluated at the fitted parameters (%.6g). The optimum is not ' ...
            'numerically reliable.'], -nll, llCheck);
    end
end

if isNIG
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
% Evaluate the GH density at points x given physical parameters p.
xm = x - p.mu;
d2 = p.delta^2 + xm.^2;
z  = sqrt(d2);

a0 = p.delta * p.gamma;
az = p.alpha  * z;

K0 = besselk(p.lambda,       a0, 1);
Kz = besselk(p.lambda - 0.5, az, 1);

logpdf = p.lambda*log(p.gamma) - p.lambda*log(p.delta) - 0.5*log(2*pi) ...
    - (p.lambda - 0.5)*log(p.alpha) - (log(K0) - a0) ...
    + (p.lambda/2 - 0.25).*log(d2) + (log(Kz) - az) + p.beta*xm;

pdf = exp(logpdf);
pdf(~isfinite(pdf)) = NaN;   % flag any density evaluation failures
end

function plot_results(data, results)
% Plot histogram of data with superimposed fitted density.
figure('Color', 'w');
histogram(data, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'k');
hold on
plot(results.FitX, results.FitPDF, 'r-', 'LineWidth', 2);
xlabel('Data Value');
ylabel('Density');
title(['Fit: ' results.Model]);
grid on;  box on;
legend('Data', 'Fitted Model', 'Location', 'best');
end
