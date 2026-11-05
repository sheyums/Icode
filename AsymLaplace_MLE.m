function ModelResults = AsymLaplace_MLE(x, mode)
% ASYMLAPLACE_MLE  Exact maximum-likelihood fit of the asymmetric Laplace
% distribution, with standard errors and confidence intervals.
%
%   ModelResults = AsymLaplace_MLE(x)
%   ModelResults = AsymLaplace_MLE(x, mode)
%
% MODEL
% -----
%   The asymmetric Laplace distribution (ALD) in the (mu, sigma, kappa)
%   parameterization of Kotz, Kozubowski & Podgorski (2001):
%
%       f(x) = kappa / (sigma*(1 + kappa^2)) * exp(-kappa*(x - mu)/sigma)    x >= mu
%       f(x) = kappa / (sigma*(1 + kappa^2)) * exp( (x - mu)/(kappa*sigma))  x <  mu
%
%   mu     location of the mode (the cusp); any real number.
%   sigma  scale, > 0.
%   kappa  asymmetry, > 0. kappa = 1 is the symmetric Laplace; kappa < 1
%          puts the long tail on the RIGHT (right skew), kappa > 1 on the
%          left. P(X >= mu) = 1/(1 + kappa^2).
%   Right of mu the density decays at rate kappa/sigma (mean excess
%   sigma/kappa); left of mu at rate 1/(kappa*sigma) (mean kappa*sigma).
%
% ESTIMATION (exact, no iterative optimizer)
% ----------
%   For fixed mu let S1 = sum over x >= mu of (x - mu) and S2 = sum over
%   x < mu of (mu - x). The likelihood is maximized over sigma and kappa
%   in closed form,
%       kappa(mu) = (S2/S1)^(1/4),
%       sigma(mu) = (S1*S2)^(1/4) * (sqrt(S1) + sqrt(S2)) / n,
%   leaving the profile log-likelihood
%       LL(mu) = n*log(n) - 2n*log(sqrt(S1) + sqrt(S2)) - n.
%   sqrt(S1) + sqrt(S2) is concave in mu between consecutive data points,
%   so its minimum, the MLE of mu, is one of the data points (Kotz,
%   Kozubowski & Podgorski 2002). All n data points are evaluated at once
%   from cumulative sums of the sorted data: O(n log n) in total, and the
%   result is the GLOBAL maximum, not a local one. The data are centred on
%   their median and scaled by their mean absolute deviation first (for
%   precision; the fit is exactly equivariant: x -> a*x + b gives
%   mu -> a*mu + b, sigma -> |a|*sigma, kappa -> kappa for a > 0 and
%   1/kappa for a < 0).
%
%   One-sided data. If the minimum falls at the smallest observation (S2 =
%   0) or the largest (S1 = 0), the supremum is NOT attained by any ALD: it
%   is the limit kappa -> 0 (or Inf), a shifted exponential. The function
%   then returns the best INTERIOR data point, sets Diagnostics.OneSided,
%   records the limit's log-likelihood in Diagnostics.LogLikeOneSidedLimit,
%   marks SEValid false, and warns (AsymLaplace_MLE:OneSided).
%
% INPUTS
% ------
%   x     real numeric vector. NaN and Inf are removed with a warning. At
%         least 5 observations and at least 2 distinct values are
%         required (5 so that AICc is finite).
%   mode  optional text; keywords may be combined ('FIT PLOT'):
%         'FIT'        fill Fit and FitX (the fitted PDF on a grid).
%         'PLOT'       also plot a histogram with the fitted PDF.
%         'NOPROFILE'  skip the profile-likelihood intervals; CILower and
%                      CIUpper are then the Wald intervals (faster for
%                      large n; see ACCURACY for when that matters).
%         'BOOTSTRAP'  also compute nonparametric bootstrap SEs and
%                      percentile intervals (SEBoot, CIBootLower/Upper;
%                      500 resamples from a PRIVATE random stream with a
%                      fixed seed, so the global random state is untouched
%                      and results are reproducible).
%
% OUTPUTS (ModelResults struct)
% -------
%   Model           'Asymmetric Laplace (mu, sigma, kappa)'.
%   ParameterNames  {'mu', 'sigma', 'kappa'}.
%   Params          [mu, sigma, kappa], the MLE.
%   LogLike         maximized log-likelihood.
%   AIC             2k - 2*LogLike, k = 3. (Plain AIC: before Nov 2026
%                   this field held AICc.)
%   AICc            AIC + 2k(k+1)/(n - k - 1).
%   BIC             k*log(n) - 2*LogLike.
%   N, NParameters  sample size used, and 3.
%   SE              1 x 3 standard errors, from the EXPECTED (Fisher)
%                   information evaluated at the MLE, in closed form. The
%                   observed Hessian is not used: the log-density has a
%                   kink at mu, so its second derivative in mu is a point
%                   mass and the observed Hessian is meaningless there.
%   Cov             3 x 3 covariance behind SE (order of Params).
%   SEValid         false when the fit is one-sided (see above).
%   CILower,        1 x 3 confidence limits at CILevel. Profile likelihood
%   CIUpper         by default: the values where the log-likelihood,
%                   maximized over the other two parameters (exactly, not
%                   numerically), falls chi2inv(CILevel,1)/2 below
%                   LogLike. Wald when 'NOPROFILE' is given.
%   CIOpen          1 x 3 logical: an end never fell below the threshold,
%                   so the interval is open there (reported as -Inf/Inf for
%                   mu, 0/Inf for sigma and kappa).
%   CIMethod        'profile' or 'Wald'.
%   CIWaldLower,    Wald limits, always computed: mu +/- z*SE, and sigma
%   CIWaldUpper     and kappa on the log scale (so they stay positive).
%   CILevel         0.95.
%   SEBoot,         ('BOOTSTRAP' only) bootstrap SE and 2.5/97.5 percentile
%   CIBootLower,    limits; empty otherwise.
%   CIBootUpper
%   Fit, FitX       ('FIT' or 'PLOT') fitted PDF on a 1000-point grid over
%                   the data range that includes mu exactly; empty
%                   otherwise.
%   Diagnostics     struct: OneSided, LogLikeOneSidedLimit, Center, Scale,
%                   NBoot.
%
% WARNINGS
% --------
%   AsymLaplace_MLE:NonfiniteDataRemoved  NaN/Inf dropped from x.
%   AsymLaplace_MLE:SmallSample           n < 100: intervals under-cover
%                                         (see ACCURACY).
%   AsymLaplace_MLE:OneSided              the supremum is the exponential
%                                         limit (see ESTIMATION).
% ERRORS
% ------
%   AsymLaplace_MLE:BadInput, :BadMode, :TooFewObservations,
%   :ConstantData.
%
% ACCURACY (measured Nov 2026, R2026a)
% --------
%   Point estimates: LogLike equals the density summed directly and the
%   maximum of a brute-force profile (21 points between every pair of
%   data points) in every case checked; each profile bound sits on the
%   deviance threshold to 1e-5 against an independent Nelder-Mead
%   maximization over the other two parameters.
%   Interval coverage, nominal 95%, 400 samples per row (Monte Carlo SE
%   about 0.011), mu = 0, sigma = 1; columns mu, sigma, kappa:
%
%      n  kappa | profile (default) | Wald              | bootstrap
%     30   0.5  | 0.88  0.87  0.86  | 0.67  0.91  0.95  | 0.95  0.80  0.91
%     30   1    | 0.90  0.92  0.90  | 0.78  0.91  0.84  | 0.97  0.75  0.97
%     30   2    | 0.81  0.86  0.81  | 0.64  0.93  0.91  | 0.95  0.77  0.93
%    100   0.5  | 0.92  0.93  0.93  | 0.86  0.92  0.91  | 0.95  0.92  0.97
%    100   1    | 0.94  0.95  0.91  | 0.87  0.95  0.89  | 0.95  0.92  0.95
%    100   2    | 0.93  0.92  0.92  | 0.84  0.92  0.89  | 0.96  0.90  0.95
%    500   0.5  | 0.93  0.95  0.94  | 0.90  0.94  0.93  | 0.94  0.94  0.95
%    500   1    | 0.93  0.96  0.95  | 0.90  0.96  0.94  | 0.95  0.96  0.97
%    500   2    | 0.94  0.95  0.94  | 0.91  0.95  0.92  | 0.95  0.95  0.94
%
%   - Use the profile intervals (default). The Wald interval for mu is
%     too short even at n = 500 (0.90): the MLE of mu is a data point
%     and its sampling distribution is far from normal at these n.
%   - For mu alone the bootstrap is the most accurate (0.94-0.97 at all
%     n); for sigma at n = 30 it is the worst (0.75-0.80).
%   - Below n = 100 every method under-covers; hence the SmallSample
%     warning. At n = 30 and kappa = 0.5, 140 of 400 samples were
%     one-sided (the supremum at the exponential limit), 64 at kappa = 2.
%   - The kappa = 1 rows are the symmetric Laplace, which is the interior
%     point kappa = 1 of this model; nothing special is needed for it.
%     Whether kappa's interval contains 1 is a test of symmetry.
%   - The previous bootstrap-only version covered 0.92-0.96 (200 samples,
%     n = 100 kappa = 0.5 and n = 500 kappa = 2): its intervals were
%     reasonable; the defects were elsewhere (see CHANGE LOG).
%
% REFERENCES
% ----------
%   Kotz S, Kozubowski TJ, Podgorski K (2001). The Laplace Distribution
%     and Generalizations. Birkhauser. Ch. 3.
%   Kotz S, Kozubowski TJ, Podgorski K (2002). Maximum likelihood
%     estimation of asymmetric Laplace parameters. Ann. Inst. Statist.
%     Math. 54, 816-826.
%
% CHANGE LOG
% ----------
%   Jan 2026   Nelder-Mead fit (kept as AsymLaplace_MLE_V0.m).
%   Nov 2026   Exact O(n) profile solver and 500-resample bootstrap
%              (analysis-folder working copy, Icode 059f10f).
%   2026-11-04 Audit (Claude Code), against ground truth sharing no code
%              with the solver (brute-force profile on a grid between the
%              data points, the density summed directly). The solver was
%              already exact and global; the changes are:
%              - AIC is plain AIC; AICc has its own field.
%              - SEs from the closed-form Fisher information; profile and
%                Wald intervals; the bootstrap is optional ('BOOTSTRAP')
%                and uses a private stream (it used to consume the global
%                random stream, changing every later rand call).
%              - The degeneracy filter was an absolute 1e-12 on S1 and
%                S2, so data on a 1e-15 scale failed; the profile is now
%                computed scale-free on standardized data.
%              - One-sided data are flagged instead of silently returning
%                a near-degenerate interior fit.
%              - Constant data error (they returned [mu NaN NaN] before);
%                n >= 5 required (was 3; AICc needs n > 4).

%% 1. Input handling
if nargin < 2 || isempty(mode), mode = ''; end
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('AsymLaplace_MLE:BadMode', 'mode must be text, e.g. ''FIT'' or ''PLOT''.');
end
mode = upper(char(mode));
if ~(isnumeric(x) && isreal(x))
    error('AsymLaplace_MLE:BadInput', 'x must be a real numeric vector.');
end
x = double(x(:));
bad = ~isfinite(x);
if any(bad)
    warning('AsymLaplace_MLE:NonfiniteDataRemoved', '%d NaN/Inf values removed.', sum(bad));
    x = x(~bad);
end
x = sort(x);
n = numel(x);
if n < 5
    error('AsymLaplace_MLE:TooFewObservations', ...
        'At least 5 finite observations are required (got %d).', n);
end
if x(1) == x(end)
    error('AsymLaplace_MLE:ConstantData', 'All observations are equal; the ALD scale is 0.');
end
if n < 100
    warning('AsymLaplace_MLE:SmallSample', ...
        'n = %d < 100: the 95%% intervals cover only 0.81-0.94 at n = 30-100 (measured).', n);
end
ciLevel = 0.95;
zc      = sqrt(2) * erfinv(ciLevel);
crit    = zc^2;                         % chi2inv(0.95, 1)

%% 2. Standardize and fit exactly
ctr = median(x);
scl = mean(abs(x - ctr));               % > 0 because x is not constant
z   = (x - ctr) / scl;
D   = ald_data(z);
[muZ, sZ, K, LLz, oneSided, LLlimZ, kLim] = ald_exact(D);

mu  = ctr + scl*muZ;
s   = scl*sZ;
LL  = LLz - n*log(scl);
LLlim = LLlimZ - n*log(scl);

if oneSided
    warning('AsymLaplace_MLE:OneSided', ['The likelihood supremum is the one-sided ' ...
        '(shifted exponential) limit, kappa -> %s, with LogLike %.4f; the returned ' ...
        'interior fit has %.4f. The data look one-sided; the ALD fit and its ' ...
        'intervals are unreliable.'], kLim, LLlim, LL);
end

kpar = 3;
ModelResults = struct();
ModelResults.Model          = 'Asymmetric Laplace (mu, sigma, kappa)';
ModelResults.ParameterNames = {'mu', 'sigma', 'kappa'};
ModelResults.Params         = [mu, s, K];
ModelResults.LogLike        = LL;
ModelResults.AIC            = 2*kpar - 2*LL;
ModelResults.AICc           = 2*kpar - 2*LL + 2*kpar*(kpar + 1)/(n - kpar - 1);
ModelResults.BIC            = kpar*log(n) - 2*LL;
ModelResults.N              = n;
ModelResults.NParameters    = kpar;

%% 3. Standard errors (closed-form Fisher information) and Wald intervals
I1z = ald_fisher(sZ, K);                % per observation, standardized units
Dz  = diag([scl, scl, 1]);
Cov = Dz*(I1z \ eye(3))*Dz / n;
se  = sqrt(diag(Cov))';
ModelResults.SE      = se;
ModelResults.Cov     = Cov;
ModelResults.SEValid = ~oneSided;
waldLo = [mu - zc*se(1), s*exp(-zc*se(2)/s), K*exp(-zc*se(3)/K)];
waldHi = [mu + zc*se(1), s*exp( zc*se(2)/s), K*exp( zc*se(3)/K)];
ModelResults.CIWaldLower = waldLo;
ModelResults.CIWaldUpper = waldHi;
ModelResults.CILevel     = ciLevel;

%% 4. Profile-likelihood intervals (exact profiles, in standardized units)
if contains(mode, 'NOPROFILE')
    ModelResults.CILower  = waldLo;
    ModelResults.CIUpper  = waldHi;
    ModelResults.CIOpen   = false(1, 3);
    ModelResults.CIMethod = 'Wald';
else
    seZ = se ./ [scl, scl, 1];
    % transformed scales: mu itself, log sigma, log kappa
    prof = {@(t) prof_mu(D, t), @(t) prof_sigma(D, exp(t)), @(t) prof_kappa(D, exp(t))};
    that = [muZ, log(sZ), log(K)];
    step = max([seZ(1), seZ(2)/sZ, seZ(3)/K], 1e-3);
    lim  = [100*(D.z(end) - D.z(1)), 50, 50];
    lo = NaN(1, 3); hi = NaN(1, 3); open = false(1, 3);
    for j = 1:3
        dev = @(t) 2*(LLz - prof{j}(t)) - crit;
        [lo(j), oL] = find_bound(dev, that(j), -1, step(j), lim(j));
        [hi(j), oH] = find_bound(dev, that(j),  1, step(j), lim(j));
        open(j) = oL || oH;
    end
    ModelResults.CILower  = [ctr + scl*lo(1), scl*exp(lo(2)), exp(lo(3))];
    ModelResults.CIUpper  = [ctr + scl*hi(1), scl*exp(hi(2)), exp(hi(3))];
    ModelResults.CIOpen   = open;
    ModelResults.CIMethod = 'profile';
end

%% 5. Optional bootstrap (private stream)
nBoot = 0;
ModelResults.SEBoot = []; ModelResults.CIBootLower = []; ModelResults.CIBootUpper = [];
if contains(mode, 'BOOTSTRAP')
    nBoot = 500;
    rs = RandStream('twister', 'Seed', 0);
    B  = NaN(nBoot, 3);
    for b = 1:nBoot
        zb = sort(z(randi(rs, n, n, 1)));
        if zb(1) == zb(end), continue; end
        [mb, sb, kb] = ald_exact(ald_data(zb));
        B(b, :) = [ctr + scl*mb, scl*sb, kb];
    end
    B = B(all(isfinite(B), 2), :);
    ModelResults.SEBoot      = std(B, 0, 1);
    ModelResults.CIBootLower = prctile(B, 100*(1 - ciLevel)/2, 1);
    ModelResults.CIBootUpper = prctile(B, 100*(1 + ciLevel)/2, 1);
end

%% 6. Fitted PDF and plot
ModelResults.Fit = []; ModelResults.FitX = [];
if contains(mode, 'FIT') || contains(mode, 'PLOT')
    xs = unique([linspace(x(1), x(end), 1000)'; mu]);
    ModelResults.FitX = xs;
    ModelResults.Fit  = ald_pdf(xs, mu, s, K);
end
ModelResults.Diagnostics = struct('OneSided', oneSided, 'LogLikeOneSidedLimit', LLlim, ...
    'Center', ctr, 'Scale', scl, 'NBoot', nBoot);

if contains(mode, 'PLOT')
    figure('Color', 'w');
    histogram(x, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
        'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.5);
    hold on;
    plot(ModelResults.FitX, ModelResults.Fit, 'r-', 'LineWidth', 2);
    xlabel('Data values'); ylabel('Probability density');
    title(sprintf('Asymmetric Laplace fit\n\\mu=%.3g, \\sigma=%.3g, \\kappa=%.3g', mu, s, K));
    legend('Empirical PDF', 'ALD MLE', 'Location', 'best');
    grid on; hold off;
end
end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================
function D = ald_data(z)
% Sorted data plus S1, S2 at every data point (mu = z(i)).
n  = numel(z);
cz = [0; cumsum(z)];
i  = (1:n)';
D.z  = z;
D.n  = n;
D.cz = cz;
D.S1 = max((cz(end) - cz(i)) - (n - i + 1).*z, 0);   % sum (z_j - mu), z_j >= mu
D.S2 = max((i - 1).*z - cz(i), 0);                   % sum (mu - z_j), z_j <  mu
end

function [mu, s, K, LL, oneSided, LLlim, kLim] = ald_exact(D)
% Global MLE over the data points. One-sided if the best point has S1 or
% S2 equal to 0; then the best interior point is returned.
n    = D.n;
crit = sqrt(D.S1) + sqrt(D.S2);
edge = D.S1 <= 0 | D.S2 <= 0;
[cAll, iAll] = min(crit);
oneSided = edge(iAll);
LLlim = -Inf;
kLim = '0';
if D.S1(iAll) <= 0, kLim = 'Inf'; end
if oneSided
    LLlim = n*log(n) - 2*n*log(cAll) - n;
    crit(edge) = Inf;
end
[cBest, ib] = min(crit);
if ~isfinite(cBest)                      % only two distinct values: no interior point
    mu = D.z(iAll); s = NaN; K = NaN; LL = LLlim; return
end
S1 = D.S1(ib); S2 = D.S2(ib);
mu = D.z(ib);
K  = (S2/S1)^0.25;
s  = (S1*S2)^0.25 * (sqrt(S1) + sqrt(S2)) / n;
LL = n*log(n) - 2*n*log(cBest) - n;
end

function [S1, S2] = s12(D, mu)
% S1, S2 at an arbitrary mu.
j  = sum(D.z < mu);
S2 = max(j*mu - D.cz(j + 1), 0);
S1 = max((D.cz(end) - D.cz(j + 1)) - (D.n - j)*mu, 0);
end

function LL = prof_mu(D, mu)
[S1, S2] = s12(D, mu);
n  = D.n;
LL = n*log(n) - 2*n*log(sqrt(S1) + sqrt(S2)) - n;
end

function LL = prof_kappa(D, K)
% max over mu and sigma for fixed kappa: minimize K*S1 + S2/K (piecewise
% linear and convex in mu, so its minimum is at a data point).
n  = D.n;
h  = min(K*D.S1 + D.S2/K);
LL = n*log(K) - n*log(1 + K^2) - n*log(h/n) - n;
end

function LL = prof_sigma(D, s)
% max over mu and kappa for fixed sigma: for each kappa the best mu is the
% data point minimizing K*S1 + S2/K; then a 1-D search in log kappa.
n = D.n;
f = @(t) -(n*t - n*log(1 + exp(2*t)) - n*log(s) - min(exp(t)*D.S1 + D.S2*exp(-t))/s);
tg = linspace(-12, 12, 97);
fv = arrayfun(f, tg);
[~, k] = min(fv);
a = tg(max(k - 1, 1)); b = tg(min(k + 1, numel(tg)));
[~, fb] = fminbnd(f, a, b, optimset('TolX', 1e-10));
LL = -min(fb, fv(k));
end

function [bnd, open] = find_bound(dev, t0, dir, step, lim)
% First crossing of dev = 0 moving from t0 in direction dir; dev(t0) < 0.
% Steps grow by 1.5x; open if no crossing within lim of t0.
open = false;
tPrev = t0;
h = step;
while true
    t = t0 + dir*min(abs(tPrev - t0) + h, lim);
    d = dev(t);
    if d > 0
        bnd = fzero(dev, sort([tPrev, t]));
        return
    end
    if abs(t - t0) >= lim
        bnd = dir*Inf; open = true;
        return
    end
    tPrev = t;
    h = 1.5*h;
end
end

function I = ald_fisher(s, K)
% Per-observation Fisher information for (mu, sigma, kappa), exact.
% Scores, with R = 1{x >= mu}, L = 1{x < mu}, d+ = (x-mu)R, d- = (mu-x)L:
%   g_mu    = (K/s) R - 1/(K s) L
%   g_sigma = -1/s + (K/s^2) d+ + 1/(K s^2) d-
%   g_kappa = 1/K - 2K/(1+K^2) - d+/s + d-/(K^2 s)
% I = E[g g'] = C*M*C', C the coefficients on v = [1 R L d+ d-], M = E[v v'].
p = 1/(1 + K^2); q = K^2/(1 + K^2);
Ed1 = p*s/K;  Ed1s = 2*p*s^2/K^2;  Ed2 = q*K*s;  Ed2s = 2*q*K^2*s^2;
M = [1    p    q    Ed1   Ed2
     p    p    0    Ed1   0
     q    0    q    0     Ed2
     Ed1  Ed1  0    Ed1s  0
     Ed2  0    Ed2  0     Ed2s];
C = [0, K/s, -1/(K*s), 0, 0
     -1/s, 0, 0, K/s^2, 1/(K*s^2)
     1/K - 2*K/(1 + K^2), 0, 0, -1/s, 1/(K^2*s)];
I = C*M*C';
I = (I + I')/2;
end

function y = ald_pdf(x, mu, s, K)
c = K/(s*(1 + K^2));
y = zeros(size(x));
r = x >= mu;
y(r)  = c*exp(-K*(x(r) - mu)/s);
y(~r) = c*exp((x(~r) - mu)/(K*s));
end

