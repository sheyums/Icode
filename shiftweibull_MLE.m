function ModelResults = shiftweibull_MLE(x, mode)
% SHIFTWEIBULL_MLE  Maximum-likelihood fit of the three-parameter (shifted)
% Weibull distribution, with standard errors and confidence intervals.
%
%   ModelResults = shiftweibull_MLE(x)
%   ModelResults = shiftweibull_MLE(x, mode)
%
% MODEL
% -----
%       f(x) = (b/a) * ((x - c)/a)^(b - 1) * exp(-((x - c)/a)^b),   x > c
%
%   a  scale, > 0.
%   b  shape, constrained to b >= 1.001 (see THE SHAPE FLOOR).
%   c  shift (threshold, location), constrained to c < min(x).
%   b = 1 is a shifted exponential, b ~ 3.6 is close to symmetric, b > 3.6
%   is left-skewed; b -> Inf with a, c -> Inf together tends to the
%   Gumbel (minimum) distribution.
%
% THE SHAPE FLOOR
% ---------------
%   For b < 1 the density is infinite at c, so the likelihood is
%   UNBOUNDED as c -> min(x): no maximum-likelihood estimate exists. The
%   function therefore restricts b >= 1.001 (kept from the Nov 2026
%   working version). This is a modelling restriction, not a result: when
%   the fit sits on the floor (data whose shape is really below ~1, e.g.
%   a true b = 0.8), ShapeAtFloor is set, a warning is issued, and the
%   shift lands just below min(x). Do not read the fitted shape as an
%   estimate there.
%
% ESTIMATION
% ----------
%   The data are standardized, z = (x - min(x))/std(x), so the fit is
%   exactly equivariant under x -> s*x + t, s > 0 (a -> s*a, c -> s*c + t,
%   b unchanged). The shift is searched as the gap g = min(z) - c_z > 0,
%   on t = log(g):
%     1. For each gap the two-parameter Weibull fit of y = z + g is exact:
%        the shape solves  sum(y^b log y)/sum(y^b) - 1/b - mean(log y) = 0,
%        which is increasing in b, so fzero finds the unique root (or the
%        floor if the root lies below it); the scale is then
%        (mean(y^b))^(1/b). This gives the profile log-likelihood in t.
%     2. That profile is evaluated on 300 points of t in [-18, 8] (gaps
%        from 1.5e-8 to 3000 standard deviations) and the best grid point
%        refined with fminbnd. The result is the GLOBAL maximum of the
%        profile to within the grid's ability to separate peaks (spacing
%        0.087 in log gap), not a local optimum from one start.
%   If the best gap is the grid's largest, the likelihood is still rising
%   toward the Gumbel limit (b -> Inf): ShiftAtLowerLimit is set and a
%   warning issued.
%
% STANDARD ERRORS AND INTERVALS
% -----------------------------
%   SEs come from the observed information: the Hessian of the exact
%   negative log-likelihood in (log a, log b, log gap) by central
%   differences with Richardson extrapolation, mapped to (a, b, c) by the
%   delta method.
%   Regularity: the threshold MLE is asymptotically normal only when
%   b > 2 (Smith 1985). For 1 < b <= 2 the information for c is not
%   finite in the limit, the Wald interval for c is unreliable, and
%   NonRegular is set. See ACCURACY for the measured consequences.
%   Profile-likelihood intervals (default) take the values where the
%   log-likelihood, maximized over the other two parameters, falls
%   chi2inv(CILevel,1)/2 below LogLike. All three profiles are computed
%   with the same exact inner solves as the fit.
%
% INPUTS
% ------
%   x     real numeric vector; NaN/Inf removed with a warning. At least 5
%         observations and 3 distinct values are required.
%   mode  optional text; keywords may be combined ('FIT PLOT'):
%         'FIT'        fill Fit and FitX (fitted PDF at the sorted data).
%         'PLOT'       also plot a histogram with the fitted PDF.
%         'NOPROFILE'  skip the profiles; CILower/CIUpper are Wald.
%         'BOOTSTRAP'  also compute nonparametric bootstrap SEs and
%                      percentile intervals (SEBoot, CIBootLower/Upper;
%                      300 resamples from a PRIVATE random stream with a
%                      fixed seed, so the global random state is untouched).
%
% OUTPUTS (ModelResults struct)
% -------
%   Model           'Shifted-Weibull (Scale, Shape, Shift)'.
%   ParameterNames  {'scale', 'shape', 'shift'}.
%   Params          [a, b, c], the MLE (same order as before).
%   LogLike         maximized log-likelihood.
%   AIC             2k - 2*LogLike, k = 3. (Plain AIC: before Nov 2026
%                   this field held AICc.)
%   AICc            AIC + 2k(k+1)/(n - k - 1).
%   BIC             k*log(n) - 2*LogLike.
%   N, NParameters  sample size used, and 3.
%   SE              1 x 3 standard errors (observed information), NaN when
%                   that matrix is not positive definite or the fit is on
%                   the shape floor.
%   Cov             3 x 3 covariance behind SE.
%   SEValid         false when SE is NaN, or when NonRegular (b <= 2) or a
%                   limit flag is set.
%   CILower,        1 x 3 confidence limits at CILevel: profile likelihood
%   CIUpper         by default, Wald with 'NOPROFILE'.
%   CIOpen          1 x 3 logical: an end never crossed the threshold. An
%                   open shape end is reported as the floor (lower) or Inf;
%                   an open shift end as -Inf; an open scale end as 0/Inf.
%   CIMethod        'profile' or 'Wald'.
%   CIWaldLower,    Wald limits, always computed: a and b on the log scale,
%   CIWaldUpper     c on log(min(x) - c) so the interval stays below min(x).
%   CILevel         0.95.
%   SEBoot,         ('BOOTSTRAP' only) bootstrap SE and percentile limits.
%   CIBootLower,
%   CIBootUpper
%   Fit, FitX       ('FIT' or 'PLOT') PDF at the sorted data; else empty.
%   Diagnostics     struct: ShapeAtFloor, ShiftAtLowerLimit, NonRegular,
%                   ShapeFloor, Gap (min(x) - c), NBoot.
%
% WARNINGS
% --------
%   shiftweibull_MLE:NonfiniteDataRemoved  NaN/Inf dropped from x.
%   shiftweibull_MLE:SmallSample           n < 50 (see ACCURACY).
%   shiftweibull_MLE:ShapeAtFloor          b on the 1.001 floor.
%   shiftweibull_MLE:ShiftAtLowerLimit     the Gumbel limit (b -> Inf).
% ERRORS
% ------
%   shiftweibull_MLE:BadInput, :BadMode, :TooFewObservations,
%   :TooFewDistinct.
%
% ACCURACY (measured Nov 2026, R2026a)
% --------
%   Point estimates: LogLike equals the density summed directly and is
%   never below a 40-start Nelder-Mead in its own coordinates; SEs equal
%   an independent physical-coordinate Hessian to 1e-3; each profile
%   bound sits on the deviance threshold against Nelder-Mead over the
%   other two parameters.
%   Interval coverage, nominal 95%, 300 samples per row (Monte Carlo SE
%   about 0.013; 100 for the bootstrap row), scale 2, shift 0; columns
%   scale, shape, shift:
%
%       n  shape | profile (default) | Wald              | bootstrap
%      30   3    | 0.94  0.94  0.93  | 0.86  0.89  0.98  |
%      30   1.5  | 0.93  0.90  0.90  | 0.79  0.85  0.85  |
%      50   3    | 0.95  0.94  0.94  | 0.86  0.90  0.99  |
%      50   1.5  | 0.94  0.92  0.91  | 0.93  0.93  0.95  |
%     100   3    | 0.96  0.95  0.96  | 0.87  0.92  0.97  |
%     100   1.5  | 0.92  0.96  0.89  | 0.91  0.94  1.00  | 0.86  0.96  0.62
%     300   5    | 0.97  0.96  0.96  | 0.94  0.94  0.98  |
%     500   1.5  | 0.94  0.96  0.96  | 0.94  0.95  1.00  |
%    1000   3    | 0.95  0.95  0.96  | 0.93  0.94  0.98  |
%
%   - Use the profile intervals (default). Wald is short for the scale at
%     n <= 100 (0.79-0.87) and over-long for the shift (0.97-1.00).
%   - Do NOT use the bootstrap for the shift: 0.62 at shape 1.5, n = 100
%     (the bootstrap is known to fail for a threshold parameter). The
%     previous bootstrap-only version measured 0.90/0.92/0.74 there and
%     0.90/0.88/0.87 at shape 3, n = 300 (100 samples each).
%   - Small samples: open intervals are common (an end never crosses the
%     threshold: 123 of 300 at n = 30, shape 3; 34 at n = 50), and at
%     shape 1.5 the fit sat on the shape floor in 43 of 300 samples at
%     n = 30 and 9 at n = 50. Hence the SmallSample warning below n = 50.
%
% REFERENCES
% ----------
%   Smith RL (1985). Maximum likelihood estimation in a class of
%     nonregular cases. Biometrika 72, 67-90.
%   Cohen AC (1975). Multi-censored sampling in the three parameter
%     Weibull distribution. Technometrics 17, 347-351 (profile equations).
%
% CHANGE LOG
% ----------
%   Jan 2026   Nested wblfit search (kept as shiftweibull_MLE_V0.m).
%   Nov 2026   Simultaneous Nelder-Mead in (log a, log(b - 1.001),
%              log gap) with a 300-resample bootstrap (analysis-folder
%              working copy, Icode 059f10f).
%   2026-11-04 Audit (Claude Code) against an independent profile-over-gap
%              fit sharing no code with the fitter. The Nelder-Mead
%              estimates were already right in all 8 cases probed; the
%              changes are:
%              - The search is the exact profile over the gap (global and
%                deterministic, no Statistics Toolbox call for the fit).
%              - AIC is plain AIC; AICc has its own field.
%              - The bootstrap called warning('on','all') afterwards,
%                switching ON every warning the caller had turned off; it
%                also consumed the global random stream. It is now
%                optional, on a private stream, and touches neither.
%              - SEs from the observed information; profile and Wald
%                intervals; flags for the shape floor, the Gumbel limit
%                and the non-regular range b <= 2.

%% 1. Input handling
if nargin < 2 || isempty(mode), mode = ''; end
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('shiftweibull_MLE:BadMode', 'mode must be text, e.g. ''FIT'' or ''PLOT''.');
end
mode = upper(char(mode));
if ~(isnumeric(x) && isreal(x))
    error('shiftweibull_MLE:BadInput', 'x must be a real numeric vector.');
end
x = double(x(:));
bad = ~isfinite(x);
if any(bad)
    warning('shiftweibull_MLE:NonfiniteDataRemoved', '%d NaN/Inf values removed.', sum(bad));
    x = x(~bad);
end
n = numel(x);
if n < 5
    error('shiftweibull_MLE:TooFewObservations', ...
        'At least 5 finite observations are required (got %d).', n);
end
if numel(unique(x)) < 3
    error('shiftweibull_MLE:TooFewDistinct', 'At least 3 distinct values are required.');
end
if n < 50
    warning('shiftweibull_MLE:SmallSample', ...
        'n = %d < 50: intervals are often open and shape-1.5 coverage is 0.90 at n = 30 (measured).', n);
end
ciLevel = 0.95;
zc      = sqrt(2) * erfinv(ciLevel);
crit    = zc^2;
FLOOR   = shape_floor();
TGRID   = [-18, 8];

%% 2. Standardize and fit
m0 = min(x);
sc = std(x);
z  = (x - m0) / sc;
[tHat, LLz, aZ, b, atLimit] = fit_core(z, TGRID);
g  = exp(tHat);
a  = sc*aZ;
c  = m0 - sc*g;
LL = LLz - n*log(sc);

atFloor  = b <= FLOOR*(1 + 1e-9);
nonReg   = b <= 2;
if atFloor
    warning('shiftweibull_MLE:ShapeAtFloor', ['Shape is on its floor (%.3f): the data ' ...
        'suggest a shape below ~1, where the likelihood is unbounded. The shift sits ' ...
        'just below min(x) and the parameters are not estimates.'], FLOOR);
end
if atLimit
    warning('shiftweibull_MLE:ShiftAtLowerLimit', ['The likelihood still rises as the ' ...
        'shift -> -Inf (Gumbel limit, shape %.3g). The data look left-skewed; the ' ...
        'parameters are not reliable.'], b);
end

kpar = 3;
ModelResults = struct();
ModelResults.Model          = 'Shifted-Weibull (Scale, Shape, Shift)';
ModelResults.ParameterNames = {'scale', 'shape', 'shift'};
ModelResults.Params         = [a, b, c];
ModelResults.LogLike        = LL;
ModelResults.AIC            = 2*kpar - 2*LL;
ModelResults.AICc           = 2*kpar - 2*LL + 2*kpar*(kpar + 1)/(n - kpar - 1);
ModelResults.BIC            = kpar*log(n) - 2*LL;
ModelResults.N              = n;
ModelResults.NParameters    = kpar;

%% 3. Observed information, SEs and Wald intervals
% theta = [log a, log b, log g] in standardized units.
th  = [log(aZ); log(b); tHat];
nll = @(v) -sw_ll(z, exp(v(1)), exp(v(2)), exp(v(3)));
seOK = false;
Cov  = NaN(3); se = NaN(1, 3); seT = NaN(1, 3);
if ~atFloor
    H = fd_hessian(nll, th);
    [R, p] = chol(H);
    if p == 0 && rcond(H) > 1e-14
        Ri  = R \ eye(3);
        CT  = Ri*Ri';                       % covariance of theta
        seT = sqrt(diag(CT))';
        J   = diag([a, b, -sc*g]);          % d(a, b, c)/d(theta), physical units
        Cov = J*CT*J';
        se  = sqrt(diag(Cov))';
        seOK = all(isfinite(se));
    end
end
ModelResults.SE      = se;
ModelResults.Cov     = Cov;
ModelResults.SEValid = seOK && ~nonReg && ~atLimit;
waldLo = [a*exp(-zc*seT(1)), b*exp(-zc*seT(2)), m0 - sc*exp(tHat + zc*seT(3))];
waldHi = [a*exp( zc*seT(1)), b*exp( zc*seT(2)), m0 - sc*exp(tHat - zc*seT(3))];
ModelResults.CIWaldLower = waldLo;
ModelResults.CIWaldUpper = waldHi;
ModelResults.CILevel     = ciLevel;

%% 4. Profile-likelihood intervals
if contains(mode, 'NOPROFILE')
    ModelResults.CILower  = waldLo;
    ModelResults.CIUpper  = waldHi;
    ModelResults.CIOpen   = false(1, 3);
    ModelResults.CIMethod = 'Wald';
else
    step = seT;
    step(~isfinite(step)) = 0.1;
    step = min(max(step, 1e-3), 1);
    dev  = @(f) @(u) 2*(LLz - f(u)) - crit;
    % scale: u = log a_z
    pa = dev(@(u) prof_scale(z, exp(u), TGRID));
    [aL, oaL] = find_bound(pa, th(1), -1, step(1), 30);
    [aH, oaH] = find_bound(pa, th(1),  1, step(1), 30);
    % shape: u = log b, lower limit at the floor
    pb = dev(@(u) prof_shape(z, exp(u), TGRID));
    [bL, obL] = find_bound(pb, th(2), -1, step(2), th(2) - log(FLOOR));
    [bH, obH] = find_bound(pb, th(2),  1, step(2), 10);
    % shift: u = log gap (upper u = lower shift)
    pc = dev(@(u) fit_given_t(z, u));
    [tL, otL] = find_bound(pc, th(3), -1, step(3), th(3) - TGRID(1));
    [tH, otH] = find_bound(pc, th(3),  1, step(3), TGRID(2) - th(3));
    if obL, bLo = FLOOR; else, bLo = exp(bL); end
    if otH, cLo = -Inf; else, cLo = m0 - sc*exp(tH); end
    aLo = sc*exp(aL); if oaL, aLo = 0; end
    aHi = sc*exp(aH); if oaH, aHi = Inf; end
    bHi = exp(bH);    if obH, bHi = Inf; end
    ModelResults.CILower  = [aLo, bLo, cLo];
    ModelResults.CIUpper  = [aHi, bHi, m0 - sc*exp(tL)];
    ModelResults.CIOpen   = [oaL || oaH, obL || obH, otL || otH];
    ModelResults.CIMethod = 'profile';
end

%% 5. Optional bootstrap (private stream)
nBoot = 0;
ModelResults.SEBoot = []; ModelResults.CIBootLower = []; ModelResults.CIBootUpper = [];
if contains(mode, 'BOOTSTRAP')
    nBoot = 300;
    rs = RandStream('twister', 'Seed', 0);
    B  = NaN(nBoot, 3);
    for k = 1:nBoot
        xb = x(randi(rs, n, n, 1));
        if numel(unique(xb)) < 3, continue; end
        mb = min(xb); sb = std(xb);
        [tb, ~, ab, bb] = fit_core((xb - mb)/sb, TGRID);
        B(k, :) = [sb*ab, bb, mb - sb*exp(tb)];
    end
    B = B(all(isfinite(B), 2), :);
    ModelResults.SEBoot      = std(B, 0, 1);
    ModelResults.CIBootLower = prctile(B, 100*(1 - ciLevel)/2, 1);
    ModelResults.CIBootUpper = prctile(B, 100*(1 + ciLevel)/2, 1);
end

%% 6. Fitted PDF and plot
ModelResults.Fit = []; ModelResults.FitX = [];
if contains(mode, 'FIT') || contains(mode, 'PLOT')
    xs = sort(x);
    u  = (xs - c)/a;
    ModelResults.FitX = xs;
    ModelResults.Fit  = (b/a) .* u.^(b - 1) .* exp(-u.^b);
end
ModelResults.Diagnostics = struct('ShapeAtFloor', atFloor, 'ShiftAtLowerLimit', atLimit, ...
    'NonRegular', nonReg, 'ShapeFloor', FLOOR, 'Gap', m0 - c, 'NBoot', nBoot);

if contains(mode, 'PLOT')
    figure('Color', 'w');
    histogram(x, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
        'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.5);
    hold on;
    plot(ModelResults.FitX, ModelResults.Fit, 'r-', 'LineWidth', 2);
    xlabel('Data values'); ylabel('Probability density');
    title(sprintf('Shifted Weibull fit\nscale=%.3g, shape=%.3g, shift=%.3g', a, b, c));
    legend('Empirical PDF', 'MLE fit', 'Location', 'best');
    grid on; hold off;
end
end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================
function f = shape_floor()
f = 1.001;
end

function [tHat, LL, a, b, atEdge] = fit_core(z, TG)
% Global maximum of the profile over t = log gap: grid, then fminbnd.
tg = linspace(TG(1), TG(2), 300);
L  = arrayfun(@(t) fit_given_t(z, t), tg);
[Lb, k] = max(L);
atEdge = k == numel(tg);
lo = tg(max(k - 1, 1)); hi = tg(min(k + 1, numel(tg)));
[tR, nR] = fminbnd(@(t) -fit_given_t(z, t), lo, hi, optimset('TolX', 1e-10));
if -nR >= Lb
    tHat = tR;
else
    tHat = tg(k);
end
[LL, a, b] = fit_given_t(z, tHat);
end

function [LL, a, b] = fit_given_t(z, t)
% Exact two-parameter Weibull fit of y = z + exp(t), shape >= floor.
y  = z + exp(t);
ly = log(y);
n  = numel(y);
lm = max(ly);
mly = mean(ly);
G = @(b) wsum_ratio(ly, lm, b) - 1/b - mly;    % increasing in b
FLOOR = shape_floor();
if G(FLOOR) >= 0
    b = FLOOR;
else
    hi = 2;
    while G(hi) < 0 && hi < 1e7, hi = 2*hi; end
    if G(hi) < 0
        b = hi;
    else
        b = fzero(G, [FLOOR, hi]);
    end
end
% a^b = mean(y^b), computed in logs
lab = lm*b + log(mean(exp(b*(ly - lm))));
a   = exp(lab/b);
LL  = n*log(b) - n*lab + (b - 1)*sum(ly) - n;
end

function r = wsum_ratio(ly, lm, b)
w = exp(b*(ly - lm));
r = sum(w.*ly)/sum(w);
end

function LL = sw_ll(z, a, b, g)
y  = (z + g)/a;
LL = numel(z)*(log(b) - log(a)) + (b - 1)*sum(log(y)) - sum(y.^b);
if ~isfinite(LL), LL = -Inf; end
end

function LL = prof_shape(z, b, TG)
% max over gap (and scale, closed form) for fixed shape b.
f  = @(t) shape_fixed(z, b, t);
LL = max_over_t(f, TG);
end

function LL = shape_fixed(z, b, t)
ly  = log(z + exp(t));
lm  = max(ly);
n   = numel(z);
lab = lm*b + log(mean(exp(b*(ly - lm))));
LL  = n*log(b) - n*lab + (b - 1)*sum(ly) - n;
end

function LL = prof_scale(z, a, TG)
% max over gap and shape for fixed scale a.
f  = @(t) scale_fixed(z, a, t);
LL = max_over_t(f, TG);
end

function LL = scale_fixed(z, a, t)
% For fixed a and gap, the shape solves
%   n/b + sum(log(y/a)) - sum((y/a)^b log(y/a)) = 0, decreasing in b.
la = log(z + exp(t)) - log(a);
n  = numel(z);
h  = @(b) n/b + sum(la) - sum(exp(b*la).*la);
FLOOR = shape_floor();
if ~(h(FLOOR) > 0)
    b = FLOOR;
else
    hi = 2;
    while h(hi) > 0 && hi < 1e7, hi = 2*hi; end
    if h(hi) > 0
        b = hi;
    else
        hb = @(bb) max(min(h(bb), realmax), -realmax);
        b  = fzero(hb, [FLOOR, hi]);
    end
end
LL = n*log(b) - n*log(a) + (b - 1)*sum(la) - sum(exp(b*la));
if ~isfinite(LL), LL = -Inf; end
end

function LL = max_over_t(f, TG)
tg = linspace(TG(1), TG(2), 80);
L  = arrayfun(f, tg);
[Lb, k] = max(L);
lo = tg(max(k - 1, 1)); hi = tg(min(k + 1, numel(tg)));
[~, nR] = fminbnd(@(t) -f(t), lo, hi, optimset('TolX', 1e-9));
LL = max(Lb, -nR);
end

function [bnd, open] = find_bound(dev, t0, dir, step, lim)
% First crossing of dev = 0 moving from t0 in direction dir; dev(t0) < 0.
open = false;
if lim <= 0
    bnd = t0; open = true; return
end
tPrev = t0;
h = step;
while true
    t = t0 + dir*min(abs(tPrev - t0) + h, lim);
    if dev(t) > 0
        bnd = fzero(dev, sort([tPrev, t]));
        return
    end
    if abs(t - t0) >= lim
        bnd = t; open = true;
        return
    end
    tPrev = t;
    h = 1.5*h;
end
end

function H = fd_hessian(f, x)
% Central-difference Hessian with one Richardson step (h, h/2).
H1 = fdh(f, x, 1e-3);
H2 = fdh(f, x, 5e-4);
H  = (4*H2 - H1)/3;
H  = (H + H')/2;
end

function H = fdh(f, x, h)
k = numel(x); H = zeros(k); f0 = f(x);
for i = 1:k
    ei = zeros(k, 1); ei(i) = h;
    H(i, i) = (f(x + ei) - 2*f0 + f(x - ei))/h^2;
    for j = i+1:k
        ej = zeros(k, 1); ej(j) = h;
        H(i, j) = (f(x+ei+ej) - f(x+ei-ej) - f(x-ei+ej) + f(x-ei-ej))/(4*h^2);
        H(j, i) = H(i, j);
    end
end
end
