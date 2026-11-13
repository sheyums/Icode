addpath('C:/Users/sheyu/Box/MATLAB scripts/Rockefeller Matlab scripts/Icode'); nFail=0; nPass=0;
% =========================================================================
%  SECTION 14: SKEWNESS OUTPUTS (Nov 2026)
% =========================================================================
fprintf('\n--- Section 14: skewness outputs ---\n');
prevW14 = warning('off', 'GeneralizedHyperbolic_MLE:NoGainOverNormal');
rs14 = RandStream('mt19937ar', 'Seed', 1401);
P14  = struct('mu', 0, 'alpha', 2, 'beta', 0.8, 'delta', 1);  P14.gamma = sqrt(P14.alpha^2 - P14.beta^2);
x14  = sample_NIG_rs(rs14, P14, 800);
rN14 = GeneralizedHyperbolic_MLE(x14, 'NIG', 'NoProfile');
rH14 = GeneralizedHyperbolic_MLE(x14, 'HYP', 'NoProfile');

% 14A: NIG moment skewness has the closed form 3*beta/(alpha*sqrt(delta*gamma)).
q = rN14.Params;  gA = sqrt(q(3)^2 - q(4)^2);
gRef = 3*q(4) / (q(3)*sqrt(q(5)*gA));
[pf,~] = run_test('14A: NIG Skewness equals 3*beta/(alpha*sqrt(delta*gamma))', ...
    @() check_condition(abs(rN14.Skewness - gRef) < 1e-9*max(1, abs(gRef)), ...
    sprintf('%.12g versus %.12g', rN14.Skewness, gRef)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14B: hyperbolic moment skewness against GIG moments integrated numerically
% (no Bessel functions, so an independent route).
gRef = skew_by_integration(rH14.Params);
[pf,~] = run_test('14B: HYP Skewness matches GIG moments by numerical integration', ...
    @() check_condition(abs(rH14.Skewness - gRef) < 1e-7*max(1, abs(gRef)), ...
    sprintf('%.12g versus %.12g', rH14.Skewness, gRef)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14C: model Bowley against quartiles of a brute-force CDF on a fine grid.
bRef = bowley_by_grid(rH14.Params);
[pf,~] = run_test('14C: SkewnessBowleyModel matches a brute-force CDF', ...
    @() check_condition(abs(rH14.SkewnessBowleyModel - bRef) < 1e-6, ...
    sprintf('%.10g versus %.10g', rH14.SkewnessBowleyModel, bRef)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14D: data Bowley equals quantile()'s on tie-free data, and the hand-worked
% mid-quantile value on tied data (Fmid = .1 .4 .7 .9 at 1 2 3 10, quartiles
% 1.5, 7/3, 4.75).
qd  = quantile(x14, [0.25 0.5 0.75]);  bQ = (qd(3) + qd(1) - 2*qd(2)) / (qd(3) - qd(1));
prevT = warning('off', 'all');
rT  = GeneralizedHyperbolic_MLE(repelem([1 2 2 3 10], 2), 'NIG', 'NoProfile');
warning(prevT);
bHand = (1.5 + 4.75 - 2*7/3) / 3.25;
[pf,~] = run_test('14D: SkewnessBowleyData = quantile() without ties, mid-quantiles with ties', ...
    @() check_condition(abs(rN14.SkewnessBowleyData - bQ) < 1e-12 && ...
    abs(rT.SkewnessBowleyData - bHand) < 1e-12, ...
    sprintf('%.12g vs %.12g; tied %.12g vs %.12g', rN14.SkewnessBowleyData, bQ, rT.SkewnessBowleyData, bHand)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14E: SkewnessSE equals the delta method through Cov (an independent path:
% Cov is in the physical parameters, the code works in theta).
gfun = @(v) skew_closed(v);
gr = zeros(5, 1);  v0 = rH14.Params;
for i = [1 3 4 5]
    e = zeros(5, 1);  e(i) = 1e-5 * max(1, abs(v0(i)));
    gr(i) = (gfun(v0 + e) - gfun(v0 - e)) / (2*e(i));
end
seRef = sqrt(gr.' * rH14.Cov * gr);
[pf,~] = run_test('14E: SkewnessSE equals the delta method through Cov', ...
    @() check_condition(abs(rH14.SkewnessSE - seRef) < 1e-4*seRef, ...
    sprintf('%.8g versus %.8g', rH14.SkewnessSE, seRef)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14F: equivariance. A shift and a rescaling by 1e-6 leave every skewness
% output unchanged; negating the data flips their signs. Also the regression
% for the Wald fix: at scale 1e-6 the SEs used to come back NaN (the
% data-unit information failed the rcond test); they must scale by 1e-6.
rS = GeneralizedHyperbolic_MLE(1e-6*x14 - 7, 'NIG', 'NoProfile');
rM = GeneralizedHyperbolic_MLE(-x14, 'NIG', 'NoProfile');
f1 = [rN14.Skewness rN14.SkewnessSE rN14.SkewnessBowleyModel rN14.SkewnessBowleyData];
fS = [rS.Skewness   rS.SkewnessSE   rS.SkewnessBowleyModel   rS.SkewnessBowleyData];
fM = [rM.Skewness   rM.SkewnessSE   rM.SkewnessBowleyModel   rM.SkewnessBowleyData];
okF = all(abs(fS - f1) < 1e-5*max(1, abs(f1))) && ...
    all(abs(fM - f1.*[-1 1 -1 -1]) < 1e-5*max(1, abs(f1))) && rS.SEValid && ...
    abs(rS.SE(1)/rN14.SE(1) - 1e-6) < 1e-9 && abs(rS.SE(5)/rN14.SE(5) - 1e-6) < 1e-9;
[pf,~] = run_test('14F: skewness outputs equivariant; Wald SEs valid at scale 1e-6', ...
    @() check_condition(okF, sprintf('orig %s scaled %s negated %s', ...
    mat2str(f1, 6), mat2str(fS, 6), mat2str(fM, 6))));
nFail = nFail+~pf; nPass = nPass+pf;
warning(prevW14);

fprintf('pass %d fail %d\n',nPass,nFail);
function [passed, msg] = check_condition(cond, failMsg)
    if cond
        passed = true;  msg = '';
    else
        passed = false;  msg = failMsg;
    end
end

function x = sample_NIG_rs(rs, p, n)
% sample_NIG drawing from RandStream rs, so Section 12 does not depend on
% how many global draws the sections before it consumed.
    gam = p.gamma;
    m   = p.delta / gam;   lam = p.delta^2;
    y   = randn(rs, n, 1).^2;
    v   = m + (m^2*y)/(2*lam) - (m/(2*lam))*sqrt(4*m*lam*y + m^2*y.^2);
    u   = rand(rs, n, 1);
    flip = u > m ./ (m + v);
    v(flip) = m^2 ./ v(flip);
    x = p.mu + p.beta*v + sqrt(v).*randn(rs, n, 1);
end

function g = skew_closed(v)
% Moment skewness from Params v = [mu; lambda; alpha; beta; delta] by the
% Bessel closed form (used only as the function differentiated in 14E).
    gam = sqrt(v(3)^2 - v(4)^2);  w = v(5)*gam;  s = v(5)/gam;  L = v(2);
    K = @(nu) besselk(nu, w);
    m1 = s*K(L+1)/K(L);  m2 = s^2*K(L+2)/K(L);  m3 = s^3*K(L+3)/K(L);
    k2 = m2 - m1^2;  k3 = m3 - 3*m1*m2 + 2*m1^3;
    g  = (3*v(4)*k2 + v(4)^3*k3) / (m1 + v(4)^2*k2)^1.5;
end

function g = skew_by_integration(v)
% Moment skewness of GH with Params v, from central moments of the GIG mixing
% variable integrated numerically in u = log V (no Bessel functions).
    gam = sqrt(v(3)^2 - v(4)^2);  L = v(2);  d = v(5);  b = v(4);
    lf  = @(u) L*u - 0.5*(d^2*exp(-u) + gam^2*exp(u));
    u0  = fminbnd(@(u) -lf(u), -50, 50);  c = lf(u0);
    % Finite limits where the log-density has fallen 700 below its peak; an
    % infinite range evaluates (exp(u) - m1)^3 * 0 = Inf*0 far in the tail.
    a   = fzero(@(u) lf(u) - c + 700, [u0 - 50, u0]);
    z   = fzero(@(u) lf(u) - c + 700, [u0, u0 + 50]);
    f   = @(u, r) exp(r*u + lf(u) - c);
    o   = {'AbsTol', 0, 'RelTol', 1e-13};
    Z   = integral(@(u) f(u, 0), a, z, o{:});
    m1  = integral(@(u) f(u, 1), a, z, o{:}) / Z;
    k2  = integral(@(u) (exp(u) - m1).^2 .* f(u, 0), a, z, o{:}) / Z;
    k3  = integral(@(u) (exp(u) - m1).^3 .* f(u, 0), a, z, o{:}) / Z;
    g   = (3*b*k2 + b^3*k3) / (m1 + b^2*k2)^1.5;
end

function b = bowley_by_grid(v)
% Bowley skewness of GH with Params v from a trapezoid CDF on a fine grid.
    mu = v(1);  L = v(2);  al = v(3);  be = v(4);  de = v(5);
    gam = sqrt(al^2 - be^2);  w = 60/(al - abs(be)) + 60*de;
    xs = linspace(mu - w, mu + w, 4e6);  q = sqrt(de^2 + (xs - mu).^2);
    l  = L*log(gam/de) - (L - 0.5)*log(al) - log(besselk(L, de*gam)) + ...
         (L - 0.5)*log(q) + log(besselk(L - 0.5, al*q, 1)) - al*q + be*(xs - mu);
    F  = cumtrapz(xs, exp(l));  F = F / F(end);
    [Fu, iu] = unique(F);
    qq = interp1(Fu, xs(iu), [0.25 0.5 0.75]);
    b  = (qq(3) + qq(1) - 2*qq(2)) / (qq(3) - qq(1));
end
function [passed,msg]=run_test(label,fn), try, [passed,msg]=fn(); catch ME, passed=false; msg=ME.message; end, if passed, fprintf("  PASS  %s\n",label); else, fprintf("  FAIL  %s  ->  %s\n",label,msg); end, end
