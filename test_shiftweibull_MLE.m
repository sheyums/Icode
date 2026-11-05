%% test_shiftweibull_MLE.m
%
% Tests for shiftweibull_MLE. Run from the folder containing it; prints
% [PASS]/[FAIL] per test and a Summary line.
%
% Ground truth shares no code with the fitter: the density summed
% directly, a 40-start Nelder-Mead in its own coordinates, Nelder-Mead
% over the other two parameters for the profile bounds, and a
% finite-difference Hessian in PHYSICAL coordinates (a, b, c).
%
%  1. Global optimum            LogLike = direct sum >= independent NM multistart
%  2. Output fields             all documented fields exist
%  3. AIC / AICc / BIC          AIC is plain AIC (Nov 2026: it held AICc)
%  4. Equivariance              s*x + t for s = 1e-6, 1e6
%  5. SEs                       = an independent physical-coordinate Hessian
%  6. Profile bounds            each on the deviance threshold (independent NM)
%  7. Wald formulas             a, b on the log scale, c on log(min(x) - c)
%  8. Shape floor               true shape 0.8: ShapeAtFloor, SE NaN, warning
%  9. Gumbel limit              ShiftAtLowerLimit, open lower shift end
% 10. Side effects              warning state and global random stream untouched
%                               (the Nov 2026 version turned ALL warnings on);
%                               bootstrap reproducible
% 11. Errors                    n = 4, 2 distinct values, bad mode, complex
% 12. NaN/Inf removal           warns, same fit as clean data
% 13. NOPROFILE                 CI = Wald, CIMethod 'Wald'
% 14. Fit output                Fit is the density at FitX; sum(log(Fit)) = LogLike
% 15. Non-regular range         shape <= 2: NonRegular set, SEValid false
% 16. SmallSample warning       at n = 30, not at n = 50
%
% Claude-written, like the fitter it tests: passing is self-consistency
% against independent computations, not external validation.

nPass = 0; nFail = 0;
fprintf('\nRunning tests for shiftweibull_MLE...\n\n');
w0 = warning;
warning('off', 'shiftweibull_MLE:SmallSample');
rs = RandStream('twister', 'Seed', 20261104);

%% 1. Global optimum
[ok, msg] = deal(true, '');
cases = [2 1.5 0 200; 2 3 5 500; 1 6 -2 300; 2 2.5 10 60];
for i = 1:size(cases, 1)
    x = wb_draw(rs, cases(i,4), cases(i,1), cases(i,2), cases(i,3));
    R = shiftweibull_MLE(x);
    llD = wb_ll(x, R.Params);
    llN = nm_max(x);
    if abs(llD - R.LogLike) > 1e-8*abs(llD) || llN > R.LogLike + 1e-6
        ok = false; msg = sprintf('case %d: LogLike %.10g direct %.10g NM %.10g', ...
            i, R.LogLike, llD, llN);
    end
end
[nPass, nFail] = report(1, 'global optimum (direct sum, 40-start Nelder-Mead)', ok, msg, nPass, nFail);

%% 2. Output fields
x = wb_draw(rs, 400, 2, 3, 5);
R = shiftweibull_MLE(x, 'FIT');
need = {'Model','ParameterNames','Params','LogLike','AIC','AICc','BIC','N','NParameters', ...
    'SE','Cov','SEValid','CILower','CIUpper','CIOpen','CIMethod','CIWaldLower', ...
    'CIWaldUpper','CILevel','SEBoot','CIBootLower','CIBootUpper','Fit','FitX','Diagnostics'};
miss = need(~isfield(R, need));
[nPass, nFail] = report(2, 'documented output fields present', isempty(miss), ...
    ['missing: ' strjoin(miss, ', ')], nPass, nFail);

%% 3. AIC / AICc / BIC
k = 3; n = R.N;
ok = abs(R.AIC - (2*k - 2*R.LogLike)) < 1e-9 && ...
     abs(R.AICc - (2*k - 2*R.LogLike + 2*k*(k+1)/(n-k-1))) < 1e-9 && ...
     abs(R.BIC - (k*log(n) - 2*R.LogLike)) < 1e-9;
[nPass, nFail] = report(3, 'AIC is plain AIC; AICc and BIC formulas', ok, ...
    sprintf('AIC %.6f AICc %.6f BIC %.6f', R.AIC, R.AICc, R.BIC), nPass, nFail);

%% 4. Equivariance
[ok, msg] = deal(true, '');
for s = [1e-6, 1e6]
    Rs = shiftweibull_MLE(s*x - 7);
    eP  = max(abs([Rs.Params(1)/(s*R.Params(1)), Rs.Params(2)/R.Params(2), ...
                   (Rs.Params(3) + 7)/(s*R.Params(3))] - 1));
    eLL = abs(Rs.LogLike + n*log(s) - R.LogLike);
    eSE = max(abs(Rs.SE ./ (R.SE .* [s 1 s]) - 1));
    eCI = max(abs([Rs.CILower(1)/s, Rs.CILower(2), (Rs.CILower(3) + 7)/s] ./ R.CILower - 1));
    if eP > 1e-6 || eLL > 1e-6 || eSE > 1e-4 || eCI > 1e-5
        ok = false; msg = sprintf('s = %g: params %.2g LL %.2g SE %.2g CI %.2g', s, eP, eLL, eSE, eCI);
    end
end
[nPass, nFail] = report(4, 'equivariant under s*x + t (s = 1e-6, 1e6)', ok, msg, nPass, nFail);

%% 5. SEs against an independent physical-coordinate Hessian
seRef = phys_se(x, R.Params);
rel = max(abs(R.SE - seRef) ./ seRef);
[nPass, nFail] = report(5, 'SEs = independent physical-coordinate Hessian (rel < 1e-3)', ...
    R.SEValid && rel < 1e-3, sprintf('max relative difference %.3g', rel), nPass, nFail);

%% 6. Profile bounds on the deviance threshold
x6 = wb_draw(rs, 150, 2, 3, 5);
R6 = shiftweibull_MLE(x6);
crit = 2*erfinv(0.95)^2;
[ok, msg] = deal(true, '');
for j = 1:3
    for bd = [R6.CILower(j), R6.CIUpper(j)]
        dv = 2*(R6.LogLike - nm_profile(x6, j, bd, R6.Params));
        if abs(dv - crit) > 1e-3
            ok = false; msg = sprintf('param %d bound %.6g: deviance %.6f', j, bd, dv);
        end
    end
end
[nPass, nFail] = report(6, 'profile bounds on the 3.8415 deviance (independent Nelder-Mead)', ok, msg, nPass, nFail);

%% 7. Wald formulas
z = sqrt(2)*erfinv(0.95); P = R.Params; S = R.SE; g = min(x) - P(3);
ok = max(abs(R.CIWaldLower - [P(1)*exp(-z*S(1)/P(1)), P(2)*exp(-z*S(2)/P(2)), min(x) - g*exp( z*S(3)/g)])) < 1e-9 && ...
     max(abs(R.CIWaldUpper - [P(1)*exp( z*S(1)/P(1)), P(2)*exp( z*S(2)/P(2)), min(x) - g*exp(-z*S(3)/g)])) < 1e-9 && ...
     R.CIWaldUpper(3) < min(x) && R.CILevel == 0.95;
[nPass, nFail] = report(7, 'Wald intervals follow their formulas; shift stays below min(x)', ok, 'formula mismatch', nPass, nFail);

%% 8. Shape floor
x8 = wb_draw(rs, 400, 3, 0.8, 1);
lastwarn(''); R8 = shiftweibull_MLE(x8); [~, wid] = lastwarn;
ok = strcmp(wid, 'shiftweibull_MLE:ShapeAtFloor') && R8.Diagnostics.ShapeAtFloor && ...
     abs(R8.Params(2) - 1.001) < 1e-9 && all(isnan(R8.SE)) && ~R8.SEValid && R8.Params(3) < min(x8);
[nPass, nFail] = report(8, 'true shape 0.8: on the floor, flagged, SEs NaN', ok, wid, nPass, nFail);

%% 9. Gumbel limit
x9 = log(-log(rand(rs, 500, 1)));                    % Gumbel (minimum)
lastwarn(''); R9 = shiftweibull_MLE(x9); [~, wid] = lastwarn;
ok = strcmp(wid, 'shiftweibull_MLE:ShiftAtLowerLimit') && R9.Diagnostics.ShiftAtLowerLimit && ...
     R9.CIOpen(3) && R9.CILower(3) == -Inf && ~R9.SEValid;
[nPass, nFail] = report(9, 'Gumbel-minimum data: limit flagged, shift open below', ok, wid, nPass, nFail);

%% 10. Side effects
rng(5); a1 = rand; rng(5);
warning('off', 'MATLAB:singularMatrix');
Rb1 = shiftweibull_MLE(x, 'BOOTSTRAP NOPROFILE');
a2 = rand;
st = warning('query', 'MATLAB:singularMatrix');
Rb2 = shiftweibull_MLE(x, 'BOOTSTRAP NOPROFILE');
ok = a1 == a2 && strcmp(st.state, 'off') && isequal(Rb1.SEBoot, Rb2.SEBoot) && ...
     all(Rb1.CIBootLower < Rb1.CIBootUpper) && Rb1.Diagnostics.NBoot == 300;
warning(w0); warning('off', 'shiftweibull_MLE:SmallSample');
[nPass, nFail] = report(10, 'warning state and global random stream untouched; bootstrap reproducible', ...
    ok, sprintf('rand equal %d, warning %s', a1 == a2, st.state), nPass, nFail);

%% 11. Errors
ids = {'shiftweibull_MLE:TooFewObservations', 'shiftweibull_MLE:TooFewDistinct', ...
       'shiftweibull_MLE:BadMode', 'shiftweibull_MLE:BadInput'};
calls = {@() shiftweibull_MLE([1 2 3 4]), @() shiftweibull_MLE([1 1 2 2 2 1]), ...
         @() shiftweibull_MLE(x, 3), @() shiftweibull_MLE(x + 1i)};
[ok, msg] = deal(true, '');
for i = 1:numel(calls)
    try
        calls{i}(); ok = false; msg = sprintf('no error for %s', ids{i});
    catch e
        if ~strcmp(e.identifier, ids{i}), ok = false; msg = e.identifier; end
    end
end
[nPass, nFail] = report(11, 'error identifiers (n = 4, 2 distinct, bad mode, complex)', ok, msg, nPass, nFail);

%% 12. NaN / Inf removal
lastwarn(''); R12 = shiftweibull_MLE([x; NaN; Inf], 'NOPROFILE'); [~, wid] = lastwarn;
R12b = shiftweibull_MLE(x, 'NOPROFILE');
ok = strcmp(wid, 'shiftweibull_MLE:NonfiniteDataRemoved') && R12.N == numel(x) && ...
     isequal(R12.Params, R12b.Params);
[nPass, nFail] = report(12, 'NaN/Inf removed with a warning, same fit', ok, wid, nPass, nFail);

%% 13. NOPROFILE
ok = strcmp(R12b.CIMethod, 'Wald') && isequal(R12b.CILower, R12b.CIWaldLower) && ...
     isequal(R12b.CIUpper, R12b.CIWaldUpper) && isequal(R12b.Params, R.Params) && ...
     strcmp(R.CIMethod, 'profile');
[nPass, nFail] = report(13, 'NOPROFILE returns the Wald intervals, same fit', ok, 'mismatch', nPass, nFail);

%% 14. Fit output
fRef = wb_pdf(R.FitX, R.Params);
ok = issorted(R.FitX) && numel(R.Fit) == R.N && max(abs(R.Fit - fRef) ./ fRef) < 1e-12 && ...
     abs(sum(log(R.Fit)) - R.LogLike) < 1e-8*abs(R.LogLike);
[nPass, nFail] = report(14, 'Fit is the density at FitX; sum(log(Fit)) = LogLike', ok, 'Fit mismatch', nPass, nFail);

%% 15. Non-regular range
x15 = wb_draw(rs, 300, 2, 1.4, 0);
R15 = shiftweibull_MLE(x15);
ok = R15.Params(2) <= 2 && R15.Diagnostics.NonRegular && ~R15.SEValid && ~R15.Diagnostics.ShapeAtFloor;
[nPass, nFail] = report(15, 'shape <= 2: NonRegular set, SEValid false', ok, ...
    sprintf('shape %.3f', R15.Params(2)), nPass, nFail);

%% 16. SmallSample warning
warning(w0);
lastwarn(''); shiftweibull_MLE(wb_draw(rs, 30, 2, 3, 5), 'NOPROFILE'); [~, w50] = lastwarn;
lastwarn(''); shiftweibull_MLE(wb_draw(rs, 50, 2, 3, 5), 'NOPROFILE'); [~, w100] = lastwarn;
ok = strcmp(w50, 'shiftweibull_MLE:SmallSample') && ~strcmp(w100, 'shiftweibull_MLE:SmallSample');
[nPass, nFail] = report(16, 'SmallSample warning at n = 30, not at n = 50', ok, ...
    sprintf('n=30 "%s", n=50 "%s"', w50, w100), nPass, nFail);

warning(w0);
fprintf('\nSummary: %d passed, %d failed, %d total\n', nPass, nFail, nPass + nFail);

%% ------------------------------------------------------------------------
function [np, nf] = report(k, name, ok, msg, np, nf)
if ok
    fprintf('[PASS] %d. %s\n', k, name); np = np + 1;
else
    fprintf('[FAIL] %d. %s: %s\n', k, name, msg); nf = nf + 1;
end
end

function x = wb_draw(rs, n, a, b, c)
x = c + a*(-log(rand(rs, n, 1))).^(1/b);
end

function f = wb_pdf(x, p)
u = (x - p(3))/p(1);
f = (p(2)/p(1)) * u.^(p(2) - 1) .* exp(-u.^p(2));
end

function ll = wb_ll(x, p)
ll = sum(log(wb_pdf(x, p)));
end

function nll = wb_nll_phys(x, q)
if q(1) <= 0 || q(2) < 1.001 || q(3) >= min(x), nll = Inf; return; end
nll = -wb_ll(x, q);
end

function llBest = nm_max(x)
% 40-start Nelder-Mead in (log a, log(b - 1.001), log(min(x) - c))
mn = min(x); s = std(x);
f = @(q) wb_nll_phys(x, [exp(q(1)), 1.001 + exp(q(2)), mn - exp(q(3))]);
o = optimset('Display', 'off', 'MaxFunEvals', 2e4, 'MaxIter', 2e4, 'TolX', 1e-10, 'TolFun', 1e-10);
best = Inf; rs = RandStream('twister', 'Seed', 9);
for k = 1:40
    q0 = [log(s*(0.3 + 3*rand(rs))), log(0.05 + 5*rand(rs)), log(s) + 3*randn(rs)];
    q = fminsearch(f, q0, o); q = fminsearch(f, q, o);
    v = f(q); if isfinite(v), best = min(best, v); end
end
llBest = -best;
end

function LL = nm_profile(x, j, v, p)
% parameter j fixed at v; maximize over the other two by multistart NM
mn = min(x); F = 1.001;
o = optimset('Display', 'off', 'MaxFunEvals', 2e4, 'MaxIter', 2e4, 'TolX', 1e-11, 'TolFun', 1e-11);
best = Inf; rs = RandStream('twister', 'Seed', 3);
for k = 1:25
    switch j
        case 1, g = @(u) wb_nll_phys(x, [v, F + exp(u(1)), mn - exp(u(2))]);
                u0 = [log(max(p(2) - F, 1e-3)) + randn(rs), log(mn - p(3)) + randn(rs)];
        case 2, g = @(u) wb_nll_phys(x, [exp(u(1)), v, mn - exp(u(2))]);
                u0 = [log(p(1)) + 0.5*randn(rs), log(mn - p(3)) + randn(rs)];
        case 3, g = @(u) wb_nll_phys(x, [exp(u(1)), F + exp(u(2)), v]);
                u0 = [log(p(1)) + 0.5*randn(rs), log(max(p(2) - F, 1e-3)) + randn(rs)];
    end
    u = fminsearch(g, u0, o); u = fminsearch(g, u, o);
    best = min(best, g(u));
end
LL = -best;
end

function se = phys_se(x, p)
% central-difference Hessian of the NLL in (a, b, c), relative steps 1e-4
h = 1e-4*abs(p(:)); q = p(:);
f = @(v) wb_nll_phys(x, v');
H = zeros(3);
for i = 1:3
    for j = i:3
        ei = zeros(3,1); ei(i) = h(i); ej = zeros(3,1); ej(j) = h(j);
        H(i,j) = (f(q+ei+ej) - f(q+ei-ej) - f(q-ei+ej) + f(q-ei-ej))/(4*h(i)*h(j));
        H(j,i) = H(i,j);
    end
end
se = sqrt(diag(inv(H)))';
end
