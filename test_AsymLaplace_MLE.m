%% test_AsymLaplace_MLE.m
%
% Tests for AsymLaplace_MLE. Run from the folder containing it; prints
% [PASS]/[FAIL] per test and a Summary line.
%
% Ground truth shares no code with the fitter: the density summed
% directly, a brute-force profile over a grid between the data points,
% Nelder-Mead over the other two parameters for the profile bounds, and a
% Monte Carlo estimate of E[score*score'] for the Fisher information.
%
%  1. Global exact optimum      LogLike = direct density sum = brute force
%  2. Output fields             all documented fields exist
%  3. AIC / AICc / BIC          AIC is plain AIC (Nov 2026: it held AICc)
%  4. Equivariance              a*x + b for a = 1e-15, 1e12, -3 (kappa -> 1/kappa)
%  5. Fisher information        closed form = Monte Carlo E[g g'] (2e6 draws)
%  6. Profile bounds            each on the deviance threshold (independent NM)
%  7. Wald formulas             mu linear, sigma and kappa on the log scale
%  8. One-sided data            warns OneSided, SEValid false, limit LL >= LL
%  9. Errors                    n = 4, constant, bad mode, complex input
% 10. NaN/Inf removal           warns, N reduced, same fit as clean data
% 11. Side effects              global random stream and warning state untouched;
%                               bootstrap reproducible
% 12. NOPROFILE                 CI = Wald, CIMethod 'Wald'
% 13. Fit output                FitX contains mu, Fit integrates to ~1
% 14. Symmetric Laplace         kappa = 1 recovered, its interval contains 1
% 15. SmallSample warning       at n = 50, not at n = 100
%
% Claude-written, like the fitter it tests: passing is self-consistency
% against independent computations, not external validation.

nPass = 0; nFail = 0;
fprintf('\nRunning tests for AsymLaplace_MLE...\n\n');
w0 = warning;
warning('off', 'AsymLaplace_MLE:SmallSample');
rs = RandStream('twister', 'Seed', 20261104);

%% 1. Global exact optimum
[ok, msg] = deal(true, '');
cases = [50 1; 200 2; 1000 0.5; 300 3];
for i = 1:size(cases, 1)
    x = ald_draw(rs, cases(i,1), 0.7, 1.5, cases(i,2));
    if i == 4, x = round(x*4)/4; end                 % ties
    R = AsymLaplace_MLE(x);
    llD = ald_ll(x, R.Params);
    llB = ald_brute(x);
    if abs(llD - R.LogLike) > 1e-8*abs(llD) || llB > R.LogLike + 1e-8*abs(llD)
        ok = false; msg = sprintf('case %d: LogLike %.10g direct %.10g brute %.10g', ...
            i, R.LogLike, llD, llB);
    end
end
[nPass, nFail] = report(1, 'global exact optimum (direct sum, brute-force profile)', ok, msg, nPass, nFail);

%% 2. Output fields
x = ald_draw(rs, 400, 0, 1, 1.4);
R = AsymLaplace_MLE(x, 'FIT');
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
for a = [1e-15, 1e12, -3]
    Ra = AsymLaplace_MLE(a*x + 5*a);
    Kw = R.Params(3); if a < 0, Kw = 1/Kw; end
    ePar = [abs((Ra.Params(1) - 5*a)/(a*R.Params(1)) - 1), ...
            abs(Ra.Params(2)/(abs(a)*R.Params(2)) - 1), abs(Ra.Params(3)/Kw - 1)];
    eLL  = abs(Ra.LogLike + n*log(abs(a)) - R.LogLike);
    seW  = R.SE .* [abs(a), abs(a), 1];
    if a < 0, seW(3) = R.SE(3)/R.Params(3)^2; end    % delta method for 1/kappa
    eSE  = max(abs(Ra.SE ./ seW - 1));
    loW  = [a*R.CILower(1) + 5*a, abs(a)*R.CILower(2), R.CILower(3)];
    hiW  = [a*R.CIUpper(1) + 5*a, abs(a)*R.CIUpper(2), R.CIUpper(3)];
    if a < 0
        [loW(1), hiW(1)] = deal(hiW(1), loW(1));
        [loW(3), hiW(3)] = deal(1/R.CIUpper(3), 1/R.CILower(3));
    end
    eCI = max(abs([Ra.CILower ./ loW, Ra.CIUpper ./ hiW] - 1));
    if max(ePar) > 1e-9 || eLL > 1e-7 || eSE > 1e-9 || eCI > 1e-6
        ok = false; msg = sprintf('a = %g: params %.2g LL %.2g SE %.2g CI %.2g', ...
            a, max(ePar), eLL, eSE, eCI);
    end
end
[nPass, nFail] = report(4, 'equivariant under a*x + b (a = 1e-15, 1e12, -3)', ok, msg, nPass, nFail);

%% 5. Fisher information
p = [0.3 1.2 0.6]; X = ald_draw(rs, 2e6, p(1), p(2), p(3));
s = p(2); K = p(3); Rt = X >= p(1); d = X - p(1);
G = [(K/s)*Rt - (~Rt)/(K*s), ...
     -1/s + (K/s^2)*max(d,0) + max(-d,0)/(K*s^2), ...
     1/K - 2*K/(1+K^2) - max(d,0)/s + max(-d,0)/(K^2*s)];
Imc = (G'*G)/size(G, 1);
Rf = AsymLaplace_MLE(ald_draw(rs, 800, p(1), p(2), p(3)), 'NOPROFILE');
Iq = inv(Rf.Cov)/Rf.N;                                % information at Rf's estimate
X2 = ald_draw(rs, 2e6, Rf.Params(1), Rf.Params(2), Rf.Params(3));
s = Rf.Params(2); K = Rf.Params(3); Rt = X2 >= Rf.Params(1); d = X2 - Rf.Params(1);
G2 = [(K/s)*Rt - (~Rt)/(K*s), -1/s + (K/s^2)*max(d,0) + max(-d,0)/(K*s^2), ...
      1/K - 2*K/(1+K^2) - max(d,0)/s + max(-d,0)/(K^2*s)];
Imc2 = (G2'*G2)/size(G2, 1);
relI = max(abs(Iq(:) - Imc2(:))) / max(abs(Imc2(:)));
[nPass, nFail] = report(5, 'closed-form Fisher information = Monte Carlo E[g g''] (2e6 draws)', ...
    relI < 0.01 && abs(Imc(1,1) - 1/p(2)^2) < 0.01, sprintf('max relative difference %.3g', relI), nPass, nFail);

%% 6. Profile bounds on the deviance threshold
x6 = ald_draw(rs, 150, 0.7, 1.5, 1.3);
R6 = AsymLaplace_MLE(x6);
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
z = sqrt(2)*erfinv(0.95); P = R6.Params; S = R6.SE;
ok = max(abs(R6.CIWaldLower - [P(1) - z*S(1), P(2)*exp(-z*S(2)/P(2)), P(3)*exp(-z*S(3)/P(3))])) < 1e-12 && ...
     max(abs(R6.CIWaldUpper - [P(1) + z*S(1), P(2)*exp( z*S(2)/P(2)), P(3)*exp( z*S(3)/P(3))])) < 1e-12 && ...
     all(R6.CIWaldLower < P & P < R6.CIWaldUpper) && R6.CILevel == 0.95;
[nPass, nFail] = report(7, 'Wald intervals follow their formulas', ok, 'formula mismatch', nPass, nFail);

%% 8. One-sided data
xe = -log(rand(rs, 500, 1));
lastwarn('');
R8 = AsymLaplace_MLE(xe);
[~, wid] = lastwarn;
ok = strcmp(wid, 'AsymLaplace_MLE:OneSided') && R8.Diagnostics.OneSided && ~R8.SEValid && ...
     R8.Diagnostics.LogLikeOneSidedLimit >= R8.LogLike;
[nPass, nFail] = report(8, 'one-sided (exponential) data flagged', ok, ...
    sprintf('warning "%s", SEValid %d', wid, R8.SEValid), nPass, nFail);

%% 9. Errors
ids = {'AsymLaplace_MLE:TooFewObservations', 'AsymLaplace_MLE:ConstantData', ...
       'AsymLaplace_MLE:BadMode', 'AsymLaplace_MLE:BadInput'};
calls = {@() AsymLaplace_MLE([1 2 3 4]), @() AsymLaplace_MLE(ones(10,1)), ...
         @() AsymLaplace_MLE(x, 3), @() AsymLaplace_MLE(x + 1i)};
[ok, msg] = deal(true, '');
for i = 1:numel(calls)
    try
        calls{i}(); ok = false; msg = sprintf('no error for %s', ids{i});
    catch e
        if ~strcmp(e.identifier, ids{i}), ok = false; msg = e.identifier; end
    end
end
[nPass, nFail] = report(9, 'error identifiers (n = 4, constant, bad mode, complex)', ok, msg, nPass, nFail);

%% 10. NaN / Inf removal
lastwarn('');
R10 = AsymLaplace_MLE([x; NaN; Inf; -Inf]);
[~, wid] = lastwarn;
ok = strcmp(wid, 'AsymLaplace_MLE:NonfiniteDataRemoved') && R10.N == numel(x) && ...
     isequal(R10.Params, R.Params);
[nPass, nFail] = report(10, 'NaN/Inf removed with a warning, same fit', ok, wid, nPass, nFail);

%% 11. Side effects
rng(5); a1 = rand; rng(5);
warning('off', 'MATLAB:singularMatrix');
Rb1 = AsymLaplace_MLE(x, 'BOOTSTRAP');
a2 = rand;
st = warning('query', 'MATLAB:singularMatrix');
Rb2 = AsymLaplace_MLE(x, 'BOOTSTRAP');
ok = a1 == a2 && strcmp(st.state, 'off') && isequal(Rb1.SEBoot, Rb2.SEBoot) && ...
     all(Rb1.CIBootLower < Rb1.CIBootUpper) && Rb1.Diagnostics.NBoot == 500;
warning(w0); warning('off', 'AsymLaplace_MLE:SmallSample');
[nPass, nFail] = report(11, 'global random stream and warning state untouched; bootstrap reproducible', ...
    ok, 'side effect or irreproducible bootstrap', nPass, nFail);

%% 12. NOPROFILE
R12 = AsymLaplace_MLE(x, 'NOPROFILE');
ok = strcmp(R12.CIMethod, 'Wald') && isequal(R12.CILower, R12.CIWaldLower) && ...
     isequal(R12.CIUpper, R12.CIWaldUpper) && isequal(R12.Params, R.Params) && ...
     strcmp(R.CIMethod, 'profile');
[nPass, nFail] = report(12, 'NOPROFILE returns the Wald intervals', ok, 'mismatch', nPass, nFail);

%% 13. Fit output
mu = R.Params(1);
xx = linspace(mu - 40*R.Params(2), mu + 40*R.Params(2), 400001)';
f  = ald_pdf_ref(xx, R.Params);
ok = any(R.FitX == mu) && issorted(R.FitX) && numel(R.Fit) == numel(R.FitX) && ...
     abs(trapz(xx, f) - 1) < 1e-6 && ...
     max(abs(R.Fit - ald_pdf_ref(R.FitX, R.Params))) < 1e-12;
[nPass, nFail] = report(13, 'FitX contains mu; Fit is the density', ok, 'Fit mismatch', nPass, nFail);

%% 14. Symmetric Laplace
xs = ald_draw(rs, 2000, 2, 0.5, 1);
R14 = AsymLaplace_MLE(xs);
ok = R14.CILower(3) <= 1 && 1 <= R14.CIUpper(3) && abs(R14.Params(3) - 1) < 0.1 && ...
     abs(R14.Params(1) - 2) < 0.1 && abs(R14.Params(2) - 0.5) < 0.05;
[nPass, nFail] = report(14, 'symmetric Laplace (kappa = 1) recovered; kappa CI contains 1', ok, ...
    mat2str([R14.Params R14.CILower(3) R14.CIUpper(3)], 4), nPass, nFail);

%% 15. SmallSample warning
warning(w0);
warning('off', 'AsymLaplace_MLE:OneSided');
lastwarn(''); AsymLaplace_MLE(ald_draw(rs, 50, 0, 1, 1)); [~, w50] = lastwarn;
lastwarn(''); AsymLaplace_MLE(ald_draw(rs, 100, 0, 1, 1)); [~, w100] = lastwarn;
ok = strcmp(w50, 'AsymLaplace_MLE:SmallSample') && ~strcmp(w100, 'AsymLaplace_MLE:SmallSample');
[nPass, nFail] = report(15, 'SmallSample warning at n = 50, not at n = 100', ok, ...
    sprintf('n=50 "%s", n=100 "%s"', w50, w100), nPass, nFail);

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

function x = ald_draw(rs, n, m, s, K)
r = rand(rs, n, 1) < 1/(1 + K^2);
x = zeros(n, 1);
x(r)  = m - s/K*log(rand(rs, sum(r), 1));
x(~r) = m + K*s*log(rand(rs, sum(~r), 1));
end

function ll = ald_ll(x, p)
f  = ald_pdf_ref(x, p);
ll = sum(log(f));
end

function f = ald_pdf_ref(x, p)
mu = p(1); s = p(2); K = p(3);
f = K/(s*(1 + K^2)) * exp(-(x >= mu).*K.*(x - mu)/s + (x < mu).*(x - mu)/(K*s));
end

function best = ald_brute(x)
% profile over a grid of mu (every data point plus 20 between each pair),
% inner sigma and kappa by the closed form
xs = sort(x); n = numel(x);
g  = unique([xs; reshape(xs(1:end-1) + diff(xs)*linspace(0, 1, 21), [], 1)]);
best = -Inf;
for mu = g'
    S1 = sum(max(x - mu, 0)); S2 = sum(max(mu - x, 0));
    if S1 <= 0 || S2 <= 0, continue; end
    k = (S2/S1)^0.25; s = (k*S1 + S2/k)/n;
    best = max(best, sum(log(ald_pdf_ref(x, [mu s k]))));
end
end

function LL = nm_profile(x, j, v, p)
% maximize the log-likelihood with parameter j fixed at v, multistart NM;
% for j = 2, 3 also over every data point as mu with the remaining
% parameter by fminbnd (the inner optimum in mu is at a data point).
o = optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxFunEvals', 1e4, 'MaxIter', 1e4);
f = @(q) -sum(log(ald_pdf_ref(x, q)));
best = Inf; rs = RandStream('twister', 'Seed', 7);
for st = 1:8
    switch j
        case 1, g = @(u) f([v exp(u(1)) exp(u(2))]); u0 = [log(p(2)) + 0.3*randn(rs), log(p(3)) + 0.3*randn(rs)];
        case 2, g = @(u) f([u(1) v exp(u(2))]);      u0 = [p(1) + 0.3*p(2)*randn(rs), log(p(3)) + 0.3*randn(rs)];
        case 3, g = @(u) f([u(1) exp(u(2)) v]);      u0 = [p(1) + 0.3*p(2)*randn(rs), log(p(2)) + 0.3*randn(rs)];
    end
    u = fminsearch(g, u0, o); u = fminsearch(g, u, o);
    best = min(best, g(u));
end
if j > 1
    for mu = sort(x)'
        if j == 2, h = @(t) f([mu v exp(t)]); else, h = @(t) f([mu exp(t) v]); end
        best = min(best, h(fminbnd(h, -15, 15, o)));
    end
end
LL = -best;
end
