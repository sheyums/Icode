%% test_multigaussian_MLE.m
%
% Tests for multigaussian_MLE. Run from the folder holding the function;
% each test prints [PASS] or [FAIL], then a summary. Ground truth is written
% here from the definitions and shares no code with the fitter: the
% closed-form K = 1 MLE, a direct mixture density, central-difference
% Hessians, and an independent Nelder-Mead multistart restricted to
% NON-SPURIOUS solutions (every sigma above 1.01x the 1e-3*IQR floor, and at
% least 5 effective observations per component). Guards marked
% (audit) failed on the pre-Nov-2026 code.
%
%  1. K = 1 is the closed-form MLE (mean, n-divisor SD)               (audit)
%  2. LogLike equals the independent density summed at Params
%  3. Fit integrates to 1 and matches the density; FitX sorted       (audit)
%  4. Equivariant under c*x + b, K = 1..3                            (audit)
%  5. AIC is plain AIC; AICc and BIC formulas                       (audit)
%  6. Sample-size guard: n >= 3K + 1                                (audit)
%  7. The caller's global random stream is untouched                (audit)
%  8. 'Select' picks K = 2 (bimodal) and K = 1 (unimodal) by BIC
%  9. Input errors and NaN removal                                  (audit)
% 10. K = 2 on three populations reaches the independent optimum    (audit)
% 11. K = 3 and K = 4 reach the independent non-degenerate optimum  (audit)
% 12. No collapsed component when a non-degenerate fit exists
% 13. SE: K = 1 analytic; K = 2 equals an independent Hessian
% 14. K = 1 profile CI ends sit on the closed-form profile threshold
% 15. Intervals contain the estimate and respect sigma > 0, w in (0,1)
% 16. SEs and CIs equivariant under c*x + b

nPassed = 0;  nFailed = 0;
fprintf('\nRunning tests for multigaussian_MLE...\n\n');
chiHalf = 0.5 * (sqrt(2)*erfinv(0.95))^2;

%% Test 1: K = 1 is the closed-form MLE
try
    rs = RandStream('twister', 'Seed', 1);
    x  = 5 + 2*randn(rs, 40, 1);
    [~, r] = evalc('multigaussian_MLE(x, 1)');
    n  = numel(x);
    LLmle = -n/2*(log(2*pi*var(x, 1)) + 1);
    assert(abs(r.Params(1) - mean(x)) < 1e-12 && ...
        abs(r.Params(2) - std(x, 1)) < 1e-12 && ...
        abs(r.LogLike - LLmle) < 1e-10, ...
        sprintf('sigma %.10g vs %.10g, LL %.10g vs %.10g', r.Params(2), ...
        std(x, 1), r.LogLike, LLmle));
    fprintf('[PASS] 1. K = 1 is the closed-form MLE\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 1. K = 1 MLE: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 2: LogLike equals the independent density sum
try
    rs = RandStream('twister', 'Seed', 2);
    x  = [randn(rs, 300, 1); 5 + 0.7*randn(rs, 200, 1)];
    [~, r] = evalc('multigaussian_MLE(x, [1 2 3])');
    for i = 1:3
        d = abs(sum(log(mix_pdf(x, r(i).Params))) - r(i).LogLike);
        assert(d < 1e-8 * abs(r(i).LogLike), sprintf('K = %d: %.3g', i, d));
        assert(abs(sum(r(i).Params(:, 3)) - 1) < 1e-12, 'weights do not sum to 1');
        assert(issorted(r(i).Params(:, 1)), 'rows not sorted by mu');
    end
    fprintf('[PASS] 2. LogLike equals the independent density sum\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 2. LogLike consistency: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 3: Fit integrates to 1 and matches the density
try
    [~, r] = evalc('multigaussian_MLE(x, 2, ''fit'')');
    Zint = integral(@(t) mix_pdf(t, r.Params), -Inf, Inf, ...
        'AbsTol', 1e-13, 'RelTol', 1e-11);
    assert(abs(Zint - 1) < 1e-8, sprintf('integral %.12g', Zint));
    assert(numel(r.Fit) == numel(x) && issorted(r.FitX), 'Fit/FitX shape');
    assert(max(abs(r.Fit - mix_pdf(r.FitX, r.Params)) ./ r.Fit) < 1e-10, ...
        'Fit disagrees with the independent density');
    fprintf('[PASS] 3. Fit integrates to 1 and matches the density\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 3. Fit: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 4: equivariance under c*x + b
% Old code: 1e-6*x moved the K = 2 LogLike by 191 beyond n*log(c).
try
    rs = RandStream('twister', 'Seed', 4);
    x  = [randn(rs, 300, 1); 6 + 0.5*randn(rs, 200, 1); -5 + 2*randn(rs, 150, 1)];
    msg = '';
    for K = 1:3
        [~, r0] = evalc('multigaussian_MLE(x, K)');
        for cb = [1e-6 0; 1e6 0; 1 1e8; -1e3 7]'
            c = cb(1);  b = cb(2);
            [~, r] = evalc('multigaussian_MLE(c*x + b, K)');
            e1 = abs(r.LogLike + numel(x)*log(abs(c)) - r0.LogLike);
            if e1 > 1e-6
                msg = [msg sprintf('K=%d c=%g b=%g: dLL %.3g; ', K, c, b, e1)]; %#ok<AGROW>
            end
        end
    end
    assert(isempty(msg), msg);
    fprintf('[PASS] 4. Equivariant under c*x + b (K = 1..3)\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 4. Equivariance: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 5: AIC / AICc / BIC
try
    rs = RandStream('twister', 'Seed', 5);
    x  = randn(rs, 60, 1);
    [~, r] = evalc('multigaussian_MLE(x, [1 2])');
    for i = 1:2
        p = 3*i - 1;  LL = r(i).LogLike;
        assert(abs(r(i).AIC - (2*p - 2*LL)) < 1e-9, 'AIC is not 2p - 2LL');
        assert(abs(r(i).AICc - (2*p - 2*LL + 2*p*(p+1)/(60 - p - 1))) < 1e-9, 'AICc');
        assert(abs(r(i).BIC - (p*log(60) - 2*LL)) < 1e-9, 'BIC');
        assert(r(i).NParameters == p && r(i).N == 60, 'N / NParameters');
    end
    fprintf('[PASS] 5. AIC, AICc and BIC formulas\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 5. Information criteria: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 6: sample-size guard
% Old code: n = 2..7 at K = 2 returned AIC from -56 to Inf.
try
    rs = RandStream('twister', 'Seed', 6);
    id = '';
    try
        evalc('multigaussian_MLE(randn(rs, 3, 1), 1)');
    catch eS
        id = eS.identifier;
    end
    assert(strcmp(id, 'multigaussian_MLE:TooFewObservations'), 'n = 3 did not error');
    x6 = randn(rs, 6, 1);
    wsAll = warning;
    warning('error', 'multigaussian_MLE:TooFewObservations');
    wid = '';
    try
        evalc('multigaussian_MLE(x6, [1 2])');
    catch eW
        wid = eW.identifier;
    end
    warning(wsAll);
    assert(strcmp(wid, 'multigaussian_MLE:TooFewObservations'), 'no warning at n = 6');
    warning('off', 'multigaussian_MLE:TooFewObservations');
    [~, r] = evalc('multigaussian_MLE(x6, [1 2])');
    warning(wsAll);
    assert(isscalar(r) && r.NParameters == 2, 'K = 2 not dropped at n = 6');
    [~, r] = evalc('multigaussian_MLE(randn(rs, 7, 1), 2)');
    assert(isfinite(r.AICc), 'n = 7, K = 2 should be allowed');
    fprintf('[PASS] 6. Sample-size guard (n >= 3K + 1)\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 6. Sample-size guard: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 7: the global random stream is untouched, results reproducible
try
    rs = RandStream('twister', 'Seed', 7);
    x  = [randn(rs, 200, 1); 3 + randn(rs, 200, 1)];
    rng(123);  s0 = rng;
    [~, ra] = evalc('multigaussian_MLE(x, 3)');
    s1 = rng;
    rng(999);
    [~, rb] = evalc('multigaussian_MLE(x, 3)');
    assert(isequal(s0.State, s1.State), 'global stream advanced');
    assert(isequal(ra.Params, rb.Params), 'result depends on the global stream');
    fprintf('[PASS] 7. Global random stream untouched; results reproducible\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 7. Random stream: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 8: 'Select' by BIC
try
    rs = RandStream('twister', 'Seed', 8);
    xb = [randn(rs, 300, 1); 6 + randn(rs, 300, 1)];
    xu = randn(rs, 500, 1);
    [~, rb] = evalc('multigaussian_MLE(xb, [1 2 3], ''select'')');
    [~, ru] = evalc('multigaussian_MLE(xu, [1 2 3], ''select'')');
    assert(isscalar(rb) && size(rb.Params, 1) == 2, ...
        sprintf('bimodal: picked K = %d', size(rb.Params, 1)));
    assert(isscalar(ru) && size(ru.Params, 1) == 1, ...
        sprintf('unimodal: picked K = %d', size(ru.Params, 1)));
    fprintf('[PASS] 8. Select picks K = 2 (bimodal) and K = 1 (unimodal)\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 8. Select: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 9: input errors and NaN removal
try
    ids = {};
    bad = {{5*ones(20, 1), 1, ''}, {randn(20, 1), 1, 7}, ...
           {randn(20, 1), 0, ''}, {'abc', 1, ''}};
    for i = 1:numel(bad)
        try
            evalc('multigaussian_MLE(bad{i}{:})');
            ids{end+1} = 'none'; %#ok<SAGROW>
        catch eIn
            ids{end+1} = eIn.identifier; %#ok<SAGROW>
        end
    end
    want = {'multigaussian_MLE:ConstantData', 'multigaussian_MLE:InvalidMode', ...
            'multigaussian_MLE:InvalidModels', 'multigaussian_MLE:InvalidData'};
    assert(isequal(ids, want), strjoin(ids, ', '));
    ws = warning('off', 'multigaussian_MLE:NonfiniteDataRemoved');
    [~, r] = evalc('multigaussian_MLE([randn(30, 1); NaN; Inf], 1)');
    warning(ws);
    assert(r.N == 30, sprintf('N = %d after removal', r.N));
    fprintf('[PASS] 9. Input errors and NaN removal\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 9. Inputs: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 10: K = 2 on three populations
% Old code: an independent search beat it by 84.7 log-likelihood units.
try
    worst = -Inf;
    for seed = 301:303
        rs = RandStream('twister', 'Seed', seed);
        x  = [randn(rs, 200, 1); 4 + 0.5*randn(rs, 150, 1); -5 + 2*randn(rs, 150, 1)];
        [~, r] = evalc('multigaussian_MLE(x, 2)');
        worst = max(worst, indep_gmm(x, 2, 30, rs) - r.LogLike);
    end
    assert(worst < 1e-6, sprintf('independent search higher by %.3g', worst));
    fprintf('[PASS] 10. K = 2 on three populations reaches the optimum (%.2g)\n', worst);
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 10. K = 2 optimum: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 11: K = 3 and K = 4 reach the independent non-degenerate optimum
try
    msg = '';
    rs = RandStream('twister', 'Seed', 1101);
    x  = [randn(rs, 150, 1); 6 + 0.6*randn(rs, 120, 1); ...
          -6 + randn(rs, 130, 1); 12 + 0.8*randn(rs, 100, 1)];
    for K = [3 4]
        [~, r] = evalc('multigaussian_MLE(x, K)');
        d = indep_gmm(x, K, 30, rs) - r.LogLike;
        if d > 1e-6
            msg = [msg sprintf('K = %d short by %.3g; ', K, d)]; %#ok<AGROW>
        end
    end
    assert(isempty(msg), msg);
    fprintf('[PASS] 11. K = 3 and K = 4 reach the independent optimum\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 11. K = 3/4 optimum: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 12: no collapsed component when a non-degenerate fit exists
try
    rs = RandStream('twister', 'Seed', 1201);
    x  = randn(rs, 500, 1) ./ sqrt(sum(randn(rs, 500, 3).^2, 2)/3);   % t(3)
    ws = warning('off', 'multigaussian_MLE:SigmaAtFloor');
    [~, r] = evalc('multigaussian_MLE(x, 4)');
    warning(ws);
    assert(~any(r.SigmaAtFloor), 'a collapsed component was returned');
    d = indep_gmm(x, 4, 30, rs) - r.LogLike;
    assert(d < 1e-6, sprintf('independent non-degenerate search higher by %.3g', d));
    fprintf('[PASS] 12. No collapsed component; non-degenerate optimum reached\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 12. Degeneracy: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 13: standard errors
try
    rs = RandStream('twister', 'Seed', 13);
    x1 = 3 + 2*randn(rs, 400, 1);
    [~, r1] = evalc('multigaussian_MLE(x1, 1)');
    s = std(x1, 1);  n = 400;
    assert(abs(r1.SE(1) - s/sqrt(n)) < 1e-6*s && abs(r1.SE(2) - s/sqrt(2*n)) < 1e-6*s, ...
        sprintf('K = 1 SE %s vs analytic %s', mat2str(r1.SE(1:2), 6), ...
        mat2str([s/sqrt(n) s/sqrt(2*n)], 6)));
    assert(isnan(r1.SE(3)), 'K = 1 weight SE should be NaN');
    x2 = [randn(rs, 300, 1); 4 + 0.6*randn(rs, 200, 1)];
    [~, r2] = evalc('multigaussian_MLE(x2, 2)');
    P = r2.Params;
    LL2 = @(q) sum(log(mix_pdf(x2, [q(1) q(3) q(5); q(2) q(4) 1 - q(5)])));
    q0 = [P(1,1) P(2,1) P(1,2) P(2,2) P(1,3)];
    seN = sqrt(diag(inv(-num_hessian(LL2, q0))))';
    seF = [r2.SE(:,1)' r2.SE(:,2)' r2.SE(1,3)];
    relErr = max(abs(seN - seF) ./ seN);
    assert(relErr < 1e-4, sprintf('K = 2 max relative SE error %.3g', relErr));
    fprintf('[PASS] 13. SEs: K = 1 analytic, K = 2 independent Hessian (%.1g)\n', relErr);
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 13. SE: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 14: K = 1 profile CI ends on the closed-form threshold
try
    [~, r] = evalc('multigaussian_MLE(x1, 1, ''profile'')');
    thr = r.LogLike - chiHalf;
    n = numel(x1);
    pMu  = @(m) -n/2*(log(2*pi*mean((x1 - m).^2)) + 1);
    pSig = @(sg) -n*log(sg) - n/2*log(2*pi) - n*var(x1, 1)/(2*sg^2);
    d = [pMu(r.CILower(1)), pMu(r.CIUpper(1)), pSig(r.CILower(2)), ...
         pSig(r.CIUpper(2))] - thr;
    assert(all(abs(d) < 1e-6) && ~any(r.CIOpen(:)), ...
        sprintf('profile - threshold: %s', mat2str(d, 3)));
    assert(startsWith(r.CIMethod, 'profile'), 'CIMethod');
    fprintf('[PASS] 14. K = 1 profile CI ends sit on the threshold\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 14. Profile CI: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 15: intervals contain the estimate and respect constraints
try
    [~, r2p] = evalc('multigaussian_MLE(x2, 2, ''profile'')');
    ok = true;
    for r = [r1, r2, r2p]
        P = r.Params;  K = size(P, 1);
        cols = 1:2;  if K > 1, cols = 1:3; end
        ok = ok && all(all(r.CILower(:, cols) < P(:, cols) & P(:, cols) < r.CIUpper(:, cols)));
        ok = ok && all(all(r.CIWaldLower(:, cols) < P(:, cols) & P(:, cols) < r.CIWaldUpper(:, cols)));
        ok = ok && all(r.CILower(:, 2) > 0) && all(r.CIWaldLower(:, 2) > 0);
        if K > 1
            ok = ok && all(r.CIWaldLower(:, 3) > 0 & r.CIWaldUpper(:, 3) < 1);
        end
        ok = ok && isequal(size(r.Cov), [3*K 3*K]) && r.CILevel == 0.95;
    end
    assert(ok, 'an interval misses the estimate or breaks a constraint');
    assert(strcmp(r2.CIMethod, 'Wald'), 'default CIMethod should be Wald');
    fprintf('[PASS] 15. Intervals contain the estimate and respect constraints\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 15. Interval sanity: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 16: SEs and CIs equivariant under c*x + b
try
    c = -1e3;  b = 7;
    [~, s2] = evalc('multigaussian_MLE(c*x2 + b, 2, ''profile'')');
    % c < 0 mirrors the data: component order reverses, mu intervals swap ends.
    eSE = max(max(abs(flipud(s2.SE) ./ [abs(c) abs(c) 1] - r2p.SE) ./ r2p.SE));
    muLo = flipud((c*r2p.CIUpper(:, 1) + b));
    muHi = flipud((c*r2p.CILower(:, 1) + b));
    eMu = max(abs([s2.CILower(:, 1) - muLo; s2.CIUpper(:, 1) - muHi]) ./ ...
        (abs(c)*[r2p.SE(:, 1); r2p.SE(:, 1)]));
    eSg = max(max(abs([flipud(s2.CILower(:, 2)), flipud(s2.CIUpper(:, 2))]/abs(c) - ...
        [r2p.CILower(:, 2), r2p.CIUpper(:, 2)]) ./ [r2p.CILower(:, 2), r2p.CIUpper(:, 2)]));
    assert(max([eSE eMu eSg]) < 1e-4, sprintf('SE %.3g, mu CI %.3g, sigma CI %.3g', ...
        eSE, eMu, eSg));
    fprintf('[PASS] 16. SEs and CIs equivariant under -1e3*x + 7\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 16. SE/CI equivariance: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Summary
fprintf('\n----------------------------------------\n');
fprintf('Results: %d passed, %d failed (of %d)\n', nPassed, nFailed, ...
    nPassed + nFailed);
if nFailed == 0
    fprintf('All tests passed.\n');
else
    fprintf('%d test(s) FAILED.\n', nFailed);
end
fprintf('----------------------------------------\n\n');

%% Local functions (independent of multigaussian_MLE)
function f = mix_pdf(t, P)
f = zeros(size(t));
for j = 1:size(P, 1)
    f = f + P(j, 3) / (P(j, 2)*sqrt(2*pi)) * exp(-0.5*((t - P(j, 1))/P(j, 2)).^2);
end
end

function best = indep_gmm(x, k, nStarts, rs)
% Nelder-Mead from random starts, restricted to non-degenerate solutions
% (every sigma above 1.01x the 1e-3*IQR floor).
xs = sort(x);
qq = interp1(linspace(0, 1, numel(xs)), xs, [0.25 0.75]);
sMin = 1e-3 * diff(qq);
best = -Inf;
o = optimset('MaxFunEvals', 4e4, 'MaxIter', 4e4, 'TolX', 1e-10, ...
    'TolFun', 1e-10, 'Display', 'off');
for t = 1:nStarts
    p0 = [sort(x(randi(rs, numel(x), 1, k)))', ...
          log(std(x)*(0.2 + rand(rs, 1, k))), randn(rs, 1, k - 1)];
    f = @(p) nd_nll(x, p, k, sMin);
    p = fminsearch(f, p0, o);
    p = fminsearch(f, p, o);
    [~, fv] = fminsearch(f, p, o);
    best = max(best, -fv);
end
end

function v = nd_nll(x, p, k, sMin)
% Non-spurious solutions only: every sigma above 1.01x the floor and every
% component carrying at least 5 effective observations.
mu = p(1:k);  s = exp(p(k+1:2*k));
if any(s <= 1.01*sMin)
    v = 1e300;  return
end
u = [p(2*k+1:end) 0];  w = exp(u - max(u));  w = w / sum(w);
R = zeros(numel(x), k);
for j = 1:k
    R(:, j) = w(j)/(s(j)*sqrt(2*pi)) * exp(-0.5*((x - mu(j))/s(j)).^2);
end
L = sum(R, 2);
Nk = sum(R ./ max(L, realmin), 1);
v = -sum(log(L));
if ~isfinite(v) || any(Nk < 5), v = 1e300; end
end

function H = num_hessian(f, q)
nq = numel(q);  h = 1e-5 * max(abs(q), 1);  H = zeros(nq);
for i = 1:nq
    for j = 1:nq
        ei = zeros(1, nq);  ej = ei;  ei(i) = h(i);  ej(j) = h(j);
        H(i, j) = (f(q+ei+ej) - f(q+ei-ej) - f(q-ei+ej) + f(q-ei-ej)) / (4*h(i)*h(j));
    end
end
H = (H + H') / 2;
end
