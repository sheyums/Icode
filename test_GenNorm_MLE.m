%% test_GenNorm_MLE.m
%
% Tests for GenNorm_MLE. Run from the folder holding GenNorm_MLE.m; each
% test prints [PASS] or [FAIL], and a summary follows. Ground truth comes
% from code written here from the definition, sharing nothing with the
% fitter: a closed-form K = 1 profile likelihood and a direct mixture
% density. Guards marked (audit) failed on the pre-Nov-2026 code.
%
%  1. K = 1 optimum matches an independent profile (alpha 0.6..4)  (audit)
%  2. LogLike equals the independent density summed at Params
%  3. Fit integrates to 1; Fit/FitX shapes                          (audit)
%  4. Equivariance under c*x + b for K = 1 and K = 2                 (audit)
%  5. AIC is plain AIC; AICc and BIC formulas                       (audit)
%  6. Sample-size guard: n = 4 errors, K too large is dropped       (audit)
%  7. AlphaAtBound flags uniform data, not Normal data              (audit)
%  8. The caller's global random stream is not touched              (audit)
%  9. The caller's warning state is restored                        (audit)
% 10. 'Select' picks K = 2 on bimodal data and K = 1 on unimodal
% 11. Input errors and NaN removal
% 12. K = 2 optimum vs an independent 30-start search               (audit)
% 13-19 Standard errors and confidence intervals (see that section)

nPassed = 0;  nFailed = 0;
fprintf('\nRunning tests for GenNorm_MLE...\n\n');

% Independent pieces
gg  = @(rs, n, a) sign(rand(rs, n, 1) - 0.5) .* gamrnd_rs(rs, 1/a, n).^(1/a);
pl  = @(x, m, a) numel(x)*(log(a) - log(2) - gammaln(1/a) - 1/a) ...
                 - numel(x)/a*log(a/numel(x)*sum(abs(x - m).^a));
mixpdf = @(t, P) mix_pdf(t, P);

%% Test 1: K = 1 optimum vs independent profile over (mu, alpha)
try
    worst = -Inf;  where = '';
    for a0 = [0.6 1 2 4]
        for seed = 1:3
            rs = RandStream('twister', 'Seed', 100*a0 + seed);
            x  = 2 + 3*gg(rs, 300, a0);
            [~, r] = evalc('GenNorm_MLE(x, 1)');
            best = indep_k1(x, pl);
            if best - r.LogLike > worst
                worst = best - r.LogLike;
                where = sprintf('alpha0 %g seed %d', a0, seed);
            end
        end
    end
    assert(worst < 1e-6, sprintf('independent optimum higher by %.3g (%s)', ...
        worst, where));
    fprintf('[PASS] 1. K = 1 optimum matches independent profile (worst %.2g)\n', worst);
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 1. K = 1 optimum: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 2: LogLike equals the independent density at Params
try
    rs = RandStream('twister', 'Seed', 2);
    x  = [randn(rs, 300, 1); 5 + 0.7*gg(rs, 200, 1.3)];
    [~, r] = evalc('GenNorm_MLE(x, [1 2])');
    for i = 1:2
        d = abs(sum(log(mixpdf(x, r(i).Params))) - r(i).LogLike);
        assert(d < 1e-8 * abs(r(i).LogLike), ...
            sprintf('K = %d: |indep - LogLike| = %.3g', i, d));
        assert(abs(sum(r(i).Params(:, 4)) - 1) < 1e-12, 'weights do not sum to 1');
    end
    fprintf('[PASS] 2. LogLike equals the independent density sum\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 2. LogLike consistency: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 3: Fit integrates to 1; Fit / FitX shapes
try
    [~, r] = evalc('GenNorm_MLE(x, 2, ''fit'')');
    Zint = integral(@(t) mixpdf(t, r.Params), -Inf, Inf, ...
        'AbsTol', 1e-13, 'RelTol', 1e-11);
    assert(abs(Zint - 1) < 1e-8, sprintf('integral %.12g', Zint));
    assert(numel(r.Fit) == numel(x) && numel(r.FitX) == numel(x), 'sizes');
    assert(issorted(r.FitX), 'FitX not sorted');
    assert(max(abs(r.Fit - mixpdf(r.FitX, r.Params)) ./ r.Fit) < 1e-10, ...
        'Fit disagrees with the independent density');
    fprintf('[PASS] 3. Fit integrates to 1 and matches the density\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 3. Fit: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 4: equivariance under c*x + b
% Old code: 1e-6*x changed the K = 2 LogLike by 10.1 beyond n*log(c).
try
    rs = RandStream('twister', 'Seed', 3);
    x1 = 2 + 3*gg(rs, 400, 1.5);
    x2 = [randn(rs, 300, 1); 6 + 0.5*randn(rs, 200, 1)];
    [~, r1] = evalc('GenNorm_MLE(x1, 1)');
    [~, r2] = evalc('GenNorm_MLE(x2, 2)');
    msg = '';
    for cb = [1e-6 0; 1e6 0; 1 1e8; -1e3 7]'
        c = cb(1);  b = cb(2);
        [~, s1] = evalc('GenNorm_MLE(c*x1 + b, 1)');
        [~, s2] = evalc('GenNorm_MLE(c*x2 + b, 2)');
        e1 = abs(s1.LogLike + 400*log(abs(c)) - r1.LogLike);
        e2 = abs(s2.LogLike + 500*log(abs(c)) - r2.LogLike);
        ea = abs(s1.Params(3) - r1.Params(3));
        if e1 > 1e-6 || e2 > 1e-6 || ea > 1e-5
            msg = [msg sprintf('c=%g b=%g: dLL %.3g / %.3g, dalpha %.3g; ', ...
                c, b, e1, e2, ea)]; %#ok<AGROW>
        end
    end
    assert(isempty(msg), msg);
    fprintf('[PASS] 4. Equivariant under c*x + b (K = 1, 2)\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 4. Equivariance: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 5: AIC / AICc / BIC
try
    rs = RandStream('twister', 'Seed', 4);
    x  = randn(rs, 60, 1);
    [~, r] = evalc('GenNorm_MLE(x, [1 2])');
    for i = 1:2
        p = 4*i - 1;  LL = r(i).LogLike;
        assert(abs(r(i).AIC - (2*p - 2*LL)) < 1e-9, 'AIC is not 2p - 2LL');
        assert(abs(r(i).AICc - (2*p - 2*LL + 2*p*(p+1)/(60 - p - 1))) < 1e-9, ...
            'AICc formula');
        assert(abs(r(i).BIC - (p*log(60) - 2*LL)) < 1e-9, 'BIC formula');
        assert(r(i).NParameters == p && r(i).N == 60, 'N / NParameters');
    end
    fprintf('[PASS] 5. AIC, AICc and BIC formulas\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 5. Information criteria: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 6: sample-size guard
% Old code: n = 2..7 at K = 2 returned finite AICc of -25 to -97; n = 8, Inf.
try
    rs = RandStream('twister', 'Seed', 6);
    id = '';
    try
        evalc('GenNorm_MLE(randn(rs, 4, 1), 1)');
    catch eSmall
        id = eSmall.identifier;
    end
    assert(strcmp(id, 'GenNorm_MLE:TooFewObservations'), 'n = 4 did not error');
    % The warning is promoted to an error and caught by identifier (lastwarn
    % is unreliable: functions called during the fit can overwrite it).
    x8 = randn(rs, 8, 1);
    wsAll = warning;
    warning('error', 'GenNorm_MLE:TooFewObservations');
    wid = '';
    try
        evalc('GenNorm_MLE(x8, [1 2])');
    catch eW
        wid = eW.identifier;
    end
    warning(wsAll);
    assert(strcmp(wid, 'GenNorm_MLE:TooFewObservations'), 'no warning at n = 8');
    warning('off', 'GenNorm_MLE:TooFewObservations');
    [~, r] = evalc('GenNorm_MLE(x8, [1 2])');
    warning(wsAll);
    assert(isscalar(r) && r.NParameters == 3, 'K = 2 not dropped at n = 8');
    [~, r] = evalc('GenNorm_MLE(randn(rs, 9, 1), 2)');
    assert(isfinite(r.AICc), 'n = 9, K = 2 should be allowed with finite AICc');
    fprintf('[PASS] 6. Sample-size guard (n >= 4K + 1)\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 6. Sample-size guard: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 7: AlphaAtBound
try
    rs = RandStream('twister', 'Seed', 7);
    ws = warning('off', 'GenNorm_MLE:AlphaAtBound');
    [~, ru] = evalc('GenNorm_MLE(rand(rs, 500, 1), 1)');
    [~, rn] = evalc('GenNorm_MLE(randn(rs, 500, 1), 1)');
    warning(ws);
    assert(isfield(ru, 'AlphaAtBound') && ru.AlphaAtBound, ...
        sprintf('uniform: alpha %.4g not flagged', ru.Params(3)));
    assert(~rn.AlphaAtBound, sprintf('Normal: alpha %.4g flagged', rn.Params(3)));
    fprintf('[PASS] 7. AlphaAtBound flags uniform, not Normal, data\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 7. AlphaAtBound: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 8: the global random stream is untouched, results reproducible
try
    rs = RandStream('twister', 'Seed', 8);
    x  = [randn(rs, 200, 1); 4 + randn(rs, 200, 1)];
    rng(123);  s0 = rng;
    [~, ra] = evalc('GenNorm_MLE(x, 2)');
    s1 = rng;
    rng(999);
    [~, rb] = evalc('GenNorm_MLE(x, 2)');
    assert(isequal(s0.State, s1.State), 'global stream advanced');
    assert(isequal(ra.Params, rb.Params), ...
        'result depends on the global random state');
    fprintf('[PASS] 8. Global random stream untouched; results reproducible\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 8. Random stream: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 9: warning state restored
% Old code forced MATLAB:nearlySingularMatrix back 'on'.
try
    ws = warning('off', 'MATLAB:nearlySingularMatrix');
    evalc('GenNorm_MLE(x, 1)');
    st = warning('query', 'MATLAB:nearlySingularMatrix');
    warning(ws);
    assert(strcmp(st.state, 'off'), 'state changed to on');
    fprintf('[PASS] 9. Warning state restored\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 9. Warning state: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 10: 'Select'
try
    rs = RandStream('twister', 'Seed', 10);
    xb = [randn(rs, 300, 1); 6 + randn(rs, 300, 1)];
    xu = randn(rs, 500, 1);
    [~, rb] = evalc('GenNorm_MLE(xb, [1 2 3], ''select'')');
    [~, ru] = evalc('GenNorm_MLE(xu, [1 2], ''select'')');
    assert(isscalar(rb) && size(rb.Params, 1) == 2, ...
        sprintf('bimodal: picked K = %d', size(rb.Params, 1)));
    assert(isscalar(ru) && size(ru.Params, 1) == 1, ...
        sprintf('unimodal: picked K = %d', size(ru.Params, 1)));
    fprintf('[PASS] 10. Select picks K = 2 (bimodal) and K = 1 (unimodal)\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 10. Select: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 11: input errors and NaN removal
try
    ids = {};
    bad = {{5*ones(20, 1), 1, ''}, {randn(20, 1), 1, 7}, ...
           {randn(20, 1), 0, ''}, {'abc', 1, ''}};
    for i = 1:numel(bad)
        try
            evalc('GenNorm_MLE(bad{i}{:})');
            ids{end+1} = 'none'; %#ok<SAGROW>
        catch eIn
            ids{end+1} = eIn.identifier; %#ok<SAGROW>
        end
    end
    want = {'GenNorm_MLE:ConstantData', 'GenNorm_MLE:InvalidMode', ...
            'GenNorm_MLE:InvalidModels', 'GenNorm_MLE:InvalidData'};
    assert(isequal(ids, want), strjoin(ids, ', '));
    ws = warning('off', 'GenNorm_MLE:NonfiniteDataRemoved');
    [~, r] = evalc('GenNorm_MLE([randn(30, 1); NaN; Inf], 1)');
    warning(ws);
    assert(r.N == 30, sprintf('N = %d after removal', r.N));
    fprintf('[PASS] 11. Input errors and NaN removal\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 11. Inputs: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 12: K = 2 optimum vs an independent 30-start search
% Three populations fitted with K = 2, so the likelihood has competing
% optima. Old code: an independent 40-start search beat it by 11.8 (seed 51)
% and 5.97 (seed 53) log-likelihood units.
try
    worst = -Inf;
    for seed = [51 53]
        rs = RandStream('twister', 'Seed', seed);
        x  = [randn(rs, 250, 1); 3 + 0.4*gg(rs, 150, 1); -4 + 2*gg(rs, 100, 4)];
        [~, r] = evalc('GenNorm_MLE(x, 2)');
        best = indep_mix(x, 2, 30, rs);
        worst = max(worst, best - r.LogLike);
    end
    assert(worst < 1e-4, sprintf('independent search higher by %.3g', worst));
    fprintf('[PASS] 12. K = 2 optimum vs independent search (worst %.2g)\n', worst);
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 12. K = 2 optimum: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Standard errors and confidence intervals (added Nov 2026)
%
% 13. K = 1: SEHessian equals an independent numerical Hessian of the
%     log-likelihood in (mu, sigma, alpha)
% 14. K = 2: SEHessian equals an independent numerical Hessian in
%     (mu1, mu2, sigma1, sigma2, alpha1, alpha2, w1), including the weights
% 15. K = 1 profile CI ends sit on the threshold of an independent profile
% 16. Intervals contain the estimate and respect sigma, alpha > 0, w in (0,1)
% 17. SEs and CIs are equivariant under c*x + b
% 18. A parameter on its bound gets NaN SE and NaN Wald interval

chiHalf = 0.5 * (sqrt(2)*erfinv(0.95))^2;
rsU = RandStream('twister', 'Seed', 1301);
xU1 = 2 + 3*gg(rsU, 400, 2.5);
xU2 = [randn(rsU, 300, 1); 6 + 0.5*gg(rsU, 200, 1.8)];
[~, rU1] = evalc('GenNorm_MLE(xU1, 1, ''profile'')');
[~, rU2] = evalc('GenNorm_MLE(xU2, 2)');

%% Test 13: K = 1 SEHessian vs an independent Hessian
try
    LL1 = @(q) sum(log(q(3)/(2*q(2)*gamma(1/q(3)))) - (abs(xU1 - q(1))/q(2)).^q(3));
    q0 = rU1.Params(1:3);
    seN = sqrt(diag(inv(-num_hessian(LL1, q0))))';
    relErr = max(abs(seN - rU1.SEHessian(1:3)) ./ seN);
    assert(relErr < 1e-3, sprintf('max relative SE error %.3g', relErr));
    assert(isnan(rU1.SE(4)), 'K = 1 weight SE should be NaN (weight fixed at 1)');
    fprintf('[PASS] 13. K = 1 SE matches an independent Hessian (%.1g)\n', relErr);
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 13. K = 1 SE: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 14: K = 2 SEHessian vs an independent Hessian (weights included)
try
    P = rU2.Params;
    LL2 = @(q) sum(log(mix_pdf(xU2, [q(1) q(3) q(5) q(7); q(2) q(4) q(6) 1 - q(7)])));
    q0 = [P(1,1) P(2,1) P(1,2) P(2,2) P(1,3) P(2,3) P(1,4)];
    C  = inv(-num_hessian(LL2, q0));
    seN = sqrt(diag(C))';
    seF = [rU2.SEHessian(:,1)' rU2.SEHessian(:,2)' rU2.SEHessian(:,3)' rU2.SEHessian(1,4)];
    relErr = max(abs(seN - seF) ./ seN);
    assert(relErr < 1e-3, sprintf('max relative SE error %.3g', relErr));
    assert(abs(rU2.SE(1,4) - rU2.SE(2,4)) < 1e-12, 'the two weight SEs differ');
    fprintf('[PASS] 14. K = 2 SE matches an independent Hessian (%.1g)\n', relErr);
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 14. K = 2 SE: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 15: K = 1 profile CI ends on the threshold
% Independent profiles: sigma closed form given (mu, alpha); with sigma
% fixed, Nelder-Mead over (mu, log alpha).
try
    thr = rU1.LogLike - chiHalf;
    n1 = numel(xU1);
    P = rU1.Params;
    pMu = @(m) -fminbnd_val(@(a) -pl(xU1, m, a), 0.5, 10);
    pAl = @(a) -fminsearch_val(@(m) -pl(xU1, m, a), P(1));
    llS = @(s, m, a) n1*(log(a) - log(2) - log(s) - gammaln(1/a)) ...
        - sum(abs(xU1 - m).^a)/s^a;
    pSg = @(s) -min(arrayfun(@(t) fminsearch_val(@(q) -llS(s, q(1), ...
        min(max(exp(q(2)), 0.5), 10)), [P(1) log(t)]), [1.5 P(3) 5]));
    d = [pMu(rU1.CILower(1)), pMu(rU1.CIUpper(1)), ...
         pSg(rU1.CILower(2)), pSg(rU1.CIUpper(2)), ...
         pAl(rU1.CILower(3)), pAl(rU1.CIUpper(3))] - thr;
    assert(all(abs(d) < 1e-4) && ~any(rU1.CIOpen(:)), ...
        sprintf('profile - threshold: %s', mat2str(d, 3)));
    assert(startsWith(rU1.CIMethod, 'profile'), 'CIMethod not profile');
    fprintf('[PASS] 15. K = 1 profile CI ends sit on the threshold\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 15. Profile CI ends: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 16: intervals contain the estimate and respect constraints
try
    ok = true;
    for r = [rU1, rU2]
        P = r.Params;  K = size(P, 1);
        cols = 1:3;  if K > 1, cols = 1:4; end
        ok = ok && all(all(r.CILower(:, cols) < P(:, cols) & P(:, cols) < r.CIUpper(:, cols)));
        ok = ok && all(all(r.CIWaldLower(:, cols) < P(:, cols) & P(:, cols) < r.CIWaldUpper(:, cols)));
        ok = ok && all(all(r.CIWaldLower(:, 2:3) > 0)) && all(all(r.CILower(:, 2:3) > 0));
        if K > 1
            ok = ok && all(r.CIWaldLower(:, 4) > 0 & r.CIWaldUpper(:, 4) < 1);
        end
        ok = ok && isequal(size(r.Cov), [4*K 4*K]) && r.CILevel == 0.95;
    end
    assert(ok, 'an interval misses the estimate or breaks a constraint');
    assert(strcmp(rU2.CIMethod, 'Wald'), 'default CIMethod should be Wald');
    fprintf('[PASS] 16. Intervals contain the estimate and respect constraints\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 16. Interval sanity: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 17: SEs and CIs equivariant under c*x + b
try
    c = -1e3;  b = 7;
    [~, s1] = evalc('GenNorm_MLE(c*xU1 + b, 1, ''profile'')');
    eSE = max(abs(s1.SE(1:3) ./ [abs(c) abs(c) 1] - rU1.SE(1:3)) ./ rU1.SE(1:3));
    % c < 0 mirrors the data: the mu interval maps with its ends swapped.
    eMu = max(abs(sort((c*[rU1.CILower(1) rU1.CIUpper(1)] + b)) - ...
        [s1.CILower(1) s1.CIUpper(1)]) / (abs(c)*rU1.SE(1)));
    eSA = max(abs([s1.CILower(2:3) s1.CIUpper(2:3)] ./ [abs(c) 1 abs(c) 1] - ...
        [rU1.CILower(2:3) rU1.CIUpper(2:3)]) ./ [rU1.CILower(2:3) rU1.CIUpper(2:3)]);
    assert(max([eSE eMu eSA]) < 1e-4, sprintf('SE %.3g, mu CI %.3g, sigma/alpha CI %.3g', ...
        eSE, eMu, eSA));
    fprintf('[PASS] 17. SEs and CIs equivariant under -1e3*x + 7\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 17. Equivariance of SE/CI: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 18: a parameter on its bound gets NaN SE and NaN Wald interval
try
    rs = RandStream('twister', 'Seed', 18);
    ws = warning('off', 'GenNorm_MLE:AlphaAtBound');
    [~, r] = evalc('GenNorm_MLE(rand(rs, 500, 1), 1)');
    warning(ws);
    assert(r.AlphaAtBound, 'uniform data did not put alpha on its bound');
    assert(isnan(r.SE(3)) && isnan(r.CIWaldLower(3)) && isnan(r.CIWaldUpper(3)), ...
        'alpha on its bound still has an SE or interval');
    assert(all(isfinite(r.SE(1:2))), 'mu/sigma SEs should still be reported');
    fprintf('[PASS] 18. Parameter on its bound: NaN SE and Wald interval\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 18. Bound handling: %s\n', e.message);  nFailed = nFailed + 1;
end

%% Test 19: cusp regime -- mu SE from OPG, CuspSE warning without 'Profile'
try
    rs = RandStream('twister', 'Seed', 19);
    xc = 2 + 3*gg(rs, 300, 0.8);
    wsAll = warning;
    warning('error', 'GenNorm_MLE:CuspSE');
    wid = '';
    try
        evalc('GenNorm_MLE(xc, 1)');
    catch eW
        wid = eW.identifier;
    end
    warning(wsAll);
    assert(strcmp(wid, 'GenNorm_MLE:CuspSE'), 'no CuspSE warning at alpha ~ 0.8');
    warning('off', 'GenNorm_MLE:CuspSE');
    [~, r] = evalc('GenNorm_MLE(xc, 1)');
    warning(wsAll);
    assert(r.Params(3) < 1.5 && r.SE(1) == r.SEOPG(1) && r.SE(2) == r.SEHessian(2), ...
        'mu SE should come from OPG and sigma SE from the Hessian');
    fprintf('[PASS] 19. Cusp regime: mu SE from OPG, CuspSE warning\n');
    nPassed = nPassed + 1;
catch e
    warning(wsAll);
    fprintf('[FAIL] 19. Cusp regime: %s\n', e.message);  nFailed = nFailed + 1;
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

%% Local functions (independent of GenNorm_MLE)
function best = indep_k1(x, pl)
% Dense grid over alpha in [0.5, 10] and mu over the data points and a grid,
% Nelder-Mead polish, then every data point as mu at the polished alpha
% (the exact optimum for alpha <= 1 lies on a data point).
as = exp(linspace(log(0.5), log(10), 120));
ms = unique([x(:)', linspace(min(x), max(x), 400)]);
best = -Inf;  bm = 0;  ba = 1;
for aa = as
    L = arrayfun(@(m) pl(x, m, aa), ms);
    [Lm, i] = max(L);
    if Lm > best, best = Lm;  bm = ms(i);  ba = aa;  end
end
f = @(q) -pl(x, q(1), min(max(exp(q(2)), 0.5), 10));
q = fminsearch(f, [bm log(ba)], optimset('TolX', 1e-12, 'TolFun', 1e-12, ...
    'MaxFunEvals', 1e4, 'MaxIter', 1e4));
best = max(best, -f(q));
aq = min(max(exp(q(2)), 0.5), 10);
[Lp, ip] = max(arrayfun(@(m) pl(x, m, aq), x));
best = max(best, Lp);
g = @(a) -pl(x, x(ip), a);
ab = fminbnd(g, 0.5, 10, optimset('TolX', 1e-12));
best = max(best, -g(ab));
end

function best = indep_mix(x, k, nStarts, rs)
% Nelder-Mead from random starts on the direct mixture likelihood, with the
% same box (sigma >= 1e-3*IQR, alpha in [0.5, 10]) imposed by clamping.
sMin = 1e-3*iqr(x);
best = -Inf;
for t = 1:nStarts
    p0 = [sort(x(randi(rs, numel(x), 1, k)))', ...
          log(std(x)*(0.2 + rand(rs, 1, k))), log(0.7 + 3*rand(rs, 1, k)), ...
          randn(rs, 1, k - 1)];
    nll = @(p) mix_nll(x, p, k, sMin);
    o = optimset('MaxFunEvals', 2e4, 'MaxIter', 2e4, 'TolX', 1e-10, ...
        'TolFun', 1e-10, 'Display', 'off');
    p = fminsearch(nll, p0, o);
    [~, fv] = fminsearch(nll, p, o);
    best = max(best, -fv);
end
end

function f = mix_pdf(t, P)
% Mixture density from the definition; P rows are [mu, sigma, alpha, w].
f = zeros(size(t));
for j = 1:size(P, 1)
    mu = P(j, 1);  s = P(j, 2);  a = P(j, 3);  w = P(j, 4);
    f = f + w * a / (2*s*gamma(1/a)) * exp(-(abs(t - mu) / s).^a);
end
end

function nll = mix_nll(x, p, k, sMin)
mu = p(1:k);  s = max(exp(p(k+1:2*k)), sMin);
a  = min(max(exp(p(2*k+1:3*k)), 0.5), 10);
v  = [p(3*k+1:end) 0];  w = exp(v - max(v));  w = w / sum(w);
L  = zeros(numel(x), 1);
for j = 1:k
    L = L + w(j)*a(j)/(2*s(j)*gamma(1/a(j))) * exp(-(abs(x - mu(j))/s(j)).^a(j));
end
nll = -sum(log(L));
if ~isfinite(nll), nll = 1e300; end
end

function x = gamrnd_rs(rs, shape, n)
% Marsaglia-Tsang gamma sampler on a given stream (shape < 1 by boosting).
boost = shape < 1;  s = shape + boost;  d = s - 1/3;  c = 1/sqrt(9*d);
x = zeros(n, 1);
for i = 1:n
    while true
        z = randn(rs);  v = (1 + c*z)^3;
        if v > 0 && log(rand(rs)) < 0.5*z^2 + d - d*v + d*log(v)
            x(i) = d*v;  break
        end
    end
end
if boost
    x = x .* rand(rs, n, 1).^(1/shape);
end
end

function H = num_hessian(f, q)
% Central-difference Hessian of a scalar function.
nq = numel(q);  h = 1e-5 * max(abs(q), 1);  H = zeros(nq);
for i = 1:nq
    for j = 1:nq
        ei = zeros(1, nq);  ej = ei;  ei(i) = h(i);  ej(j) = h(j);
        H(i, j) = (f(q+ei+ej) - f(q+ei-ej) - f(q-ei+ej) + f(q-ei-ej)) / (4*h(i)*h(j));
    end
end
H = (H + H') / 2;
end

function v = fminsearch_val(f, q0)
[~, v] = fminsearch(f, q0, optimset('TolX', 1e-12, 'TolFun', 1e-12, ...
    'MaxFunEvals', 1e4, 'MaxIter', 1e4));
end

function v = fminbnd_val(f, a, b)
[~, v] = fminbnd(f, a, b, optimset('TolX', 1e-12));
end
