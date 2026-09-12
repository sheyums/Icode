function test_FitExponentiatedWeibullMLE()
%TEST_FITEXPONENTIATEDWEIBULLMLE Smoke tests for FitExponentiatedWeibullMLE.
%
% Run from the directory containing FitExponentiatedWeibullMLE.m:
%
%     test_FitExponentiatedWeibullMLE
%
% Each test prints [PASS] or [FAIL]; a summary line follows.
%
% Every test that is not mode-specific by nature runs in BOTH continuous
% and discrete mode, and reports which mode failed. The mode-specific ones
% are marked below.
%
% Tests covered
% -------------
%  1.  Law recovery                    - both modes: logL >= truth's
%                                        and survival curve matches
%                                        (parameters are only weakly
%                                        identified; see the test)
%  2.  Required fields present         - both modes
%  3.  Likelihood reimplemented        - both modes, vs independent
%                                        recomputation from the CDF
%  4.  Normalization                   - both modes: density integrates,
%                                        pmf sums, to 1
%  5.  PointwiseLogLik                 - both modes, sums to LogLik
%  6.  FixAlpha gives nested Weibull   - both modes
%  7.  Hazard shape, decreasing        - both modes
%  8.  Hazard shape, increasing        - both modes
%  9.  AIC / AICc formulas             - both modes x FixAlpha
% 10.  Alpha identifiability guard     - both modes, fires deep in tail
% 11.  No false-positive alpha guard   - both modes, xmin in the bulk
% 12.  Standard errors calibrated      - both modes, vs Monte Carlo SD
% 13.  Soft return and SE conventions  - both modes
% 14.  Mode convergence as dt -> 0     - discrete logL - n*log(dt) tends
%                                        to the continuous logL
% 15.  Unit invariance                 - discrete logL unit-free;
%                                        continuous logL shifts n*log(c)
% 16.  Comparable with hyperexp        - both modes, Vuong computable
% 17.  Discrete n_min clamp            - DISCRETE ONLY
% 18.  Discrete truncation consistency - DISCRETE ONLY
% 19.  Grid mismatch guard             - DISCRETE warns, CONTINUOUS silent
% 20.  Warning state not leaked        - both modes
% 21.  VERBOSE mode                    - both modes
% 22.  SamplingInterval is required   - omitting it errors, as does a
%                                       non-positive value
%
% Written as a function file so the local simulators at the bottom are in
% scope for every test.

fitfun = @FitExponentiatedWeibullMLE;

FAST = {'nStartsBase', 8, 'nStartsPerParameter', 6, 'maxStarts', 26, ...
        'Verbose', false};
MODES = {'continuous', 'discrete'};

nPassed = 0;
nFailed = 0;

fprintf('\nRunning tests for FitExponentiatedWeibullMLE...\n\n');

%% Test 1: Parameter recovery, both modes
try
    % The exponentiated Weibull's two shape parameters are strongly
    % correlated, so the LAW is identified far better than the individual
    % parameters: triples differing by ~10% describe survival curves that
    % agree to <0.01 everywhere, and the likelihood genuinely prefers a
    % point near but not at the truth. Recovery is therefore tested as
    % (a) the fit reaches at least the truth's log-likelihood, and (b) the
    % fitted truncated survival function matches the true one. A
    % per-parameter tolerance would be testing the identifiability of the
    % model, not the correctness of this code.
    lamT = 800; kT = 0.8; aT = 0.6; xmin = 100; bad = {}; dev = zeros(1,2);
    for m = 1:2
        rng(100+m);
        d = simEW(lamT, kT, aT, xmin, 3000, MODES{m}, 1);
        H = fitfun(d, xmin, 'DistributionType', MODES{m}, ...
            'SamplingInterval', 1, 'nStartsBase', 24, 'nStartsPerParameter', 16, ...
            'maxStarts', 72, 'Verbose', false, 'RandomSeed', 11);
        llTrue = recomputeLogLik(d, xmin, 1, lamT, kT, aT, m == 2);
        ds = sort(d);
        tg = linspace(xmin, ds(max(1, round(0.999*numel(ds)))), 2000)';
        dev(m) = max(abs(StruncEW(tg, lamT, kT, aT, xmin) ...
                       - StruncEW(tg, H.Lambda, H.K, H.Alpha, xmin)));
        gross = abs([H.Lambda H.K H.Alpha] - [lamT kT aT]) ./ [lamT kT aT];
        if H.LogLik < llTrue - 1e-6
            bad{end+1} = sprintf('%s: logL %.4f below truth %.4f', MODES{m}, ...
                H.LogLik, llTrue);
        elseif dev(m) > 0.03
            bad{end+1} = sprintf('%s: max|dS| = %.4f', MODES{m}, dev(m));
        elseif any(gross > 1.0)
            bad{end+1} = sprintf('%s: gross parameter error (%.0f %.3f %.3f)', ...
                MODES{m}, H.Lambda, H.K, H.Alpha);
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 1, ...
        sprintf('law recovered: logL >= truth, max|dS| = %.4f cont, %.4f disc', ...
        dev(1), dev(2)), strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 1, '', err.message);
end

%% Test 2: Required fields present, both modes
try
    need = {'Lambda','K','Alpha','LambdaSE','KSE','AlphaSE','HazardShape', ...
            'k','n','LogLik','PointwiseLogLik','SurvivalHandle','AIC','AICc','CovValid', ...
            'Success','Converged','ExitFlag','BestParamVector','FixAlpha', ...
            'DistributionType','xmin','Failed','FailureReason','Diagnostics'};
    needD = {'nNonFinite','nBelowXmin','nAtXmin','XminEqualsDataMin', ...
             'LooksGridded','GridSpacing','GridMismatch','nDistinctData', ...
             'nDistinctOnGrid','ExpMinusUmin','AlphaIdentifiable','Recommendation'};
    bad = {};
    for m = 1:2
        rng(110+m);
        d = simEW(800, 0.8, 0.6, 100, 400, MODES{m}, 1);
        H = fitfun(d, 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:});
        miss = [need(~isfield(H, need)), needD(~isfield(H.Diagnostics, needD))];
        if ~isempty(miss), bad{end+1} = [MODES{m} ': ' strjoin(miss, ',')]; end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 2, ...
        'all documented fields present in both modes', strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 2, '', err.message);
end

%% Test 3: LogLik matches an independent recomputation, both modes
try
    rng(103); xmin = 100; dt = 1; bad = {};
    d = round(simEW(800, 0.8, 0.6, xmin, 900, 'discrete', dt));
    for m = 1:2
        H = fitfun(d, xmin, 'DistributionType', MODES{m}, ...
            'SamplingInterval', dt, FAST{:});
        ll = recomputeLogLik(d, xmin, dt, H.Lambda, H.K, H.Alpha, m == 2);
        if abs(ll - H.LogLik) > 1e-6
            bad{end+1} = sprintf('%s off by %.3e', MODES{m}, ll - H.LogLik);
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 3, ...
        'LogLik matches independent recomputation in both modes', strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 3, '', err.message);
end

%% Test 4: normalization -- density integrates to 1, pmf sums to 1
try
    rng(104); bad = {};
    % discrete: sum the fitted pmf over the whole support
    xmin = 20; dt = 1;
    d = simEW(200, 0.9, 0.7, xmin, 600, 'discrete', dt);
    Hd = fitfun(d, xmin, 'DistributionType', 'discrete', ...
        'SamplingInterval', dt, FAST{:});
    nmin = max(1, round(xmin/dt));
    nn = (nmin:200000)';
    lFn = logEWCDF(nn*dt, Hd.Lambda, Hd.K, Hd.Alpha);
    lFp = logEWCDF((nn-1)*dt, Hd.Lambda, Hd.K, Hd.Alpha);
    totD = sum(exp(lFn) - exp(lFp)) / (1 - exp(lFp(1)));
    if abs(totD - 1) > 1e-6, bad{end+1} = sprintf('discrete sums to %.10f', totD); end
    % continuous: trapezoidal integral of the fitted density, checked
    % against the CDF increment over the same finite range
    xminC = 100;
    dc = simEW(800, 0.8, 0.6, xminC, 600, 'continuous', 1);
    Hc = fitfun(dc, xminC, 'SamplingInterval', 1, 'DistributionType', 'continuous', FAST{:});
    U = Hc.Lambda * 400;
    tg = linspace(xminC, U, 400001)';
    S = 1 - exp(logEWCDF(xminC, Hc.Lambda, Hc.K, Hc.Alpha));
    fg = ewPDF(tg, Hc.Lambda, Hc.K, Hc.Alpha) / S;
    totC = trapz(tg, fg);
    cdfC = (exp(logEWCDF(U, Hc.Lambda, Hc.K, Hc.Alpha)) ...
            - exp(logEWCDF(xminC, Hc.Lambda, Hc.K, Hc.Alpha))) / S;
    if abs(totC - cdfC) > 1e-5 || abs(totC - 1) > 1e-4
        bad{end+1} = sprintf('continuous integral %.8f vs CDF %.8f', totC, cdfC);
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 4, ...
        sprintf('pmf sums to %.8f, density integrates to %.6f', totD, totC), ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 4, '', err.message);
end

%% Test 5: PointwiseLogLik sums to LogLik, both modes
try
    bad = {};
    for m = 1:2
        rng(150+m);
        d = simEW(800, 0.8, 0.6, 100, 500, MODES{m}, 1);
        H = fitfun(d, 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:});
        if numel(H.PointwiseLogLik) ~= H.n || abs(sum(H.PointwiseLogLik) - H.LogLik) > 1e-8
            bad{end+1} = sprintf('%s: numel=%d n=%d resid=%.3e', MODES{m}, ...
                numel(H.PointwiseLogLik), H.n, sum(H.PointwiseLogLik) - H.LogLik);
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 5, ...
        'PointwiseLogLik is n x 1 and sums to LogLik in both modes', strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 5, '', err.message);
end

%% Test 6: FixAlpha fits the nested Weibull, both modes
try
    bad = {}; lrs = zeros(1,2);
    for m = 1:2
        rng(160+m);
        d = simEW(800, 0.8, 0.6, 100, 1500, MODES{m}, 1);
        Hew = fitfun(d, 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:}, 'RandomSeed', 5);
        Hwb = fitfun(d, 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, 'FixAlpha', true, ...
            FAST{:}, 'RandomSeed', 5);
        lrs(m) = 2*(Hew.LogLik - Hwb.LogLik);
        if ~(Hwb.Alpha == 1 && Hwb.k == 2 && Hew.k == 3 && Hwb.AlphaSE == 0 ...
             && Hwb.LogLik <= Hew.LogLik + 1e-6)
            bad{end+1} = sprintf('%s: alpha=%g k=%d/%d LR=%.3f', MODES{m}, ...
                Hwb.Alpha, Hwb.k, Hew.k, lrs(m));
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 6, ...
        sprintf('FixAlpha nests Weibull (LR = %.2f cont, %.2f disc, 1 df)', lrs(1), lrs(2)), ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 6, '', err.message);
end

%% Test 7: decreasing-hazard regime, both modes
try
    bad = {};
    for m = 1:2
        rng(170+m);
        d = simEW(400, 0.7, 0.5, 40, 2500, MODES{m}, 1);   % alpha*k = 0.35
        H = fitfun(d, 40, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:}, 'RandomSeed', 7);
        if ~strcmp(H.HazardShape, 'decreasing')
            bad{end+1} = sprintf('%s: "%s" (k=%.3f ak=%.3f)', MODES{m}, ...
                H.HazardShape, H.K, H.Alpha*H.K);
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 7, ...
        'decreasing hazard classified correctly in both modes', strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 7, '', err.message);
end

%% Test 8: increasing-hazard regime, both modes
try
    bad = {};
    for m = 1:2
        rng(180+m);
        d = simEW(400, 2.0, 1.5, 40, 2500, MODES{m}, 1);   % alpha*k = 3
        H = fitfun(d, 40, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:}, 'RandomSeed', 8);
        if ~strcmp(H.HazardShape, 'increasing')
            bad{end+1} = sprintf('%s: "%s" (k=%.3f ak=%.3f)', MODES{m}, ...
                H.HazardShape, H.K, H.Alpha*H.K);
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 8, ...
        'increasing hazard classified correctly in both modes', strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 8, '', err.message);
end

%% Test 9: AIC / AICc arithmetic, both modes x FixAlpha
try
    bad = {};
    for m = 1:2
        rng(190+m);
        d = simEW(800, 0.8, 0.6, 100, 400, MODES{m}, 1);
        for fa = [false true]
            H = fitfun(d, 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, 'FixAlpha', fa, FAST{:});
            p = 3 - double(fa);
            aic = 2*p - 2*H.LogLik;
            aicc = aic + 2*p*(p+1)/(H.n - p - 1);
            if H.k ~= p || abs(H.AIC - aic) > 1e-9 || abs(H.AICc - aicc) > 1e-9
                bad{end+1} = sprintf('%s FixAlpha=%d', MODES{m}, fa);
            end
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 9, ...
        'AIC/AICc arithmetic exact in both modes, both parameter counts', ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 9, '', err.message);
end

%% Test 10: alpha identifiability guard fires deep in the tail, both modes
try
    bad = {}; vals = zeros(1,2);
    for m = 1:2
        rng(200+m);
        d = simEW(100, 1.2, 0.8, 900, 700, MODES{m}, 1);
        prevWarn = warning('on', 'FitExponentiatedWeibullMLE:AlphaNotIdentified');
        lastwarn('');
        fprintf('--- begin expected warning (%s) ---\n', MODES{m});
        H = fitfun(d, 900, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:}, 'RandomSeed', 10);
        fprintf('--- end expected warning ---\n');
        [~, wid] = lastwarn();
        warning(prevWarn);
        vals(m) = H.Diagnostics.ExpMinusUmin;
        if ~(strcmp(wid, 'FitExponentiatedWeibullMLE:AlphaNotIdentified') ...
             && ~H.Diagnostics.AlphaIdentifiable && vals(m) < 1e-3)
            bad{end+1} = sprintf('%s: wid=%s flagged=%d val=%.3g', MODES{m}, ...
                wid, ~H.Diagnostics.AlphaIdentifiable, vals(m));
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 10, ...
        sprintf('alpha guard fired in both modes (exp(-u_min) = %.2g, %.2g)', ...
        vals(1), vals(2)), strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 10, '', err.message);
end

%% Test 11: no false-positive alpha guard when xmin is in the bulk
try
    bad = {}; vals = zeros(1,2);
    for m = 1:2
        rng(210+m);
        d = simEW(800, 0.8, 0.6, 100, 2000, MODES{m}, 1);
        H = fitfun(d, 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:}, 'RandomSeed', 11);
        vals(m) = H.Diagnostics.ExpMinusUmin;
        if ~(H.Diagnostics.AlphaIdentifiable && vals(m) > 1e-3)
            bad{end+1} = sprintf('%s: flagged=%d val=%.3g', MODES{m}, ...
                H.Diagnostics.AlphaIdentifiable, vals(m));
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 11, ...
        sprintf('no alpha warning with xmin in the bulk (%.3f, %.3f)', vals(1), vals(2)), ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 11, '', err.message);
end

%% Test 12: standard errors calibrated against Monte-Carlo spread
try
    bad = {}; msg = {};
    for m = 1:2
        rng(220+m);
        lamT = 600; kT = 0.9; aT = 0.7; xmin = 80; R = 15;
        est = []; ses = [];
        for r = 1:R
            d = simEW(lamT, kT, aT, xmin, 500, MODES{m}, 1);
            H = fitfun(d, xmin, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:}, ...
                'RandomSeed', 400+r);
            if H.CovValid && isfinite(H.KSE)
                est(end+1,:) = [H.Lambda H.K];
                ses(end+1,:) = [H.LambdaSE H.KSE];
            end
        end
        if size(est,1) < 8
            bad{end+1} = sprintf('%s: only %d usable replicates', MODES{m}, size(est,1));
            continue
        end
        ratio = median(ses,1) ./ std(est,0,1);
        msg{end+1} = sprintf('%s lambda %.2f k %.2f', MODES{m}, ratio(1), ratio(2));
        if any(ratio < 0.5) || any(ratio > 2.0)
            bad{end+1} = sprintf('%s: SE/MC-SD = %s', MODES{m}, mat2str(round(ratio*100)/100));
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 12, ...
        ['SE/MC-SD within 2x: ' strjoin(msg, ', ')], strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 12, '', err.message);
end

%% Test 13: soft return and SE conventions, both modes
try
    bad = {};
    for m = 1:2
        threw = false;
        try
            fitfun([110; 120; 130], 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:});
        catch err
            threw = strcmp(err.identifier, 'FitExponentiatedWeibullMLE:TooFewData');
        end
        H = fitfun([110; 120; 130], 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, ...
            'ErrorOnNoValidFit', false, FAST{:});
        okSoft = threw && H.Failed && isnan(H.Lambda) && ~H.Success ...
                 && ~isempty(strfind(H.FailureReason, 'TooFewData')) ...
                 && isnan(H.LambdaSE) && isnan(H.KSE);
        rng(230+m);
        d = simEW(800, 0.8, 0.6, 100, 400, MODES{m}, 1);
        Hf = fitfun(d, 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, 'FixAlpha', true, FAST{:});
        okSE = Hf.AlphaSE == 0 && (Hf.CovValid || isnan(Hf.LambdaSE));
        if ~(okSoft && okSE)
            bad{end+1} = sprintf('%s: soft=%d se=%d', MODES{m}, okSoft, okSE);
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 13, ...
        'TooFewData throws then soft-returns; SE conventions hold in both modes', ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 13, '', err.message);
end

%% Test 14: the two modes converge as dt -> 0
try
    rng(240);
    xmin = 100;
    d = simEW(800, 0.8, 0.6, xmin, 1500, 'continuous', 1);
    Hc = fitfun(d, xmin, 'SamplingInterval', 1, 'DistributionType', 'continuous', FAST{:}, 'RandomSeed', 3);
    % A pmf is a density times a bin width, so logL_disc - n*log(dt) should
    % approach logL_cont from below as dt shrinks.
    % These fits deliberately bin genuinely continuous data, so the
    % grid-spacing guard fires each time and is correct to: rounding really
    % does discard resolution here, and the message says so without
    % inventing an acquisition interval. Marked expected so it does not
    % read as a failure.
    gaps = [];
    fprintf('--- begin expected warnings (continuous data binned on purpose) ---\n');
    for dt = [1, 0.1, 0.01]
        Hd = fitfun(d, xmin, 'DistributionType', 'discrete', ...
            'SamplingInterval', dt, FAST{:}, 'RandomSeed', 3);
        gaps(end+1) = (Hd.LogLik - Hd.n*log(dt)) - Hc.LogLik;
    end
    fprintf('--- end expected warnings ---\n');
    shrinking = abs(gaps(3)) < abs(gaps(1)) && abs(gaps(3)) < 0.05*numel(d);
    [nPassed, nFailed] = report(shrinking, nPassed, nFailed, 14, ...
        sprintf('discrete -> continuous as dt->0 (gap %.2f, %.2f, %.3f nats)', ...
        gaps(1), gaps(2), gaps(3)), sprintf('gaps %s', mat2str(round(gaps*100)/100)));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 14, '', err.message);
end

%% Test 15: unit invariance in both modes (with the right invariant each)
try
    rng(250);
    dsec = ceil(simEW(800, 0.9, 0.7, 120, 1500, 'discrete', 1));
    c = 60; bad = {};
    % discrete: a pmf is unit-free, so logL is identical
    Hs = fitfun(dsec, 120, 'DistributionType', 'discrete', ...
        'SamplingInterval', 1, FAST{:}, 'RandomSeed', 3);
    Hm = fitfun(dsec/c, 120/c, 'DistributionType', 'discrete', ...
        'SamplingInterval', 1/c, FAST{:}, 'RandomSeed', 3);
    if abs(Hs.LogLik - Hm.LogLik) > 1e-6
        bad{end+1} = sprintf('discrete logL %.6f vs %.6f', Hs.LogLik, Hm.LogLik);
    end
    if abs(Hm.Lambda*c - Hs.Lambda)/Hs.Lambda > 1e-3 || abs(Hm.K - Hs.K) > 1e-4 ...
            || abs(Hm.Alpha - Hs.Alpha) > 1e-4
        bad{end+1} = 'discrete parameters not equivariant';
    end
    % continuous: a density per unit time scales, so logL shifts by n*log(c)
    Cs = fitfun(dsec, 120, 'SamplingInterval', 1, 'DistributionType', 'continuous', FAST{:}, 'RandomSeed', 3);
    Cm = fitfun(dsec/c, 120/c, 'SamplingInterval', 1, 'DistributionType', 'continuous', FAST{:}, 'RandomSeed', 3);
    if abs((Cm.LogLik - Cs.LogLik) - Cs.n*log(c)) > 1e-3
        bad{end+1} = sprintf('continuous shift %.4f vs n*log(c)=%.4f', ...
            Cm.LogLik - Cs.LogLik, Cs.n*log(c));
    end
    if abs(Cm.Lambda*c - Cs.Lambda)/Cs.Lambda > 1e-3
        bad{end+1} = 'continuous lambda not equivariant';
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 15, ...
        'discrete logL unit-free; continuous logL shifts by exactly n*log(c)', ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 15, '', err.message);
end

%% Test 16: comparable with FitHyperexponentialMLE in both modes
try
    bad = {}; vs = zeros(1,2);
    rng(260);
    d = 300 + ceil(-700*log(rand(1200,1)));     % single exponential, gridded
    for m = 1:2
        Hew = fitfun(d, 300, 'DistributionType', MODES{m}, ...
            'SamplingInterval', 1, FAST{:}, 'RandomSeed', 2);
        Hhx = FitHyperexponentialMLE(d, 300, 'MaxComponents', 1, ...
            'DistributionType', MODES{m}, 'SamplingInterval', 1, 'Verbose', false);
        D = Hhx.Selected.PointwiseLogLik - Hew.PointwiseLogLik;
        vs(m) = sqrt(numel(D))*mean(D)/std(D,1);
        ok = Hew.n == Hhx.n && numel(D) == Hew.n && isfinite(vs(m));
        % In discrete mode both are probabilities, so both logL must be <= 0.
        if m == 2, ok = ok && Hew.LogLik <= 0 && Hhx.Selected.LogLik <= 0; end
        if ~ok, bad{end+1} = sprintf('%s: V=%g', MODES{m}, vs(m)); end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 16, ...
        sprintf('same observations and measure; Vuong V = %.2f cont, %.2f disc', ...
        vs(1), vs(2)), strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 16, '', err.message);
end

%% Test 17: discrete n_min clamp (DISCRETE ONLY)
try
    rng(270);
    d = ceil(simEW(30, 1.0, 1.0, 1, 400, 'discrete', 1));
    H = fitfun(d, 0.4, 'DistributionType', 'discrete', 'SamplingInterval', 1, FAST{:});
    [nPassed, nFailed] = report(H.LogLik <= 0, nPassed, nFailed, 17, ...
        sprintf('xmin<dt/2 clamped to n_min=1, logL = %.2f <= 0', H.LogLik), ...
        sprintf('positive logL %.4f', H.LogLik));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 17, '', err.message);
end

%% Test 18: discrete truncation keeps values rounding into support (DISCRETE ONLY)
try
    rng(280);
    d = [repmat(1.6, 20, 1); ceil(simEW(40, 1.0, 1.0, 2, 300, 'discrete', 1))];
    H = fitfun(d, 2, 'DistributionType', 'discrete', 'SamplingInterval', 1, FAST{:});
    [nPassed, nFailed] = report(H.n == numel(d), nPassed, nFailed, 18, ...
        sprintf('values rounding into the support retained (n=%d)', H.n), ...
        sprintf('n=%d of %d', H.n, numel(d)));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 18, '', err.message);
end

%% Test 19: grid mismatch guard -- warns in discrete, silent in continuous
try
    rng(290);
    dsec = ceil(simEW(800, 0.8, 0.6, 100, 600, 'discrete', 1));
    prevWarn = warning('on', 'FitExponentiatedWeibullMLE:GridMismatch');
    lastwarn('');
    fprintf('--- begin expected warning ---\n');
    Hbad = fitfun(dsec/60, 100/60, 'SamplingInterval', 1, 'DistributionType', 'discrete', FAST{:});
    fprintf('--- end expected warning ---\n');
    [~, wid] = lastwarn();
    Hok = fitfun(dsec, 100, 'DistributionType', 'discrete', ...
        'SamplingInterval', 1, FAST{:});
    % continuous mode: the mismatch is flagged but must not warn, since no
    % rounding happens there
    lastwarn('');
    Hc = fitfun(dsec/60, 100/60, 'SamplingInterval', 1, 'DistributionType', 'continuous', FAST{:});
    [~, wid2] = lastwarn();
    warning(prevWarn);
    ok = strcmp(wid, 'FitExponentiatedWeibullMLE:GridMismatch') ...
         && Hbad.Diagnostics.GridMismatch && ~Hok.Diagnostics.GridMismatch ...
         && Hc.Diagnostics.GridMismatch ...
         && ~strcmp(wid2, 'FitExponentiatedWeibullMLE:GridMismatch');
    [nPassed, nFailed] = report(ok, nPassed, nFailed, 19, ...
        sprintf('discrete warns (%d -> %d distinct), continuous only flags', ...
        Hbad.Diagnostics.nDistinctData, Hbad.Diagnostics.nDistinctOnGrid), ...
        sprintf('wid=%s bad=%d ok=%d cont=%d wid2=%s', wid, ...
        Hbad.Diagnostics.GridMismatch, Hok.Diagnostics.GridMismatch, ...
        Hc.Diagnostics.GridMismatch, wid2));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 19, '', err.message);
end

%% Test 20: warning state not leaked, both modes
try
    warning('on', 'MATLAB:singularMatrix');
    before = warning('query', 'MATLAB:singularMatrix');
    bad = {};
    for m = 1:2
        rng(300+m);
        d = simEW(800, 0.8, 0.6, 100, 300, MODES{m}, 1);
        H = fitfun(d, 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, FAST{:}); %#ok<NASGU>
        after = warning('query', 'MATLAB:singularMatrix');
        if ~strcmp(before.state, after.state)
            bad{end+1} = sprintf('%s: %s -> %s', MODES{m}, before.state, after.state);
        end
    end
    [nPassed, nFailed] = report(isempty(bad), nPassed, nFailed, 20, ...
        sprintf('singularMatrix warning state preserved (%s)', before.state), ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 20, '', err.message);
end

%% Test 21: VERBOSE mode does not error, both modes
try
    for m = 1:2
        rng(310+m);
        d = simEW(800, 0.8, 0.6, 100, 400, MODES{m}, 1);
        fprintf('--- begin expected verbose output (%s) ---\n', MODES{m});
        H = fitfun(d, 100, 'SamplingInterval', 1, 'DistributionType', MODES{m}, 'nStartsBase', 8, ...
            'nStartsPerParameter', 6, 'maxStarts', 26, 'Verbose', true); %#ok<NASGU>
        fprintf('--- end expected verbose output ---\n');
    end
    [nPassed, nFailed] = report(true, nPassed, nFailed, 21, ...
        'verbose mode runs in both modes without error', '');
catch err
    [nPassed, nFailed] = report(false, nPassed, nFailed, 21, '', err.message);
end

%% Test 22: SamplingInterval is required
try
    rng(320);
    d = simEW(800, 0.8, 0.6, 100, 200, 'discrete', 1);
    missingId = '';
    try
        fitfun(d, 100, 'Verbose', false);
    catch err
        missingId = err.identifier;
    end
    badId = '';
    try
        fitfun(d, 100, 'SamplingInterval', -1, 'Verbose', false);
    catch err
        badId = err.identifier;
    end
    H = fitfun(d, 100, 'SamplingInterval', 1, FAST{:});
    okMissing = strcmp(missingId, 'FitExponentiatedWeibullMLE:SamplingIntervalRequired');
    okBad     = strcmp(badId, 'FitExponentiatedWeibullMLE:InvalidSamplingInterval');
    okFits    = H.Success && ~H.Failed && isfinite(H.LogLik);
    if okMissing && okBad && okFits
        fprintf('[PASS] Test 22: SamplingInterval required (omitted and non-positive both error)\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 22: missing="%s" nonpositive="%s" fitsWhenGiven=%d\n', ...
            missingId, badId, okFits);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 22: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 23: SurvivalHandle is consistent with the likelihood
try
    rng(231);
    dtS = 1; dS = simEW(800, 0.8, 0.6, 100, 400, 'discrete', dtS);
    HS = fitfun(dS, 100, 'DistributionType', 'discrete', 'SamplingInterval', dtS, FAST{:});
    nn = round(dS/dtS);
    p  = HS.SurvivalHandle((nn-1)*dtS) - HS.SurvivalHandle(nn*dtS);
    % A bin's probability must equal exp of that observation's pointwise
    % log-likelihood, or the goodness-of-fit expected counts and the plotted
    % curve would describe a different model from the one that was fitted.
    okBin  = max(abs(p - exp(HS.PointwiseLogLik))) < 1e-10;
    okOne  = abs(HS.SurvivalHandle((max(1,round(100/dtS))-1)*dtS) - 1) < 1e-12;
    tg     = (100:20:100*20)';
    okMono = all(diff(HS.SurvivalHandle(tg)) <= 1e-15);
    if okBin && okOne && okMono
        fprintf(['[PASS] Test 23: SurvivalHandle matches the likelihood ' ...
            '(max dev %.1e), is 1 below xmin, and is monotone\n'], ...
            max(abs(p - exp(HS.PointwiseLogLik))));
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 23: bins=%d one=%d mono=%d\n', okBin, okOne, okMono);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 23: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Summary
fprintf('\n----------------------------------------\n');
fprintf('Summary: %d passed, %d failed, %d total\n', nPassed, nFailed, nPassed + nFailed);
fprintf('----------------------------------------\n\n');

end

%% Local functions

function [nP, nF] = report(ok, nP, nF, idx, passMsg, failMsg)
if ok
    fprintf('[PASS] Test %d: %s\n', idx, passMsg); nP = nP + 1;
else
    fprintf('[FAIL] Test %d: %s\n', idx, failMsg); nF = nF + 1;
end
end

function d = simEW(lambda, k, alpha, xmin, n, mode, dt)
% Draw n exponentiated-Weibull durations conditioned on T >= xmin, by
% inverting F(t) = [1-exp(-(t/lambda)^k)]^alpha on (F(xmin), 1). In
% discrete mode the draws are placed on the dt grid the way the model
% assumes: n = ceil(t/dt), so the observed value is n*dt exactly.
Fmin = (1 - exp(-(xmin/lambda)^k))^alpha;
if strcmp(mode, 'discrete')
    nmin = max(1, round(xmin/dt));
    out = zeros(n,1); filled = 0;
    while filled < n
        m = 4*(n - filled) + 500;
        p = Fmin + (1 - Fmin)*rand(m,1);
        u = -log(max(1 - p.^(1/alpha), realmin));
        t = lambda * u.^(1/k);
        nn = ceil(t/dt); nn = nn(nn >= nmin);
        take = min(numel(nn), n - filled);
        out(filled+1:filled+take) = nn(1:take)*dt;
        filled = filled + take;
    end
    d = out;
else
    p = Fmin + (1 - Fmin)*rand(n, 1);
    u = -log(max(1 - p.^(1/alpha), realmin));
    d = lambda * u.^(1/k);
end
end

function S = StruncEW(t, lambda, k, alpha, xmin)
% Left-truncated survival function S(t)/S(xmin) of the exponentiated
% Weibull, used to compare two parameter triples as distributions rather
% than coordinate-by-coordinate.
S = -expm1(logEWCDF(t, lambda, k, alpha)) ...
    ./ -expm1(logEWCDF(xmin, lambda, k, alpha));
end

function lf = logEWCDF(t, lambda, k, alpha)
% log F(t), with log(1-exp(-u)) taken on its accurate branch.
t = t(:);
lf = -inf(size(t));
pos = t > 0;
u = exp(k*(log(t(pos)) - log(lambda)));
y = zeros(size(u));
sm = u < log(2);
y(sm) = log(-expm1(-u(sm)));
y(~sm) = log1p(-exp(-u(~sm)));
lf(pos) = alpha * y;
end

function f = ewPDF(t, lambda, k, alpha)
% Untruncated exponentiated-Weibull density.
t = t(:);
u = (t/lambda).^k;
f = (alpha*k/lambda) * (t/lambda).^(k-1) .* exp(-u) .* (1 - exp(-u)).^(alpha-1);
end

function ll = recomputeLogLik(d, xmin, dt, lambda, k, alpha, isDiscrete)
% Independent recomputation of the left-truncated log-likelihood straight
% from the CDF, to check the fitter's internal version.
d = d(:);
if isDiscrete
    nn = round(d/dt);
    nmin = max(1, round(xmin/dt));
    P = exp(logEWCDF(nn*dt, lambda, k, alpha)) ...
        - exp(logEWCDF((nn-1)*dt, lambda, k, alpha));
    S = 1 - exp(logEWCDF((nmin-1)*dt, lambda, k, alpha));
    ll = sum(log(P)) - numel(nn)*log(S);
else
    u = (d/lambda).^k;
    logf = log(alpha) + log(k) - log(lambda) + (k-1)*log(d/lambda) ...
        - u + (alpha-1)*log(1 - exp(-u));
    S = 1 - exp(logEWCDF(xmin, lambda, k, alpha));
    ll = sum(logf) - numel(d)*log(S);
end
end
