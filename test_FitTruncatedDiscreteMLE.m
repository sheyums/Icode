function test_FitTruncatedDiscreteMLE()
%TEST_FITTRUNCATEDDISCRETEMLE Smoke tests for the shared engine and the
%thin per-model wrappers built on it.
%
%   Run from the directory containing FitTruncatedDiscreteMLE.m:
%
%       test_FitTruncatedDiscreteMLE
%
% Sample sizes and multistart counts are deliberately small. These are
% correctness tests -- does each family fit, is the likelihood the one the
% header claims, do the guards fire -- not an assessment of estimator
% precision. Octave's gammainc is interpreted and roughly two orders of
% magnitude slower than MATLAB's, so a full-size sweep is impractical
% there; keeping n small keeps the suite runnable in both.
%
% Tests covered
% -------------
%  1.  Every registry model fits      - finite logL, right parameter count
%  2.  Required fields present        - engine output and diagnostics
%  3.  Likelihood reimplemented       - vs an independent CDF recomputation
%  4.  Truncated pmf normalization    - sums to 1 over the support
%  5.  PointwiseLogLik                - sums to LogLik for every model
%  6.  AIC / AICc / BIC formulas      - exact arithmetic
%  7.  SamplingInterval is required   - omitted and non-positive both error
%  8.  Beta requires UpperBound       - and SupportCoverage is right
%  9.  Bad inputs error               - unknown model, missing Shape
% 10.  Recovery: gamma                - shape and scale near truth
% 11.  Recovery: weibull              - scale and shape near truth
% 12.  Recovery: powerlaw             - alpha near truth
% 13.  Pearson III location guard     - fires when it runs up to xmin
% 14.  powerlaw_cutoff IS gamma       - identical logL, AICc, alpha=1-shape
% 15.  Erlang integer-shape sweep     - integer shape, NaN SE, k=2, sweep
% 16.  Unit invariance                - discrete logL is unit-free
% 17.  Mode convergence as dt -> 0    - discrete logL - n*log(dt) -> cont.
% 18.  Soft return                    - TooFewData throws, then Failed=true
% 19.  Grid mismatch guard            - wrong SamplingInterval units warn
% 20.  Toolbox cross-check            - vs gamcdf/wblcdf, SKIPPED if absent
% 21.  VERBOSE mode                   - does not error
% 22.  SurvivalHandle consistency     - matches the fitted likelihood
% 23.  logL is never positive         - the truncated-probability invariant

FAST = {'SamplingInterval', 1, 'nStarts', 3, 'Verbose', false};
ALL  = {'gamma','chisquared','pearson3','weibull','beta','powerlaw'};

nPassed = 0; nFailed = 0;
fprintf('\nRunning tests for FitTruncatedDiscreteMLE...\n\n');

%% Test 1: every registry model fits
try
    rng(1); xmin = 100;
    d = simGammaDisc(0.6, 900, xmin, 1, 250);
    bad = {};
    expect = struct('gamma',2,'chisquared',1,'pearson3',3,'weibull',2,'beta',2,'powerlaw',1);
    for m = ALL
        H = fitModel(d, xmin, m{1}, FAST);
        if ~isfinite(H.LogLik) || H.LogLik > 0 || H.k ~= expect.(m{1}) ...
                || numel(H.Params) ~= H.k || numel(H.ParamNames) ~= H.k
            bad{end+1} = sprintf('%s (logL=%g k=%d)', m{1}, H.LogLik, H.k);
        end
    end
    Hf = FitTruncatedDiscreteMLE(d, xmin, "gamma_fixedshape", FAST{:}, 'Shape', 2);
    if Hf.k ~= 1, bad{end+1} = 'gamma_fixedshape k~=1'; end
    [nPassed,nFailed] = rep(isempty(bad),nPassed,nFailed,1, ...
        'all 7 registry entries fit with finite logL and the right k', strjoin(bad,'; '));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,1,'',err.message);
end

%% Test 2: required fields present
try
    rng(2); d = simGammaDisc(0.6, 900, 100, 1, 200);
    H = fitModel(d, 100, 'gamma', FAST);
    need = {'Model','ParamNames','Params','ParamSE','SurvivalHandle','CovValid','k','n', ...
            'LogLik','PointwiseLogLik','AIC','AICc','BIC','Success', ...
            'Converged','ExitFlag','BestParamVector','DistributionType', ...
            'SamplingInterval','xmin','Failed','FailureReason','Diagnostics'};
    needD = {'nNonFinite','nBelowXmin','nAtXmin','XminEqualsDataMin', ...
             'LooksGridded','GridSpacing','GridMismatch','nDistinctData', ...
             'nDistinctOnGrid','LocationGap','SupportCoverage','GuardOK', ...
             'Recommendation'};
    miss = [need(~isfield(H,need)), needD(~isfield(H.Diagnostics,needD))];
    [nPassed,nFailed] = rep(isempty(miss),nPassed,nFailed,2, ...
        'all documented fields present', strjoin(miss,','));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,2,'',err.message);
end

%% Test 3: LogLik matches an independent recomputation from the CDF
try
    rng(3); xmin = 100; dt = 1; bad = {};
    d = simGammaDisc(0.6, 900, xmin, dt, 250);
    for m = {'gamma','weibull','powerlaw','beta'}
        H = fitModel(d, xmin, m{1}, FAST);
        ll = recomputeLogLik(d, xmin, dt, m{1}, H.Params, 86400);
        if abs(ll - H.LogLik) > 1e-8
            bad{end+1} = sprintf('%s off by %.3e', m{1}, ll - H.LogLik);
        end
    end
    [nPassed,nFailed] = rep(isempty(bad),nPassed,nFailed,3, ...
        'LogLik matches an independent CDF recomputation', strjoin(bad,'; '));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,3,'',err.message);
end

%% Test 4: the fitted truncated pmf sums to 1
try
    rng(4); xmin = 50; dt = 1; bad = {};
    d = simWeibullDisc(300, 0.8, xmin, dt, 250);
    for m = {'gamma','weibull','powerlaw'}
        H = fitModel(d, xmin, m{1}, FAST);
        nmin = max(1, round(xmin/dt));
        N = 400000;
        nn = (nmin:N)';
        s0 = (nmin-1)*dt;
        S0 = modelSF((nmin-1)*dt, m{1}, H.Params, 86400, s0);
        p = (modelSF((nn-1)*dt, m{1}, H.Params, 86400, s0) - ...
             modelSF(nn*dt, m{1}, H.Params, 86400, s0)) / S0;
        % The power law's tail decays too slowly to sum to convergence --
        % truncating at N left 0.7% of its mass out -- so the remaining
        % tail is added analytically as S(N*dt)/S(xmin). That makes the
        % check exact for every family rather than only the light-tailed
        % ones.
        tailRemainder = modelSF(N*dt, m{1}, H.Params, 86400, s0) / S0;
        total = sum(p) + tailRemainder;
        if abs(total - 1) > 1e-6
            bad{end+1} = sprintf('%s sums to %.10f', m{1}, total);
        end
    end
    [nPassed,nFailed] = rep(isempty(bad),nPassed,nFailed,4, ...
        'truncated pmf sums to 1 for gamma, weibull and powerlaw', strjoin(bad,'; '));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,4,'',err.message);
end

%% Test 5: PointwiseLogLik sums to LogLik
try
    rng(5); d = simGammaDisc(0.6, 900, 100, 1, 200); bad = {};
    for m = ALL
        H = fitModel(d, 100, m{1}, FAST);
        if numel(H.PointwiseLogLik) ~= H.n || abs(sum(H.PointwiseLogLik)-H.LogLik) > 1e-8
            bad{end+1} = m{1};
        end
    end
    [nPassed,nFailed] = rep(isempty(bad),nPassed,nFailed,5, ...
        'PointwiseLogLik is n x 1 and sums to LogLik for every model', strjoin(bad,', '));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,5,'',err.message);
end

%% Test 6: AIC / AICc / BIC arithmetic
try
    rng(6); d = simGammaDisc(0.6, 900, 100, 1, 200); bad = {};
    for m = ALL
        H = fitModel(d, 100, m{1}, FAST);
        aic = 2*H.k - 2*H.LogLik;
        aicc = aic + 2*H.k*(H.k+1)/(H.n - H.k - 1);
        bic = H.k*log(H.n) - 2*H.LogLik;
        if abs(H.AIC-aic) > 1e-9 || abs(H.AICc-aicc) > 1e-9 || abs(H.BIC-bic) > 1e-9
            bad{end+1} = m{1};
        end
    end
    [nPassed,nFailed] = rep(isempty(bad),nPassed,nFailed,6, ...
        'AIC, AICc and BIC exact for every model', strjoin(bad,', '));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,6,'',err.message);
end

%% Test 7: SamplingInterval is required
try
    rng(7); d = simGammaDisc(0.6, 900, 100, 1, 200);
    a = ''; try, FitTruncatedDiscreteMLE(d,100,"gamma",'Verbose',false); catch e, a = e.identifier; end
    b = ''; try, FitTruncatedDiscreteMLE(d,100,"gamma",'SamplingInterval',0,'Verbose',false); catch e, b = e.identifier; end
    H = fitModel(d, 100, 'gamma', FAST);
    ok = strcmp(a,'FitTruncatedDiscreteMLE:SamplingIntervalRequired') && ...
         strcmp(b,'FitTruncatedDiscreteMLE:InvalidSamplingInterval') && ...
         H.Success && isfinite(H.LogLik);
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,7, ...
        'omitted and non-positive SamplingInterval both error; a valid one fits', ...
        sprintf('a="%s" b="%s"',a,b));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,7,'',err.message);
end

%% Test 8: Beta requires UpperBound, and SupportCoverage is right
try
    rng(8); d = simGammaDisc(0.6, 900, 100, 1, 200);
    a = ''; try, FitTruncatedDiscreteMLE(d,100,"beta",FAST{:}); catch e, a = e.identifier; end
    Hb = FitTruncatedDiscreteMLE(d,100,"beta",FAST{:},'UpperBound',86400);
    Hg = fitModel(d, 100, 'gamma', FAST);
    ok = strcmp(a,'FitTruncatedDiscreteMLE:UpperBoundRequired') && ...
         abs(Hb.Diagnostics.SupportCoverage - max(d)/86400) < 1e-12 && ...
         isnan(Hg.Diagnostics.SupportCoverage);
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,8, ...
        sprintf('UpperBound required; SupportCoverage = %.4f, NaN without a bound', ...
        Hb.Diagnostics.SupportCoverage), sprintf('id="%s"',a));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,8,'',err.message);
end

%% Test 9: bad inputs error
try
    rng(9); d = simGammaDisc(0.6, 900, 100, 1, 200);
    a = ''; try, FitTruncatedDiscreteMLE(d,100,"nosuchmodel",FAST{:}); catch e, a = e.identifier; end
    b = ''; try, FitTruncatedDiscreteMLE(d,100,"gamma_fixedshape",FAST{:}); catch e, b = e.identifier; end
    ok = strcmp(a,'FitTruncatedDiscreteMLE:UnknownModel') && ...
         strcmp(b,'FitTruncatedDiscreteMLE:ShapeRequired');
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,9, ...
        'unknown model and missing Shape both error', sprintf('a="%s" b="%s"',a,b));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,9,'',err.message);
end

%% Test 10: recovery, gamma
try
    rng(10); kT = 0.7; thT = 800; xmin = 60;
    d = simGammaDisc(kT, thT, xmin, 1, 4000);
    H = FitGammaMLE(d, xmin, 'SamplingInterval', 1, 'nStarts', 6, 'Verbose', false);
    rel = abs(H.Params - [kT thT]) ./ [kT thT];
    [nPassed,nFailed] = rep(all(rel < 0.35),nPassed,nFailed,10, ...
        sprintf('gamma recovered shape=%.3f scale=%.0f (truth %.2f %.0f)', ...
        H.Params(1),H.Params(2),kT,thT), sprintf('rel err %s',mat2str(round(rel*100)/100)));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,10,'',err.message);
end

%% Test 11: recovery, weibull
try
    rng(11); lamT = 500; kT = 0.8; xmin = 60;
    d = simWeibullDisc(lamT, kT, xmin, 1, 4000);
    H = FitWeibullMLE(d, xmin, 'SamplingInterval', 1, 'nStarts', 6, 'Verbose', false);
    rel = abs(H.Params - [lamT kT]) ./ [lamT kT];
    [nPassed,nFailed] = rep(all(rel < 0.35),nPassed,nFailed,11, ...
        sprintf('weibull recovered scale=%.0f shape=%.3f (truth %.0f %.2f)', ...
        H.Params(1),H.Params(2),lamT,kT), sprintf('rel err %s',mat2str(round(rel*100)/100)));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,11,'',err.message);
end

%% Test 12: recovery, powerlaw
try
    rng(12); aT = 1.4; xmin = 100; dt = 1;
    d = simParetoDisc(xmin-dt, aT, xmin, dt, 4000);
    H = FitPowerLawMLE(d, xmin, 'SamplingInterval', dt, 'nStarts', 6, 'Verbose', false);
    rel = abs(H.Params(1) - aT)/aT;
    [nPassed,nFailed] = rep(rel < 0.15,nPassed,nFailed,12, ...
        sprintf('powerlaw recovered alpha=%.3f (truth %.2f)',H.Params(1),aT), ...
        sprintf('rel err %.3f',rel));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,12,'',err.message);
end

%% Test 13: Pearson III location guard
try
    rng(13); xmin = 300;
    % This is a CONTINUOUS-mode pathology, and the test has to be built for
    % it. With shape < 1 the gamma density diverges as t -> location, so
    % observations sitting exactly at the smallest possible duration drag
    % the location onto them. In DISCRETE mode no such thing can happen:
    % the pmf is a difference of CDFs and is bounded by 1, exactly as ties
    % at xmin are fatal to a continuous mixture but harmless to a discrete
    % one. An earlier version of this test used tightly clustered discrete
    % data and the guard (correctly) stayed silent -- the optimizer instead
    % pushed the location far NEGATIVE with a large shape, which fits a
    % narrow cluster perfectly well.
    d = [repmat(xmin, 40, 1); xmin + (-400*log(rand(360,1)))];
    prev = warning('on','FitTruncatedDiscreteMLE:LocationAtBoundary');
    lastwarn('');
    fprintf('--- begin expected warning ---\n');
    H = FitPearson3MLE(d, xmin, 'SamplingInterval', 1, ...
        'DistributionType', 'continuous', 'nStarts', 6, 'Verbose', false);
    fprintf('--- end expected warning ---\n');
    [~, wid] = lastwarn(); warning(prev);
    fired = strcmp(wid,'FitTruncatedDiscreteMLE:LocationAtBoundary');
    ok = fired && ~H.Diagnostics.GuardOK && H.Diagnostics.LocationGap < 1e-3 ...
         && H.Params(3) < xmin;
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,13, ...
        sprintf('location guard fired (gap = %.2e, location %.4f < xmin)', ...
        H.Diagnostics.LocationGap, H.Params(3)), ...
        sprintf('fired=%d GuardOK=%d gap=%.3g',fired,H.Diagnostics.GuardOK,H.Diagnostics.LocationGap));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,13,'',err.message);
end

%% Test 14: powerlaw_cutoff is the gamma, re-parametrized
try
    rng(14); d = simGammaDisc(0.6, 900, 100, 1, 300);
    C = {'SamplingInterval',1,'nStarts',4,'Verbose',false,'RandomSeed',9};
    Hg = FitGammaMLE(d, 100, C{:});
    Hc = FitPowerLawCutoffMLE(d, 100, C{:});
    ok = abs(Hc.LogLik-Hg.LogLik) < 1e-12 && abs(Hc.AICc-Hg.AICc) < 1e-12 ...
         && Hc.k == Hg.k && abs(Hc.Params(1)-(1-Hg.Params(1))) < 1e-12 ...
         && abs(Hc.Params(2)-Hg.Params(2)) < 1e-12 && ~isempty(Hc.EquivalentTo);
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,14, ...
        sprintf('cutoff alias identical to gamma (alpha = 1-shape = %.4f)',Hc.Params(1)), ...
        sprintf('dlogL=%.3e dAICc=%.3e',Hc.LogLik-Hg.LogLik,Hc.AICc-Hg.AICc));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,14,'',err.message);
end

%% Test 15: Erlang integer-shape sweep
try
    rng(15); d = simGammaDisc(2, 400, 100, 1, 400);
    H = FitErlangMLE(d, 100, 'SamplingInterval', 1, 'nStarts', 3, ...
        'Verbose', false, 'MaxShape', 4);
    ok = H.Params(1) == round(H.Params(1)) && H.Params(1) >= 1 && H.k == 2 ...
         && isnan(H.ParamSE(1)) && isfinite(H.ParamSE(2)) ...
         && numel(H.ShapeSweep) == 4 && strcmp(H.Model,'erlang') ...
         && abs(H.LogLik - max(H.ShapeSweep)) < 1e-9;
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,15, ...
        sprintf('integer shape %d, NaN shape SE, k=2, sweep reported and consistent', ...
        H.Params(1)), sprintf('shape=%g k=%d se=%s',H.Params(1),H.k,mat2str(H.ParamSE)));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,15,'',err.message);
end

%% Test 16: unit invariance -- a pmf is unit-free
try
    rng(16); c = 60; xmin = 120;
    d = simGammaDisc(0.7, 900, xmin, 1, 600);
    Hs = FitGammaMLE(d,      xmin,   'SamplingInterval', 1,   'nStarts',4,'Verbose',false,'RandomSeed',2);
    Hm = FitGammaMLE(d/c,    xmin/c, 'SamplingInterval', 1/c, 'nStarts',4,'Verbose',false,'RandomSeed',2);
    ok = abs(Hs.LogLik - Hm.LogLik) < 1e-6 ...
         && abs(Hm.Params(1) - Hs.Params(1)) < 1e-4 ...
         && abs(Hm.Params(2)*c - Hs.Params(2))/Hs.Params(2) < 1e-3;
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,16, ...
        'discrete logL identical under a unit change; shape invariant, scale scales', ...
        sprintf('dlogL=%.3e',Hs.LogLik-Hm.LogLik));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,16,'',err.message);
end

%% Test 17: discrete -> continuous as dt -> 0
try
    rng(17); xmin = 100;
    d = 100 + (-600*log(rand(500,1)));      % genuinely continuous
    Hc = FitGammaMLE(d, xmin, 'SamplingInterval', 1, ...
        'DistributionType','continuous','nStarts',4,'Verbose',false,'RandomSeed',2);
    fprintf('--- begin expected warnings (continuous data binned on purpose) ---\n');
    gaps = [];
    for dt = [1, 0.1]
        Hd = FitGammaMLE(d, xmin, 'SamplingInterval', dt, 'nStarts',4, ...
            'Verbose',false,'RandomSeed',2);
        gaps(end+1) = (Hd.LogLik - Hd.n*log(dt)) - Hc.LogLik;
    end
    fprintf('--- end expected warnings ---\n');
    ok = abs(gaps(2)) < abs(gaps(1)) && abs(gaps(2)) < 0.05*numel(d);
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,17, ...
        sprintf('gap shrinks with dt: %.3f then %.4f nats',gaps(1),gaps(2)), ...
        sprintf('gaps %s',mat2str(round(gaps*1000)/1000)));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,17,'',err.message);
end

%% Test 18: soft return
try
    a = ''; try, FitGammaMLE([110;120;130],100,FAST{:}); catch e, a = e.identifier; end
    H = FitGammaMLE([110;120;130],100,FAST{:},'ErrorOnNoValidFit',false);
    ok = strcmp(a,'FitTruncatedDiscreteMLE:TooFewData') && H.Failed && ~H.Success ...
         && all(isnan(H.Params)) && ~isempty(strfind(H.FailureReason,'TooFewData'));
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,18, ...
        'TooFewData throws by default and soft-returns when asked', ...
        sprintf('id="%s" Failed=%d',a,H.Failed));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,18,'',err.message);
end

%% Test 19: grid mismatch guard
try
    rng(19); dsec = simGammaDisc(0.7, 900, 120, 1, 400);
    prev = warning('on','FitTruncatedDiscreteMLE:GridMismatch');
    lastwarn('');
    fprintf('--- begin expected warning ---\n');
    Hbad = FitGammaMLE(dsec/60, 2, 'SamplingInterval', 1, 'nStarts',3,'Verbose',false);
    fprintf('--- end expected warning ---\n');
    [~, wid] = lastwarn(); warning(prev);
    Hok = FitGammaMLE(dsec, 120, 'SamplingInterval', 1, 'nStarts',3,'Verbose',false);
    ok = strcmp(wid,'FitTruncatedDiscreteMLE:GridMismatch') && ...
         Hbad.Diagnostics.GridMismatch && ~Hok.Diagnostics.GridMismatch;
    [nPassed,nFailed] = rep(ok,nPassed,nFailed,19, ...
        sprintf('wrong-unit SamplingInterval caught (%d -> %d distinct)', ...
        Hbad.Diagnostics.nDistinctData,Hbad.Diagnostics.nDistinctOnGrid), ...
        sprintf('wid="%s"',wid));
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,19,'',err.message);
end

%% Test 20: cross-check against the Statistics Toolbox, if present
try
    if exist('gamcdf','file') ~= 2 || exist('wblcdf','file') ~= 2
        fprintf(['[PASS] Test 20: skipped, Statistics Toolbox not available ' ...
                 '(gamcdf/wblcdf absent)\n']);
        nPassed = nPassed + 1;
    else
        rng(20); xmin = 100; dt = 1; bad = {};
        d = simGammaDisc(0.6, 900, xmin, dt, 300);
        nn = round(d/dt); nmin = max(1,round(xmin/dt));
        Hg = FitGammaMLE(d, xmin, FAST{:});
        p = gamcdf(nn*dt,Hg.Params(1),Hg.Params(2)) - gamcdf((nn-1)*dt,Hg.Params(1),Hg.Params(2));
        S = 1 - gamcdf((nmin-1)*dt,Hg.Params(1),Hg.Params(2));
        if abs(sum(log(p)) - numel(nn)*log(S) - Hg.LogLik) > 1e-6
            bad{end+1} = 'gamma vs gamcdf';
        end
        Hw = FitWeibullMLE(d, xmin, FAST{:});
        pw = wblcdf(nn*dt,Hw.Params(1),Hw.Params(2)) - wblcdf((nn-1)*dt,Hw.Params(1),Hw.Params(2));
        Sw = 1 - wblcdf((nmin-1)*dt,Hw.Params(1),Hw.Params(2));
        if abs(sum(log(pw)) - numel(nn)*log(Sw) - Hw.LogLik) > 1e-6
            bad{end+1} = 'weibull vs wblcdf';
        end
        [nPassed,nFailed] = rep(isempty(bad),nPassed,nFailed,20, ...
            'engine logL matches gamcdf and wblcdf independently', strjoin(bad,'; '));
    end
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,20,'',err.message);
end

%% Test 21: VERBOSE mode
try
    rng(21); d = simGammaDisc(0.6, 900, 100, 1, 200);
    fprintf('--- begin expected verbose output ---\n');
    H = FitGammaMLE(d, 100, 'SamplingInterval', 1, 'nStarts', 3); %#ok<NASGU>
    fprintf('--- end expected verbose output ---\n');
    [nPassed,nFailed] = rep(true,nPassed,nFailed,21,'verbose mode runs without error','');
catch err
    [nPassed,nFailed] = rep(false,nPassed,nFailed,21,'',err.message);
end

%% Test 22: SurvivalHandle is consistent with the likelihood
try
    rng(221);
    dtS = 1; dS = simGammaDisc(0.6, 900, 100, dtS, 250);
    HS = FitGammaMLE(dS, 100, 'SamplingInterval', dtS, 'nStarts', 3, 'Verbose', false);
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
        fprintf(['[PASS] Test 22: SurvivalHandle matches the likelihood ' ...
            '(max dev %.1e), is 1 below xmin, and is monotone\n'], ...
            max(abs(p - exp(HS.PointwiseLogLik))));
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 22: bins=%d one=%d mono=%d\n', okBin, okOne, okMono);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 22: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 23: the log-likelihood can never be positive
%  A left-truncated log-likelihood is a sum of log CONDITIONAL
%  probabilities, so it is at most 0 by construction. An earlier version
%  floored the bin probability p at realmin without touching the
%  normalizer S(xmin); where S(xmin) had decayed into the subnormals
%  (4.9e-324, which still passes S>0) the floored ratio p/S reached 4.5e15
%  and each observation contributed +36 instead of a negative number. A
%  Weibull with scale 10.7 fitted to data starting at 300 reported
%  logL = +39972 and beat every honest model in a comparison table.
%
%  A large xmin/scale ratio is what drives the multistarts through that
%  region, so this fits every family at xmin=300 with dt=30 -- the sleep
%  protocol's own numbers -- and checks the invariant. It also checks the
%  fit is a real one: the manufactured optimum sat at a scale thirty times
%  below xmin, so an implausibly small scale is the fingerprint.
try
    rng(23); xminF = 300; dtF = 30;
    dF = simGammaDisc(0.5, 1200, xminF, dtF, 300);
    FAR = {'SamplingInterval', dtF, 'nStarts', 6, 'Verbose', false};
    bad = {};
    for m = ALL
        H = fitModel(dF, xminF, m{1}, FAR);
        if ~(H.LogLik < 0)
            bad{end+1} = sprintf('%s: logL=%+.2f', m{1}, H.LogLik); %#ok<AGROW>
            continue
        end
        if ~all(H.PointwiseLogLik <= 0)
            bad{end+1} = sprintf('%s: %d positive pointwise terms', m{1}, ...
                nnz(H.PointwiseLogLik > 0)); %#ok<AGROW>
        end
    end
    % and specifically the family that exhibited it, with its scale checked
    Hw = fitModel(dF, xminF, 'weibull', FAR);
    if Hw.Params(1) < xminF/10
        bad{end+1} = sprintf('weibull scale %.3g is far below xmin', Hw.Params(1));
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 23, ...
        sprintf('every family keeps logL < 0 at xmin/dt=%d (weibull scale %.4g)', ...
            round(xminF/dtF), Hw.Params(1)), ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 23, '', err.message);
end

%% Summary
fprintf('\n----------------------------------------\n');
fprintf('Summary: %d passed, %d failed, %d total\n', nPassed, nFailed, nPassed+nFailed);
fprintf('----------------------------------------\n\n');

end

%% Local functions

function [nP,nF] = rep(ok,nP,nF,idx,passMsg,failMsg)
if ok
    fprintf('[PASS] Test %d: %s\n', idx, passMsg); nP = nP + 1;
else
    fprintf('[FAIL] Test %d: %s\n', idx, failMsg); nF = nF + 1;
end
end

function H = fitModel(d, xmin, name, FAST)
% Beta needs its bound; everything else does not take one.
if strcmp(name,'beta')
    H = FitTruncatedDiscreteMLE(d, xmin, name, FAST{:}, 'UpperBound', 86400);
else
    H = FitTruncatedDiscreteMLE(d, xmin, name, FAST{:});
end
end

function x = gammaDraw(k, theta, n)
% Marsaglia-Tsang, with the standard boost for shape < 1. Used instead of
% gamrnd so the suite does not require the Statistics Toolbox.
x = zeros(n,1);
for i = 1:n
    a = k; boost = 1;
    if a < 1
        boost = rand^(1/a); a = a + 1;
    end
    dd = a - 1/3; c = 1/sqrt(9*dd);
    while true
        z = randn; v = (1 + c*z)^3;
        if v > 0 && log(rand) < 0.5*z^2 + dd - dd*v + dd*log(v)
            x(i) = dd*v*boost; break
        end
    end
end
x = x * theta;
end

function d = simGammaDisc(k, theta, xmin, dt, n)
% Gamma draws placed on the dt grid the discrete model assumes
% (n = ceil(t/dt)) and conditioned on clearing xmin.
nmin = max(1, round(xmin/dt)); out = zeros(n,1); filled = 0;
while filled < n
    t = gammaDraw(k, theta, 4*(n-filled) + 200);
    nv = ceil(t/dt); nv = nv(nv >= nmin);
    take = min(numel(nv), n-filled);
    out(filled+1:filled+take) = nv(1:take)*dt; filled = filled + take;
end
d = out;
end

function d = simWeibullDisc(lambda, k, xmin, dt, n)
nmin = max(1, round(xmin/dt)); out = zeros(n,1); filled = 0;
while filled < n
    m = 4*(n-filled) + 200;
    t = lambda * (-log(rand(m,1))).^(1/k);
    nv = ceil(t/dt); nv = nv(nv >= nmin);
    take = min(numel(nv), n-filled);
    out(filled+1:filled+take) = nv(1:take)*dt; filled = filled + take;
end
d = out;
end

function d = simParetoDisc(s0, alpha, xmin, dt, n)
nmin = max(1, round(xmin/dt)); out = zeros(n,1); filled = 0;
while filled < n
    m = 4*(n-filled) + 200;
    t = s0 * rand(m,1).^(-1/alpha);
    nv = ceil(t/dt); nv = nv(nv >= nmin);
    take = min(numel(nv), n-filled);
    out(filled+1:filled+take) = nv(1:take)*dt; filled = filled + take;
end
d = out;
end

function S = modelSF(t, name, th, B, s0)
% Survival function for the registry families, mirroring the engine.
t = max(t, 0);
switch name
    case 'gamma',      S = gammainc(t/th(2), th(1), 'upper');
    case 'chisquared', S = gammainc(t/2, th(1)/2, 'upper');
    case 'pearson3',   S = gammainc(max(t-th(3),0)/th(2), th(1), 'upper');
    case 'weibull',    S = exp(-(t/th(1)).^th(2));
    case 'beta',       S = betainc(min(t/B,1), th(1), th(2), 'upper');
    case 'powerlaw',   S = min((s0./max(t,s0)).^th(1), 1);
    otherwise,         error('modelSF: %s', name);
end
end

function ll = recomputeLogLik(d, xmin, dt, name, th, B)
% Independent recomputation straight from the CDF, to check the engine.
d = d(:); nn = round(d/dt); nmin = max(1, round(xmin/dt));
switch name
    case 'gamma'
        F = @(t) gammainc(max(t,0)/th(2), th(1), 'lower');
    case 'weibull'
        F = @(t) -expm1(-(max(t,0)/th(1)).^th(2));
    case 'beta'
        F = @(t) betainc(min(max(t,0)/B,1), th(1), th(2), 'lower');
    case 'powerlaw'
        s0 = (nmin-1)*dt;
        F = @(t) 1 - min((s0./max(t,s0)).^th(1), 1);
    otherwise
        error('recomputeLogLik: %s', name);
end
p = F(nn*dt) - F((nn-1)*dt);
S = 1 - F((nmin-1)*dt);
ll = sum(log(p)) - numel(nn)*log(S);
end
