function test_FitHyperexponentialMLE()
%TEST_FITHYPEREXPONENTIALMLE Smoke tests for FitHyperexponentialMLE.
%
% Run from the directory containing FitHyperexponentialMLE.m:
%
%     test_FitHyperexponentialMLE
%
% Each test prints [PASS] or [FAIL] with a short description.
% A summary line is printed at the end.
%
% Tests covered
% -------------
%  1.  Parameter recovery              - tau and observed weights near truth
%  2.  Required fields present         - all documented fields exist
%  3.  Component pairing and ordering  - Tau ascending, vectors aligned
%  4.  WeightsObserved formula         - matches w_j*exp(-lam_j*xmin)/S(xmin)
%  5.  w vs q differ under truncation  - the reporting trap is real
%  6.  Unbounded likelihood detected   - ties at xmin are caught and warned
%  7.  xmin = min(data) flagged        - XminEqualsDataMin diagnostic
%  8.  Degeneracy gating               - unidentified K excluded from selection
%  9.  AIC / AICc formulas             - exact arithmetic check
% 10.  Single exponential -> K=1       - no spurious extra components
% 11.  MaxRate ceiling respected       - no tau below SamplingInterval
% 12.  Discrete n_min clamp            - xmin < dt/2 cannot inflate logL
% 13.  Discrete truncation consistency - rounds into support, not dropped
% 14.  SE fields are NaN, not 0        - when CovValid is false and K>1
% 15.  Warning state not leaked        - singularMatrix state preserved
% 16.  Insufficient data               - correct error identifier
% 17.  Discrete-mode recovery          - tau and q near truth on a grid
% 18.  VERBOSE mode                    - does not error
% 19.  Grid mismatch detected           - wrong SamplingInterval units warn
% 20.  No false-positive grid warning   - correct units, and continuous data
% 21.  Unit invariance                  - seconds vs minutes give same fit
% 22.  Soft return on NoValidFit        - ErrorOnNoValidFit=false
% 23.  Soft return on TooFewData        - ErrorOnNoValidFit=false
% 24.  Batch struct-array compatibility - mixed successes and failures

% Written as a function file rather than a script so that the local
% simulator at the bottom is in scope for every test: MATLAB requires a
% script's local functions to follow all script code, and a function file
% satisfies that without constraining where the tests sit.

fitfun = @FitHyperexponentialMLE;

% Light multistart settings throughout: these are correctness smoke tests,
% not an assessment of how well the multistart explores a hard likelihood.
FAST = {'nStartsBase', 4, 'nStartsPerComponent', 3, 'maxStarts', 12, ...
        'Verbose', false};

nPassed = 0;
nFailed = 0;

fprintf('\nRunning tests for FitHyperexponentialMLE...\n\n');

%% Test 1: Parameter recovery
try
    rng(42);
    wTrue = [0.8 0.2]; tauTrue = [150 1500]; xmin = 300;
    d = simTrunc(wTrue, tauTrue, xmin, 6000, 1e-6);
    H = fitfun(d, xmin, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 3, FAST{:});
    s = H.Selected;
    qTrue = wTrue .* exp(-xmin ./ tauTrue); qTrue = qTrue / sum(qTrue);
    okTau = H.SelectedK == 2 && all(abs(s.Tau - tauTrue) ./ tauTrue < 0.25);
    okQ   = H.SelectedK == 2 && all(abs(s.WeightsObserved - qTrue) < 0.08);
    if okTau && okQ
        fprintf('[PASS] Test 1: parameter recovery (K=%d, tau=[%.0f %.0f])\n', ...
            H.SelectedK, s.Tau(1), s.Tau(2));
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 1: parameter recovery (K=%d, tau=%s, q=%s)\n', ...
            H.SelectedK, mat2str(round(s.Tau)), mat2str(round(s.WeightsObserved*1000)/1000));
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 1: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 2: Required fields present
try
    rng(1);
    d = simTrunc([0.8 0.2], [150 1500], 300, 400, 1e-6);
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 2, FAST{:});
    needTop = {'SelectedK','AllFits','Selected','DistributionType','xmin','n','Diagnostics'};
    needFit = {'K','k','n','WeightsObserved','WeightsUntruncated','Rates','Tau','RateSE', ...
               'TauSE','WeightsObservedSE','WeightsUntruncatedSE','CovValid','LogLik','AIC', ...
               'AICc','Success','Converged','Degenerate','DegenerateReason', ...
               'AtRateBound','ExitFlag','BestParamVector'};
    missTop = needTop(~isfield(H, needTop));
    missFit = needFit(~isfield(H.Selected, needFit));
    if isempty(missTop) && isempty(missFit)
        fprintf('[PASS] Test 2: all documented fields present\n'); nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 2: missing %s %s\n', strjoin(missTop,','), strjoin(missFit,','));
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 2: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 3: Component pairing and ordering
try
    rng(7);
    d = simTrunc([0.7 0.3], [100 2000], 300, 1500, 1e-6);
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 3, FAST{:});
    s = H.Selected; K = H.SelectedK;
    sorted   = all(diff(s.Tau) > 0);
    lengthsOK = all(cellfun(@(f) numel(s.(f)) == K, ...
        {'WeightsUntruncated','WeightsObserved','Rates','Tau','RateSE','TauSE','WeightsUntruncatedSE','WeightsObservedSE'}));
    reciprocal = all(abs(s.Tau .* s.Rates - 1) < 1e-10);
    sumsOne  = abs(sum(s.WeightsObserved) - 1) < 1e-8 && abs(sum(s.WeightsUntruncated) - 1) < 1e-8;
    if sorted && lengthsOK && reciprocal && sumsOne
        fprintf('[PASS] Test 3: components sorted by tau, all vectors aligned and normalized\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 3: sorted=%d lengths=%d recip=%d sums=%d\n', ...
            sorted, lengthsOK, reciprocal, sumsOne);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 3: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 4: WeightsObserved matches its closed form
try
    rng(3);
    xmin = 300;
    d = simTrunc([0.7 0.3], [120 1800], xmin, 1200, 1e-6);
    H = fitfun(d, xmin, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 2, FAST{:});
    s = H.AllFits(2);
    qExpected = s.WeightsUntruncated .* exp(-s.Rates * xmin);
    qExpected = qExpected / sum(qExpected);
    if max(abs(s.WeightsObserved - qExpected)) < 1e-10
        fprintf('[PASS] Test 4: WeightsObserved equals w_j*exp(-lam_j*xmin)/S(xmin)\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 4: max deviation %.3e\n', max(abs(s.WeightsObserved - qExpected)));
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 4: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 5: untruncated and observed weights genuinely differ
try
    rng(5);
    xmin = 300;
    d = simTrunc([0.85 0.15], [90 1600], xmin, 2500, 1e-6);
    H = fitfun(d, xmin, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 2, FAST{:});
    s = H.AllFits(2);
    gap = max(abs(s.WeightsUntruncated - s.WeightsObserved));
    if gap > 0.2
        fprintf('[PASS] Test 5: w and q differ by %.2f under truncation (report q)\n', gap);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 5: w and q differ by only %.3f, expected a large gap\n', gap);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 5: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 6: unbounded continuous likelihood is detected and warned about
try
    rng(11);
    d = round(simTrunc([0.8 0.2], [150 1500], 300, 800, 1e-6));
    d(1:6) = 300;                      % force ties exactly at xmin
    % round() also maps naturally-drawn values in (300, 300.5) onto 300, so
    % the tie count is stream-dependent and must not be hard-coded: the
    % invariant is that the diagnostic equals the true number of ties.
    nTies = nnz(d == 300);
    % The warning must stay ENABLED for lastwarn to record it, so the
    % warning text below is expected output for this test.
    prevWarn = warning('on', 'FitHyperexponentialMLE:UnboundedLikelihood');
    lastwarn('');
    fprintf('--- begin expected warning ---\n');
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 3, FAST{:});
    fprintf('--- end expected warning ---\n');
    [~, wid] = lastwarn();
    warning(prevWarn);
    flagged = H.Diagnostics.UnboundedContinuousLikelihood;
    warned  = strcmp(wid, 'FitHyperexponentialMLE:UnboundedLikelihood');
    counted = H.Diagnostics.nAtXmin == nTies;
    if flagged && warned && counted
        fprintf('[PASS] Test 6: all %d ties at xmin detected, unbounded likelihood warned\n', nTies);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 6: flagged=%d warned=%d nAtXmin=%d (true ties %d)\n', ...
            flagged, warned, H.Diagnostics.nAtXmin, nTies);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 6: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 7: xmin = min(data) is flagged
try
    rng(13);
    d = round(simTrunc([0.8 0.2], [150 1500], 300, 600, 1e-6));
    prevWarn = warning('off', 'FitHyperexponentialMLE:UnboundedLikelihood');
    H = fitfun(d, min(d), 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 2, FAST{:});
    warning(prevWarn);
    if H.Diagnostics.XminEqualsDataMin && H.Diagnostics.nAtXmin >= 1
        fprintf('[PASS] Test 7: xmin = min(data) flagged (nAtXmin=%d)\n', H.Diagnostics.nAtXmin);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 7: XminEqualsDataMin=%d nAtXmin=%d\n', ...
            H.Diagnostics.XminEqualsDataMin, H.Diagnostics.nAtXmin);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 7: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 8: degeneracy gating excludes unidentified model orders
try
    rng(17);
    % Genuinely one exponential, so K=3 has nothing real to fit: it must be
    % excluded rather than selected on an unidentified likelihood gain.
    d = 300 - 700*log(rand(500,1));
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 3, FAST{:});
    excludedSomething = ~isempty(H.Diagnostics.ExcludedK);
    reasonsGiven = true;
    for K = H.Diagnostics.ExcludedK
        if isempty(H.AllFits(K).DegenerateReason) || ~H.AllFits(K).Degenerate
            reasonsGiven = false;
        end
    end
    selectedNotDegenerate = ~H.Selected.Degenerate;
    if excludedSomething && reasonsGiven && selectedNotDegenerate
        fprintf('[PASS] Test 8: excluded K=%s with reasons; selected K=%d is identified\n', ...
            mat2str(H.Diagnostics.ExcludedK), H.SelectedK);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 8: excluded=%s reasons=%d selNotDegen=%d\n', ...
            mat2str(H.Diagnostics.ExcludedK), reasonsGiven, selectedNotDegenerate);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 8: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 9: AIC and AICc formulas
try
    rng(19);
    d = simTrunc([0.8 0.2], [150 1500], 300, 500, 1e-6);
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 2, FAST{:});
    ok = true;
    for K = 1:2
        f = H.AllFits(K);
        if ~f.Success, continue; end
        aic  = 2*f.k - 2*f.LogLik;
        aicc = aic + (2*f.k*(f.k+1)) / (f.n - f.k - 1);
        if abs(f.AIC - aic) > 1e-9 || abs(f.AICc - aicc) > 1e-9 || f.k ~= 2*K-1
            ok = false;
        end
    end
    if ok
        fprintf('[PASS] Test 9: AIC/AICc arithmetic and k=2K-1 exact\n'); nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 9: AIC/AICc mismatch\n'); nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 9: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 10: single exponential selects K=1
try
    rng(23);
    d = 300 - 800*log(rand(3000,1));
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 3, FAST{:});
    tauOK = abs(H.Selected.Tau(1) - 800)/800 < 0.15;
    if H.SelectedK == 1 && tauOK
        fprintf('[PASS] Test 10: single exponential -> K=1, tau=%.0f\n', H.Selected.Tau(1));
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 10: SelectedK=%d tau=%s\n', H.SelectedK, mat2str(round(H.Selected.Tau)));
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 10: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 11: MaxRate ceiling respected
try
    rng(29);
    d = round(simTrunc([0.8 0.2], [150 1500], 300, 800, 1e-6));
    d(1:8) = 300;                      % bait the divergence
    prevWarn = warning('off', 'FitHyperexponentialMLE:UnboundedLikelihood');
    H = fitfun(d, 300, 'MaxComponents', 4, 'DistributionType', 'continuous', 'SamplingInterval', 1, FAST{:});
    warning(prevWarn);
    minTau = Inf;
    for K = 1:4
        if H.AllFits(K).Success
            minTau = min(minTau, min(H.AllFits(K).Tau));
        end
    end
    if minTau >= 1 - 1e-6
        fprintf('[PASS] Test 11: no tau below SamplingInterval (min tau = %.4g)\n', minTau);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 11: min tau = %.4g, below the ceiling\n', minTau);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 11: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 12: discrete n_min clamp keeps the normalizer <= 1
try
    rng(31);
    % xmin below dt/2 would give n_min = 0 and a survival above 1, which
    % would show up as a positive log-likelihood contribution.
    d = ceil(-30*log(rand(400,1)));
    H = fitfun(d, 0.4, 'MaxComponents', 2, 'DistributionType', 'discrete', ...
        'SamplingInterval', 1, FAST{:});
    ok = true;
    for K = 1:2
        if H.AllFits(K).Success && H.AllFits(K).LogLik > 0
            ok = false;
        end
    end
    if ok
        fprintf('[PASS] Test 12: xmin<dt/2 clamped to n_min=1, logL stays <= 0\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 12: positive log-likelihood, normalizer exceeded 1\n');
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 12: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 13: discrete truncation keeps values that round into the support
try
    rng(37);
    dt = 1; xmin = 2;
    d = [repmat(1.6, 20, 1); ceil(-20*log(rand(300,1))) + 2];
    H = fitfun(d, xmin, 'MaxComponents', 1, 'DistributionType', 'discrete', ...
        'SamplingInterval', dt, FAST{:});
    % 1.6 rounds to n=2 = n_min, so it is inside the support and must be kept.
    if H.n == numel(d)
        fprintf('[PASS] Test 13: values rounding into the support are retained (n=%d)\n', H.n);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 13: n=%d of %d, in-support values were dropped\n', H.n, numel(d));
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 13: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 14: SE fields are NaN (not 0) when CovValid is false
try
    rng(41);
    % Whether any given fit ends up non-identified depends on the RNG
    % stream, so scan several adversarial datasets rather than relying on
    % one of them to degenerate.
    cases = {[310; 320; 340; 400; 410; 450; 900], ...
             [301; 302; 303; 304; 305; 306], ...
             [305; 305; 306; 306; 307; 307; 308; 308], ...
             300 + [1; 2; 2; 3; 3; 4; 5; 800; 1600]};
    ok = true; checked = 0;
    for ci = 1:numel(cases)
        try
            H = fitfun(cases{ci}, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 3, FAST{:});
        catch innerErr
            % NoValidFit is the correct outcome when every model order is
            % degenerate for a dataset this small; move on to the next.
            if strcmp(innerErr.identifier, 'FitHyperexponentialMLE:NoValidFit')
                continue
            end
            rethrow(innerErr);
        end
        for K = 2:3
            f = H.AllFits(K);
            if f.Success && ~f.CovValid
                checked = checked + 1;
                if ~all(isnan(f.WeightsUntruncatedSE)) || ~all(isnan(f.WeightsObservedSE)) ...
                        || ~all(isnan(f.TauSE)) || ~all(isnan(f.RateSE))
                    ok = false;
                end
            end
        end
    end
    if ok && checked > 0
        fprintf('[PASS] Test 14: all SE fields NaN when CovValid is false (%d case(s) exercised)\n', ...
            checked);
        nPassed = nPassed + 1;
    elseif checked == 0
        fprintf('[FAIL] Test 14: no CovValid=false case arose, assertion never exercised\n');
        nFailed = nFailed + 1;
    else
        fprintf('[FAIL] Test 14: an SE field was 0 rather than NaN with CovValid=false\n');
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 14: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 15: warning state is not leaked
try
    rng(43);
    warning('on', 'MATLAB:singularMatrix');
    before = warning('query', 'MATLAB:singularMatrix');
    d = simTrunc([0.8 0.2], [150 1500], 300, 300, 1e-6);
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 3, FAST{:}); %#ok<NASGU>
    after = warning('query', 'MATLAB:singularMatrix');
    if strcmp(before.state, after.state)
        fprintf('[PASS] Test 15: singularMatrix warning state preserved (%s)\n', after.state);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 15: warning state leaked (%s -> %s)\n', before.state, after.state);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 15: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 16: insufficient data raises the documented error
try
    caught = '';
    try
        fitfun([310; 320; 330], 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 1, FAST{:});
    catch err
        caught = err.identifier;
    end
    if strcmp(caught, 'FitHyperexponentialMLE:TooFewData')
        fprintf('[PASS] Test 16: TooFewData error identifier correct\n'); nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 16: got identifier "%s"\n', caught); nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 16: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 17: discrete-mode recovery on gridded data
try
    rng(47);
    wTrue = [0.75 0.25]; tauTrue = [4 60]; dt = 1; xmin = 2;
    pTrue = 1 - exp(-dt ./ tauTrue);
    nWant = 8000; d = zeros(nWant,1); filled = 0;
    while filled < nWant
        m = 4*(nWant - filled) + 2000;
        pv = pTrue(:);
        j = 1 + (rand(m,1) > wTrue(1));
        nn = ceil(log(rand(m,1)) ./ log(1 - pv(j)));
        nn = nn(nn >= round(xmin/dt));
        take = min(numel(nn), nWant - filled);
        d(filled+1:filled+take) = nn(1:take) * dt;
        filled = filled + take;
    end
    H = fitfun(d, xmin, 'MaxComponents', 3, 'DistributionType', 'discrete', ...
        'SamplingInterval', dt, FAST{:});
    s = H.Selected;
    qTrue = wTrue .* (1-pTrue).^(round(xmin/dt)-1); qTrue = qTrue/sum(qTrue);
    okK = H.SelectedK == 2;
    okT = okK && all(abs(s.Tau - tauTrue)./tauTrue < 0.25);
    okQ = okK && all(abs(s.WeightsObserved - qTrue) < 0.08);
    if okK && okT && okQ
        fprintf('[PASS] Test 17: discrete recovery (K=2, tau=[%.1f %.1f])\n', s.Tau(1), s.Tau(2));
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 17: K=%d tau=%s q=%s\n', H.SelectedK, ...
            mat2str(round(s.Tau*10)/10), mat2str(round(s.WeightsObserved*1000)/1000));
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 17: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 18: VERBOSE mode does not error
try
    rng(53);
    d = simTrunc([0.8 0.2], [150 1500], 300, 400, 1e-6);
    fprintf('--- begin expected verbose output ---\n');
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 2, 'nStartsBase', 4, ...
        'nStartsPerComponent', 3, 'maxStarts', 12, 'Verbose', true); %#ok<NASGU>
    fprintf('--- end expected verbose output ---\n');
    fprintf('[PASS] Test 18: verbose mode runs without error\n'); nPassed = nPassed + 1;
catch err
    fprintf('[FAIL] Test 18: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 19: grid mismatch detected when SamplingInterval is in wrong units
try
    rng(59);
    % Integer-second durations expressed in MINUTES, with SamplingInterval
    % left at its default 1: round(t/1) rounds every duration to a whole
    % minute, which is the silent failure this guard exists to catch.
    dsec = 300 + ceil(-600*log(rand(500,1)));
    dmin = dsec / 60;
    prevWarn = warning('on', 'FitHyperexponentialMLE:GridMismatch');
    lastwarn('');
    fprintf('--- begin expected warning ---\n');
    H = fitfun(dmin, 5, 'SamplingInterval', 1, 'MaxComponents', 1, 'DistributionType', 'discrete', FAST{:});
    fprintf('--- end expected warning ---\n');
    [~, wid] = lastwarn();
    warning(prevWarn);
    warned   = strcmp(wid, 'FitHyperexponentialMLE:GridMismatch');
    spacing  = abs(H.Diagnostics.GridSpacing - 1/60) < 1e-9;
    mism     = H.Diagnostics.GridMismatch;
    collapsed = H.Diagnostics.nDistinctOnGrid < H.Diagnostics.nDistinctData;
    if warned && spacing && mism && collapsed
        fprintf('[PASS] Test 19: wrong-unit SamplingInterval caught (%d -> %d distinct)\n', ...
            H.Diagnostics.nDistinctData, H.Diagnostics.nDistinctOnGrid);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 19: warned=%d spacing=%d mismatch=%d collapsed=%d\n', ...
            warned, spacing, mism, collapsed);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 19: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 20: no false-positive grid warning
try
    rng(61);
    dsec = 300 + ceil(-600*log(rand(500,1)));
    % (a) correct units: minutes with SamplingInterval = 1/60
    Ha = fitfun(dsec/60, 5, 'MaxComponents', 1, 'DistributionType', 'discrete', ...
        'SamplingInterval', 1/60, FAST{:});
    % (b) correct units: seconds with SamplingInterval = 1
    Hb = fitfun(dsec, 300, 'MaxComponents', 1, 'DistributionType', 'discrete', ...
        'SamplingInterval', 1, FAST{:});
    % (c) genuinely continuous durations: the smallest gap is arbitrary and
    %     must NOT be reported as a grid
    dcont = 300 - 900*log(rand(500,1));
    Hc = fitfun(dcont, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 1, FAST{:});
    okA = ~Ha.Diagnostics.GridMismatch && ...
          Ha.Diagnostics.nDistinctOnGrid == Ha.Diagnostics.nDistinctData;
    okB = ~Hb.Diagnostics.GridMismatch && ...
          Hb.Diagnostics.nDistinctOnGrid == Hb.Diagnostics.nDistinctData;
    okC = ~Hc.Diagnostics.GridMismatch && ~Hc.Diagnostics.LooksGridded;
    if okA && okB && okC
        fprintf('[PASS] Test 20: no grid warning for correct units or continuous data\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 20: minutes=%d seconds=%d continuous=%d\n', okA, okB, okC);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 20: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 21: unit invariance -- seconds and minutes must give the same fit
try
    rng(67);
    wT = [0.7 0.3]; tauT = [120 1200]; xs = 300;
    pT = 1 - exp(-1 ./ tauT);
    nWant = 2500; dsec = zeros(nWant,1); filled = 0;
    while filled < nWant
        mm = 4*(nWant - filled) + 2000;
        j = 1 + (rand(mm,1) > wT(1));
        pv = pT(:);
        nn = ceil(log(rand(mm,1)) ./ log(1 - pv(j)));
        nn = nn(nn >= xs);
        take = min(numel(nn), nWant - filled);
        dsec(filled+1:filled+take) = nn(1:take);
        filled = filled + take;
    end
    Hs = fitfun(dsec, xs, 'MaxComponents', 2, 'DistributionType', 'discrete', ...
        'SamplingInterval', 1, 'RandomSeed', 3, FAST{:});
    Hm = fitfun(dsec/60, xs/60, 'MaxComponents', 2, 'DistributionType', 'discrete', ...
        'SamplingInterval', 1/60, 'RandomSeed', 3, FAST{:});
    sameK = Hs.SelectedK == Hm.SelectedK;
    % A pmf is unit-free, so the discrete log-likelihood must be identical.
    sameL = abs(Hs.Selected.LogLik - Hm.Selected.LogLik) < 1e-6;
    sameT = sameK && max(abs(Hm.Selected.Tau*60 - Hs.Selected.Tau) ./ Hs.Selected.Tau) < 1e-3;
    sameQ = sameK && max(abs(Hm.Selected.WeightsObserved - Hs.Selected.WeightsObserved)) < 1e-4;
    if sameK && sameL && sameT && sameQ
        fprintf('[PASS] Test 21: seconds and minutes agree (K=%d, logL identical, tau ratio 60)\n', ...
            Hs.SelectedK);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 21: sameK=%d sameLogL=%d sameTau=%d sameQ=%d\n', ...
            sameK, sameL, sameT, sameQ);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 21: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 22: soft return when no model order is identifiable
try
    rng(71);
    % Six observations within 6 s of xmin: no model order survives the
    % identifiability gating, which by default raises NoValidFit.
    d = [301; 302; 303; 304; 305; 306];
    threw = false;
    try
        fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 3, FAST{:});
    catch err
        threw = strcmp(err.identifier, 'FitHyperexponentialMLE:NoValidFit');
    end
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 3, 'ErrorOnNoValidFit', false, FAST{:});
    okFail   = H.Failed && isnan(H.SelectedK);
    okReason = ~isempty(strfind(H.FailureReason, 'NoValidFit'));
    okFits   = numel(H.AllFits) == 3;
    okSel    = isempty(H.Selected.Tau);
    if threw && okFail && okReason && okFits && okSel
        fprintf('[PASS] Test 22: NoValidFit throws by default, returns Failed=true when asked\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 22: threw=%d failed=%d reason=%d fits=%d emptySel=%d\n', ...
            threw, okFail, okReason, okFits, okSel);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 22: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 23: soft return when there are too few observations
try
    d = [310; 320; 330];
    threw = false;
    try
        fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 2, FAST{:});
    catch err
        threw = strcmp(err.identifier, 'FitHyperexponentialMLE:TooFewData');
    end
    H = fitfun(d, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 2, 'ErrorOnNoValidFit', false, FAST{:});
    okFail   = H.Failed && isnan(H.SelectedK) && H.n == 3;
    okReason = ~isempty(strfind(H.FailureReason, 'TooFewData'));
    if threw && okFail && okReason
        fprintf('[PASS] Test 23: TooFewData throws by default, returns Failed=true when asked\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 23: threw=%d failed=%d reason=%d\n', threw, okFail, okReason);
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 23: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Test 24: successes and failures coexist in one struct array
try
    rng(73);
    % The actual batch-loop pattern: per-animal fits collected into a
    % struct array, some of which failed. MATLAB rejects this assignment if
    % the field names or their order differ between the two outputs.
    good = simTrunc([0.8 0.2], [150 1500], 300, 400, 1e-6);
    samples = {good, [301; 302; 303; 304; 305; 306], [310; 320; 330], good};
    clear results
    for ii = 1:numel(samples)
        results(ii) = fitfun(samples{ii}, 300, 'DistributionType', 'continuous', 'SamplingInterval', 1, 'MaxComponents', 2, ...
            'ErrorOnNoValidFit', false, FAST{:});
    end
    failedFlags = [results.Failed];
    nOK = sum(~failedFlags);
    sel  = [results(~failedFlags).Selected];
    taus = [sel.Tau];
    if numel(results) == 4 && nOK == 2 && all(failedFlags == [false true true false]) ...
            && numel(taus) == 4 && all(isfinite(taus))
        fprintf('[PASS] Test 24: mixed batch assigns into one struct array (%d/%d usable)\n', ...
            nOK, numel(results));
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 24: n=%d nOK=%d flags=%s\n', numel(results), nOK, ...
            mat2str(failedFlags));
        nFailed = nFailed + 1;
    end
catch err
    fprintf('[FAIL] Test 24: errored (%s)\n', err.message); nFailed = nFailed + 1;
end

%% Summary
fprintf('\n----------------------------------------\n');
fprintf('Summary: %d passed, %d failed, %d total\n', nPassed, nFailed, nPassed + nFailed);
fprintf('----------------------------------------\n\n');

end

% Simulate a left-truncated hyperexponential. xmin is a genuine threshold
% here (draws below it are discarded), and jitter keeps observations off the
% truncation point so the continuous likelihood stays bounded.
function d = simTrunc(w, tau, xmin, n, jitter)
    d = zeros(n, 1);
    filled = 0;
    while filled < n
        m = max(2000, 4*(n - filled));
        tauv = tau(:);
        j = 1 + sum(rand(m,1) > cumsum(w(:)') , 2);
        t = -tauv(j) .* log(rand(m,1));
        t = t(t > xmin + jitter);
        take = min(numel(t), n - filled);
        d(filled+1 : filled+take) = t(1:take);
        filled = filled + take;
    end
end
