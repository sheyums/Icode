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
% 14.  WeightSE is NaN, not 0          - when CovValid is false and K>1
% 15.  Warning state not leaked        - singularMatrix state preserved
% 16.  Insufficient data               - correct error identifier
% 17.  Discrete-mode recovery          - tau and q near truth on a grid
% 18.  VERBOSE mode                    - does not error

% MATLAB runs the real function. Octave cannot parse an "arguments" block,
% so under Octave these tests run against the auto-generated twin
% FitHyperexponentialMLE_oct.m (identical body, name-value parsing shim).
%
% Written as a function file rather than a script so that the local
% simulator below is visible in both MATLAB (local functions must follow
% all script code) and Octave (a script's functions must be defined
% before first use). A function file satisfies both.

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    fitfun = @FitHyperexponentialMLE_oct;
else
    fitfun = @FitHyperexponentialMLE;
end

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
    H = fitfun(d, xmin, 'MaxComponents', 3, FAST{:});
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
    H = fitfun(d, 300, 'MaxComponents', 2, FAST{:});
    needTop = {'SelectedK','AllFits','Selected','DistributionType','xmin','n','Diagnostics'};
    needFit = {'K','k','n','Weights','WeightsObserved','Rates','Tau','RateSE', ...
               'TauSE','WeightSE','WeightObservedSE','CovValid','LogLik','AIC', ...
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
    H = fitfun(d, 300, 'MaxComponents', 3, FAST{:});
    s = H.Selected; K = H.SelectedK;
    sorted   = all(diff(s.Tau) > 0);
    lengthsOK = all(cellfun(@(f) numel(s.(f)) == K, ...
        {'Weights','WeightsObserved','Rates','Tau','RateSE','TauSE','WeightSE','WeightObservedSE'}));
    reciprocal = all(abs(s.Tau .* s.Rates - 1) < 1e-10);
    sumsOne  = abs(sum(s.WeightsObserved) - 1) < 1e-8 && abs(sum(s.Weights) - 1) < 1e-8;
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
    H = fitfun(d, xmin, 'MaxComponents', 2, FAST{:});
    s = H.AllFits(2);
    qExpected = s.Weights .* exp(-s.Rates * xmin);
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
    H = fitfun(d, xmin, 'MaxComponents', 2, FAST{:});
    s = H.AllFits(2);
    gap = max(abs(s.Weights - s.WeightsObserved));
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
    % The warning must stay ENABLED for lastwarn to record it, so the
    % warning text below is expected output for this test.
    prevWarn = warning('on', 'FitHyperexponentialMLE:UnboundedLikelihood');
    lastwarn('');
    fprintf('--- begin expected warning ---\n');
    H = fitfun(d, 300, 'MaxComponents', 3, FAST{:});
    fprintf('--- end expected warning ---\n');
    [~, wid] = lastwarn();
    warning(prevWarn);
    flagged = H.Diagnostics.UnboundedContinuousLikelihood;
    warned  = strcmp(wid, 'FitHyperexponentialMLE:UnboundedLikelihood');
    counted = H.Diagnostics.nAtXmin == 6;
    if flagged && warned && counted
        fprintf('[PASS] Test 6: %d ties at xmin detected, unbounded likelihood warned\n', ...
            H.Diagnostics.nAtXmin);
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] Test 6: flagged=%d warned=%d nAtXmin=%d\n', ...
            flagged, warned, H.Diagnostics.nAtXmin);
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
    H = fitfun(d, min(d), 'MaxComponents', 2, FAST{:});
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
    H = fitfun(d, 300, 'MaxComponents', 3, FAST{:});
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
    H = fitfun(d, 300, 'MaxComponents', 2, FAST{:});
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
    H = fitfun(d, 300, 'MaxComponents', 3, FAST{:});
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
    H = fitfun(d, 300, 'MaxComponents', 4, 'SamplingInterval', 1, FAST{:});
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

%% Test 14: WeightSE is NaN (not 0) when CovValid is false
try
    rng(41);
    d = [310; 320; 340; 400; 410; 450; 900];
    H = fitfun(d, 300, 'MaxComponents', 3, FAST{:});
    ok = true; checked = false;
    for K = 2:3
        f = H.AllFits(K);
        if f.Success && ~f.CovValid
            checked = true;
            if ~all(isnan(f.WeightSE)) || ~all(isnan(f.WeightObservedSE)) ...
                    || ~all(isnan(f.TauSE)) || ~all(isnan(f.RateSE))
                ok = false;
            end
        end
    end
    if ok && checked
        fprintf('[PASS] Test 14: all SE fields NaN when CovValid is false\n'); nPassed = nPassed + 1;
    elseif ~checked
        fprintf('[PASS] Test 14: skipped, every fit was identified (no CovValid=false case)\n');
        nPassed = nPassed + 1;
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
    H = fitfun(d, 300, 'MaxComponents', 3, FAST{:}); %#ok<NASGU>
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
        fitfun([310; 320; 330], 300, 'MaxComponents', 1, FAST{:});
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
    H = fitfun(d, 300, 'MaxComponents', 2, 'nStartsBase', 4, ...
        'nStartsPerComponent', 3, 'maxStarts', 12, 'Verbose', true); %#ok<NASGU>
    fprintf('--- end expected verbose output ---\n');
    fprintf('[PASS] Test 18: verbose mode runs without error\n'); nPassed = nPassed + 1;
catch err
    fprintf('[FAIL] Test 18: errored (%s)\n', err.message); nFailed = nFailed + 1;
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
