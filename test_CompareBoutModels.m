function test_CompareBoutModels()
%TEST_COMPAREBOUTMODELS Smoke tests for the model-comparison wrapper.
%
%   Run from the directory containing CompareBoutModels.m:
%
%       test_CompareBoutModels
%
% These test the WRAPPER, not the fitters -- each fitter has its own suite.
% What is checked here is everything the wrapper itself is responsible for:
% that every candidate is fitted on the same data with the same measure, so
% the information criteria are comparable; that the table reports all
% parameters with the right names and keeps k separate from the number
% reported; that ranking, Akaike weights and the walk-down behave; and that
% the expected counts the G test is built from agree with each fit's own
% survival handle to full precision.
%
% Sample sizes and multistart counts are deliberately small: Octave's
% gammainc is interpreted and roughly two orders of magnitude slower than
% MATLAB's, and the library fits many models. Where a test only needs one
% family it restricts Models= to that family.
%
% Tests covered
% -------------
%  1.  SamplingInterval is required   - and non-positive is rejected
%  2.  End-to-end run                 - every documented field present
%  3.  Params holds ALL parameters    - hyperexp K reports 2K, k stays 2K-1
%  4.  Comparability enforced         - one n, one dt, one mode, one xmin
%  5.  n matches the fitters' own     - wrapper's retention rule is theirs
%  6.  powerlaw_cutoff IS gamma       - identical logL and AICc in-table
%  7.  Criterion arithmetic           - AIC, AICc, BIC recomputed from k, n
%  8.  Akaike weights                 - sum to 1, best row has dAICc = 0
%  9.  Ranking follows RankBy         - AICc and BIC orders both monotone
% 10.  Expected counts sum to n       - survival handle and bin convention
% 11.  A true model passes the G test - gamma data, gamma fit
% 12.  A wrong model is rejected      - two-scale mixture read as one gamma
% 13.  Walk-down picks the first pass  - Selected agrees with SelectionPath
% 14.  Degenerate rows cannot win     - flagged, ranked last, not Selected
% 15.  Models filter                  - restricts, and errors when empty
% 16.  Beta needs UpperBound          - absent from the library without one
% 17.  Erlang's integer shape         - k=2, integer shape, NaN SE
% 18.  Continuous mode runs           - and gives a different logL scale
% 19.  Parametric bootstrap runs      - finite p, respects a forced B
% 20.  RandomSeed reproducibility     - identical log-likelihoods twice
% 21.  ParamText                      - non-empty and names every parameter
% 22.  Plot=false leaves no figure    - R.Figure empty, nothing opened
% 23.  Plot=true actually draws      - plotFit's only coverage; SKIPPED
%                                      where there is no graphics toolkit
% 24.  Order LRT wiring              - off/on, and what comes back

if exist('OCTAVE_VERSION', 'builtin')
    cbm = @CompareBoutModels_oct;
    hyp = @FitHyperexponentialMLE_oct;
    eng = @FitTruncatedDiscreteMLE_oct;
else
    cbm = @CompareBoutModels;
    hyp = @FitHyperexponentialMLE;
    eng = @FitTruncatedDiscreteMLE;
end

DT = 10; XMIN = 100;
% OrderLRT off by default here: at B=999 an auto-triggered order test
% would dominate the suite's runtime, and it is exercised deliberately in
% test 24 with a small B rather than incidentally wherever two orders
% happen to land close.
BASE = {'SamplingInterval', DT, 'Plot', false, 'Verbose', false, ...
        'OrderLRT', false};
SMALL = {'hyperexponential', 'gamma', 'weibull', 'powerlaw'};

nPassed = 0; nFailed = 0;
fprintf('\nRunning tests for CompareBoutModels...\n\n');

%% Test 1: SamplingInterval is required
try
    rng(1);
    d = simGammaDisc(0.7, 400, XMIN, DT, 300);
    e1 = ''; e2 = '';
    try
        cbm(d, XMIN, 'Plot', false, 'Verbose', false);
    catch err
        e1 = err.identifier;
    end
    try
        cbm(d, XMIN, 'SamplingInterval', 0, 'Plot', false, 'Verbose', false);
    catch err
        e2 = err.identifier;
    end
    ok = strcmp(e1, 'CompareBoutModels:SamplingIntervalRequired') && ~isempty(e2);
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 1, ...
        'SamplingInterval is required and must be positive', ...
        sprintf('got "%s" and "%s"', e1, e2));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 1, '', err.message);
end

%% Test 2: end-to-end, every documented field present
try
    rng(2);
    d = simHyperDisc([150 1200], [0.65 0.35], XMIN, DT, 600);
    R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 2, 'Models', SMALL, ...
        'GoFBootstrap', 0);
    need = {'Table','Fits','Selected','SelectionPath','GoF','Nesting', ...
            'Figure','n','xmin','SamplingInterval','DistributionType','RankBy'};
    missing = need(~isfield(R, need));
    ok = isempty(missing) && R.n == numel(d) && R.xmin == XMIN && ...
        R.SamplingInterval == DT && numel(R.Fits) == numel(tableRows(R.Table));
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 2, ...
        sprintf('end-to-end run, %d rows, all fields present', numel(R.Fits)), ...
        sprintf('missing: %s', strjoin(missing, ', ')));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 2, '', err.message);
end

%% Test 3: Params holds ALL parameters, k counts only the free ones
try
    rng(3);
    d = simHyperDisc([150 1200], [0.65 0.35], XMIN, DT, 600);
    R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 3, ...
        'Models', {'hyperexponential'}, 'GoFBootstrap', 0);
    T = tableRows(R.Table);
    bad = {};
    for i = 1:numel(T)
        K = sscanf(T(i).Model, 'hyperexp K=%d');
        pn = T(i).ParamNames; pv = T(i).Params; pse = T(i).ParamSE;
        wantNames = [arrayfun(@(j) sprintf('tau%d', j), 1:K, 'UniformOutput', false), ...
                     arrayfun(@(j) sprintf('q%d',   j), 1:K, 'UniformOutput', false)];
        if T(i).nReported ~= 2*K || T(i).k ~= 2*K-1 || numel(pv) ~= 2*K || ...
                numel(pse) ~= 2*K || ~isequal(pn, wantNames)
            bad{end+1} = sprintf('K=%d: nReported=%d k=%d numel=%d', ...
                K, T(i).nReported, T(i).k, numel(pv)); %#ok<AGROW>
            continue
        end
        if ~T(i).Degenerate
            % all K weights are reported even though sum(q)=1 makes one
            % of them redundant -- that is the point of the field
            q = pv(K+1:end);
            if abs(sum(q) - 1) > 1e-8 || any(pv(1:K) <= 0)
                bad{end+1} = sprintf('K=%d: sum(q)=%.12g', K, sum(q)); %#ok<AGROW>
            end
        end
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 3, ...
        'Params reports 2K numbers at order K while k stays 2K-1, sum(q)=1', ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 3, '', err.message);
end

%% Test 4: comparability -- one n, one dt, one mode, one xmin for all rows
try
    rng(4);
    d = simGammaDisc(0.7, 500, XMIN, DT, 500);
    R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 2, 'GoFBootstrap', 0, ...
        'UpperBound', 86400);
    T = tableRows(R.Table);
    ns = [T.n];
    ok = all(ns == ns(1)) && ns(1) == R.n;
    % and the fits themselves agree about the measure they were taken under
    for i = 1:numel(R.Fits)
        F = R.Fits(i).Full;
        if isfield(F, 'DistributionType') && ~strcmp(char(F.DistributionType), 'discrete')
            ok = false;
        end
        if isfield(F, 'SamplingInterval') && F.SamplingInterval ~= DT
            ok = false;
        end
    end
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 4, ...
        sprintf('all %d rows share n=%d, dt=%g, discrete mode', numel(T), R.n, DT), ...
        sprintf('n ranged over %s', mat2str(unique(ns))));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 4, '', err.message);
end

%% Test 5: the wrapper's retention rule is the fitters' own
try
    rng(5);
    % put mass in [(n_min-1/2)*dt, xmin), which round(t/dt) >= n_min keeps
    % but a naive t >= xmin filter would discard
    d = simGammaDisc(0.7, 500, XMIN, DT, 400);
    d(1:20) = XMIN - DT/2 + 1e-9;
    R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 1, 'Models', {'gamma'}, ...
        'GoFBootstrap', 0);
    H = eng(d, XMIN, "gamma", 'SamplingInterval', DT, 'nStarts', 3, 'Verbose', false);
    ok = R.n == H.n && R.n > nnz(d >= XMIN);
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 5, ...
        sprintf('n=%d matches the fitter and exceeds the naive count %d', ...
            R.n, nnz(d >= XMIN)), ...
        sprintf('wrapper n=%d, fitter n=%d', R.n, H.n));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 5, '', err.message);
end

%% Test 6: powerlaw_cutoff is the gamma relabelled, in the table
try
    rng(6);
    d = simGammaDisc(0.6, 600, XMIN, DT, 400);
    R = cbm(d, XMIN, BASE{:}, 'Models', {'gamma','powerlaw_cutoff'}, ...
        'GoFBootstrap', 0);
    T = tableRows(R.Table);
    ig = find(strcmp({T.Model}, 'gamma'));
    ic = find(strcmp({T.Model}, 'powerlaw_cutoff'));
    dl = abs(T(ig).LogLik - T(ic).LogLik);
    da = abs(T(ig).AICc - T(ic).AICc);
    shape = T(ig).Params(1); alpha = T(ic).Params(1);
    ok = dl < 1e-6 && da < 1e-6 && abs(alpha - (1 - shape)) < 1e-4 && ...
        T(ig).k == T(ic).k;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 6, ...
        sprintf('identical logL (d=%.2g) and alpha=1-shape (%.4g vs %.4g)', ...
            dl, alpha, 1-shape), ...
        sprintf('dlogL=%.3g dAICc=%.3g alpha=%.4g 1-shape=%.4g', dl, da, alpha, 1-shape));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 6, '', err.message);
end

%% Test 7: criterion arithmetic, recomputed from k and n
try
    rng(7);
    d = simHyperDisc([150 1200], [0.6 0.4], XMIN, DT, 500);
    R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 3, 'Models', SMALL, ...
        'GoFBootstrap', 0);
    T = tableRows(R.Table);
    worst = 0;
    for i = 1:numel(T)
        if T(i).Degenerate || ~isfinite(T(i).LogLik), continue; end
        k = T(i).k; n = T(i).n; L = T(i).LogLik;
        aic = 2*k - 2*L;
        if (n - k - 1) > 0, aicc = aic + 2*k*(k+1)/(n-k-1); else, aicc = Inf; end
        bic = k*log(n) - 2*L;
        worst = max([worst, abs(aic - T(i).AIC), abs(aicc - T(i).AICc), ...
                     abs(bic - T(i).BIC)]);
    end
    ok = worst < 1e-9;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 7, ...
        sprintf('AIC, AICc and BIC reproduce from k and n (max err %.2g)', worst), ...
        sprintf('max discrepancy %.3g', worst));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 7, '', err.message);
end

%% Test 8: Akaike weights sum to 1 and the best row has dAICc = 0
try
    rng(8);
    d = simHyperDisc([150 1200], [0.6 0.4], XMIN, DT, 500);
    R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 3, 'Models', SMALL, ...
        'GoFBootstrap', 0);
    T = tableRows(R.Table);
    w = [T.AkaikeWeight]; dA = [T.dAICc]; dB = [T.dBIC];
    live = ~[T.Degenerate] & isfinite([T.AICc]);
    ok = abs(sum(w(live)) - 1) < 1e-12 && all(w(~live) == 0) && ...
        min(dA(live)) == 0 && min(dB(live)) == 0 && all(dA(live) >= -1e-12);
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 8, ...
        sprintf('weights over %d live rows sum to 1, best dAICc = 0', nnz(live)), ...
        sprintf('sum=%.12g, min dAICc=%.3g', sum(w(live)), min(dA(live))));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 8, '', err.message);
end

%% Test 9: ranking follows RankBy
try
    rng(9);
    d = simHyperDisc([150 1200], [0.6 0.4], XMIN, DT, 500);
    bad = {};
    for key = {'AICc', 'BIC'}
        R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 3, 'Models', SMALL, ...
            'GoFBootstrap', 0, 'RankBy', key{1});
        T = tableRows(R.Table);
        live = ~[T.Degenerate] & isfinite([T.AICc]);
        v = [T.(key{1})];
        v = v(live);
        if any(diff(v) < -1e-9)
            bad{end+1} = sprintf('%s not monotone', key{1}); %#ok<AGROW>
        end
        if ~all(live(1:nnz(live)))
            bad{end+1} = sprintf('%s: degenerate row not last', key{1}); %#ok<AGROW>
        end
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 9, ...
        'rows sort by the chosen criterion, degenerate rows last', ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 9, '', err.message);
end

%% Test 10: expected counts sum to n -- the handle and the bin convention
%  agree. This is the one test that would catch an off-by-one grid slip
%  between a fitter's SurvivalHandle and the bin edges the G test uses.
try
    rng(10);
    bad = {};
    for fam = {'gamma','weibull','powerlaw','hyperexponential'}
        d = simGammaDisc(0.7, 500, XMIN, DT, 400);
        R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 2, 'Models', fam, ...
            'GoFBootstrap', 0);
        g = R.GoF;
        if abs(g.ProbabilityMass - 1) > 1e-9 || ...
                abs(sum(g.Expected) - R.n) > 1e-6 || sum(g.Observed) ~= R.n
            bad{end+1} = sprintf('%s: mass=%.12g sumE=%.6f sumO=%d n=%d', ...
                fam{1}, g.ProbabilityMass, sum(g.Expected), sum(g.Observed), R.n); %#ok<AGROW>
        end
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 10, ...
        'bin probabilities sum to 1 and expected counts to n for every family', ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 10, '', err.message);
end

%% Test 11: a true model passes the G test
try
    rng(11);
    d = simGammaDisc(0.7, 500, XMIN, DT, 800);
    R = cbm(d, XMIN, BASE{:}, 'Models', {'gamma'}, 'GoFBootstrap', 0);
    ok = ~isempty(R.Selected) && strcmp(R.Selected.Model, 'gamma') && R.GoF.Passed;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 11, ...
        sprintf('gamma data accepted as gamma: G=%.2f, p in [%.3f, %.3f]', ...
            R.GoF.G, R.GoF.pLower, R.GoF.pUpper), ...
        sprintf('G=%.2f p=[%.4g %.4g] passed=%d', R.GoF.G, R.GoF.pLower, ...
            R.GoF.pUpper, R.GoF.Passed));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 11, '', err.message);
end

%% Test 12: a wrong model is rejected
%  Two widely separated exponential scales are not a power law; at n=1500
%  the G test should say so. This is the test that the machinery can fail
%  a model at all -- test 11 alone would pass with a broken statistic.
try
    rng(12);
    d = simHyperDisc([120 4000], [0.7 0.3], XMIN, DT, 1500);
    R = cbm(d, XMIN, BASE{:}, 'Models', {'powerlaw'}, 'GoFBootstrap', 0);
    ok = isempty(R.Selected) && ~R.GoF.Passed && R.GoF.G > 0;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 12, ...
        sprintf('two-scale mixture rejected as a power law: G=%.1f, p=%.3g', ...
            R.GoF.G, R.GoF.pUpper), ...
        sprintf('G=%.2f pUpper=%.4g passed=%d selected empty=%d', ...
            R.GoF.G, R.GoF.pUpper, R.GoF.Passed, isempty(R.Selected)));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 12, '', err.message);
end

%% Test 13: the walk-down selects the first row that passes
try
    rng(13);
    d = simHyperDisc([120 4000], [0.7 0.3], XMIN, DT, 900);
    R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 2, 'Models', SMALL, ...
        'GoFBootstrap', 0);
    T = tableRows(R.Table);
    acc = find(~cellfun(@isempty, strfind(R.SelectionPath, 'accepted'))); %#ok<STRCL1>
    if isempty(R.Selected)
        ok = isempty(acc) && numel(R.SelectionPath) == numel(T);
        msg = 'nothing passed, path walked the whole table';
    else
        % exactly one acceptance, it is the last path entry, and it names
        % the selected model
        ok = numel(acc) == 1 && acc == numel(R.SelectionPath) && ...
            ~isempty(strfind(R.SelectionPath{acc}, R.Selected.Model)) && ...
            strcmp(R.SelectedFit.Model, R.Selected.Model);
        msg = sprintf('walked %d row(s), accepted %s', numel(R.SelectionPath), ...
            R.Selected.Model);
    end
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 13, msg, ...
        sprintf('path=%d acceptances=%d', numel(R.SelectionPath), numel(acc)));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 13, '', err.message);
end

%% Test 14: degenerate rows are flagged, ranked last, and cannot win
%  Asking for many components on a small sample makes the high orders
%  non-identifiable, which is exactly what the mixture fitter's gating is
%  for. They must still appear, with a reason.
try
    rng(14);
    d = simHyperDisc([200 1500], [0.6 0.4], XMIN, DT, 200);
    R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 6, ...
        'Models', {'hyperexponential'}, 'GoFBootstrap', 0);
    T = tableRows(R.Table);
    dg = [T.Degenerate];
    ok = any(dg) && all(diff(double(dg)) >= 0) && ...
        all(~cellfun(@isempty, {T(dg).Reason}));
    if ~isempty(R.Selected), ok = ok && ~R.Selected.Degenerate; end
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 14, ...
        sprintf('%d of %d orders excluded with reasons, all ranked last', ...
            nnz(dg), numel(T)), ...
        sprintf('degenerate mask %s', mat2str(dg)));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 14, '', err.message);
end

%% Test 15: the Models filter
try
    rng(15);
    d = simGammaDisc(0.7, 500, XMIN, DT, 300);
    R = cbm(d, XMIN, BASE{:}, 'Models', {'gamma','weibull'}, 'GoFBootstrap', 0);
    T = tableRows(R.Table);
    eid = '';
    try
        cbm(d, XMIN, BASE{:}, 'Models', {'lognormal'});
    catch err
        eid = err.identifier;
    end
    ok = numel(T) == 2 && isempty(setdiff({T.Model}, {'gamma','weibull'})) && ...
        strcmp(eid, 'CompareBoutModels:NoModelsSelected');
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 15, ...
        'Models restricts the library and an unknown name errors', ...
        sprintf('%d rows, error id "%s"', numel(T), eid));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 15, '', err.message);
end

%% Test 16: beta enters the library only with an UpperBound
try
    rng(16);
    d = simGammaDisc(0.7, 500, XMIN, DT, 300);
    R1 = cbm(d, XMIN, BASE{:}, 'MaxComponents', 1, 'GoFBootstrap', 0);
    R2 = cbm(d, XMIN, BASE{:}, 'MaxComponents', 1, 'GoFBootstrap', 0, ...
        'UpperBound', 86400);
    m1 = {tableRows(R1.Table).Model}; m2 = {tableRows(R2.Table).Model};
    ok = ~any(strcmp(m1, 'beta')) && any(strcmp(m2, 'beta'));
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 16, ...
        'beta appears only when UpperBound is supplied', ...
        sprintf('without: %d beta rows; with: %d', nnz(strcmp(m1,'beta')), ...
            nnz(strcmp(m2,'beta'))));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 16, '', err.message);
end

%% Test 17: Erlang's integer shape, as the wrapper reports it
try
    rng(17);
    d = simGammaDisc(1.8, 400, XMIN, DT, 300);
    R = cbm(d, XMIN, BASE{:}, 'Models', {'erlang'}, 'GoFBootstrap', 0);
    T = tableRows(R.Table);
    pv = T(1).Params; pse = T(1).ParamSE;
    ok = T(1).k == 2 && T(1).nReported == 2 && pv(1) == round(pv(1)) && ...
        pv(1) >= 1 && isnan(pse(1)) && isfinite(pse(2)) && ...
        isequal(T(1).ParamNames, {'shape','scale'});
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 17, ...
        sprintf('erlang reports integer shape %g with NaN SE, k=2', pv(1)), ...
        sprintf('shape=%g SE=%g k=%d', pv(1), pse(1), T(1).k));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 17, '', err.message);
end

%% Test 18: continuous mode runs, and sits on a different measure
try
    rng(18);
    d = simGammaDisc(0.7, 500, XMIN, DT, 400);
    Rd = cbm(d, XMIN, BASE{:}, 'Models', {'gamma'}, 'GoFBootstrap', 0);
    Rc = cbm(d, XMIN, BASE{:}, 'Models', {'gamma'}, 'GoFBootstrap', 0, ...
        'DistributionType', 'continuous');
    Td = tableRows(Rd.Table); Tc = tableRows(Rc.Table);
    % the shift is the Jacobian of binning, about n*log(dt); check it is of
    % that order rather than equal, since the estimates differ slightly too
    shift = Td(1).LogLik - Tc(1).LogLik;
    ok = isfinite(Tc(1).LogLik) && strcmp(char(Rc.DistributionType), 'continuous') && ...
        abs(shift - Rd.n*log(DT)) < 0.05 * abs(Rd.n*log(DT));
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 18, ...
        sprintf('continuous mode runs; logL shifts by %.1f vs n*log(dt)=%.1f', ...
            shift, Rd.n*log(DT)), ...
        sprintf('shift=%.3f expected about %.3f', shift, Rd.n*log(DT)));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 18, '', err.message);
end

%% Test 19: the parametric bootstrap runs and respects a forced B
try
    rng(19);
    d = simGammaDisc(0.7, 500, XMIN, DT, 300);
    R = cbm(d, XMIN, BASE{:}, 'Models', {'gamma'}, 'GoFBootstrap', 30);
    g = R.GoF;
    ok = isfinite(g.pBootstrap) && g.pBootstrap > 0 && g.pBootstrap <= 1 && ...
        g.B >= 20 && g.B <= 30 && strcmp(g.Basis, 'bootstrap');
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 19, ...
        sprintf('bootstrap p=%.3f from %d valid replicates', g.pBootstrap, g.B), ...
        sprintf('p=%g B=%d basis=%s', g.pBootstrap, g.B, g.Basis));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 19, '', err.message);
end

%% Test 20: RandomSeed makes the whole comparison reproducible
try
    rng(20);
    d = simHyperDisc([150 1200], [0.6 0.4], XMIN, DT, 400);
    A = cbm(d, XMIN, BASE{:}, 'MaxComponents', 2, 'Models', SMALL, ...
        'GoFBootstrap', 0, 'RandomSeed', 99);
    B = cbm(d, XMIN, BASE{:}, 'MaxComponents', 2, 'Models', SMALL, ...
        'GoFBootstrap', 0, 'RandomSeed', 99);
    Ta = tableRows(A.Table); Tb = tableRows(B.Table);
    ok = isequal({Ta.Model}, {Tb.Model}) && ...
        max(abs([Ta.LogLik] - [Tb.LogLik])) < 1e-12;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 20, ...
        'the same RandomSeed reproduces the ranking and every logL', ...
        sprintf('max logL difference %.3g', max(abs([Ta.LogLik] - [Tb.LogLik]))));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 20, '', err.message);
end

%% Test 21: ParamText names every parameter of every fitted row
try
    rng(21);
    d = simHyperDisc([150 1200], [0.6 0.4], XMIN, DT, 400);
    R = cbm(d, XMIN, BASE{:}, 'MaxComponents', 2, 'Models', SMALL, ...
        'GoFBootstrap', 0);
    T = tableRows(R.Table);
    bad = {};
    for i = 1:numel(T)
        if T(i).Degenerate, continue; end
        s = T(i).ParamText;
        if isempty(s) || strcmp(s, '(no fit)')
            bad{end+1} = T(i).Model; continue %#ok<AGROW>
        end
        if strncmp(T(i).Model, 'hyperexp', 8)
            if isempty(strfind(s, 'tau=[')) || isempty(strfind(s, 'q=['))
                bad{end+1} = T(i).Model; %#ok<AGROW>
            end
        else
            for j = 1:numel(T(i).ParamNames)
                if isempty(strfind(s, T(i).ParamNames{j}))
                    bad{end+1} = sprintf('%s/%s', T(i).Model, ...
                        T(i).ParamNames{j}); %#ok<AGROW>
                end
            end
        end
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 21, ...
        'ParamText names every parameter it reports', strjoin(bad, ', '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 21, '', err.message);
end

%% Test 22: Plot=false opens nothing
try
    rng(22);
    d = simGammaDisc(0.7, 500, XMIN, DT, 300);
    before = numel(findall(0, 'Type', 'figure'));
    R = cbm(d, XMIN, BASE{:}, 'Models', {'gamma'}, 'GoFBootstrap', 0);
    after = numel(findall(0, 'Type', 'figure'));
    ok = isempty(R.Figure) && after == before;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 22, ...
        'Plot=false leaves R.Figure empty and opens no figure', ...
        sprintf('figures %d -> %d, R.Figure empty=%d', before, after, isempty(R.Figure)));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 22, '', err.message);
end

%% Test 23: the plot path actually runs
%  Test 22 only checks that Plot=false opens NOTHING, which would pass
%  just as happily if plotFit were broken. This one calls it. plotFit had
%  no coverage at all until this test: the container it was written in has
%  no graphics toolkit, so in MATLAB this is its only exercise.
try
    rng(23);
    d = simGammaDisc(0.7, 500, XMIN, DT, 300);
    hasGraphics = true;
    if exist('OCTAVE_VERSION', 'builtin')
        hasGraphics = ~isempty(available_graphics_toolkits());
    end
    if ~hasGraphics
        fprintf(['[PASS] Test 23: skipped, no graphics toolkit available ' ...
            '(plotFit NOT exercised)\n']);
        nPassed = nPassed + 1;
    else
        R = cbm(d, XMIN, 'SamplingInterval', DT, 'Models', {'gamma'}, ...
            'GoFBootstrap', 0, 'Plot', true, 'Verbose', false);
        % two panels: the survival curve and the Pearson residuals
        nAx = numel(findall(R.Figure, 'Type', 'axes'));
        ok = ~isempty(R.Figure) && ishandle(R.Figure) && nAx >= 2;
        if ~isempty(R.Figure) && ishandle(R.Figure), close(R.Figure); end
        [nPassed, nFailed] = rep(ok, nPassed, nFailed, 23, ...
            sprintf('plotFit drew a figure with %d axes', nAx), ...
            sprintf('R.Figure empty=%d, axes=%d', isempty(R.Figure), nAx));
    end
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 23, '', err.message);
end

%% Test 24: the order LRT wiring
%  Three things: OrderLRT=false leaves R.OrderLRT empty, OrderLRT=true runs
%  it whatever the gap, and what comes back is the LRT's own output for the
%  two best admissible orders. B is tiny here -- this tests the WIRING, not
%  the test itself, which has its own suite.
try
    rng(24);
    d = simHyperDisc([150 1600], [0.6 0.4], XMIN, DT, 500);
    args = {'SamplingInterval', DT, 'Plot', false, 'Verbose', false, ...
            'MaxComponents', 3, 'Models', {'hyperexponential'}, ...
            'GoFBootstrap', 0};
    Roff = cbm(d, XMIN, args{:}, 'OrderLRT', false);
    Ron  = cbm(d, XMIN, args{:}, 'OrderLRT', true, 'OrderLRTReplicates', 29);
    L = Ron.OrderLRT;
    ok = isempty(Roff.OrderLRT) && ~isempty(L) && ...
         L.K0 < L.K1 && isfinite(L.LR) && ...
         L.pValue > 0 && L.pValue <= 1 && L.B == 29 && ...
         numel(L.LRNull) == L.BValid;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 24, ...
        sprintf('OrderLRT off leaves it empty; on runs K=%d vs K=%d, p=%.3g', ...
            L.K0, L.K1, L.pValue), ...
        sprintf('off empty=%d, on empty=%d', isempty(Roff.OrderLRT), isempty(L)));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 24, '', err.message);
end

fprintf('\nSummary: %d passed, %d failed, %d total\n\n', ...
    nPassed, nFailed, nPassed + nFailed);
end

% ======================================================================
%                              helpers
% ======================================================================
function [nP, nF] = rep(ok, nP, nF, idx, passMsg, failMsg)
if ok
    fprintf('[PASS] Test %d: %s\n', idx, passMsg); nP = nP + 1;
else
    fprintf('[FAIL] Test %d: %s\n', idx, failMsg); nF = nF + 1;
end
end

function T = tableRows(Tin)
% R.Table is a table in MATLAB and the underlying struct array in Octave.
% Both support the field/variable access these tests use once turned back
% into a struct array.
if isstruct(Tin)
    T = Tin;
else
    % table2struct unboxes the cell columns rowsToTable boxed, so the
    % parameter fields come back in the same shape Octave's struct array
    % already has: a cellstr of names and a numeric row of values.
    T = table2struct(Tin);
end
T = T(:).';
end

function d = simHyperDisc(tau, w, xmin, dt, n)
% Mixture of exponentials on the dt grid the discrete model assumes, then
% conditioned on clearing xmin. Note ceil, matching the other suites'
% simulators: the fitters map observations with round, so this deliberately
% does NOT assume the two agree.
nmin = max(1, round(xmin/dt)); out = zeros(n,1); filled = 0;
w = w(:).' / sum(w);
cw = cumsum(w);
while filled < n
    m = 4*(n - filled) + 200;
    u = rand(m,1);
    comp = ones(m,1);
    for j = 1:numel(cw)-1, comp = comp + (u > cw(j)); end
    t = -tau(comp)' .* log(rand(m,1));
    nv = ceil(t(:)/dt); nv = nv(nv >= nmin);
    take = min(numel(nv), n - filled);
    out(filled+1:filled+take) = nv(1:take)*dt; filled = filled + take;
end
d = out;
end

function x = gammaDraw(k, theta, n)
% Marsaglia-Tsang, with the standard boost for shape < 1, so the suite does
% not require the Statistics Toolbox.
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
nmin = max(1, round(xmin/dt)); out = zeros(n,1); filled = 0;
while filled < n
    t = gammaDraw(k, theta, 4*(n - filled) + 200);
    nv = ceil(t/dt); nv = nv(nv >= nmin);
    take = min(numel(nv), n - filled);
    out(filled+1:filled+take) = nv(1:take)*dt; filled = filled + take;
end
d = out;
end
