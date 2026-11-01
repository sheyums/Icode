function test_sqw_stage1()
%TEST_SQW_STAGE1  Stage 1 of the sqw_ref port: the observable core, against the golden files and by stress tests.
%
%   Run by name from the sqw folder (Icode CLAUDE.md, Testing):
%
%       cd <Icode>/sqw; test_sqw_stage1
%
%   Golden tests (hmm_audit/golden, made by sqw_ref.py; 0-BASED, converted here):
%    1. G1  sleepBouts, bridged and unbridged, every bridging case
%    2. G2  sqwSymbols: y, refs, bouts, number of wake events
%    3. G3  hmmForward log-likelihood and one hmmEmStep (relative 1e-10)
%    4. G9  hmmEmStep with unknown-wake minutes (relative 1e-10)
%    5. G4  routingPredicted, stationary_dist, modelEnergetics (1e-12; Inf kept)
%    6. G5  routeEvents, routeCounts, pair counts, xyCounts, EPR bound, J, cObsMinutes
%    7. G7  observablePairCounts on the G2 symbols (exact)
%    8. G10 onsagerFit, eprJModels (1e-12); onsagerBootstrap median within 5% of L
%    9. G11 energyLandscape (1e-12)
%   10. G12 chainFitPhased (loglik 1e-6, P 1e-9), routeCountsObserved, observableQuantities (1e-9)
%   Stress tests (no golden; checked against brute force or an identity):
%   11. sleepBouts against a literal loop over the bridging rule, 300 random traces
%   12. observablePairCounts / cObsMinutes / routeCountsObserved against loops, random sequences with gaps
%   13. chainFitPhased on complete data = count MLE; with gaps the likelihood never decreases
%   14. modelEnergetics: sign(J) = sign(A) = sign(c), EPR >= 0, detailed balance gives J = EPR = 0,
%       J equal on all three edges, one-way edge gives Inf
%   15. jCycleCount = epr_count_based's J_cycle; EPR bound >= 0 on random counts
%   16. onsagerFit recovers L exactly on noiseless points; tQuantile = tinv (if present); quantile7 = numpy
%   17. readCounts: round trip, and refuses a header, a negative and a fraction
%   18. sqwSymbols: S minutes = >= 30 s asleep; "unknown" codes 4 exactly where "missing" codes 0
%   19. Argument validation: wrong lengths, bad symbols, non-stochastic P are errors

here = fileparts(mfilename('fullpath'));
assert(isequal(fileparts(which('chainFitPhased')), here), ...
    'chainFitPhased resolves to %s: run from the sqw folder.', which('chainFitPhased'));
assert(isequal(fileparts(which('epr_count_based')), here), ...
    'epr_count_based resolves to %s, not the sqw copy.', which('epr_count_based'));
gd = fullfile(here, '..', 'hmm_audit', 'golden');
G = @(f) jsondecode(fileread(fullfile(gd, f)));
nP = 0; nF = 0;
fprintf('\nRunning stage 1 tests for the sqw port...\n\n');

%% 1. G1 sleep bouts
try
    g = G('G1_sleep_bouts.json'); ok = true; bad = '';
    for f = fieldnames(g)'
        c = g.(f{1});
        e1 = c.bouts_bridged; if ~isempty(e1), e1(:, 1) = e1(:, 1) + 1; else, e1 = zeros(0, 2); end
        e0 = c.bouts_unbridged; if ~isempty(e0), e0(:, 1) = e0(:, 1) + 1; else, e0 = zeros(0, 2); end
        b1 = sleepBouts(c.x); b0 = sleepBouts(c.x, Bridge=false);
        if ~isequal(b1, reshape(e1, [], 2)) || ~isequal(b0, reshape(e0, [], 2))
            ok = false; bad = [bad ' ' f{1}]; %#ok<AGROW>
        end
    end
    [nP, nF] = rep(ok, nP, nF, 1, 'G1 sleep bouts, all 8 bridging cases', ['failed:' bad]);
catch err, [nP, nF] = rep(false, nP, nF, 1, '', err.message); end

%% 2. G2 symbols
try
    g = G('G2_symbols.json');
    [y, refs, info] = sqwSymbols(g.x, Theta=g.theta);
    yE = pyToSym(g.y); bE = g.bouts; bE(:, 1) = bE(:, 1) + 1;
    ok = isequal(y, yE) && abs(refs.light - g.refs.light) < 1e-15 && ...
        abs(refs.dark - g.refs.dark) < 1e-15 && isequal(info.bouts, bE) && ...
        numel(info.events) == g.n_wake_events;
    [nP, nF] = rep(ok, nP, nF, 2, sprintf('G2 symbols: %d minutes, %d bouts, %d events, refs %.4f / %.4f', ...
        numel(y), size(bE, 1), g.n_wake_events, g.refs.light, g.refs.dark), ...
        sprintf('y mismatches %d, bouts equal %d, events %d vs %d', nnz(y ~= yE), ...
        isequal(info.bouts, bE), numel(info.events), g.n_wake_events));
catch err, [nP, nF] = rep(false, nP, nF, 2, '', err.message); end

%% 3. G3 forward and one EM step
try
    g = G('G3_hmm.json');
    [y, ph, P] = deal(pyToSym(g.y), g.ph + 1, pyP(g.P));
    [~, ~, ll] = hmmForward(g.pi, P, g.B, y, ph);
    [p1, P1, B1, llIn] = hmmEmStep(g.pi, P, g.B, y, ph);
    e = g.em_step;
    r = [relErr(ll, g.loglik), relErr(llIn, e.loglik_input), maxRel(p1, e.pi), ...
        maxRel(P1, pyP(e.P)), maxRel(B1, e.B)];
    [nP, nF] = rep(all(r < 1e-10), nP, nF, 3, sprintf('G3 loglik %.6f and one EM step (max rel err %.1e)', ...
        ll, max(r)), sprintf('relative errors %s', mat2str(r, 3)));
catch err, [nP, nF] = rep(false, nP, nF, 3, '', err.message); end

%% 4. G9 unknown-wake minutes
try
    g = G('G9_wake_unknown.json'); g3 = G('G3_hmm.json');
    y = pyToSym(g.y); ph = g.ph + 1;
    [~, ~, ll] = hmmForward(g3.pi, pyP(g3.P), g3.B, y, ph);
    [p1, P1, B1, llIn] = hmmEmStep(g3.pi, pyP(g3.P), g3.B, y, ph);
    e = g.em_step;
    r = [relErr(ll, g.loglik), relErr(llIn, e.loglik_input), maxRel(p1, e.pi), ...
        maxRel(P1, pyP(e.P)), maxRel(B1, e.B)];
    [nP, nF] = rep(nnz(y == 4) > 0 && all(r < 1e-10), nP, nF, 4, ...
        sprintf('G9 unknown-wake minutes (%d coded 4): loglik and EM step, max rel err %.1e', nnz(y == 4), max(r)), ...
        sprintf('relative errors %s', mat2str(r, 3)));
catch err, [nP, nF] = rep(false, nP, nF, 4, '', err.message); end

%% 5. G4 matrices
try
    g = G('G4_matrices.json'); ok = true; msg = '';
    for f = fieldnames(g.matrices)'
        m = g.matrices.(f{1});
        [xn, yn, c] = routingPredicted(m.P);
        E = modelEnergetics(m.P);
        pi0 = stationary_dist(m.P);
        d = [xn - m.xn, yn - m.yn, c - m.c, max(abs(pi0 - m.pi)), E.J - num(m.J), ...
            infDiff(E.EPR, num(m.EPR)), infDiff(E.A, num(m.A))];
        if any(abs(d) > 1e-12), ok = false; msg = sprintf('%s %s %s;', msg, f{1}, mat2str(d, 3)); end
    end
    [nP, nF] = rep(ok, nP, nF, 5, 'G4 routing, stationary distribution, J / EPR / A (one-way edge: EPR Inf, A -Inf)', msg);
catch err, [nP, nF] = rep(false, nP, nF, 5, '', err.message); end

%% 6. G5 counting
try
    g = G('G5_counts.json');
    z = g.z + 1; ph = g.ph + 1;
    [r, t] = routeEvents(z);
    names = {'N312', 'N32', 'N213', 'N23'};
    eN = cellfun(@(c) c{1}, g.route_events, 'UniformOutput', false);
    eT = cellfun(@(c) c{2}, g.route_events) + 1;
    ok1 = isequal(names(r)', eN) && isequal(t, eT);
    ok2 = isequal(routeCounts(z, ph, 2), g.route_counts);
    ok3 = isequal(observablePairCounts(z, ph, 2), pyP(g.pair_counts));
    [X, Y] = xyCounts([5 3 2 6; 1 1 4 0]);          % make_golden.py line 62
    ok4 = isequaln(X, num(g.xy_counts{1}(:))) && isequaln(Y, num(g.xy_counts{2}(:)));
    ok5 = true; msg5 = '';
    for f = fieldnames(g.epr_bound)'
        c = g.epr_bound.(f{1});
        [e, Jc] = epr_count_based(c.N, c.hours * 3600, 'bound');
        if isnan(e), e = 0; end                       % all-zero counts: no edge contributed
        Jj = jCycleCount(c.N, c.hours);
        if abs(e - c.epr) > 1e-12 || abs(Jj - c.J) > 1e-12 || abs(Jc - c.J) > 1e-12
            ok5 = false; msg5 = sprintf('%s %s (EPR %.15g vs %.15g)', msg5, f{1}, e, c.epr);
        end
    end
    co = cObsMinutes(pyToSym(g.c_obs_minutes.y), g.c_obs_minutes.ph + 1, 2);
    ok6 = max(abs(co - g.c_obs_minutes.value)) < 1e-15;
    [nP, nF] = rep(ok1 && ok2 && ok3 && ok4 && ok5 && ok6, nP, nF, 6, ...
        'G5 route events / counts, pair counts, Xn Yn, EPR bound and J (4 cases), c_obs', ...
        sprintf('events %d counts %d pairs %d xy %d epr %d%s cobs %d', ok1, ok2, ok3, ok4, ok5, msg5, ok6));
catch err, [nP, nF] = rep(false, nP, nF, 6, '', err.message); end

%% 7. G7 observable pair counts on the G2 symbols
try
    g = G('G7_roundtrips.json'); g2 = G('G2_symbols.json');
    y = pyToSym(g2.y);
    N = observablePairCounts(y, minutePhases(numel(y)), 2);
    [nP, nF] = rep(isequal(N, pyP(g.observable_pair_counts.value)), nP, nF, 7, ...
        sprintf('G7 observable pair counts on G2 symbols (%d changes)', sum(N(:))), 'counts differ');
catch err, [nP, nF] = rep(false, nP, nF, 7, '', err.message); end

%% 8. G10 Onsager and EPR-J models
try
    g = G('G10_onsager.json');
    F = onsagerFit(g.J, g.EPR); e = g.onsager_fit;
    d1 = [F.s - e.s, F.se_s - e.se_s, F.L - e.L, F.se_L - e.se_L, F.ci_L - e.ci_L', F.n - e.n, F.r2 - e.r2];
    M = eprJModels(g.J, g.EPR); e = g.epr_j_models;
    d2 = [M.a - e.a, M.b_linear - e.b_linear, M.rss_linear - e.rss_linear, M.c - e.c, ...
        M.b_quadratic - e.b_quadratic, M.rss_quadratic - e.rss_quadratic, M.exponent - e.exponent, ...
        M.n_points - e.n_points, M.n_exponent - e.n_exponent];
    Lb = onsagerBootstrap((1:12)', g.J, g.EPR, NBoot=2000, RandomSeed=20);
    d3 = abs(median(Lb) / F.L - 1);
    [nP, nF] = rep(all(abs(d1) < 1e-12) && all(abs(d2) < 1e-12) && d3 < 0.05, nP, nF, 8, ...
        sprintf('G10 Onsager L %.6f, CI [%.4f %.4f]; models; bootstrap median off by %.2f%%', ...
        F.L, F.ci_L, 100 * d3), sprintf('fit %s | models %s | boot %.3f', mat2str(d1, 2), mat2str(d2, 2), d3));
catch err, [nP, nF] = rep(false, nP, nF, 8, '', err.message); end

%% 9. G11 landscape
try
    g = G('G11_route_duration.json'); e = g.landscape;
    L = energyLandscape(e.P);
    fl = {'E2_minus_E1', 'E3_minus_E1', 'b12', 'b13', 'E31c', 'E21c', 'E32', 'E23', 'F_SW', 'F_WS'};
    d = cellfun(@(f) L.(f) - e.(f), fl);
    [xn, yn] = routingPredicted(e.P);
    [nP, nF] = rep(all(abs(d) < 1e-12) && abs(xn - L.F_SW) < 1e-12 && abs(yn - L.F_WS) < 1e-12, nP, nF, 9, ...
        'G11 landscape energies, barriers, F_SW = Xn_pred, F_WS = Yn_pred', mat2str(d, 3));
catch err, [nP, nF] = rep(false, nP, nF, 9, '', err.message); end

%% 10. G12 chain fit and observable quantities
try
    g = G('G12_planA.json');
    y = pyToSym(g.y); ph = g.ph + 1;
    [ll, ~, P] = chainFitPhased(y, ph, 2);
    okA = abs(ll - g.mc_fit.loglik) < 1e-6 && max(abs(P - pyP(g.mc_fit.P)), [], 'all') < 1e-9;
    okB = isequal(routeCountsObserved(y, ph, 2), g.routes_observed);
    Q = observableQuantities(y, ph, 2, pyP(g.mc_fit.P));
    fl = {'c_obs', 'J', 'EPR', 'rule_of_three', 'J_chain', 'A_chain', 'EPR_chain', 'E21', 'E31', ...
        'F_SW', 'F_WS', 'Xn_cnt', 'Yn_cnt'};
    d = zeros(2, numel(fl));
    for p = 1:2
        for k = 1:numel(fl)
            d(p, k) = infDiff(double(Q(p).(fl{k})), num(g.observable_quantities(p).(fl{k})));
        end
    end
    okC = all(abs(d(:)) < 1e-9);
    [nP, nF] = rep(okA && okB && okC, nP, nF, 10, ...
        sprintf('G12 chain loglik %.6f (golden %.6f), routes, 13 observable quantities x 2 phases', ll, g.mc_fit.loglik), ...
        sprintf('chain %d routes %d quantities %d (max diff %.2g)', okA, okB, okC, max(abs(d(:)))));
catch err, [nP, nF] = rep(false, nP, nF, 10, '', err.message); end

%% 11. sleepBouts against a literal loop
try
    rs = RandStream('twister', 'Seed', 11); ok = true;
    for it = 1:300
        x = makeTrace(rs);
        ok = ok && isequal(sleepBouts(x), slowBouts(x, true, 300)) && ...
            isequal(sleepBouts(x, Bridge=false), slowBouts(x, false, 300)) && ...
            isequal(sleepBouts(x, MinLen=60), slowBouts(x, true, 60));
        if ~ok, break, end
    end
    [nP, nF] = rep(ok, nP, nF, 11, 'sleepBouts = literal loop over the bridging rule (300 random traces, 3 settings)', ...
        sprintf('first mismatch at trace %d', it));
catch err, [nP, nF] = rep(false, nP, nF, 11, '', err.message); end

%% 12. Counting against loops on sequences with gaps
try
    rs = RandStream('twister', 'Seed', 12); ok = true;
    for it = 1:200
        T = 50 + randi(rs, 400); np = randi(rs, 4);
        y = randsample3(rs, T); ph = sort(randi(rs, np, T, 1));
        [N1, c1, R1] = slowCounts(y, ph, np);
        ok = ok && isequal(observablePairCounts(y, ph, np), N1) && ...
            isequaln(cObsMinutes(y, ph, np), c1) && isequal(routeCountsObserved(y, ph, np), R1);
        yc = y; yc(yc == 0 | yc == 4) = 2;
        ok = ok && isequal(routeCountsObserved(yc, ph, np), routeCounts(yc, ph, np));
        if ~ok, break, end
    end
    [nP, nF] = rep(ok, nP, nF, 12, 'pair counts, c_obs, observed routes = loops (200 sequences with gaps); = routeCounts when complete', ...
        sprintf('mismatch at sequence %d', it));
catch err, [nP, nF] = rep(false, nP, nF, 12, '', err.message); end

%% 13. Chain fit: count MLE on complete data, monotone EM with gaps
try
    rs = RandStream('twister', 'Seed', 13);
    Ptrue = cat(3, [0.8 0.15 0.05; 0.06 0.9 0.04; 0.02 0.03 0.95], [0.7 0.2 0.1; 0.1 0.85 0.05; 0.03 0.05 0.92]);
    ph = minutePhases(2880); y = simChain(rs, Ptrue, ph);
    [~, ~, P, info] = chainFitPhased(y, ph, 2);
    t = (2:numel(y))';
    Nc = accumarray([y(t - 1), y(t), ph(t)], 1, [3 3 2]); Pmle = Nc ./ sum(Nc, 2);
    okA = max(abs(P - Pmle), [], 'all') < 1e-10;
    yg = y; yg(500:530) = 0; yg(2000:2003) = 0;
    idx = (1200:1210)'; yg(idx(y(idx) ~= 3)) = 4;      % wake minutes of unknown type; S stays observed
    [pi1, P1] = deal(ones(3, 1) / 3, info.startP);
    lls = zeros(30, 1);
    for k = 1:30
        [pi1, P1, ~, lls(k)] = hmmEmStep(pi1, P1, eye(3), yg, ph);
    end
    okB = all(diff(lls) >= -1e-9 * abs(lls(1:end - 1)));
    [ll2, ~, ~, info2] = chainFitPhased(yg, ph, 2);
    okC = info2.converged && ll2 >= lls(end) - 1e-9;
    [nP, nF] = rep(okA && okB && okC, nP, nF, 13, sprintf(['chain fit = count MLE on complete data (%d it); ' ...
        'EM monotone with gaps and unknown wake; converged in %d it'], info.iterations, info2.iterations), ...
        sprintf('MLE %d monotone %d converged %d', okA, okB, okC));
catch err, [nP, nF] = rep(false, nP, nF, 13, '', err.message); end

%% 14. Energetics identities
try
    rs = RandStream('twister', 'Seed', 14); ok = true; nNo = 0;
    for it = 1:2000
        P = rand(rs, 3) .^ 3; P = P ./ sum(P, 2);
        E = modelEnergetics(P); [~, ~, c] = routingPredicted(P); p = E.pi;
        J12 = 60 * (p(2) * P(2, 1) - p(1) * P(1, 2));
        J13 = 60 * (p(1) * P(1, 3) - p(3) * P(3, 1));
        J32 = 60 * (p(3) * P(3, 2) - p(2) * P(2, 3));
        if abs(E.J) < 1e-12, nNo = nNo + 1; continue, end
        if E.oneWay, ok = ok && isinf(E.EPR); continue, end
        ok = ok && sign(E.J) == sign(E.A) && sign(E.J) == sign(c) && E.EPR >= 0 && ...
            abs(J12 - J13) < 1e-12 && abs(J13 - J32) < 1e-12 && ...
            E.EPR >= E.J * E.A - 1e-9;                 % EPR >= J A for a single cycle (equality)
        ok = ok && abs(E.EPR - E.J * E.A) < 1e-9 * max(1, E.EPR);
    end
    % detailed balance: P from symmetric weights
    Wt = rand(rs, 3); Wt = Wt + Wt'; Pdb = Wt ./ sum(Wt, 2); Edb = modelEnergetics(Pdb);
    okDB = abs(Edb.J) < 1e-13 && abs(Edb.EPR) < 1e-12 && abs(Edb.A) < 1e-12;
    Pow = [0.7 0.3 0; 0.1 0.85 0.05; 0.05 0 0.95]; Eow = modelEnergetics(Pow);
    okOW = isinf(Eow.EPR) && Eow.EPR > 0 && Eow.oneWay;
    [nP, nF] = rep(ok && okDB && okOW, nP, nF, 14, ...
        'energetics on 2000 random matrices: sign(J)=sign(A)=sign(c), J equal on all edges, EPR = J A >= 0; detailed balance 0; one-way Inf', ...
        sprintf('identities %d detailed balance %d one-way %d', ok, okDB, okOW));
catch err, [nP, nF] = rep(false, nP, nF, 14, '', err.message); end

%% 15. Count J and the EPR bound on random counts
try
    rs = RandStream('twister', 'Seed', 15); ok = true;
    for it = 1:5000
        N = randi(rs, [0 12], 3) .* (rand(rs, 3) > 0.3); N(1:4:9) = 0; h = 1 + 20 * rand(rs);
        [e, Jc] = epr_count_based(N, h * 3600, 'bound');
        ok = ok && abs(jCycleCount(N, h) - Jc) < 1e-12 && (isnan(e) || e >= 0);
    end
    [nP, nF] = rep(ok, nP, nF, 15, 'jCycleCount = epr_count_based J_cycle, EPR bound >= 0 (5000 random count matrices)', 'violation');
catch err, [nP, nF] = rep(false, nP, nF, 15, '', err.message); end

%% 16. Onsager exactness, t quantile, type-7 quantile
try
    J = [-0.3 -0.1 0.05 0.2 0.4 0.6]'; F = onsagerFit(J, J .^ 2 / 0.25);
    okA = abs(F.L - 0.25) < 1e-14 && F.se_s < 1e-12 && abs(F.r2 - 1) < 1e-14;
    nu = [1 2 3 5 10 30 100]; pp = [0.6 0.9 0.975 0.995];
    T7 = tQuantile(pp', nu);
    okB = true;
    if exist('tinv', 'file') == 2
        [PP, NU] = ndgrid(pp, nu); okB = max(abs(T7 - tinv(PP, NU)), [], 'all') < 1e-10;
    end
    okB = okB && abs(tQuantile(0.975, 1) - tan(pi * 0.475)) < 1e-10 && ...      % nu = 1: Cauchy
        abs(tQuantile(0.9, 2) - 0.8 / sqrt(2 * 0.9 * 0.1)) < 1e-12;               % nu = 2: (2p-1)/sqrt(2p(1-p))
    okC = abs(quantile7((1:10)', 0.25) - 3.25) < 1e-15 && abs(quantile7([3 1 2]', 0.5) - 2) < 1e-15 && ...
        isequal(quantile7([5 1]', [0 1]), [1 5]);
    [nP, nF] = rep(okA && okB && okC, nP, nF, 16, 'Onsager exact on noiseless points; tQuantile = tinv / scipy; quantile7 = numpy', ...
        sprintf('onsager %d tq %d q7 %d', okA, okB, okC));
catch err, [nP, nF] = rep(false, nP, nF, 16, '', err.message); end

%% 17. readCounts
try
    f = [tempname '.txt']; c1 = onCleanup(@() delete(f));
    x = [0 3 0 0 12 1]'; writematrix(x, f); okA = isequal(readCounts(f), x);
    okB = threw(@() readCounts(writeText(f, "count\n1\n2\n")), 'readCounts:parse') && ...
        threw(@() readCounts(writeText(f, "1\n-2\n")), 'readCounts:values') && ...
        threw(@() readCounts(writeText(f, "1\n2.5\n")), 'readCounts:values');
    [nP, nF] = rep(okA && okB, nP, nF, 17, 'readCounts round trip; header, negative and fraction refused', ...
        sprintf('roundtrip %d refusals %d', okA, okB));
catch err, [nP, nF] = rep(false, nP, nF, 17, '', err.message); end

%% 18. sqwSymbols invariants
try
    g = G('G2_symbols.json');
    [y, ~, info] = sqwSymbols(g.x, FullDays=false);
    okA = isequal(y == 3, info.sleepSeconds >= 30);
    rs = RandStream('twister', 'Seed', 18);
    x = double(rand(rs, 2 * 86400, 1) < 0.15);          % restless: hardly any sleep, no references
    x(30000:40000) = 0; x(60000:70000) = 0;              % two long bouts only
    [ym, rm] = sqwSymbols(x); yu = sqwSymbols(x, UnknownWake="unknown");
    okB = isequal(ym == 0, yu == 4) && isequal(ym(ym ~= 0), yu(yu ~= 4)) && isnan(rm.light) && any(ym == 0);
    [nP, nF] = rep(okA && okB, nP, nF, 18, 'sqwSymbols: S = >= 30 s asleep; "unknown" codes 4 exactly where "missing" codes 0', ...
        sprintf('S rule %d unknown %d', okA, okB));
catch err, [nP, nF] = rep(false, nP, nF, 18, '', err.message); end

%% 19. Validation
try
    ok = threw(@() observablePairCounts([1 2 3]', [1 1]', 1), 'sqw:length') && ...
        threw(@() observablePairCounts([1 5 3]', [1 1 1]', 1), 'sqw:symbols') && ...
        threw(@() modelEnergetics([0.5 0.5 0; 0 1 0; 0 0 0.9]), 'sqw:notStochastic') && ...
        threw(@() onsagerFit([1 NaN]', [1 1]'), 'onsagerFit:finite') && ...
        threw(@() routeEvents([1 0 2]'), 'routeEvents:states');
    [nP, nF] = rep(ok, nP, nF, 19, 'bad lengths, symbols, matrices and NaN points are refused', 'a bad input passed');
catch err, [nP, nF] = rep(false, nP, nF, 19, '', err.message); end

%% 20. Drivers: sqwObservableFly = the functions it chains; group table, summary, CSV
try
    g = G('G2_symbols.json');
    F = sqwObservableFly(g.x, SixHour=true);
    y = pyToSym(g.y); ph = minutePhases(numel(y));
    [~, ~, P] = chainFitPhased(y, ph, 2); Q = observableQuantities(y, ph, 2, P);
    okA = isequal(F.y, y) && isequal(F.chain.P, P) && isequaln(rmfield(F.Q, {'refOK', 'hasCycle'}), Q) && ...
        numel(F.six.Q) == 4 && all(F.refOK);
    d = [tempname '_sqw']; mkdir(d); c2 = onCleanup(@() rmdir(d, 's'));
    sub = fullfile(d, 'sub'); mkdir(sub);
    writematrix(g.x, fullfile(d, 'b.txt')); writematrix(g.x, fullfile(sub, 'a.txt'));
    writematrix(g.x(end:-1:1), fullfile(d, 'c.txt'));
    Gr = sqwObservableGroup(string(d), "test", OutDir=string(d), NBoot=200);
    okB = numel(Gr.files) == 3 && height(Gr.perFly) == 6 && isequal(Gr.summary.n, [3; 3]) && ...
        isfile(fullfile(d, 'test_observable_perfly.csv')) && isfile(fullfile(d, 'test_observable_summary.csv'));
    C = readtable(fullfile(d, 'test_observable_perfly.csv'));
    okC = height(C) == 6 && all(ismember({'J', 'EPR', 'N_SW', 'N_WS', 'N312'}, C.Properties.VariableNames)) && ...
        max(abs(C.J - Gr.perFly.J)) < 1e-12;
    okD = abs(Gr.summary.J_med(1) - median(Gr.perFly.J(Gr.perFly.phase == 1))) < 1e-15;
    [nP, nF] = rep(okA && okB && okC && okD, nP, nF, 20, ...
        'drivers: fly = chained functions; group finds files recursively, summary medians, CSV round trip', ...
        sprintf('fly %d group %d csv %d median %d', okA, okB, okC, okD));
catch err, [nP, nF] = rep(false, nP, nF, 20, '', err.message); end

%% 21. Onsager group: Wald / BH / bootstrap difference on synthetic genotypes with known L
try
    rs = RandStream('twister', 'Seed', 21); Gs = struct('genotype', {}, 'perFly', {});
    Ls = [0.25 0.25 0.15];
    for k = 1:3
        nfl = 40; J = 0.05 + 0.4 * rand(rs, 2 * nfl, 1);
        T = table(repelem((1:nfl)', 2), J, J .^ 2 / Ls(k) .* (1 + 0.05 * randn(rs, 2 * nfl, 1)), ...
            'VariableNames', {'fly', 'J', 'EPR'});
        T.J_chain = T.J; T.EPR_chain = T.EPR; T.refOK = true(height(T), 1); T.hasCycle = T.refOK;
        T.file = "f" + T.fly; T.phaseName = repmat(["light"; "dark"], nfl, 1);
        Gs(k).genotype = "g" + k; Gs(k).perFly = T;
    end
    O = sqwOnsagerGroup(Gs, Reference="g1", NBoot=500);
    ok = all(abs(O.L - Ls') ./ Ls' < 0.03) && O.diff_p(2) > 0.05 && O.diff_p(3) < 0.01 && ...
        O.wald_p(3) < 0.01 && isnan(O.wald_p(1)) && all(O.L_lo < O.L & O.L < O.L_hi);
    [nP, nF] = rep(ok, nP, nF, 21, sprintf('Onsager group recovers L %s; equal L not different, lower L detected', ...
        mat2str(round(O.L', 3))), sprintf('L %s diff p %s', mat2str(O.L', 3), mat2str(O.diff_p', 3)));
catch err, [nP, nF] = rep(false, nP, nF, 21, '', err.message); end

fprintf('\n%d passed, %d failed\n', nP, nF);
if nF > 0, error('test_sqw_stage1:failed', '%d test(s) failed', nF); end
end

% ----------------------------------------------------------------------- helpers
function [nP, nF] = rep(ok, nP, nF, k, what, why)
if ok
    nP = nP + 1; fprintf('[PASS] %2d  %s\n', k, what);
else
    nF = nF + 1; fprintf('[FAIL] %2d  %s\n        %s\n', k, what, why);
end
end

function y = pyToSym(v)
% Python Q 0, W 1, S 2, missing -1, unknown -2  ->  MATLAB 1, 2, 3, 0, 4
y = v(:) + 1; y(v(:) == -1) = 0; y(v(:) == -2) = 4;
end

function P = pyP(v)
% jsondecode gives Python's (nP x K x K) as a [nP K K] array; MATLAB uses K x K x nP.
P = permute(v, [2 3 1]);
end

function v = num(v)
% golden non-finite values arrive as the strings "NaN", "Inf", "-Inf"
if iscell(v), v = cellfun(@num, v); elseif ischar(v) || isstring(v), v = str2double(v); end
v = double(v);
end

function d = infDiff(a, b)
if isequaln(a, b), d = 0; else, d = a - b; end
if isnan(d), d = Inf; end
end

function r = relErr(a, b), r = abs(a - b) / max(abs(b), realmin); end
function r = maxRel(a, b), r = max(abs(a(:) - b(:)) ./ max(abs(b(:)), 1e-300 + max(abs(b(:))) * 1e-12)); end

function tf = threw(f, id)
tf = false;
try, f(); catch e, tf = strcmp(e.identifier, id); end
end

function f = writeText(f, s)
fid = fopen(f, 'w'); fprintf(fid, s); fclose(fid);
end

function x = makeTrace(rs)
% Random trace built from inactive runs near the 300-s threshold and short active spells.
x = zeros(0, 1);
while numel(x) < 4000
    if rand(rs) < 0.5
        L = randi(rs, [250 360]);
    else
        L = randi(rs, 60);
    end
    x = [x; zeros(L, 1)]; %#ok<AGROW>
    k = randi(rs, 3);
    x = [x; randi(rs, 5, k, 1)]; %#ok<AGROW>
end
if rand(rs) < 0.5, x = x(randi(rs, 50):end); end
end

function b = slowBouts(x, bridge, minLen)
a = x > 0; n = numel(a);
r0 = zeros(0, 3); s = 1;                      % runs of the UNBRIDGED trace: start, length, isRest
while s <= n
    e = s; while e < n && a(e + 1) == a(s), e = e + 1; end
    r0(end + 1, :) = [s, e - s + 1, ~a(s)]; %#ok<AGROW>
    s = e + 1;
end
if bridge
    for i = 2:size(r0, 1) - 1
        if ~r0(i, 3) && r0(i, 2) == 1 && r0(i - 1, 3) && r0(i + 1, 3) && ...
                r0(i - 1, 2) >= minLen && r0(i + 1, 2) >= minLen
            a(r0(i, 1)) = false;
        end
    end
end
b = zeros(0, 2); s = 1;
while s <= n
    e = s; while e < n && a(e + 1) == a(s), e = e + 1; end
    if ~a(s) && e - s + 1 >= minLen, b(end + 1, :) = [s, e - s + 1]; end %#ok<AGROW>
    s = e + 1;
end
end

function y = randsample3(rs, T)
% Sticky random symbols with gaps (0) and unknown-wake minutes (4).
y = zeros(T, 1); y(1) = randi(rs, 3);
for t = 2:T
    if rand(rs) < 0.7, y(t) = y(t - 1); else, y(t) = randi(rs, 3); end
    if y(t) == 0 || y(t) == 4, y(t) = randi(rs, 3); end
end
g = rand(rs, T, 1) < 0.04; y(g) = 0;
u = rand(rs, T, 1) < 0.03 & y ~= 3; y(u) = 4;
end

function [N, c, R] = slowCounts(y, ph, np)
N = zeros(3, 3, np); on = zeros(np, 2); wk = zeros(np, 2);
for t = 2:numel(y)
    a = y(t - 1); b = y(t);
    if a >= 1 && a <= 3 && b >= 1 && b <= 3 && a ~= b, N(a, b, ph(t)) = N(a, b, ph(t)) + 1; end
    if (a == 1 || a == 2) && b == 3, on(ph(t), :) = on(ph(t), :) + [a == 1, 1]; end
    if a == 3 && (b == 1 || b == 2), wk(ph(t), :) = wk(ph(t), :) + [b == 1, 1]; end
end
c = NaN(np, 1);
for p = 1:np
    if on(p, 2) > 0 && wk(p, 2) > 0, c(p) = on(p, 1) / on(p, 2) - wk(p, 1) / wk(p, 2); end
end
R = zeros(np, 4); t = 1; T = numel(y);
while t <= T                                  % each run of observed minutes on its own
    if y(t) < 1 || y(t) > 3, t = t + 1; continue, end
    u = t; while u < T && y(u + 1) >= 1 && y(u + 1) <= 3, u = u + 1; end
    seg = y(t:u); st = []; at = [];
    for k = 1:numel(seg)
        if k == 1 || seg(k) ~= seg(k - 1), st(end + 1) = seg(k); at(end + 1) = t + k - 1; end %#ok<AGROW>
    end
    for k = 2:numel(st)
        r = 0;
        if st(k - 1) == 3 && st(k) == 2, r = 2;
        elseif st(k - 1) == 2 && st(k) == 3, r = 4;
        elseif k >= 3 && st(k - 1) == 1 && st(k - 2) == 3 && st(k) == 2, r = 1;
        elseif k >= 3 && st(k - 1) == 1 && st(k - 2) == 2 && st(k) == 3, r = 3;
        end
        if r, R(ph(at(k)), r) = R(ph(at(k)), r) + 1; end
    end
    t = u + 1;
end
end

function y = simChain(rs, P, ph)
T = numel(ph); y = zeros(T, 1); y(1) = randi(rs, 3);
for t = 2:T
    y(t) = find(rand(rs) <= cumsum(P(y(t - 1), :, ph(t))), 1);
end
end
