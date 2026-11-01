function test_sqw_stage3()
%TEST_SQW_STAGE3  Stage 3 of the sqw_ref port: the model-free probes, against golden G12 and by stress tests.
%
%   Run by name from the sqw folder:   cd <Icode>/sqw; test_sqw_stage3
%
%   Golden (hmm_audit/golden/G12_planA.json):
%    1. chain2FitPhased log-likelihood (1e-6)
%    2. blockKL (all symbols and wake-only) and blockKLStats (1e-9)
%    3. turEstimate: J, D, B, n for tau = 180, 120, 360, including the NaN at 360 (1e-9)
%    4. turExact, turExactCT and the chain EPR on the true light matrix; turExact on hmmObservedChain (1e-9)
%    5. fillEmptyRows example (exact)
%   Stress:
%    6. chain2FitPhased: complete record = triple-count MLE; structural zeros kept with gaps;
%       a first-order record gives P(c | a, b) independent of a (within sampling error)
%    7. blockKL = a literal loop over blocks (random symbols, gaps, unknown wake, 3 phases, n <= 5)
%    8. blockKL on a long first-order chain: d2 = the chain's exact EPR, d3 = d4 = d2 (bias small)
%    9. turEstimate = a literal loop over blocks (random records with gaps and phase changes)
%   10. turExact: J = sum pi_i T_ij d_ij (exact), = 3 x the cycle current; long record (1e6 min)
%       recovers J and D; continuous-time TUR B <= sigma over 1000 random chains
%   11. hmmObservedChain: stochastic; B = I reproduces the chain; J = the symbol current of a
%       simulated HMM record (statistical)
%   12. fillEmptyRows leaves non-empty rows alone; pair chain simulated through hmmSimulatePhased
%       respects the pair structure
%   13. null calibration: Delta3 p-values (20p null A) on records FROM a first-order chain are uniform

here = fileparts(mfilename('fullpath'));
assert(isequal(fileparts(which('blockKL')), here), 'Run from the sqw folder.');
g = jsondecode(fileread(fullfile(here, '..', 'hmm_audit', 'golden', 'G12_planA.json')));
y = pyToSym(g.y); ph = g.ph + 1;
Ptrue = [0.80 0.15 0.05; 0.06 0.90 0.04; 0.02 0.03 0.95];   % make_golden_planA.py, light
nP = 0; nF = 0;
fprintf('\nRunning stage 3 tests for the sqw port...\n\n');

%% 1. second-order chain
try
    ll = chain2FitPhased(y, ph, 2);
    [nP, nF] = rep(abs(ll - g.mc2_fit.loglik) < 1e-6, nP, nF, 1, ...
        sprintf('G12 chain2FitPhased loglik %.6f (golden %.6f)', ll, g.mc2_fit.loglik), sprintf('diff %.3g', ll - g.mc2_fit.loglik));
catch err, [nP, nF] = rep(false, nP, nF, 1, '', err.message); end

%% 2. block KL
try
    D = blockKL(y, ph, 2); Dq = blockKL(y, ph, 2, QW=true); S = blockKLStats(y, ph, 2);
    fl = {'d2', 'd3', 'd4', 'D3', 'D4', 'D3QW'}; d = zeros(2, 6);
    for p = 1:2, for k = 1:6, d(p, k) = S(p).(fl{k}) - g.block_kl_stats(p).(fl{k}); end, end
    e = [max(abs(D - g.block_kl), [], 'all'), max(abs(Dq - g.block_kl_qw), [], 'all'), max(abs(d(:)))];
    [nP, nF] = rep(all(e < 1e-9), nP, nF, 2, sprintf('G12 blockKL, QW, blockKLStats (max diff %.1g)', max(e)), mat2str(e, 3));
catch err, [nP, nF] = rep(false, nP, nF, 2, '', err.message); end

%% 3. TUR estimate
try
    R = turEstimate(y, ph, 2); e = 0;
    for p = 1:2
        G = g.tur_estimate(p); e = max(e, infDiff(R(p).J, num(G.J)));
        for q = 1:3
            tau = R(p).Taus(q);
            e = max([e, infDiff(R(p).D(q), num(G.(sprintf('D%d', tau)))), infDiff(R(p).B(q), num(G.(sprintf('B%d', tau)))), ...
                infDiff(R(p).n(q), num(G.(sprintf('n%d', tau))))]);
        end
    end
    [nP, nF] = rep(e < 1e-9, nP, nF, 3, sprintf('G12 turEstimate J, D, B, n x 3 taus x 2 phases (max diff %.1g)', e), sprintf('%.3g', e));
catch err, [nP, nF] = rep(false, nP, nF, 3, '', err.message); end

%% 4. exact TUR
try
    [J, Dd, B] = turExact(Ptrue); [Jc, Dc, Bc] = turExactCT(Ptrue);
    en = modelEnergetics(Ptrue);
    [Th, dh] = hmmObservedChain(Ptrue, [0.9 0.1 0; 0.05 0.95 0; 0 0 1]); [Jh, Dh, Bh] = turExact(Th, dh);
    e = [[J Dd B] - g.tur_exact_true_light', [Jc Dc Bc] - g.tur_exact_ct_true_light', en.EPR - g.sigma_true_light, ...
        [Jh Dh Bh] - g.hmm_observed_chain_check'];
    % D (and B = J^2/D) is a second central difference of a log Perron root: a 1e-15 change in T
    % moves it by ~2e-8, so the golden 1e-9 is unattainable for D in a second implementation
    % (QUESTIONS_lClaude Q3). J and sigma at 1e-9; D and B at 1e-6 relative.
    gold = [g.tur_exact_true_light', g.tur_exact_ct_true_light', g.sigma_true_light, g.hmm_observed_chain_check'];
    tolv = 1e-9 * ones(1, 10); tolv([2 3 5 6 9 10]) = 1e-6 * abs(gold([2 3 5 6 9 10]));
    [nP, nF] = rep(all(abs(e) <= tolv), nP, nF, 4, sprintf(['G12 turExact, turExactCT, sigma, hmmObservedChain (J, sigma 1e-9; D, B 1e-6 rel; max rel %.1g); ' ...
        'B/sigma discrete %.4f, continuous %.4f'], max(abs(e ./ max(abs(gold), 1))), B / en.EPR, Bc / en.EPR), mat2str(e, 3));
catch err, [nP, nF] = rep(false, nP, nF, 4, '', err.message); end

%% 5. fillEmptyRows example
try
    out = fillEmptyRows(pyP(g.sim_rows_example.input));
    [nP, nF] = rep(max(abs(out - pyP(g.sim_rows_example.output)), [], 'all') < 1e-15, nP, nF, 5, ...
        'G12 fillEmptyRows example', mat2str(out, 4));
catch err, [nP, nF] = rep(false, nP, nF, 5, '', err.message); end

%% 6. second-order chain stress
try
    rs = RandStream('twister', 'Seed', 6);
    P2 = randDir(rs, 3, 9)';                               % P(c | a, b), rows k = 3(a-1)+b
    T = 6000; yy = zeros(T, 1); yy(1:2) = [1; 2];
    for t = 3:T, yy(t) = find(rand(rs) <= cumsum(P2(3 * (yy(t - 2) - 1) + yy(t - 1), :)), 1); end
    php = ones(T, 1);
    [~, ~, Pc] = chain2FitPhased(yy, php, 1);
    C = zeros(9, 3);
    for t = 3:T, k = 3 * (yy(t - 2) - 1) + yy(t - 1); C(k, yy(t)) = C(k, yy(t)) + 1; end
    Pmle = C ./ sum(C, 2);
    Pfit = zeros(9, 3);
    for k = 1:9, b = mod(k - 1, 3) + 1; Pfit(k, :) = Pc(k, 3 * (b - 1) + (1:3)); end
    % EM also counts the first transition, from the unobserved pair (y0, y1): one fractional
    % count lands in the rows whose second symbol is y1. Other rows must equal the triple MLE.
    nRow = sum(C, 2); r1 = mod((1:9)' - 1, 3) + 1 == yy(1);
    dA = max(abs(Pfit - Pmle), [], 2);
    okA = all(dA(~r1) < 1e-6) && all(dA(r1) <= 1 ./ nRow(r1));
    yg = yy; yg(1000:1040) = 0; yg(3000:3002) = 0; yg(4000) = 4;
    [llg, ~, Pg, info] = chain2FitPhased(yg, php, 1);
    [~, Allow] = pairChainStructure();
    okB = all(Pg(~Allow) == 0) && info.converged && isfinite(llg);
    % first-order data: the a-dependence of P(c | a, b) is sampling noise only
    yf = simChain(rs, Ptrue, ones(20000, 1));
    [~, ~, Pf] = chain2FitPhased(yf, ones(20000, 1), 1);
    dev = 0;
    for b = 1:3
        rows = 3 * ((1:3) - 1) + b; cols = 3 * (b - 1) + (1:3);
        n = arrayfun(@(a) nnz(yf(1:end - 2) == a & yf(2:end - 1) == b), (1:3)');
        M = Pf(rows, cols); se = sqrt(Ptrue(b, :) .* (1 - Ptrue(b, :)) ./ n);
        dev = max(dev, max(abs(M - Ptrue(b, :)) ./ se, [], 'all'));
    end
    okC = dev < 4.5;                                       % max |z| over 27 cells
    [nP, nF] = rep(okA && okB && okC, nP, nF, 6, sprintf(['chain2FitPhased: complete = triple MLE (%.1g), structure kept ' ...
        'with gaps, first-order data max |z| of P(c|a,b) - P(c|b) %.2f'], max(dA(~r1)), dev), ...
        sprintf('MLE %d gaps %d first-order %d', okA, okB, okC));
catch err, [nP, nF] = rep(false, nP, nF, 6, '', err.message); end

%% 7. blockKL against a literal loop
try
    rs = RandStream('twister', 'Seed', 7); e = 0;
    for rep_ = 1:6
        T = 700; yr = randi(rs, 3, T, 1);
        yr(rand(rs, T, 1) < 0.04) = 0; yr(rand(rs, T, 1) < 0.02) = 4;
        phr = ones(T, 1); phr(200:449) = 2; phr(450:end) = 3;
        for qw = [false true]
            D = blockKL(yr, phr, 3, QW=qw, NMax=5);
            e = max(e, max(abs(D - slowBlockKL(yr, phr, 3, qw, 5)), [], 'all'));
        end
    end
    [nP, nF] = rep(e < 1e-12, nP, nF, 7, sprintf('blockKL = literal loop (6 records x 2 alphabets x n <= 5, max diff %.1g)', e), sprintf('%.3g', e));
catch err, [nP, nF] = rep(false, nP, nF, 7, '', err.message); end

%% 8. blockKL on a long first-order chain
try
    rs = RandStream('twister', 'Seed', 8);
    % a strongly driven chain: with Ptrue (EPR 0.02 nats/h) the sampling error of d2 is ~20%
    Pd = [0.80 0.04 0.16; 0.16 0.80 0.04; 0.04 0.16 0.80];   % forward W->Q, Q->S, S->W = 0.16
    T = 400000; yl = simChain(rs, Pd, ones(T, 1));
    S = blockKLStats(yl, ones(T, 1), 1); en = modelEnergetics(Pd);
    ok = abs(S.d2 - en.EPR) < 0.1 * en.EPR && abs(S.D3) < 0.1 * en.EPR && abs(S.D4) < 0.1 * en.EPR;
    [nP, nF] = rep(ok, nP, nF, 8, sprintf('long first-order chain (4e5 min): d2 %.4f vs exact EPR %.4f; D3 %+.4f D4 %+.4f nats/h', ...
        S.d2, en.EPR, S.D3, S.D4), 'd2 differs from the chain EPR or D3/D4 not ~0');
catch err, [nP, nF] = rep(false, nP, nF, 8, '', err.message); end

%% 9. turEstimate against a literal loop
try
    rs = RandStream('twister', 'Seed', 9); e = 0;
    for rep_ = 1:6
        T = 3 * 1440; yr = simChain(rs, Ptrue, ones(T, 1));
        yr(rand(rs, T, 1) < 0.003) = 0; yr(2000:2030) = 0; yr(rand(rs, T, 1) < 0.002) = 4;
        phr = minutePhases(T);
        taus = [60 180 120 360 7];
        R = turEstimate(yr, phr, 2, Taus=taus); Rs = slowTur(yr, phr, 2, taus);
        for p = 1:2
            e = max([e, infDiff(R(p).J, Rs(p).J), max(arrayfun(@infDiff, R(p).D, Rs(p).D)), ...
                max(arrayfun(@infDiff, R(p).B, Rs(p).B)), max(abs(R(p).n - Rs(p).n))]);
        end
    end
    [nP, nF] = rep(e < 1e-9, nP, nF, 9, sprintf('turEstimate = literal loop (6 records, 5 taus, gaps, max diff %.1g)', e), sprintf('%.3g', e));
catch err, [nP, nF] = rep(false, nP, nF, 9, '', err.message); end

%% 10. exact TUR
try
    rs = RandStream('twister', 'Seed', 10);
    [J, Dd] = turExact(Ptrue); pis = stationary_dist(Ptrue); pis = pis(:);
    Jdirect = sum(pis .* Ptrue .* cycleIncrements(), 'all') * 60;
    en = modelEnergetics(Ptrue);
    okA = abs(J - Jdirect) < 1e-6 && abs(J - 3 * en.J) < 1e-6;
    T = 1e6; yl = simChain(rs, Ptrue, ones(T, 1)); R = turEstimate(yl, ones(T, 1), 1);
    seJ = sqrt(2 * Dd / (T / 60)); seD = Dd * sqrt(2 / (R.n(1) - 1));
    okB = abs(R.J - J) < 4 * seJ && abs(R.D(1) - Dd) < 4 * seD;
    worst = 0; worstD = 0;
    for k = 1:1000
        P = randDir(rs, 3, 3)'; P = 0.5 * P + 0.5 * eye(3) * 1.9 * rand(rs); P = P ./ sum(P, 2);
        s = modelEnergetics(P).EPR;
        if isfinite(s) && s > 1e-9
            [~, ~, Bct] = turExactCT(P); [~, ~, Bd] = turExact(P);
            worst = max(worst, Bct / s); worstD = max(worstD, Bd / s);
        end
    end
    okC = worst <= 1 + 1e-6;
    [nP, nF] = rep(okA && okB && okC, nP, nF, 10, sprintf(['turExact J = sum pi T d = 3 x cycle J (%.2g); 1e6-min record J %.4f ' ...
        '(exact %.4f), D %.4f (exact %.4f); max B/sigma over 1000 chains: continuous %.6f, discrete %.4f'], ...
        abs(J - Jdirect), R.J, J, R.D(1), Dd, worst, worstD), sprintf('direct %d long %d TUR %d', okA, okB, okC));
catch err, [nP, nF] = rep(false, nP, nF, 10, '', err.message); end

%% 11. hmmObservedChain
try
    rs = RandStream('twister', 'Seed', 11);
    [Ti, di] = hmmObservedChain(Ptrue, eye(3)); [a1, a2, a3] = turExact(Ti, di); [b1, b2, b3] = turExact(Ptrue);
    okA = max(abs([a1 a2 a3] - [b1 b2 b3])) < 1e-6 * max(abs([b1 b2 b3]));
    Bm = [0.85 0.10 0.05; 0.10 0.85 0.05; 0.02 0.03 0.95];
    [Th, dh] = hmmObservedChain(Ptrue, Bm); okB = max(abs(sum(Th, 2) - 1)) < 1e-12;
    Jh = turExact(Th, dh);
    T = 1e6; [~, ys] = hmmSimulatePhased([1 0 0]', Ptrue, Bm, ones(T, 1), rs);
    R = turEstimate(ys, ones(T, 1), 1);
    [~, Dh] = turExact(Th, dh); seJ = sqrt(2 * Dh / (T / 60));
    okC = abs(R.J - Jh) < 4 * seJ;
    [nP, nF] = rep(okA && okB && okC, nP, nF, 11, sprintf(['hmmObservedChain: B = I reproduces the chain, stochastic; ' ...
        'symbol current exact %.4f vs 1e6-min simulation %.4f (z %.2f)'], Jh, R.J, (R.J - Jh) / seJ), sprintf('B=I %d rows %d sim %d', okA, okB, okC));
catch err, [nP, nF] = rep(false, nP, nF, 11, '', err.message); end

%% 12. fillEmptyRows and pair-chain simulation
try
    rs = RandStream('twister', 'Seed', 12);
    P = cat(3, Ptrue, [0.2 0.8 0; 0 0 0; 0 0.5 0.5]);
    F = fillEmptyRows(P);
    okA = isequal(F(:, :, 1), Ptrue) && isequal(F([1 3], :, 2), P([1 3], :, 2)) && isequal(F(2, :, 2), [1 1 1] / 3);
    [B2, Allow] = pairChainStructure();
    P9 = zeros(9); P9(Allow) = rand(rs, nnz(Allow), 1); P9 = P9 ./ sum(P9, 2); P9(5, :) = 0;
    F9 = fillEmptyRows(P9, Allow);
    okB = isequal(F9(5, :), double(Allow(5, :)) / 3) && isequal(F9([1:4 6:9], :), P9([1:4 6:9], :));
    [z, ys] = hmmSimulatePhased(ones(9, 1) / 9, F9, B2, ones(5000, 1), rs);
    okC = all(mod(z(1:end - 1) - 1, 3) + 1 == floor((z(2:end) - 1) / 3) + 1) && isequal(ys, mod(z - 1, 3) + 1);
    [nP, nF] = rep(okA && okB && okC, nP, nF, 12, 'fillEmptyRows (3-state and pair chain); simulated pair path consistent', ...
        sprintf('3-state %d pair %d simulation %d', okA, okB, okC));
catch err, [nP, nF] = rep(false, nP, nF, 12, '', err.message); end

%% 13. null calibration of the 20p p-value
try
    rs = RandStream('twister', 'Seed', 13);
    nData = 100; nNull = 100; pv = zeros(nData, 1);
    T = 2880; phc = ones(T, 1);
    for k = 1:nData
        yd = simChain(rs, Ptrue, phc); miss = rand(rs, T, 1) < 0.01; yd(miss) = 0;
        [~, pc, Pc] = chain2Or1(yd, phc);
        o = blockKLStats(yd, phc, 1).D3; s = zeros(nNull, 1);
        for r = 1:nNull
            [~, ys] = hmmSimulatePhased(pc, Pc, eye(3), phc, rs); ys(miss) = 0;
            s(r) = blockKLStats(ys, phc, 1).D3;
        end
        pv(k) = (1 + sum(s >= o)) / (1 + nNull);
    end
    f05 = mean(pv < 0.05); f50 = mean(pv < 0.5);
    ok = f05 <= 0.12 && f50 > 0.35 && f50 < 0.65;
    [nP, nF] = rep(ok, nP, nF, 13, sprintf(['null calibration (100 first-order records, 100 null sims each): ' ...
        'p < 0.05 in %.2f, p < 0.5 in %.2f'], f05, f50), 'p-values not uniform under the null');
catch err, [nP, nF] = rep(false, nP, nF, 13, '', err.message); end

fprintf('\n%d passed, %d failed\n', nP, nF);
end

%% ------------------------------------------------------------------ helpers
function [nP, nF] = rep(ok, nP, nF, k, what, why)
if ok
    nP = nP + 1; fprintf('[PASS] %2d  %s\n', k, what);
else
    nF = nF + 1; fprintf('[FAIL] %2d  %s\n        %s\n', k, what, why);
end
end

function y = pyToSym(v), y = v(:) + 1; y(v(:) == -1) = 0; y(v(:) == -2) = 4; end
function P = pyP(v), P = permute(v, [2 3 1]); end
function X = randDir(rs, k, n), g = -log(rand(rs, k, n)); X = g ./ sum(g, 1); end

function v = num(v)
if iscell(v), v = cellfun(@num, v); elseif ischar(v) || isstring(v), v = str2double(v); end
v = double(v);
end

function d = infDiff(a, b)
if isequaln(a, b), d = 0; else, d = abs(a - b); end
if isnan(d), d = Inf; end
end

function y = simChain(rs, P, ph)
T = numel(ph); y = zeros(T, 1); y(1) = randi(rs, 3); u = rand(rs, T, 1); cP = cumsum(P, 2);
for t = 2:T, y(t) = find(u(t) <= cP(y(t - 1), :, ph(t)), 1); end
end

function [ll, pi1, P] = chain2Or1(y, ph)
% the 20p null A chain: first-order fit, empty rows filled for simulation
[ll, pi1, P] = chainFitPhased(y, ph, max(ph)); P = fillEmptyRows(P);
end

function D = slowBlockKL(y, ph, nph, qw, nmax)
% literal transcription of the definition, one block at a time
D = zeros(nph, nmax);
for n = 2:nmax
    nb = 3^n; C = zeros(nb, nph);
    for i = 1:numel(y) - n + 1
        b = y(i:i + n - 1); pp = ph(i:i + n - 1);
        if any(pp ~= pp(1)), continue, end
        if qw, if ~all(b == 1 | b == 2), continue, end
        elseif ~all(b >= 1 & b <= 3), continue, end
        code = sum((b(:)' - 1) .* 3.^(n - 1:-1:0)) + 1;
        C(code, pp(1)) = C(code, pp(1)) + 1;
    end
    for p = 1:nph
        pr = zeros(nb, 1); s = 0;
        for c = 1:nb
            dg = mod(floor((c - 1) ./ 3.^(n - 1:-1:0)), 3);
            if ~qw || all(dg < 2), pr(c) = C(c, p) + 0.5; s = s + pr(c); end
        end
        pr = pr / s; v = 0;
        for c = 1:nb
            if pr(c) == 0, continue, end
            dg = mod(floor((c - 1) ./ 3.^(n - 1:-1:0)), 3);
            rc = sum(fliplr(dg) .* 3.^(n - 1:-1:0)) + 1;
            v = v + pr(c) * log(pr(c) / pr(rc));
        end
        D(p, n) = v;
    end
end
end

function R = slowTur(y, ph, nph, taus)
% literal transcription of sqw_ref.tur_estimate with 1-based indices
T = numel(y); dC = cycleIncrements();
val = false(T, 1); inc = zeros(T, 1);
for t = 2:T
    if y(t) >= 1 && y(t) <= 3 && y(t - 1) >= 1 && y(t - 1) <= 3
        val(t) = true; inc(t) = dC(y(t - 1), y(t));
    end
end
st = 1; en = [];
for t = 2:T, if ph(t) ~= ph(t - 1), en(end + 1) = t - 1; st(end + 1) = t; end, end %#ok<AGROW>
en(end + 1) = T;
for p = 1:nph
    m = ph == p & val; h = nnz(m) / 60; R(p).J = sum(inc(m)) / h; %#ok<AGROW>
    for q = 1:numel(taus)
        tau = taus(q); Nb = [];
        for s = 1:numel(st)
            if ph(st(s)) ~= p, continue, end
            for b0 = st(s):tau:en(s) - tau + 1
                if b0 == 1 || ~all(val(b0:b0 + tau - 1)), continue, end
                Nb(end + 1) = sum(inc(b0:b0 + tau - 1)); %#ok<AGROW>
            end
        end
        R(p).n(q) = numel(Nb);
        if numel(Nb) >= 4, R(p).D(q) = var(Nb) / (2 * tau / 60); else, R(p).D(q) = NaN; end
        if isfinite(R(p).D(q)) && R(p).D(q) > 0, R(p).B(q) = R(p).J^2 / R(p).D(q); else, R(p).B(q) = NaN; end
    end
end
end
