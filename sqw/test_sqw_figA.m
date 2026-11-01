function test_sqw_figA()
%TEST_SQW_FIGA  Plan A figure exports (FIGURE_RECIPES K, Figs A1-A4): stress tests on synthetic inputs.
%
%   Run by name from the sqw folder:   cd <Icode>/sqw; test_sqw_figA
%
%    1. planAStyle: genotype order rule, colours exactly the recipe's hex, light before dark,
%       one block per DD genotype
%    2. planAFlyHMM with B = I and no missing minute: every posterior path IS the record, so
%       J_hmm and EPR_hmm equal the counted J and EPR exactly
%    3. planAFlyHMM with a noisy B: J_hmm (mean over paths) = the HMM's expected count current
%       given the data (forward-backward two-minute marginals), within 4 Monte Carlo SE
%    4. exportFigA1: summary median = median of the exported points per block; fly-phases
%       without a reference absent; HMM median over non-degenerate flies only; wrong genotype
%       order refused; figure and every axes white (theme pinned)
%    5. exportFigA2: noiseless EPR = J^2 / L0 recovers L0 per genotype to 1e-12; only
%       reference + cycle + finite points enter; curves = s J^2
%    6. exportFigA3: only cycle fly-phases; summary medians = point medians
%    7. exportFigA4: medians from unclipped values (a value of 40 moves the median, is drawn
%       at 5); use_180 and hasRef filters; CS DD read from the halves file; filled flags as
%       the recipe; dBIC share above 0 computed, not hard-coded

here = fileparts(mfilename('fullpath'));
assert(isequal(fileparts(which('exportFigA1')), here), 'Run from the sqw folder.');
nP = 0; nF = 0; tmp = tempname; mkdir(tmp); c0 = onCleanup(@() rmdir(tmp, 's'));
fprintf('\nRunning Plan A figure-export tests...\n\n');

%% 1. style
try
    S = planAStyle();
    hex = ["2a78d6" "6250d6" "1baf7a" "eda100" "eb6834" "e87ba4"];
    okC = all(arrayfun(@(k) isequal(round(S.col(k, :) * 255), sscanf(hex(k), '%2x%2x%2x')'), 1:6));
    okO = isequal(S.order, ["CS" "CSDD" "per0LD" "per0DD" "fumin" "sss"]);
    B = S.blocks; okB = height(B) == 10 && isequal(B.g', [1 1 2 3 3 4 5 5 6 6]) && ...
        isequal(B.phase', [1 2 1 1 2 1 1 2 1 2]) && isequal(B.filled', logical([0 1 1 0 1 1 0 1 0 1]));
    [nP, nF] = rep(okC && okO && okB, nP, nF, 1, 'planAStyle: order rule, recipe colours, light before dark, DD one block', ...
        sprintf('colours %d order %d blocks %d', okC, okO, okB));
catch err, [nP, nF] = rep(false, nP, nF, 1, '', err.message); end

%% 2. B = I: posterior paths are the record
try
    rs = RandStream('twister', 'Seed', 2);
    T = 2880; ph = minutePhases(T); P = repmat([0.9 0.06 0.04; 0.05 0.9 0.05; 0.03 0.04 0.93], 1, 1, 2);
    P(:, :, 2) = [0.85 0.05 0.10; 0.10 0.85 0.05; 0.04 0.06 0.90];
    [~, y] = hmmSimulatePhased([1 1 1]' / 3, P, eye(3), ph, rs);
    f = struct('pi', [1 1 1]' / 3, 'P', P, 'B', eye(3), 'deg', false);
    H = planAFlyHMM(y, ph, 2, f, NPaths=20, RandomSeed=3);
    Q = observableQuantities(y, ph, 2, P);
    d = max(abs([[H.J_hmm] - [Q.J], [H.EPR_hmm] - [Q.EPR]]));
    [nP, nF] = rep(d < 1e-12, nP, nF, 2, sprintf('planAFlyHMM, B = I: J_hmm, EPR_hmm = counted values (max diff %.2g)', d), 'differ');
catch err, [nP, nF] = rep(false, nP, nF, 2, '', err.message); end

%% 3. noisy B: path-mean J = exact posterior expected current
try
    Bn = [0.9 0.1 0; 0.15 0.85 0; 0 0 1];
    [~, y] = hmmSimulatePhased([1 1 1]' / 3, P, Bn, ph, rs);
    f = struct('pi', [1 1 1]' / 3, 'P', P, 'B', Bn, 'deg', false);
    nPath = 400; H = planAFlyHMM(y, ph, 2, f, NPaths=nPath, RandomSeed=4);
    Z = hmmFFBS(f.pi, P, Bn, y, ph, NPaths=nPath, RandomSeed=4);
    hrs = accumarray(ph, 1, [2 1]) / 60; Jp = zeros(nPath, 2);
    for k = 1:nPath, Jp(k, :) = jCycleCount(pairCounts(Z(k, :)', ph, 2), hrs)'; end
    Jex = exactCurrent(f.pi, P, Bn, y, ph, hrs);
    z = ([H.J_hmm] - Jex') ./ (std(Jp) / sqrt(nPath));
    [nP, nF] = rep(all(abs(z) < 4), nP, nF, 3, sprintf('planAFlyHMM, noisy B: J_hmm vs exact posterior current, z = %s', mat2str(z, 3)), 'outside 4 SE');
catch err, [nP, nF] = rep(false, nP, nF, 3, '', err.message); end

%% fake six-genotype data for 4-6
D = fakeD(RandStream('twister', 'Seed', 5));

%% 4. Fig A1
try
    X = exportFigA1(D, fullfile(tmp, 'a1'));
    ok = true; why = "";
    for b = 1:height(S.blocks)
        lab = erase(erase(erase(S.blocks.label(b), "\it"), "\rm"), "^");
        pts = X.points(X.points.block == lab, :);
        for q = ["c_obs" "J" "EPR"]
            r = X.summary(X.summary.quantity == q & X.summary.label == lab, :);
            if abs(r.median - median(pts.(q), 'omitnan')) > 1e-12, ok = false; why = why + lab + " " + q + "; "; end
        end
        T = D.G(S.blocks.g(b)).perFly; T = T(T.phase == S.blocks.phase(b), :);
        if height(pts) ~= nnz(T.refOK), ok = false; why = why + lab + " refOK; "; end
        nd = T.refOK & ~T.hmm_deg; r = X.summary(X.summary.quantity == "J" & X.summary.label == lab, :);
        if any(nd) && abs(r.hmmMedian - median(T.J_hmm(nd))) > 1e-12, ok = false; why = why + lab + " hmm; "; end
    end
    th = X.fig.Color; ax = findobj(X.fig, 'Type', 'axes');
    okTheme = isequal(th, [1 1 1]) && all(arrayfun(@(a) isequal(a.Color, [1 1 1]), ax));
    close(X.fig);
    D2 = D; D2.order = D.order([2 1 3:6]); refused = false;
    try exportFigA1(D2, fullfile(tmp, 'a1b'), Plot=false); catch, refused = true; end
    [nP, nF] = rep(ok && okTheme && refused, nP, nF, 4, 'exportFigA1: medians = point medians, refOK, HMM over non-degenerate, order, light theme', ...
        sprintf('%s theme %d order refused %d', why, okTheme, refused));
catch err, [nP, nF] = rep(false, nP, nF, 4, '', err.message); end

%% 5. Fig A2
try
    X = exportFigA2(D, fullfile(tmp, 'a2'), NBoot=50); close(X.fig);
    L0 = [0.27 0.26 0.29 0.35 0.24 0.15];
    dL = max(abs(X.L.L' - L0));
    okPts = true;
    for k = 1:6
        T = D.G(k).perFly; m = T.refOK & T.hasCycle & isfinite(T.J) & isfinite(T.EPR);
        okPts = okPts && nnz(X.points.genotype == S.order(k)) == nnz(m);
        c = X.curves(X.curves.genotype == S.order(k), :);
        okPts = okPts && max(abs(c.EPR_fit - c.J .^ 2 / L0(k))) < 1e-9;
    end
    [nP, nF] = rep(dL < 1e-12 && okPts, nP, nF, 5, sprintf('exportFigA2: L recovered (max err %.2g), point filter, curves', dL), 'mismatch');
catch err, [nP, nF] = rep(false, nP, nF, 5, '', err.message); end

%% 6. Fig A3
try
    X = exportFigA3(D, fullfile(tmp, 'a3')); close(X.fig);
    ok = true;
    for b = 1:height(S.blocks)
        lab = erase(erase(erase(S.blocks.label(b), "\it"), "\rm"), "^");
        T = D.G(S.blocks.g(b)).perFly; T = T(T.phase == S.blocks.phase(b) & T.refOK & T.hasCycle, :);
        pts = X.points(X.points.block == lab, :);
        r = X.summary(X.summary.quantity == "E31" & X.summary.label == lab, :);
        ok = ok && height(pts) == height(T) && abs(r.median - median(T.E31)) < 1e-12;
    end
    [nP, nF] = rep(ok, nP, nF, 6, 'exportFigA3: cycle fly-phases only, medians = point medians', 'mismatch');
catch err, [nP, nF] = rep(false, nP, nF, 6, '', err.message); end

%% 7. Fig A4
try
    rs = RandStream('twister', 'Seed', 7); pf = fullfile(tmp, 'pf'); mkdir(pf);
    [L, truth] = fakeA4(rs, pf);
    X = exportFigA4(L, pf, fullfile(tmp, 'a4')); close(X.fig);
    r = X.rows;
    okB = abs(r.median(r.panel == "b" & r.group == "CS" & r.what == "R_mem") - truth.RmemMedCS) < 1e-12;
    okA = abs(r.fracAbove0or1(r.panel == "a" & r.group == "fumin") - truth.fracFumin) < 1e-12;
    okC = abs(r.median(r.panel == "c" & r.group == "CS DD, ZT 12-24") - truth.csddHalf2) < 1e-12 && ...
        r.n(r.panel == "c" & r.group == "CS light") == truth.nCSlight;
    okD = abs(r.median(r.panel == "d" & r.group == "per0 DD") - truth.per0DDratio) < 1e-12 && ...
        r.n(r.panel == "d" & r.group == "per0 DD") == truth.nPer0DD;
    fc = r.group(r.panel == "c" & r.filled)'; fd = r.group(r.panel == "d" & r.filled)';
    okF = isequal(sort(fc), sort(["per0 DD" "fmn dark"])) && isequal(sort(fd), sort(["CS DD" "per0 DD"]));
    [nP, nF] = rep(okA && okB && okC && okD && okF, nP, nF, 7, ...
        'exportFigA4: unclipped medians, dBIC share, halves file, use_180 / hasRef, filled flags', ...
        sprintf('a %d b %d c %d d %d filled %d', okA, okB, okC, okD, okF));
catch err, [nP, nF] = rep(false, nP, nF, 7, '', err.message); end

fprintf('\n%d passed, %d failed\n', nP, nF);
end

% ------------------------------------------------------------------------
function J = exactCurrent(pi0, P, B, y, ph, hrs)
% Posterior expected J from forward-backward two-minute marginals (not from sampled paths).
T = numel(y); E = hmmEmission(B, y);
A = zeros(T, 3); c = zeros(T, 1); a = pi0(:)' .* E(1, :); c(1) = sum(a); A(1, :) = a / c(1);
for t = 2:T, a = (A(t - 1, :) * P(:, :, ph(t))) .* E(t, :); c(t) = sum(a); A(t, :) = a / c(t); end
Bk = ones(T, 3);
for t = T - 1:-1:1, Bk(t, :) = (P(:, :, ph(t + 1)) * (E(t + 1, :) .* Bk(t + 1, :))')' / c(t + 1); end
N = zeros(3, 3, 2);
for t = 2:T
    X = (A(t - 1, :)' * (E(t, :) .* Bk(t, :))) .* P(:, :, ph(t)) / c(t);
    X(1:4:end) = 0; N(:, :, ph(t)) = N(:, :, ph(t)) + X;
end
J = jCycleCount(N, hrs);
end

function D = fakeD(rs)
S = planAStyle(); D.order = S.order; L0 = [0.27 0.26 0.29 0.35 0.24 0.15];
for k = 1:6
    nFly = 8 + k; nPh = 1 + ~S.dd(k); rows = cell(nFly, 1);
    for i = 1:nFly
        p = (1:nPh)'; n = nPh;
        J = 0.05 + 0.3 * rand(rs, n, 1);
        nm = ["light"; "dark"];
        t = table(repmat(S.order(k), n, 1), repmat(i, n, 1), repmat("f" + i + ".txt", n, 1), p, ...
            nm(1:n), 'VariableNames', {'genotype', 'fly', 'file', 'phase', 'phaseName'});
        if S.dd(k), t.phaseName(:) = "whole"; end
        t.refOK = rand(rs, n, 1) > 0.15; t.hasCycle = rand(rs, n, 1) > 0.1;
        t.c_obs = randn(rs, n, 1) * 0.1; t.J = J; t.EPR = J .^ 2 / L0(k); t.rule_of_three = rand(rs, n, 1) > 0.5;
        t.hmm_deg = rand(rs, n, 1) > 0.8; t.J_hmm = J + 0.02 * randn(rs, n, 1); t.EPR_hmm = t.EPR * 1.5;
        t.E21 = -0.5 + 0.2 * randn(rs, n, 1); t.E31 = -1 + 0.5 * randn(rs, n, 1);
        t.F_SW = 2 * rand(rs, n, 1) - 1; t.Xn_cnt = t.F_SW + 0.05 * randn(rs, n, 1); t.F_WS = 2 * rand(rs, n, 1) - 1;
        t.Yn_cnt = t.F_WS; t.Yn_cnt(1) = NaN;
        rows{i} = t;
    end
    G = struct('genotype', S.order(k), 'perFly', vertcat(rows{:}));
    D.G(k) = G;
end
end

function [L, truth] = fakeA4(rs, pf)
g = [repmat("CS", 12, 1); repmat("fumin", 15, 1)]; n = numel(g);
L = table(g, (0:n - 1)', "f" + (1:n)' + ".txt", repmat(5760, n, 1), false(n, 1), 'VariableNames', {'genotype', 'fly0', 'file', 'T', 'deg'});
L.MC1_12 = -3000 + 10 * randn(rs, n, 1); L.HMM_12 = L.MC1_12 + 20 * rand(rs, n, 1); L.MC2_12 = L.MC1_12 + 30 * rand(rs, n, 1);
L.MC1_6 = L.MC1_12 + 10 * rand(rs, n, 1); L.HMM_6 = L.MC1_6 + 15 * rand(rs, n, 1);
L.MC2_12(2) = L.MC1_12(2) + 40 * (L.HMM_12(2) - L.MC1_12(2));        % R_mem = 40 for fly 2: drawn at 5
L.HMM_12(1) = L.MC1_12(1) + 50;
lnT = log(L.T); L.dBIC_T = (-2 * L.MC1_12 + 14 * lnT) - (-2 * L.HMM_12 + 20 * lnT);
L.gain = L.HMM_12 - L.MC1_12; L.R_mem = (L.MC2_12 - L.MC1_12) ./ L.gain; L.R_drift = 1 - (L.HMM_6 - L.MC1_6) ./ L.gain;
m = L.genotype == "CS" & L.gain > 3; truth.RmemMedCS = median(L.R_mem(m));
truth.fracFumin = mean(L.dBIC_T(L.genotype == "fumin") > 0);
for gg = ["CS" "per0LD" "per0DD" "fumin" "sss" "CSDD"]
    dd = ismember(gg, ["per0DD" "CSDD"]); nf = 10;
    if dd, phs = "whole"; else, phs = ["light"; "dark"]; end
    [fl, pp] = ndgrid(1:nf, 1:numel(phs)); fl = fl(:); ph = phs(pp(:));
    M = table(repmat(gg, numel(fl), 1), "f" + fl + ".txt", fl, ph, double(rand(rs, numel(fl), 1) > 0.2), randn(rs, numel(fl), 1) * 0.05, ...
        'VariableNames', {'genotype', 'file', 'fly', 'phase', 'hasRef', 'E_A'});
    writetable(M, fullfile(pf, gg + "_memory_perfly.csv"));
    if gg == "CS", truth.nCSlight = nnz(M.phase == "light" & M.hasRef == 1); end
    Q = table(repmat(gg, numel(fl), 1), "f" + fl + ".txt", fl, ph, 0.5 + rand(rs, numel(fl), 1), 0.1 + rand(rs, numel(fl), 1) * 0.1, ...
        double(rand(rs, numel(fl), 1) > 0.3), 'VariableNames', {'genotype', 'file', 'fly', 'phase', 'B_c_180', 'sig_mc', 'use_180'});
    Q.B_c_180 = Q.B_c_180 .* Q.sig_mc;
    writetable(Q, fullfile(pf, gg + "_tur_perfly.csv"));
    if gg == "per0DD"
        u = Q.use_180 == 1; truth.per0DDratio = median(Q.B_c_180(u) ./ Q.sig_mc(u)); truth.nPer0DD = nnz(u);
    end
    if gg == "CSDD"
        hp = ["ZT0-12"; "ZT12-24"]; [fl2, pp2] = ndgrid(1:nf, 1:2); fl2 = fl2(:);
        Hh = table(repmat(gg, numel(fl2), 1), "f" + fl2 + ".txt", fl2, hp(pp2(:)), double(rand(rs, numel(fl2), 1) > 0.2), ...
            randn(rs, numel(fl2), 1) * 0.05, 'VariableNames', {'genotype', 'file', 'fly', 'phase', 'hasRef', 'E_A2'});
        writetable(Hh, fullfile(pf, "CSDD_memory_halves_perfly.csv"));
        k = Hh.phase == "ZT12-24" & Hh.hasRef == 1; truth.csddHalf2 = median(Hh.E_A2(k));
    end
end
end

function [nP, nF] = rep(ok, nP, nF, k, what, why)
if ok
    nP = nP + 1; fprintf('[PASS] %2d  %s\n', k, what);
else
    nF = nF + 1; fprintf('[FAIL] %2d  %s\n        %s\n', k, what, why);
end
end
