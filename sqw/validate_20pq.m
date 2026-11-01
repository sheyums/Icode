function V = validate_20pq(dataRoot, jsonlFile, outDir, options)
%VALIDATE_20PQ  Stage-3 check of the port: reproduce the PRESPEC 20p and 20q ledger lines from raw data.
%
%   V = validate_20pq(dataRoot, jsonlFile, outDir)
%   V = validate_20pq(..., NBootMC=200, Workers=[])
%
%   dataRoot   folders CS, per0LD, per0DD, fumin, sss (all 57), CSDD
%   jsonlFile  hmm_audit/sqw20b_fit2_theta0.5.jsonl (cClaude's labelled 20b fits; CS DD is
%              not in it and is fitted here as 20p did: K = 3, 20 starts, cached in outDir)
%   outDir     results: <g>_memory_*.csv, <g>_tur_*.csv, validate_20pq.mat, compare CSVs
%
%   Every ledger number is compared with the port's:
%     exact   deterministic given the data (d_n, Delta3_obs, n, J, D, B_obs, sigma_MC, ...):
%             must agree to the ledger's printed precision
%     mc      depends on simulated records (E, fractions, S1/S2, predictions): z =
%             (port - ledger) / sqrt(2 SE^2 + rounding^2 / 3), SE the port's Monte Carlo
%             SE from resampling each fly's own simulated records (NBootMC times).
%             Python's SE is the same (same numbers of records), hence the 2.
%     p       p-values (Wilcoxon, binomial, Holm): reported beside the ledger, not judged
%             one by one (they inherit Monte Carlo error and, for Wilcoxon with ties,
%             MATLAB's exact method differs from scipy's); the pre-specified VERDICTS
%             (20p QUALIFIES, 20q S1/S2 support) must agree.
%
%   See also SQWMEMORYGROUP, SQWTURGROUP, CHECK20BEXACT.

arguments
    dataRoot (1,1) string {mustBeFolder}
    jsonlFile (1,1) string {mustBeFile}
    outDir (1,1) string
    options.NBootMC (1,1) double {mustBePositive, mustBeInteger} = 200
end
if ~isfolder(outDir), mkdir(outDir); end
geno = ["CS" "per0LD" "per0DD" "fumin" "sss" "CSDD"]; ld = [true true false true true false];
S = readStored20bFits(jsonlFile);

M = struct(); Q = struct();
for k = 1:numel(geno)
    g = geno(k); files = sqwListFiles(fullfile(dataRoot, g));
    fits = hmmFitsFor(g, files, ld(k), S, outDir);
    fprintf('%s: %d flies ... ', g, numel(files)); t0 = tic;
    M.(g) = sqwMemoryGroup(files, g, LightDark=ld(k), HMMFits=fits, OutDir=outDir, Plot=false, RandomSeed=15000 + 1000 * k);
    Q.(g) = sqwTURGroup(files, g, LightDark=ld(k), HMMFits=fits, OutDir=outDir, Plot=false, RandomSeed=16000 + 1000 * k);
    fprintf('%.0f s\n', toc(t0));
end

% ------------------------------------------------------------ 20p
blocks = ["CS" "light"; "CS" "dark"; "per0LD" "light"; "per0LD" "dark"; "fumin" "light"; "fumin" "dark"; ...
    "sss" "light"; "sss" "dark"; "per0DD" "whole"; "CSDD" "whole"];
P = table();
for b = 1:size(blocks, 1)
    s = M.(blocks(b, 1)).summary; r = s(s.phase == blocks(b, 2), :);
    P = [P; [table(blocks(b, 1), 'VariableNames', {'g'}), r]]; %#ok<AGROW>
end
P.holm_A = holmAdjust(P.pW_A); P.holm_B = holmAdjust(P.pW_B);
P.qualifies = P.holm_A < 0.05 & P.E_A > 0 & P.pbin_A < 0.05 & P.holm_B < 0.05 & P.E_B > 0 & P.pbin_B < 0.05;
se = mcSE20p(M, blocks, options.NBootMC);
H = table();
for g = ["per0DD" "CSDD"]
    H = [H; [table(repmat(g, 2, 1), 'VariableNames', {'g'}), M.(g).halvesSummary]]; %#ok<AGROW>
end
H.holm_A = holmAdjust(H.pW_A); H.holm_B = holmAdjust(H.pW_B);
H.qualifies = H.holm_A < 0.05 & H.E_A > 0 & H.pbin_A < 0.05 & H.holm_B < 0.05 & H.E_B > 0 & H.pbin_B < 0.05;

% ------------------------------------------------------------ 20q
qb = ["per0DD" "whole"; "CSDD" "whole"; "CS" "light"; "CS" "dark"; "per0LD" "light"; "per0LD" "dark"; ...
    "fumin" "light"; "fumin" "dark"; "sss" "light"; "sss" "dark"];
R = table();
for b = 1:size(qb, 1)
    s = Q.(qb(b, 1)).summary; r = s(s.phase == qb(b, 2), :);
    R = [R; [table(qb(b, 1), 'VariableNames', {'g'}), r]]; %#ok<AGROW>
end
pr = 1:2;
R.S1_holm = NaN(height(R), 1); R.S2_holm = NaN(height(R), 1);
R.S1_holm(pr) = holmAdjust(R.S1_p(pr)); R.S2_holm(pr) = holmAdjust(R.S2_p(pr));
R.S1_supported = R.S1_holm < 0.05 & R.S1_med > 1 & R.S1_bin < 0.05;
R.S2_supported = R.S2_holm < 0.05 & R.S2_med > 1;
seq = mcSE20q(Q, qb, options.NBootMC);

% ------------------------------------------------------------ compare with the ledger
C = [compare(P, ledger20p(), se, "20p"); compare(H, ledger20pHalves(), [], "20p halves"); compare(R, ledger20q(), seq, "20q")];
writetable(P, fullfile(outDir, 'port_20p_blocks.csv')); writetable(H, fullfile(outDir, 'port_20p_halves.csv'));
writetable(R, fullfile(outDir, 'port_20q_blocks.csv')); writetable(C, fullfile(outDir, 'compare_20pq.csv'));
V = struct('p20', P, 'halves', H, 'q20', R, 'compare', C);
save(fullfile(outDir, 'validate_20pq.mat'), 'V', 'M', 'Q', '-v7.3');

ex = C(C.kind == "exact", :); mc = C(C.kind == "mc" & isfinite(C.z), :);
fprintf('\nEXACT (to printed precision): %d of %d agree\n', nnz(ex.ok), height(ex));
disp(ex(~ex.ok, {'analysis', 'g', 'phase', 'stat', 'ledger', 'port', 'diff'}));
fprintf('MONTE CARLO: %d statistics, |z| > 2 in %.3f (expect ~0.05), |z| > 3 in %d, max |z| %.2f\n', ...
    height(mc), mean(abs(mc.z) > 2), nnz(abs(mc.z) > 3), max(abs(mc.z)));
disp(mc(abs(mc.z) > 2.5, {'analysis', 'g', 'phase', 'stat', 'ledger', 'port', 'se', 'z'}));
fprintf('VERDICTS: 20p qualifies (port / ledger):\n');
disp([P(:, {'g', 'phase'}), table(P.qualifies, ledger20p().qualifies, 'VariableNames', {'port', 'ledger'})]);
fprintf('20q S1, S2 supported, primary blocks (ledger: none):\n'); disp(R(pr, {'g', 'phase', 'S1_supported', 'S2_supported'}));
end

% ===================================================================== HMM fits
function fits = hmmFitsFor(g, files, lightDark, S, outDir)
if g ~= "CSDD"
    fits = repmat(struct('pi', [], 'P', [], 'B', [], 'deg', false), numel(files), 1);
    for k = 1:numel(files)
        s = S([S.g] == g & [S.fly0] == k - 1);
        if numel(s) ~= 1, error('validate:fits', 'No stored fit for %s fly %d.', g, k - 1); end
        fits(k) = struct('pi', s.pi, 'P', s.P, 'B', s.B, 'deg', s.deg);
    end
    return
end
cache = fullfile(outDir, 'CSDD_hmm_fits.mat');
if isfile(cache), L = load(cache, 'fits'); fits = L.fits; return, end
n = numel(files); c = cell(n, 1);
parfor k = 1:n                               % 20p: K = 3, 20 starts, seed 20 + 1000*5 + fly (0-based)
    [y, ph, nP] = sqwFlyInput(files(k), lightDark);
    f = hmmFitPhased(y, ph, nP, 3, Starts=20, RandomSeed=20 + 5000 + k - 1);
    [P, B, pi0] = labelStates(f.P, f.B, f.pi);
    c{k} = struct('pi', pi0, 'P', P, 'B', B, 'deg', isDegenerateH1(B), 'll', f.ll);
end
fits = vertcat(c{:});
save(cache, 'fits');
end

% ===================================================================== Monte Carlo SE
function se = mcSE20p(M, blocks, nBoot)
% resample each fly's own simulated records; recompute the group statistics
stats = ["E_A" "frac_A" "E_B" "frac_B" "E_D4" "E_D3QW" "pred_HMM" "in_HMM" "pred_MC2" "in_MC2"];
se = table();
rs = RandStream('twister', 'Seed', 99);
for b = 1:size(blocks, 1)
    G = M.(blocks(b, 1)); p = find(unique(G.perFly.phase, 'stable') == blocks(b, 2));
    keep = G.perFly.hasRef(G.perFly.phase == blocks(b, 2));
    F = G.flies(keep); nf = numel(F); V = NaN(nBoot, numel(stats));
    for r = 1:nBoot
        x = NaN(nf, numel(stats));
        for f = 1:nf
            o = F(f).obs(p); s = F(f).sims;
            ia = randi(rs, size(s.A_D3, 1), size(s.A_D3, 1), 1); ib = randi(rs, size(s.B_D3, 1), size(s.B_D3, 1), 1);
            A3 = s.A_D3(ia, p); B3 = s.B_D3(ib, p);
            x(f, 1) = o.D3 - median(A3); x(f, 2) = (1 + sum(A3 >= o.D3)) / (1 + numel(A3)) < 0.05;
            x(f, 3) = o.D3 - median(B3); x(f, 4) = (1 + sum(B3 >= o.D3)) / (1 + numel(B3)) < 0.05;
            x(f, 5) = o.D4 - median(s.A_D4(ia, p)); x(f, 6) = o.D3QW - median(s.A_D3QW(ia, p));
            [x(f, 7), x(f, 8)] = pred(s.HMM_D3(:, p), o.D3, rs); [x(f, 9), x(f, 10)] = pred(s.MC2_D3(:, p), o.D3, rs);
        end
        V(r, :) = [median(x(:, 1)), mean(x(:, 2)), median(x(:, 3)), mean(x(:, 4)), median(x(:, 5)), median(x(:, 6)), ...
            median(x(:, 7), 'omitnan'), mean(x(:, 8), 'omitnan'), median(x(:, 9), 'omitnan'), mean(x(:, 10), 'omitnan')];
    end
    for j = 1:numel(stats)
        se = [se; table(blocks(b, 1), blocks(b, 2), stats(j), std(V(:, j)), 'VariableNames', {'g', 'phase', 'stat', 'se'})]; %#ok<AGROW>
    end
end
end

function [m, in] = pred(s, o, rs)
if isempty(s), m = NaN; in = NaN; return, end
s = s(randi(rs, numel(s), numel(s), 1)); m = median(s); q = quantile7(s, [0.025 0.975]); in = q(1) <= o && o <= q(2);
end

function se = mcSE20q(Q, qb, nBoot)
stats = ["S1_med" "S1_frac" "S2_med" "S2_frac" "Q_obs" "HMM_enough" "HMM_inside" "Q_hmm"];
se = table(); rs = RandStream('twister', 'Seed', 98);
for b = 1:size(qb, 1)
    G = Q.(qb(b, 1)); p = find(unique(G.perFly.phase, 'stable') == qb(b, 2));
    T = G.perFly(G.perFly.phase == qb(b, 2), :); use = T.use_180;
    F = G.flies(use); T = T(use, :); nf = numel(F); V = NaN(nBoot, numel(stats));
    for r = 1:nBoot
        x = NaN(nf, 7);
        for f = 1:nf
            Bo = F(f).obs(p).B(1); Bs = F(f).simMC(:, p, 1); Bs = Bs(isfinite(Bs)); Bs = Bs(randi(rs, numel(Bs), numel(Bs), 1));
            Bc = Bo * T.B_mc(f) / median(Bs) / T.kappa(f);
            Bh = F(f).simH(:, p, 1); Bh = Bh(isfinite(Bh)); inH = NaN;
            if numel(Bh) >= 50 && T.deg(f) ~= 1
                Bh = Bh(randi(rs, numel(Bh), numel(Bh), 1)); q = quantile7(Bh, [0.025 0.975]); inH = q(1) <= Bo && Bo <= q(2);
            end
            sh = T.sig_hmm(f); qh = 0; if isfinite(sh) && sh > 0, qh = Bc / sh; end
            x(f, :) = [Bo / median(Bs), (1 + sum(Bs >= Bo)) / (1 + numel(Bs)) < 0.05, Bc / T.sig_mc(f), sh >= Bc, inH, qh, T.deg(f) == 1];
        end
        nd = x(:, 7) == 0;
        V(r, :) = [median(x(:, 1)), mean(x(:, 2)), median(x(:, 3)), mean(x(:, 3) > 1), median(x(:, 3)), ...
            mean(x(nd, 4)), mean(x(nd, 5), 'omitnan'), median(x(nd, 6))];
    end
    for j = 1:numel(stats)
        se = [se; table(qb(b, 1), qb(b, 2), stats(j), std(V(:, j)), 'VariableNames', {'g', 'phase', 'stat', 'se'})]; %#ok<AGROW>
    end
end
end

% ===================================================================== comparison
function C = compare(port, L, se, name)
% L: ledger table with g, phase and one column per statistic; L.Properties.UserData has kind and decimals
spec = L.Properties.UserData;
C = table();
for i = 1:height(L)
    r = port(port.g == L.g(i) & port.phase == L.phase(i), :);
    for j = 1:height(spec)
        st = spec.stat(j); if ~ismember(st, L.Properties.VariableNames), continue, end
        lv = L.(st)(i); if isnan(lv), continue, end
        pv = r.(st); d = pv - lv; rnd = 0.5 * 10^(-spec.dec(j));
        c = table(name, L.g(i), L.phase(i), st, spec.kind(j), lv, pv, d, NaN, NaN, false, ...
            'VariableNames', {'analysis', 'g', 'phase', 'stat', 'kind', 'ledger', 'port', 'diff', 'se', 'z', 'ok'});
        switch spec.kind(j)
            case "exact"
                c.ok = abs(d) <= rnd * (1 + 1e-9) + 1e-12;
            case "mc"
                s = NaN;
                if ~isempty(se)
                    k = se.g == L.g(i) & se.phase == L.phase(i) & se.stat == st;
                    if any(k), s = se.se(k); end
                end
                c.se = s; c.z = d / sqrt(2 * s^2 + rnd^2 / 3); c.ok = abs(c.z) <= 3;
            otherwise
                c.ok = true;
        end
        C = [C; c]; %#ok<AGROW>
    end
end
end

function L = mkLedger(cols, kinds, decs, rows)
L = cell2table(rows, 'VariableNames', cols);
L.g = string(L.g); L.phase = string(L.phase);
L.Properties.UserData = table(string(cols(3:end))', string(kinds)', decs(:), 'VariableNames', {'stat', 'kind', 'dec'});
end

function L = ledger20p()
% RESULTS_LEDGER.md, PRESPEC 20p (code d0cdfca, report ced44cc)
cols = {'g', 'phase', 'n', 'd2', 'd3', 'd4', 'D3', 'E_A', 'holm_A', 'frac_A', 'pbin_A', 'E_B', 'holm_B', 'frac_B', 'pbin_B', ...
    'E_D4', 'pW_D4', 'E_D3QW', 'pW_D3QW', 'pred_HMM', 'in_HMM', 'pred_MC2', 'in_MC2', 'ratio_d3_d2', 'qualifies'};
kinds = {'exact', 'exact', 'exact', 'exact', 'exact', 'mc', 'p', 'mc', 'p', 'mc', 'p', 'mc', 'p', ...
    'mc', 'p', 'mc', 'p', 'mc', 'mc', 'mc', 'mc', 'exact', 'p'};
decs = [0 3 3 3 3 3 2 2 2 3 2 2 2 3 2 3 2 3 2 3 2 2 0];
rows = {
 'CS', 'light', 25, 0.030, 0.108, 0.362, 0.071, 0.005, 0.57, 0.08, 0.36, 0.002, 0.76, 0.08, 0.36, -0.060, 0.11, -0.002, 0.96, 0.084, 0.84, 0.153, 1.00, 3.37, 0
 'CS', 'dark', 25, 0.118, 0.222, 0.442, 0.134, 0.053, 0.01, 0.04, 0.72, 0.061, 0.0073, 0.08, 0.36, 0.031, 0.51, 0.065, 0.00049, 0.096, 0.88, 0.163, 1.00, 2.28, 0
 'per0LD', 'light', 27, 0.011, 0.133, 0.376, 0.089, 0.011, 0.57, 0.11, 0.15, 0.008, 0.76, 0.07, 0.39, 0.003, 0.82, 0.011, 0.0096, 0.080, 0.93, 0.133, 1.00, 6.84, 0
 'per0LD', 'dark', 27, 0.072, 0.177, 0.453, 0.036, -0.032, 0.57, 0.07, 0.39, -0.040, 0.76, 0.07, 0.39, -0.006, 0.92, 0.015, 0.14, 0.089, 0.81, 0.100, 0.96, 2.39, 0
 'fumin', 'light', 49, 0.129, 0.246, 0.493, 0.080, 0.010, 0.57, 0.04, 0.71, 0.010, 0.75, 0.04, 0.71, 0.072, 0.072, 0.001, 0.003, 0.122, 0.90, 0.160, 0.98, 1.97, 0
 'fumin', 'dark', 53, 0.090, 0.254, 0.534, 0.097, 0.040, 0.00024, 0.13, 0.016, 0.045, 0.00064, 0.11, 0.048, 0.080, 0.0014, 0.002, 4.6e-05, 0.072, 0.91, 0.137, 1.00, 2.19, 1
 'sss', 'light', 50, 0.032, 0.118, 0.329, 0.053, 0.000, 0.45, 0.02, 0.92, 0.001, 0.3, 0.04, 0.72, 0.040, 0.028, 0.000, 0.031, 0.068, 0.96, 0.091, 1.00, 1.86, 0
 'sss', 'dark', 53, 0.055, 0.139, 0.393, 0.067, 0.017, 0.093, 0.08, 0.27, 0.010, 0.13, 0.06, 0.5, 0.007, 0.073, 0.000, 0.056, 0.068, 0.92, 0.105, 0.98, 2.49, 0
 'per0DD', 'whole', 25, 0.109, 0.177, 0.329, 0.104, 0.059, 0.013, 0.20, 0.0072, 0.051, 0.01, 0.20, 0.0072, -0.047, 0.075, 0.068, 2.2e-05, 0.058, 0.80, 0.113, 1.00, 2.31, 1
 'CSDD', 'whole', 25, 0.052, 0.177, 0.312, 0.085, 0.044, 0.012, 0.28, 0.00017, 0.043, 0.026, 0.28, 0.00017, 0.007, 0.81, 0.015, 2e-06, 0.049, 0.76, 0.116, 0.96, 2.41, 1};
L = mkLedger(cols, kinds, decs, rows);
end

function L = ledger20pHalves()
% RESULTS_LEDGER.md, PRESPEC 20p amendment A (DD in 12-h clock halves). MC SE not bootstrapped here.
cols = {'g', 'phase', 'n', 'd2', 'd3', 'd4', 'D3', 'E_A', 'holm_A', 'frac_A', 'pbin_A', 'E_B', 'holm_B', 'frac_B', 'pbin_B', 'qualifies'};
kinds = {'exact', 'exact', 'exact', 'exact', 'exact', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p'};
decs = [0 3 3 3 3 3 2 2 2 3 2 2 2 0];
rows = {
 'per0DD', 'ZT0-12', 25, 0.116, 0.253, 0.542, 0.139, 0.062, 0.077, 0.12, 0.13, 0.052, 0.069, 0.12, 0.13, 0
 'per0DD', 'ZT12-24', 25, 0.087, 0.268, 0.570, 0.105, 0.032, 0.1, 0.16, 0.034, 0.026, 0.14, 0.16, 0.034, 0
 'CSDD', 'ZT0-12', 25, 0.083, 0.267, 0.443, 0.121, 0.045, 0.064, 0.24, 0.0012, 0.039, 0.069, 0.24, 0.0012, 0
 'CSDD', 'ZT12-24', 25, 0.064, 0.209, 0.396, 0.095, 0.021, 0.1, 0.08, 0.36, 0.014, 0.14, 0.08, 0.36, 0};
L = mkLedger(cols, kinds, decs, rows);
end

function L = ledger20q()
% RESULTS_LEDGER.md, PRESPEC 20q (report c463b63) with amendment E rows for fumin and sss (a0e5ca3)
cols = {'g', 'phase', 'n', 'J', 'D', 'blocks', 'B_obs', 'sig_mc', 'S1_med', 'S1_p', 'S1_frac', 'S1_bin', 'S2_med', 'S2_p', 'S2_frac', ...
    'HMM_enough', 'HMM_inside', 'Q_mc', 'Q_obs', 'Q_hmm', 'kappa_max', 'S1_120', 'S1_360', 'S2_120', 'S2_360'};
kinds = {'exact', 'exact', 'exact', 'exact', 'exact', 'exact', 'mc', 'p', 'mc', 'p', 'mc', 'p', 'mc', ...
    'mc', 'mc', 'exact', 'mc', 'mc', 'exact', 'p', 'p', 'p', 'p'};
decs = [0 3 3 0 3 3 2 2 2 2 2 2 2 2 2 2 2 2 4 2 2 2 2];
rows = {
 'per0DD', 'whole', 24, 0.514, 3.237, 31, 0.107, 0.111, 0.92, 0.14, 0.00, 1, 0.90, 0.11, 0.42, 0.87, 0.96, 0.99, 0.90, 0.66, 1.0045, 0.91, 0.94, 0.90, 0.92
 'CSDD', 'whole', 25, 0.323, 2.399, 31, 0.055, 0.053, 0.75, 0.0028, 0.00, 1, 0.75, 0.0023, 0.36, 0.88, 0.96, 0.99, 0.75, 0.45, 1.0031, 0.98, 0.80, 0.97, 0.78
 'CS', 'light', 24, 0.177, 2.944, 15, 0.024, 0.037, 0.70, 0.00065, 0.00, 1, 0.70, 0.00043, 0.21, 0.67, 0.96, 1.00, 0.70, 0.73, 1.0023, 0.83, 0.77, 0.83, 0.77
 'CS', 'dark', 25, 0.521, 1.842, 16, 0.117, 0.115, 0.96, 0.89, 0.00, 1, 0.95, 0.77, 0.48, 0.76, 0.96, 0.99, 0.95, 0.43, 1.0080, 1.07, 1.03, 1.05, 0.95
 'per0LD', 'light', 26, 0.135, 1.985, 15, 0.020, 0.018, 0.52, 0.00017, 0.00, 1, 0.51, 0.00013, 0.19, 0.88, 1.00, 1.00, 0.51, 0.22, 1.0051, 0.76, 0.46, 0.74, 0.44
 'per0LD', 'dark', 27, 0.417, 2.860, 16, 0.099, 0.072, 0.74, 0.0055, 0.00, 1, 0.73, 0.0039, 0.22, 0.81, 0.89, 0.99, 0.73, 0.45, 1.0035, 0.90, 0.62, 0.89, 0.60
 'fumin', 'light', 46, 0.334, 1.063, 7, 0.107, 0.149, 0.98, 0.35, 0.02, 0.91, 0.84, 0.13, 0.43, 0.76, 0.91, 0.95, 0.84, 0.44, 1.0096, 0.91, NaN, 0.88, NaN
 'fumin', 'dark', 49, 0.292, 0.640, 8, 0.130, 0.093, 0.92, 0.67, 0.04, 0.71, 0.81, 0.25, 0.41, 0.82, 0.98, 0.92, 0.81, 0.16, 1.0090, 1.00, 1.08, 0.87, 0.99
 'sss', 'light', 36, -0.125, 0.421, 7, 0.088, 0.035, 0.93, 0.4, 0.00, 1, 0.80, 0.23, 0.36, 0.78, 0.97, 0.97, 0.80, 0.29, 1.0043, 0.74, NaN, 0.69, NaN
 'sss', 'dark', 48, -0.042, 0.546, 8, 0.064, 0.019, 0.80, 0.062, 0.00, 1, 0.73, 0.029, 0.42, 0.75, 1.00, 0.99, 0.73, 0.01, 1.0041, 0.93, 1.00, 0.84, 0.97};
L = mkLedger(cols, kinds, decs, rows);
end
