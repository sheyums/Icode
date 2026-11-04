function R = windowDirectionTests(dataRoot, options)
%WINDOWDIRECTIONTESTS  PRESPEC 25: direction of measurable cycles in 6-h / 4-h single-day windows within flies (sqw25.py).
%
%   R = windowDirectionTests(dataRoot)
%   R = windowDirectionTests(dataRoot, NSim=500, Sizes=[360 240], RandomSeed=25, Cache="windows.mat")
%
%   dataRoot  folder with the count-file subfolders per0LD, per0DD, CS, CSDD, sss (all 57)
%   Windows of Sizes minutes from ZT0 (record start in DD). Per window N_SW (S -> W) and N_WS (W -> S) over
%   consecutive minute pairs with BOTH minutes inside the window. Measurable if
%   min(1, 2 binocdf(min(a, b), a + b, 0.5)) < 0.05 (a = N_SW, b = N_WS): forward (+1) if a > b, backward
%   (-1) if b > a, else 0; a window with any missing minute is excluded (class 9; Amendment A, sss only).
%   Nulls (NSim records per fly, starting in the first scored state): clock-matched chain (one matrix per
%   window position of the day; an empty row takes the 12-h matrix of the window's first minute;
%   primary) and the 12-h chain (secondary), classified the same way. Per block: B and F (backward /
%   forward windows summed over flies), null means and 95% ranges, p_B = (1 + #{B_null >= B}) / (NSim + 1),
%   Holm over the 4 blocks of {per0 LD, per0 DD} x sizes and of {CS LD, CS DD} x sizes; Wilcoxon of
%   B - mean B_null per fly. Counts and classes are exact; null quantities are Monte Carlo.
%   R  struct: blocks (table: key, genotype, size_min, n, windowsPerFly, F, B, flies_anyB, and per null:
%      B_null_mean, B_null_lo, B_null_hi, p_B, holm_B, F_null_mean, p_F, wilcoxon_p), flies (cell per block:
%      table of fly, nSW, nWS (per window), cls)
%   Expected (port_ref/s25_windows.json): 4 measurable forward and 0 backward windows in total.
%
%   Example (from the Icode\sqw folder):
%       R = windowDirectionTests("D:\flies", Cache="windows.mat"); disp(R.blocks(:, ["key" "F" "B"]))
%
%   See also LOOPPERSISTENCETESTS, WILCOXONP, HOLMADJUST.

arguments
    dataRoot (1,1) string {mustBeFolder}
    options.NSim (1,1) double = 500
    options.Sizes (1,:) double = [360 240]
    options.RandomSeed (1,1) double = 25
    options.Cache (1,1) string = ""
end
if strlength(options.Cache) > 0 && isfile(options.Cache)
    C = load(options.Cache, 'R'); R = C.R; fprintf('windowDirectionTests: loaded %s\n', options.Cache); return
end
rs = RandStream('twister', 'Seed', options.RandomSeed); NR = options.NSim;
G = ["per0LD" "per0DD" "CS" "CSDD" "sss"]; dd = ["per0DD" "CSDD"]; B = table(); FL = {};
for g = G
    files = sqwListFiles(fullfile(dataRoot, g)); Y = cell(numel(files), 1); PH = Y;
    for i = 1:numel(files), [Y{i}, PH{i}] = sqwFlyInput(files(i), ~ismember(g, dd)); end
    for L = options.Sizes
        nF = numel(files); fl = table(); Bn = struct('clock', zeros(NR, 1), 'phase', zeros(NR, 1)); Fn = Bn; per = struct('clock', zeros(nF, 1), 'phase', zeros(nF, 1));
        for i = 1:nF
            y = Y{i}; y(y == 4) = 0; ph = PH{i}; T = numel(y); miss = y == 0;
            Pph = chainMatrices(y, ph, max(ph), []);
            wpos = floor(mod((0:T - 1)', 1440) / L) + 1; nwp = 1440 / L;
            fb = zeros(3, 3, nwp); for m = 1:nwp, fb(:, :, m) = Pph(:, :, ph(find(wpos == m, 1))); end
            Pclk = chainMatrices(y, wpos, nwp, fb);
            [a, b] = windowCounts(y', L); wbad = any(reshape(miss, L, []), 1); cls = classify(a, b); cls(wbad) = 9;
            fl = [fl; table(i, {a}, {b}, {cls}, nnz(cls == -1), nnz(cls == 1), 'VariableNames', {'fly', 'nSW', 'nWS', 'cls', 'B', 'F'})]; %#ok<AGROW>
            y0 = y(find(~miss, 1)); if isempty(y0), y0 = 1; end
            nulls = {"clock", Pclk, wpos; "phase", Pph, ph};
            for k = 1:2
                X = simulateChain(nulls{k, 2}, nulls{k, 3}, y0, T, NR, rs); [sa, sb] = windowCounts(X, L); sc = classify(sa, sb); sc(:, wbad) = 9;
                Bs = sum(sc == -1, 2); Fs = sum(sc == 1, 2); nm = nulls{k, 1};
                Bn.(nm) = Bn.(nm) + Bs; Fn.(nm) = Fn.(nm) + Fs; per.(nm)(i) = mean(Bs);
            end
        end
        Bo = sum(fl.B); Fo = sum(fl.F);
        row = table(g + "|" + L / 60 + "h", g, L, nF, numel(fl.cls{1}), Fo, Bo, nnz(fl.B > 0), 'VariableNames', ...
            {'key', 'genotype', 'size_min', 'n', 'windowsPerFly', 'F', 'B', 'flies_anyB'});
        for nm = ["clock" "phase"]
            q = quantile7(Bn.(nm), [0.025 0.975]); d = fl.B - per.(nm); pw = 1; if any(d ~= 0), pw = wilcoxonP(d); end
            row.("B_null_mean_" + nm) = mean(Bn.(nm)); row.("B_null_lo_" + nm) = q(1); row.("B_null_hi_" + nm) = q(2);
            row.("p_B_" + nm) = (1 + nnz(Bn.(nm) >= Bo)) / (NR + 1); row.("holm_B_" + nm) = NaN;
            row.("F_null_mean_" + nm) = mean(Fn.(nm)); row.("p_F_" + nm) = (1 + nnz(Fn.(nm) >= Fo)) / (NR + 1); row.("wilcoxon_p_" + nm) = pw;
        end
        B = [B; row]; FL{end + 1, 1} = fl; %#ok<AGROW>
        fprintf('windowDirectionTests: %s done\n', row.key);
    end
end
for fam = {["per0LD" "per0DD"], ["CS" "CSDD"]}
    k = ismember(B.genotype, fam{1});
    for nm = ["clock" "phase"], B.("holm_B_" + nm)(k) = holmAdjust(B.("p_B_" + nm)(k)); end
end
R.blocks = B; R.flies = FL;
if strlength(options.Cache) > 0, save(options.Cache, 'R'); end
end

function [a, b] = windowCounts(Y, L)
% Y (R x T): per window, S -> W and W -> S counts over consecutive pairs inside the window
[R, T] = size(Y); nw = T / L; t = 1:T - 1;                 % pair (t, t + 1), 1-based minutes
inside = mod(t, L) ~= 0; w = floor((t - 1) / L) + 1;
sw = Y(:, 1:end - 1) == 3 & Y(:, 2:end) == 2; ws = Y(:, 1:end - 1) == 2 & Y(:, 2:end) == 3;
a = zeros(R, nw); b = a;
for k = 1:nw, m = inside & w == k; a(:, k) = sum(sw(:, m), 2); b(:, k) = sum(ws(:, m), 2); end
end

function c = classify(a, b)
n = a + b; k = min(a, b); cdf = ones(size(n)); s = k < n;
cdf(s) = betainc(0.5, n(s) - k(s), k(s) + 1);              % P(X <= k), X ~ Binomial(n, 1/2)
p = min(1, 2 * cdf); meas = p < 0.05 & n > 0;
c = zeros(size(n)); c(meas & a > b) = 1; c(meas & b > a) = -1;
end

function P = chainMatrices(y, idx, nm, fb)
% as LOOPPERSISTENCEFLY: transition matrix per index of the arrival minute, consecutive scored minutes
C = zeros(3, 3, nm); a = y(1:end - 1); b = y(2:end); m = idx(2:end); ok = a > 0 & b > 0;
C = C + reshape(accumarray(sub2ind([3 3 nm], a(ok), b(ok), m(ok)), 1, [9 * nm 1]), 3, 3, nm);
P = zeros(3, 3, nm);
for k = 1:nm
    for s = 1:3
        rr = sum(C(s, :, k));
        if rr > 0, P(s, :, k) = C(s, :, k) / rr; elseif isempty(fb), P(s, :, k) = 1 / 3; else, P(s, :, k) = fb(s, :, k); end
    end
end
end

function X = simulateChain(P, idx, y0, T, R, rs)
Cc = cumsum(P, 2); Cc(:, 3, :) = 1; X = zeros(R, T); cur = repmat(y0, R, 1); X(:, 1) = cur;
for t = 2:T
    C = Cc(:, :, idx(t)); u = rand(rs, R, 1);
    cur = 1 + (u > C(cur, 1)) + (u > C(cur, 2)); X(:, t) = cur;
end
end
