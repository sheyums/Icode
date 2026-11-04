function R = loopPersistenceTests(dataRoot, options)
%LOOPPERSISTENCETESTS  PRESPEC 26: is loop direction more persistent than a clock-matched chain predicts? (sqw26.py)
%
%   R = loopPersistenceTests(dataRoot)
%   R = loopPersistenceTests(dataRoot, NSim=500, MinPairs=10, Cache="loopseq.mat")
%
%   dataRoot  folder with the count-file subfolders per0DD, per0LD, CS, CSDD, fumin, sss (all 57 sss files)
%   Per fly LOOPPERSISTENCEFLY (seed 2600 + 1000 x genotype index + fly index, as sqw26.py). Flies with at
%   least MinPairs consecutive loop pairs enter the tests:
%     Wilcoxon signed-rank (two-sided, WILCOXONP) of E against 0, clock null (primary) and 12-h null;
%     Holm within {per0 DD, per0 LD} and within {CS LD, CS DD}; fumin and sss descriptive;
%     Mann-Whitney (two-sided, normal approximation with tie and continuity correction, as scipy's
%     default for these sizes) of E_clock, per0 DD vs CS DD and per0 LD vs CS LD.
%   R  struct: flies (table, one row per fly: genotype, fly (1-based), file, nLoops, nPairs, fracFwd,
%      S_obs, S_null_clock, E_clock, S_null_phase, E_phase, S_shuffle), loops (cell: each fly's loopT,
%      loopDir, runs), genotypes (table: n flies, excluded, loops, medians, E med / mean, share E > 0,
%      p and Holm per null), secondary (table: comparison, median difference, p)
%   Expected (port_ref/s26_loopseq.json): no Holm < 0.05 (per0 DD clock Holm 0.068); ranksum p 0.30 / 0.65;
%   loops per0 DD 1472, CS LD 1282 (exact).
%
%   Example (from the Icode\sqw folder):
%       R = loopPersistenceTests("D:\flies", Cache="loopseq.mat"); disp(R.genotypes)
%       X = exportLoopPersistenceFigure(R, "figures");
%
%   See also LOOPPERSISTENCEFLY, EXPORTLOOPPERSISTENCEFIGURE, WILCOXONP, HOLMADJUST.

arguments
    dataRoot (1,1) string {mustBeFolder}
    options.NSim (1,1) double = 500
    options.MinPairs (1,1) double = 10
    options.Cache (1,1) string = ""
end
if strlength(options.Cache) > 0 && isfile(options.Cache)
    C = load(options.Cache, 'R'); R = C.R; fprintf('loopPersistenceTests: loaded %s\n', options.Cache); return
end
G = ["per0DD" "per0LD" "CS" "CSDD" "fumin" "sss"]; dd = ["per0DD" "CSDD"];
F = table(); L = {};
for gi = 1:numel(G)
    g = G(gi); files = sqwListFiles(fullfile(dataRoot, g)); rows = cell(numel(files), 1); lp = cell(numel(files), 1); o = options;
    parfor i = 1:numel(files)
        [y, ph] = sqwFlyInput(files(i), ~ismember(g, dd));
        r = loopPersistenceFly(y, ph, NSim=o.NSim, RandomSeed=2600 + 1000 * (gi - 1) + (i - 1));
        [~, n, e] = fileparts(files(i));
        rows{i} = table(g, i, string(n) + e, r.nLoops, r.nPairs, r.fracFwd, r.S_obs, r.S_null_clock, r.E_clock, r.S_null_phase, ...
            r.E_phase, r.S_shuffle, 'VariableNames', {'genotype', 'fly', 'file', 'nLoops', 'nPairs', 'fracFwd', 'S_obs', ...
            'S_null_clock', 'E_clock', 'S_null_phase', 'E_phase', 'S_shuffle'});
        lp{i} = struct('loopT', r.loopT, 'loopDir', r.loopDir, 'runs', r.runs_obs, 'runsNullClock', r.runs_null_clock);
    end
    F = [F; vertcat(rows{:})]; L = [L; lp]; %#ok<AGROW>
    fprintf('loopPersistenceTests: %s done\n', g);
end
B = table();
for g = G
    U = F(F.genotype == g & F.nPairs >= options.MinPairs & isfinite(F.E_clock), :);
    row = table(g, height(U), nnz(F.genotype == g) - height(U), sum(U.nLoops), median(U.S_obs), median(U.S_shuffle), median(U.fracFwd), ...
        'VariableNames', {'genotype', 'n', 'excluded', 'loops', 'S_obs_med', 'S_shuffle_med', 'fracFwd_med'});
    for nm = ["clock" "phase"]
        e = U.("E_" + nm);
        row.("S_null_med_" + nm) = median(U.("S_null_" + nm)); row.("E_med_" + nm) = median(e); row.("E_mean_" + nm) = mean(e);
        row.("fracEpos_" + nm) = mean(e > 0); row.("p_" + nm) = wilcoxonP(e); row.("holm_" + nm) = NaN;
    end
    B = [B; row]; %#ok<AGROW>
end
for fam = {["per0DD" "per0LD"], ["CS" "CSDD"]}
    k = ismember(B.genotype, fam{1});
    for nm = ["clock" "phase"], B.("holm_" + nm)(k) = holmAdjust(B.("p_" + nm)(k)); end
end
S = table();
for pr = {["per0DD" "CSDD"], ["per0LD" "CS"]}
    a = F.E_clock(F.genotype == pr{1}(1) & F.nPairs >= options.MinPairs); b = F.E_clock(F.genotype == pr{1}(2) & F.nPairs >= options.MinPairs);
    S = [S; table(pr{1}(1) + "_vs_" + pr{1}(2), median(a) - median(b), mannWhitneyP(a, b), 'VariableNames', {'comparison', 'diff_med', 'p'})]; %#ok<AGROW>
end
R.flies = F; R.loops = L; R.genotypes = B; R.secondary = S;
if strlength(options.Cache) > 0, save(options.Cache, 'R'); end
end

function p = mannWhitneyP(a, b)
% two-sided Mann-Whitney U, normal approximation with tie correction and continuity correction (scipy default)
a = a(isfinite(a)); b = b(isfinite(b)); n1 = numel(a); n2 = numel(b); x = [a; b]; n = n1 + n2;
[s, i] = sort(x); r = zeros(n, 1); k = 1; t = [];
while k <= n
    j = k; while j < n && s(j + 1) == s(k), j = j + 1; end
    r(i(k:j)) = (k + j) / 2; t(end + 1) = j - k + 1; k = j + 1; %#ok<AGROW>
end
U = sum(r(1:n1)) - n1 * (n1 + 1) / 2; mu = n1 * n2 / 2;
sd = sqrt(n1 * n2 / 12 * ((n + 1) - sum(t .^ 3 - t) / (n * (n - 1))));
z = (abs(U - mu) - 0.5) / sd; p = min(1, erfc(z / sqrt(2)));
end
