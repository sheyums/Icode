function G = sqwMemoryGroup(files, genotype, options)
%SQWMEMORYGROUP  Block irreversibility (PRESPEC 20p) for a group of flies: per fly-phase table, phase summary, CSV, plot.
%
%   G = sqwMemoryGroup(files, genotype)
%   G = sqwMemoryGroup(files, genotype, LightDark=true, HMMFits=[], Halves=~LightDark, ...
%                      NNullA=200, NNullB=200, NModel=50, RandomSeed=15000, OutDir="", Plot=true)
%
%   files     count files (1 Hz, starting at ZT0), one per fly; a folder is
%             searched recursively for *.txt. Flies are taken in sorted order.
%   genotype  label for tables and file names
%   LightDark true: LD record, light/dark phases; false: DD record, one phase
%   HMMFits   struct array aligned with the sorted files, fields pi, P, B (labelled
%             K = 3 fit) and deg; [] skips the HMM predictions
%   Halves    DD only (amendment A): also the two 12-h clock halves, null A2
%   RandomSeed  fly k uses RandomSeed + k
%
%   G  struct:
%     perFly   one row per fly x phase: genotype, file, fly, phase, hasRef (Q/W
%              reference exists), deg, and SQWMEMORYFLY's phase columns
%     summary  one row per phase, over flies WITH a Q/W reference (20p amendment B):
%              n; medians of d2, d3, d4, D3; E_A median, its Wilcoxon signed-rank p
%              (two-sided), fraction of flies with p_A < 0.05 and its one-sided
%              Binomial(n, 0.05) p; the same for null B; E_D4, E_D3QW medians and
%              Wilcoxon p; HMM and MC2 predicted D3 (median) and fraction inside the
%              central 95%; median d3/d2 (d2 > 0). Holm across blocks is the
%              caller's (HOLMADJUST), since the family spans genotypes.
%     halves, halvesSummary  (Halves) the same for the 12-h halves (null A2, B)
%     flies    the SQWMEMORYFLY structs
%   Writes <genotype>_memory_perfly.csv and _summary.csv (and _halves) to OutDir.
%
%   Only the excess E over the fly's own null is interpretable; raw d3/d2 is
%   estimator bias (FINDINGS / ledger 20p).
%
%   Example (from the Icode\sqw folder):
%       G = sqwMemoryGroup("D:\flies\CS", "CS", OutDir="perfly");   % writes perfly\CS_memory_perfly.csv
%       disp(G.summary)
%
%   See also SQWMEMORYFLY, SQWTURGROUP, HOLMADJUST.

arguments
    files string
    genotype (1,1) string
    options.LightDark (1,1) logical = true
    options.HMMFits = []
    options.Halves = []
    options.NNullA (1,1) double {mustBeNonnegative, mustBeInteger} = 200
    options.NNullB (1,1) double {mustBeNonnegative, mustBeInteger} = 200
    options.NModel (1,1) double {mustBeNonnegative, mustBeInteger} = 50
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 15000
    options.OutDir (1,1) string = ""
    options.Plot (1,1) logical = true
end
files = sqwListFiles(files);
halves = ~options.LightDark; if ~isempty(options.Halves), halves = logical(options.Halves) && ~options.LightDark; end
n = numel(files); out = cell(n, 1); fits = options.HMMFits; o = options;
if ~isempty(fits) && numel(fits) ~= n, error('sqw:fits', 'HMMFits must have one entry per file (%d).', n); end
parfor k = 1:n
    [y, ph, nP, refs, names] = sqwFlyInput(files(k), o.LightDark);
    if isempty(fits), h = []; deg = NaN; else, h = fits(k); deg = h.deg; end
    F = sqwMemoryFly(y, ph, nP, HMM=h, NNullA=o.NNullA, NNullB=o.NNullB, NModel=o.NModel, ...
        RandomSeed=o.RandomSeed + k, Halves=halves);
    t = F.phase;
    t = [table(repmat(genotype, nP, 1), repmat(fileName(files(k)), nP, 1), repmat(k, nP, 1), names, ...
        arrayfun(@(nm) isfinite(refs.(nm)), names), repmat(deg, nP, 1), ...
        'VariableNames', {'genotype', 'file', 'fly', 'phase', 'hasRef', 'deg'}), t];
    F.rows = t;
    if halves
        F.halfRows = [table(repmat(genotype, 2, 1), repmat(fileName(files(k)), 2, 1), repmat(k, 2, 1), ["ZT0-12"; "ZT12-24"], ...
            repmat(isfinite(refs.whole), 2, 1), 'VariableNames', {'genotype', 'file', 'fly', 'phase', 'hasRef'}), F.halves];
    end
    out{k} = F;
end
G.flies = vertcat(out{:});
G.perFly = vertcat(G.flies.rows);
G.summary = summarise(G.perFly, true);
if halves
    G.halves = vertcat(G.flies.halfRows);
    G.halvesSummary = summarise(renamevars(G.halves, {'E_A2', 'p_A2'}, {'E_A', 'p_A'}), false);
end
if strlength(options.OutDir) > 0
    if ~isfolder(options.OutDir), mkdir(options.OutDir); end
    writetable(G.perFly, fullfile(options.OutDir, genotype + "_memory_perfly.csv"));
    writetable(G.summary, fullfile(options.OutDir, genotype + "_memory_summary.csv"));
    if halves
        writetable(G.halves, fullfile(options.OutDir, genotype + "_memory_halves_perfly.csv"));
        writetable(G.halvesSummary, fullfile(options.OutDir, genotype + "_memory_halves_summary.csv"));
    end
end
if options.Plot, plotMemory(G, genotype); end
end

function S = summarise(T, full)
S = table();
for ph = unique(T.phase, 'stable')'
    D = T(T.phase == ph & T.hasRef, :); r = table(); r.phase = ph; r.n = height(D);
    for k = ["d2" "d3" "d4" "D3"], r.(k) = median(D.(k)); end
    for nul = ["A" "B"]
        E = D.("E_" + nul); pv = D.("p_" + nul);
        r.("E_" + nul) = median(E); r.("pW_" + nul) = wilcoxonP(E);
        r.("frac_" + nul) = mean(pv < 0.05); r.("pbin_" + nul) = binomGreater(sum(pv < 0.05), numel(pv), 0.05);
    end
    if full
        for k = ["D4" "D3QW"], r.("E_" + k) = median(D.("E_" + k)); r.("pW_" + k) = wilcoxonP(D.("E_" + k)); end
        for m = ["HMM" "MC2"]
            ok = ~isnan(D.("in_" + m));
            r.("pred_" + m) = median(D.("pred_" + m)(ok)); r.("in_" + m) = mean(D.("in_" + m)(ok));
        end
        r.ratio_d3_d2 = median(D.d3(D.d2 > 0) ./ D.d2(D.d2 > 0));
    end
    S = [S; r]; %#ok<AGROW>
end
end

function plotMemory(G, genotype)
T = G.perFly(G.perFly.hasRef, :); phs = unique(T.phase, 'stable');
fh = figure('Name', genotype + " block irreversibility (20p)", 'Color', 'w');
try theme(fh, 'light'); catch, end
tl = tiledlayout(fh, 1, numel(phs), 'TileSpacing', 'compact');
for k = 1:numel(phs)
    D = T(T.phase == phs(k), :); ax = nexttile(tl); hold(ax, 'on');
    nullMed = D.D3 - D.E_A;
    scatter(ax, nullMed, D.D3, 30, D.p_A < 0.05, 'filled');
    lim = [min([nullMed; D.D3]) max([nullMed; D.D3])]; plot(ax, lim, lim, 'k-');
    colormap(ax, [0.4 0.4 0.4; 0.85 0.2 0.1]);
    xlabel(ax, 'median \Delta_3 under null A (nats/h)'); ylabel(ax, 'observed \Delta_3 (nats/h)');
    title(ax, sprintf('%s %s: median E %+.3f, p_A<.05 in %.2f', genotype, phs(k), median(D.E_A), mean(D.p_A < 0.05)));
    box(ax, 'on');
end
title(tl, 'Each point a fly; red: p_A < 0.05; on the line: no excess over the fly''s own first-order chain');
end

function s = fileName(f)
[~, a, b] = fileparts(f); s = a + b;
end
