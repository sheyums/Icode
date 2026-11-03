function X = exportSleepMemoryFigure(dataRoot, jsonlFile, outDir, options)
%EXPORTSLEEPMEMORYFIGURE  Fig C7, sleep-length memory at the observed level (FIGURE_RECIPES L note 2026-11-01; fig_R1_obs.py; PRESPEC 20k).
%
%   X = exportSleepMemoryFigure(dataRoot, jsonlFile, outDir)
%   X = exportSleepMemoryFigure(dataRoot, jsonlFile, outDir, NPerm=2000, NBoot=2000, RandomSeed=20, Cache="")
%
%   dataRoot   count files per genotype: CS, fumin, per0LD, per0DD, sss (all 57; the 20b fits
%              are indexed in each sorted list)
%   jsonlFile  hmm_audit/sqw20b_fit2_theta0.5.jsonl: 20k used the flies whose 20b fit is NOT
%              degenerate (sqw10.py), so the same flies are used here
%   Cache      .mat for the per-fly route records (loaded if it exists, written otherwise)
%
%   Recomputes the observed half of 20k for R1 (exit route vs sleep-bout length) from the raw
%   records: per fly x phase D = ROUTEDURATIOND; per block the mean over flies with a finite D,
%   the within-fly label-permutation null (NPerm; sqw10_report.py), p = 2 min(P(null >= D),
%   P(null <= D)) capped at 1, its 2.5-97.5% range, Holm over the 9 blocks in 20k's order
%   (CS, fumin, per0 LD, per0 DD, sss). X.blocks carries the 20k ledger values beside them.
%
%   (a) per-fly D per block (open light, filled dark / DD; faded where not significant), clipped
%       at +-1.5, bar = MEAN over flies, grey = the null's 95% range (this run), * = Holm < 0.05.
%       Significance (star, fading) is the 20k LEDGER's Holm by default (SigFrom="ledger"), as
%       Fig C12 copies its flags: per0 LD light sits at Holm 8p with p ~ 0.004, so at 2000
%       permutations the star turns on Monte Carlo luck (this port's seed 20 drew p 0.009, Holm
%       0.072; 20000 permutations, two seeds, gave p 0.0041 / 0.0043, Holm 0.033 / 0.034, as
%       the ledger). SigFrom="port" uses this run's Holm. X.blocks carries both.
%   (b) CS dark: P(exit through quiet wake) per sleep-length bin [5 10 20 40 80 Inf) min, mean of
%       per-fly proportions (bins with >= 10 flies; ROUTE_BINS.fly_equal_bins), fly-bootstrap 95%
%       interval; dashed: the same flies' overall P(quiet exit), which a Markov chain predicts at
%       every length (recipe: 0.202).
%   Writes figC7_perfly_D.csv, figC7_blocks.csv, figC7_bins.csv, figC7.mat and
%   figC7_sleep_memory_observed.png / .pdf.
%
%   Example (from the Icode\sqw folder):
%       dataRoot = "D:\flies";          % one subfolder of 1-s count files per genotype
%       jsonl = fullfile("..", "hmm_audit", "sqw20b_fit2_theta0.5.jsonl");
%       X = exportSleepMemoryFigure(dataRoot, jsonl, "figures");   % writes figures\figC7.mat
%
%   See also ROUTEDURATIONOBSERVED, ROUTEDURATIOND, HOLMADJUST.

arguments
    dataRoot (1,1) string {mustBeFolder}
    jsonlFile (1,1) string {mustBeFile}
    outDir (1,1) string
    options.NPerm (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
    options.Cache (1,1) string = ""
    options.SigFrom (1,1) string {mustBeMember(options.SigFrom, ["ledger" "port"])} = "ledger"
end
if ~isfolder(outDir), mkdir(outDir); end
G = ["CS" "fumin" "per0LD" "per0DD" "sss"];
if strlength(options.Cache) > 0 && isfile(options.Cache)
    C = load(options.Cache, 'flies'); flies = C.flies;
else
    fits = readStoredHmmFits(jsonlFile); flies = table();
    for g = G
        files = sqwListFiles(fullfile(dataRoot, g)); Fg = fits([fits.g] == g);
        assert(numel(Fg) == numel(files), 'exportSleepMemoryFigure:fits', '%s: %d fits, %d files.', g, numel(Fg), numel(files));
        deg = arrayfun(@(i) Fg([Fg.fly0] == i - 1).deg, (1:numel(files))');
        use = find(~deg); rec = cell(numel(use), 1); ld = g ~= "per0DD";
        parfor j = 1:numel(use)
            rec{j} = routeDurationObserved(readCounts(files(use(j))), LightDark=ld);
        end
        [~, nm, ex] = fileparts(files(use));
        flies = [flies; table(repmat(g, numel(use), 1), use - 1, nm + ex, rec, 'VariableNames', {'genotype', 'fly0', 'file', 'records'})]; %#ok<AGROW>
        fprintf('exportSleepMemoryFigure: %-6s %d of %d flies (non-degenerate 20b fit)\n', g, numel(use), numel(files));
    end
    if strlength(options.Cache) > 0, save(options.Cache, 'flies'); end
end

% ---- 20k observed R1 statistics, blocks in 20k's order (Holm family)
rs = RandStream('twister', 'Seed', options.RandomSeed);
blk = table(["CS"; "CS"; "fumin"; "fumin"; "per0LD"; "per0LD"; "per0DD"; "sss"; "sss"], [1; 2; 1; 2; 1; 2; 1; 1; 2], ...
    ["light"; "dark"; "light"; "dark"; "light"; "dark"; "whole"; "light"; "dark"], 'VariableNames', {'genotype', 'phase', 'phaseName'});
% ledger (20k, observed, R1): n, D, null mean, p, Holm
led = [25 0.087 0.001 0.125 0.75; 24 -0.293 0.002 0 0; 25 -0.002 -0.000 0.994 1; 15 -0.061 0.001 0.349 1; ...
       26 -0.228 0.000 0.005 0.04; 26 -0.121 0.002 0.091 0.637; 23 -0.028 -0.001 0.589 1; 32 0.010 0.000 0.848 1; 31 0.008 0.000 0.873 1];
pts = table();
for b = 1:height(blk)
    F = flies(flies.genotype == blk.genotype(b), :);
    D = cellfun(@(r) routeDurationD(r, "R1_exit_vs_sleep", blk.phase(b)), F.records);
    ok = isfinite(D); null = zeros(options.NPerm, 1);
    T = cell(nnz(ok), 1); Q = T; k = find(ok);
    for j = 1:numel(k)
        r = F.records{k(j)}; q = [r.q_out]'; m = [r.phase]' == blk.phase(b) & ~isnan(q);
        T{j} = log([r(m).sleep_min]'); Q{j} = q(m);
    end
    for i = 1:options.NPerm
        s = 0;
        for j = 1:numel(k)
            lab = Q{j}(randperm(rs, numel(Q{j})));
            s = s + mean(T{j}(lab == 1)) - mean(T{j}(lab == 0));
        end
        null(i) = s / numel(k);
    end
    real = mean(D(ok)); qn = quantile7(null, [0.025 0.975]);
    blk.n(b) = nnz(ok); blk.D(b) = real; blk.nullMean(b) = mean(null); blk.null_lo(b) = qn(1); blk.null_hi(b) = qn(2);
    blk.p(b) = min(1, 2 * min(mean(null >= real), mean(null <= real)));
    pts = [pts; table(repmat(blk.genotype(b), nnz(ok), 1), repmat(blk.phaseName(b), nnz(ok), 1), F.file(ok), F.fly0(ok), D(ok), ...
        'VariableNames', {'genotype', 'phase', 'file', 'fly0', 'D'})]; %#ok<AGROW>
end
blk.holm = holmAdjust(blk.p);
blk.ledger_n = led(:, 1); blk.ledger_D = led(:, 2); blk.ledger_nullMean = led(:, 3); blk.ledger_p = led(:, 4); blk.ledger_holm = led(:, 5);
if options.SigFrom == "ledger", blk.sig = blk.ledger_holm < 0.05; else, blk.sig = blk.holm < 0.05; end
blk.sigFrom = repmat(options.SigFrom, height(blk), 1);

% ---- panel (b): CS dark, fly-equal bins
F = flies(flies.genotype == "CS", :); edges = [5 10 20 40 80 Inf]; bins = table(); allp = NaN(height(F), 1);
for j = 1:height(F)
    r = F.records{j}; q = [r.q_out]'; m = [r.phase]' == 2 & ~isnan(q);
    if any(m), allp(j) = mean(q(m)); end
end
flat = mean(allp, 'omitnan');
for a = 1:numel(edges) - 1
    props = [];
    for j = 1:height(F)
        r = F.records{j}; q = [r.q_out]'; L = [r.sleep_min]';
        m = [r.phase]' == 2 & ~isnan(q) & L >= edges(a) & L < edges(a + 1);
        if any(m), props(end + 1, 1) = mean(q(m)); end %#ok<AGROW>
    end
    if numel(props) < 10, continue, end
    bs = zeros(options.NBoot, 1);
    for i = 1:options.NBoot, bs(i) = mean(props(randi(rs, numel(props), numel(props), 1))); end
    ci = quantile7(bs, [0.025 0.975]);
    bins = [bins; table(edges(a), edges(a + 1), sqrt(edges(a) * min(edges(a + 1), 2 * edges(a))), mean(props), ci(1), ci(2), ...
        numel(props), 'VariableNames', {'from_min', 'to_min', 'x_min', 'Pquiet', 'lo', 'hi', 'nFlies'})]; %#ok<AGROW>
end
X.flies = flies; X.points = pts; X.blocks = blk; X.bins = bins; X.flat = flat;

% ---- figure
S = paperFigureStyle(); fh = newPaperFigure("Fig C7", 6.2, 2.9); X.fig = fh;
tl = tiledlayout(fh, 1, 9, 'TileSpacing', 'loose', 'Padding', 'compact');
ax = nexttile(tl, [1 5]); hold(ax, 'on');
order = [1 2 5 6 7 3 4 8 9];       % figure order: CS, per0 LD, per0 DD, fumin, sss (light before dark)
nameS = containers.Map(["CS" "per0LD" "per0DD" "fumin" "sss"], ["CS" "\itper^0\rm LD" "\itper^0\rm DD" "\itfmn\rm" "\itsss\rm"]);
labs = strings(1, 9);
for x = 1:9
    b = order(x); g = blk.genotype(b); c = S.col(S.order == g, :); xx = x - 1;
    v = pts.D(pts.genotype == g & pts.phase == blk.phaseName(b)); v = min(max(v, -1.5), 1.5);
    patch(ax, xx + [-0.33 0.33 0.33 -0.33], [blk.null_lo(b) blk.null_lo(b) blk.null_hi(b) blk.null_hi(b)], [0.86 0.86 0.86], 'EdgeColor', 'none');
    if blk.phaseName(b) == "light", mf = 'w'; else, mf = c; end
    sc = scatter(ax, xx + rand(rs, numel(v), 1) * 0.36 - 0.18, v, 7, 'MarkerFaceColor', mf, 'MarkerEdgeColor', c, 'LineWidth', 0.6);
    if ~blk.sig(b), sc.MarkerFaceAlpha = 0.55; sc.MarkerEdgeAlpha = 0.55; end
    plot(ax, xx + [-0.28 0.28], mean(pts.D(pts.genotype == g & pts.phase == blk.phaseName(b))) * [1 1], '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.6);
    if blk.sig(b), text(ax, xx, 1.38, '*', 'HorizontalAlignment', 'center', 'FontSize', 11); end
    if blk.phaseName(b) == "whole", labs(x) = nameS(g); else, labs(x) = nameS(g) + " " + extractBefore(upper(blk.phaseName(b)), 2); end
end
yline(ax, 0, 'Color', [0.53 0.53 0.53], 'LineWidth', 0.5);
set(ax, 'XTick', 0:8, 'XTickLabel', labs, 'XTickLabelRotation', 60, 'YLim', [-1.55 1.55], 'XLim', [-0.6 8.6], 'TickLabelInterpreter', 'tex');
ax.XAxis.FontSize = 6.5;
ylabel(ax, 'D = \langleln t\rangle_{quiet exit} - \langleln t\rangle_{active exit}');
title(ax, '(a) does the waking route depend on time asleep?', 'FontSize', 8, 'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02]);
ax = nexttile(tl, [1 4]); hold(ax, 'on'); c = S.col(1, :);
h1 = errorbar(ax, bins.x_min, bins.Pquiet, bins.Pquiet - bins.lo, bins.hi - bins.Pquiet, 'o-', 'Color', c, 'MarkerFaceColor', c, ...
    'MarkerSize', 4, 'LineWidth', 1.0, 'CapSize', 3);
h2 = yline(ax, flat, '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 0.9);
set(ax, 'XScale', 'log', 'XLim', [4.5 160], 'XTick', [5 10 20 40 80], 'XTickLabel', {'5', '10', '20', '40', '80'}, 'YLim', [0 0.5], 'XMinorTick', 'off');
xlabel(ax, 'sleep-bout length (min)'); ylabel(ax, 'P(wake through quiet wake)');
title(ax, '(b) CS, night', 'FontSize', 8, 'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02]);
legend(ax, [h1 h2], {'observed', 'no memory (Markov): same at every length'}, 'Box', 'off', 'FontSize', 6, 'Location', 'northeast');
xlabel(tl, {['(a) Points: flies (open light, filled dark/DD; faded where not significant); bar: mean; grey: 95% range of the mean ' ...
    'if the route ignored time asleep'], ['(within-fly shuffle); *: Holm p < 0.05; |D| clipped at 1.5. (b) Mean over flies of ' ...
    'each fly''s proportion per bin (>= 10 flies), bootstrap 95% intervals.']}, 'FontSize', 5.6, 'Color', [0.2 0.2 0.2], 'FontName', S.font);
exportgraphics(fh, fullfile(outDir, "figC7_sleep_memory_observed.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "figC7_sleep_memory_observed.pdf"), 'ContentType', 'vector');
writetable(pts, fullfile(outDir, "figC7_perfly_D.csv"));
writetable(blk, fullfile(outDir, "figC7_blocks.csv"));
writetable(bins, fullfile(outDir, "figC7_bins.csv"));
Xs = rmfield(X, {'fig', 'flies'}); save(fullfile(outDir, "figC7.mat"), '-struct', 'Xs');
end
