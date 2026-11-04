function X = exportLoopPersistenceFigure(R, outDir)
%EXPORTLOOPPERSISTENCEFIGURE  SI Fig S11a,b: loop barcodes and per-fly excess persistence E (PRESPEC 26; fig26_si.py -> fig26_SI_ab).
%
%   X = exportLoopPersistenceFigure(R, outDir)
%
%   R       LOOPPERSISTENCETESTS output
%   outDir  folder for fig26_SI_ab.png / .pdf and fig26_SI_points.csv
%
%   (a) barcodes: every completed loop of one example fly at its completion time (blue forward, orange
%       backward), for CS LD, CS DD and per0 DD; the example is the fly (with >= 10 loop pairs) whose loop
%       count is closest to its genotype's median; grey: dark phase (CS LD).
%   (b) per-fly excess persistence E (clock null) by genotype, median with fly-bootstrap 95% CI (2000), and
%       the Wilcoxon p with Holm (primary and reference families) or "(descr.)".
%   X  struct: examples (table), points (table), fig
%
%   Example (from the Icode\sqw folder):
%       R = loopPersistenceTests("D:\flies", Cache="loopseq.mat");
%       X = exportLoopPersistenceFigure(R, "figures");
%
%   See also LOOPPERSISTENCETESTS, LOOPPERSISTENCEFLY.

arguments
    R (1,1) struct
    outDir (1,1) string
end
if ~isfolder(outDir), mkdir(outDir); end
S = paperFigureStyle(); rs = RandStream('twister', 'Seed', 26); F = R.flies;
G = ["CS" "CSDD" "per0LD" "per0DD" "fumin" "sss"];
SH = ["CS LD" "CS DD" "\itper^0\rm LD" "\itper^0\rm DD" "\itfumin\rm" "\itsss\rm"];
fwdC = [42 120 214] / 255; bwdC = [235 104 52] / 255;
fh = newPaperFigure("Fig S11ab", 7.2, 4.6); X.fig = fh;
tl = tiledlayout(fh, 5, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
ax = nexttile(tl, [2 1]); hold(ax, 'on'); ex = table();
for g = ["CS" "CSDD" "per0DD"]
    k = find(F.genotype == g & F.nPairs >= 10); m = median(F.nLoops(k)); [~, j] = min(abs(F.nLoops(k) - m)); ex = [ex; F(k(j), :)]; %#ok<AGROW>
end
for row = 1:height(ex)
    yb = height(ex) - row; j = find(F.genotype == ex.genotype(row) & F.fly == ex.fly(row)); Lp = R.loops{j};
    if ex.genotype(row) == "CS"
        for d = 0:3, patch(ax, 24 * d + [12 24 24 12], yb + [-0.4 -0.4 0.4 0.4], [0.9 0.9 0.9], 'EdgeColor', 'none'); end
    end
    t = (Lp.loopT - 1) / 60;
    for s = [1 -1]
        tt = t(Lp.loopDir == s); if s == 1, cc = fwdC; else, cc = bwdC; end
        plot(ax, [tt'; tt'], repmat(yb + [-0.35; 0.35], 1, numel(tt)), '-', 'Color', cc, 'LineWidth', 1.1);
    end
    text(ax, -1, yb, sprintf('%s (fly %d)', SH(G == ex.genotype(row)), ex.fly(row)), 'HorizontalAlignment', 'right', 'FontSize', 7);
end
set(ax, 'YTick', [], 'XLim', [-16 96], 'YLim', [-0.6 height(ex) - 0.4], 'XTick', 0:24:96, 'YColor', 'none'); ax.XRuler.FirstCrossoverValue = -0.6; xlabel(ax, 'time (h)');
title(ax, 'each tick: one completed loop (blue forward, orange backward); grey: dark phase', 'FontSize', 7, 'FontWeight', 'normal', ...
    'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02]);
text(ax, -0.12, 1.08, 'a', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
ax = nexttile(tl, [3 1]); hold(ax, 'on'); pts = table(); B = R.genotypes;
for k = 1:numel(G)
    g = G(k); e = F.E_clock(F.genotype == g & F.nPairs >= 10); e = e(isfinite(e)); cc = S.col(S.order == g, :);
    scatter(ax, k + rand(rs, numel(e), 1) * 0.4 - 0.2, e, 9, 'MarkerFaceColor', cc, 'MarkerEdgeColor', 'w', 'LineWidth', 0.3, 'MarkerFaceAlpha', 0.85);
    bs = zeros(2000, 1); for b = 1:2000, bs(b) = median(e(randi(rs, numel(e), numel(e), 1))); end
    ci = quantile7(bs, [0.025 0.975]); m = median(e);
    plot(ax, k + [-0.3 0.3], [m m], 'k-', 'LineWidth', 1.3); plot(ax, [k k], ci, 'k-', 'LineWidth', 0.8);
    q = B(B.genotype == g, :); lab = sprintf('p %.2f', q.p_clock);
    if isfinite(q.holm_clock), lab = {lab, sprintf('Holm %.2f', q.holm_clock)}; else, lab = {lab, '(descr.)'}; end
    text(ax, k, 1.0, lab, 'Units', 'data', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 6);
    pts = [pts; table(repmat(g, numel(e), 1), e, 'VariableNames', {'genotype', 'E_clock'})]; %#ok<AGROW>
end
yl = ylim(ax); for c = ax.Children', if isa(c, 'matlab.graphics.primitive.Text'), c.Position(2) = yl(2); end, end
yline(ax, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
set(ax, 'XTick', 1:numel(G), 'XTickLabel', SH, 'XLim', [0.4 numel(G) + 0.6], 'TickLabelInterpreter', 'tex');
ylabel(ax, {'excess persistence E', '(same-direction fraction, data - chain)'});
xlabel(ax, 'above 0: runs of same-direction loops longer than the chain predicts (episodes); below 0: more alternation', 'FontSize', 6.5);
text(ax, -0.12, 1.12, 'b', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
X.examples = ex; X.points = pts;
exportgraphics(fh, fullfile(outDir, "fig26_SI_ab.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "fig26_SI_ab.pdf"), 'ContentType', 'vector');
writetable(pts, fullfile(outDir, "fig26_SI_points.csv"));
end
