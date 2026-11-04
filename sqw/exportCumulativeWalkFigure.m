function X = exportCumulativeWalkFigure(D, outDir)
%EXPORTCUMULATIVEWALKFIGURE  Act 1 opening figure: each CS fly's cumulative net forward cycles, LD and DD (fig_act1_walk.py).
%
%   X = exportCumulativeWalkFigure(D, outDir)
%
%   D       BUILDPERFLYDATA output (uses CS and CSDD: every fly's symbols, and the loop rates of
%           COMPLETEDLOOPRATES for the inset text)
%   outDir  folder for figAct1_walk_mean.csv (time, mean, SD, n per panel), figAct1_walk_ends.csv
%           (each fly's final value), figAct1.mat and figAct1_cumulative_walk.png / .pdf
%
%   Per fly the walk CYCLEWALK (+1 per S -> W, -1 per W -> S over consecutive scored minutes; equal to
%   the net completed cycles within one). Thin lines: flies; black: mean over flies, truncated to the
%   shortest record; band: mean +- 1 SD across flies (population SD, normalised by n). Bar on top: light
%   (open) and dark (filled) halves; grey for constant darkness. Inset: completed loops per hour,
%   hours per net cycle and median loop duration (COMPLETEDLOOPRATES). Corner: share of flies ending
%   above 0 and the mean end. Descriptive, no test.
%   Check (fig_act1_walk.py): CS LD 92% of 25 flies end forward, mean +9.7 in 96 h; CS DD 80%, +10.0.
%
%   X  struct: walks (1 x 2 cell, flies x minutes), t (h), mean, sd, ends, loopRates (B of
%      COMPLETEDLOOPRATES), summary (table: panel, n, hours, fracForward, meanEnd), fig
%
%   Example (from the Icode\sqw folder):
%       D = buildPerFlyData("D:\flies", fullfile("..", "hmm_audit", "sqw20b_fit2_theta0.5.jsonl"), Cache="perfly_cache.mat");
%       X = exportCumulativeWalkFigure(D, "figures"); disp(X.summary)
%
%   See also CYCLEWALK, COMPLETEDLOOPS, COMPLETEDLOOPRATES, NEWPAPERFIGURE.

arguments
    D (1,1) struct
    outDir (1,1) string
end
if ~isfolder(outDir), mkdir(outDir); end
S = paperFigureStyle(); B = completedLoopRates(D);
P = {"CS", 'CS, light–dark', true; "CSDD", 'CS, constant darkness', false};
fh = newPaperFigure("Act 1 walk", 7.0, 3.3); X.fig = fh;
tl = tiledlayout(fh, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
X.walks = cell(1, 2); X.t = cell(1, 2); X.mean = cell(1, 2); X.sd = cell(1, 2); X.ends = cell(1, 2);
summ = table(); meanTab = table(); endTab = table(); ax0 = []; axs = gobjects(1, 2);
for c = 1:2
    g = P{c, 1}; ld = P{c, 3}; G = D.G(S.order == g); col = S.col(S.order == g, :);
    W = arrayfun(@(f) cycleWalk(f.y), G.flies, 'UniformOutput', false);
    T = min(cellfun(@numel, W)); M = cell2mat(cellfun(@(w) w(1:T)', W, 'UniformOutput', false));
    t = (0:T - 1) / 60; mu = mean(M, 1); sd = std(M, 1, 1); e = M(:, end);
    ax = nexttile(tl); hold(ax, 'on'); axs(c) = ax; if isempty(ax0), ax0 = ax; end
    fill(ax, [t fliplr(t)], [mu - sd fliplr(mu + sd)], col, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(ax, t, M', '-', 'Color', [col 0.35], 'LineWidth', 0.6);
    plot(ax, t, mu, 'k-', 'LineWidth', 1.8); yline(ax, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    set(ax, 'XLim', [0 t(end)], 'XTick', 0:24:t(end), 'Clipping', 'off'); xlabel(ax, 'time (h)');
    title(ax, P{c, 2}, 'FontSize', 9, 'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.06]);
    X.walks{c} = M; X.t{c} = t; X.mean{c} = mu; X.sd{c} = sd; X.ends{c} = e;
    meanTab = [meanTab; table(repmat(g, T, 1), t', mu', sd', repmat(numel(e), T, 1), 'VariableNames', {'panel', 't_h', 'mean', 'sd', 'n'})]; %#ok<AGROW>
    endTab = [endTab; table(repmat(g, numel(e), 1), arrayfun(@(f) string(f.file), G.flies), e, 'VariableNames', {'panel', 'file', 'end'})]; %#ok<AGROW>
    summ = [summ; table(g, numel(e), t(end), mean(e > 0), mean(e), 'VariableNames', {'panel', 'n', 'hours', 'fracForward', 'meanEnd'})]; %#ok<AGROW>
    text(ax, 0.98, 0.02, {sprintf('%.0f%% of %d flies end net forward', 100 * mean(e > 0), numel(e)), ...
        sprintf('mean %+.1f cycles in %.0f h', mean(e), t(end))}, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', 'FontSize', 6.5);
    b = B(B.genotype == g, :);
    if ld
        d = b(b.phase == 1, :); n = b(b.phase == 2, :);
        txt = {'completed loops per hour', sprintf('  day:   %.2f fwd, %.2f back', d.fwd_per_h, d.back_per_h), ...
            sprintf('  night: %.2f fwd, %.2f back', n.fwd_per_h, n.back_per_h), ...
            sprintf('one net cycle per ~%.0f h by day, ~%.0f h at night', d.h_per_net_cycle, n.h_per_net_cycle), ...
            sprintf('one loop lasts ~%.0f–%.0f min (median)', d.dur_fwd_median_min, n.dur_fwd_median_min)};
    else
        txt = {'completed loops per hour', sprintf('  %.2f fwd, %.2f back', b.fwd_per_h, b.back_per_h), ...
            sprintf('one net cycle per ~%.0f h', b.h_per_net_cycle), sprintf('one loop lasts ~%.0f min (median)', b.dur_fwd_median_min)};
    end
    text(ax, 0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 6.5, ...
        'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);
    text(ax, -0.10, 1.06, char('a' + c - 1), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
end
linkaxes(axs, 'y'); yl = [min(cellfun(@(m) min(m(:)), X.walks)) max(cellfun(@(m) max(m(:)), X.walks))];
yl = yl + [-0.03 0.08] .* diff(yl); ylim(axs(1), yl); set(axs(2), 'YTickLabel', []);
for c = 1:2                                         % phase bar above each panel, drawn in data coordinates
    ax = axs(c); t = X.t{c}; ld = P{c, 3}; y0 = yl(2) + 0.01 * diff(yl); hb = 0.035 * diff(yl);
    for d = 0:floor(t(end) / 24)
        if 24 * d >= t(end), break, end
        if ld
            rectangle(ax, 'Position', [24 * d y0 12 hb], 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 0.5, 'Clipping', 'off');
            rectangle(ax, 'Position', [24 * d + 12 y0 12 hb], 'FaceColor', 'k', 'EdgeColor', 'k', 'LineWidth', 0.5, 'Clipping', 'off');
        else
            rectangle(ax, 'Position', [24 * d y0 24 hb], 'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none', 'Clipping', 'off');
        end
    end
end
ylabel(ax0, {'net forward cycles completed', '(W\rightarrowQ\rightarrowS\rightarrowW minus reverse)'});
xlabel(tl, {['Thin lines: one fly each. Black: mean; band: \pm1 SD across flies. Bar on top: light (open) and dark (filled) ' ...
    'phases; grey bar: constant darkness. The slope of the black line is the cycle flux J;'], ...
    'single flies wander forward and back, and the forward drift stands clear of their spread only over days.'}, 'FontSize', 6.5, 'FontName', S.font);
X.loopRates = B; X.summary = summ;
exportgraphics(fh, fullfile(outDir, "figAct1_cumulative_walk.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "figAct1_cumulative_walk.pdf"), 'ContentType', 'vector');
writetable(meanTab, fullfile(outDir, "figAct1_walk_mean.csv")); writetable(endTab, fullfile(outDir, "figAct1_walk_ends.csv"));
writetable(B, fullfile(outDir, "loop_rates.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "figAct1.mat"), '-struct', 'Xs');
end
