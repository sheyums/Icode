out = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\homeo_run';
load(fullfile(out, 'homeo_results.mat'), 'res');
fh = figure('Color', 'w', 'Position', [40 40 1500 850]); try theme(fh, 'light'); catch, end
tl = tiledlayout(fh, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% A1: per0 r_pre vs r_post per fly
T = res.A.per0DD.PerFly;
ax = nexttile(tl); hold(ax, 'on');
scatter(ax, T.r_post, T.r_pre, 28, [0.3 0.45 0.7], 'filled');
l = [-0.3 0.5]; plot(ax, l, l, 'k-'); xline(ax, 0, ':'); yline(ax, 0, ':');
xlabel(ax, 'Sleep vs the FOLLOWING wake bout (r_{post})'); ylabel(ax, 'Sleep vs the PRECEDING wake bout (r_{pre})');
title(ax, {'A  per0 DD, 48 flies: flies sit on the diagonal', ['mean r_{pre} ' sprintf('%.3f', mean(T.r_pre)) ', r_{post} ' sprintf('%.3f', mean(T.r_post)) ', Delta ' sprintf('%.3f (p %.2f)', mean(T.Delta), res.A.per0DD.Group.MeanP(3))]}, 'FontWeight', 'normal');
axis(ax, 'square'); box(ax, 'on');
% A2: group means per dataset for the statistics
ax = nexttile(tl, [1 2]); hold(ax, 'on');
st = ["r_pre","r_post","Delta","r_pre2","q_pre","q_post","DeltaW"];
lab = {'sleep vs wake before','sleep vs wake after','\Delta (before - after)','sleep vs wake 2 back', ...
       'wake vs sleep before','wake vs sleep after','\Delta wake'};
ds = ["per0DD","CSLD","CSDD"]; dl = ["per0 DD (48)","CS LD (40)","CS DD (25)"]; cols = [0.3 0.45 0.7; 0.85 0.45 0.2; 0.4 0.65 0.35];
for i = 1:3
    G = res.A.(ds(i)).Group; T = res.A.(ds(i)).PerFly;
    x = (1:7) + (i - 2) * 0.22;
    for k = 1:7
        v = T.(st(k)); v = v(isfinite(v));
        se = std(v) / sqrt(numel(v));
        plot(ax, x(k) * [1 1], G.Mean(k) + 1.96 * se * [-1 1], '-', 'Color', cols(i, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    plot(ax, x, G.Mean, 'o', 'MarkerFaceColor', cols(i, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'DisplayName', dl(i));
end
yline(ax, 0, ':', 'HandleVisibility', 'off');
set(ax, 'XTick', 1:7, 'XTickLabel', lab, 'XLim', [0.5 7.5]); xtickangle(ax, 20);
ylabel(ax, 'Mean Spearman correlation over flies (\pm 1.96 SE)'); legend(ax, 'Location', 'northeast'); box(ax, 'on');
title(ax, 'B  Adjacent sleep and wake bouts co-vary in every dataset, symmetrically in time, and only with the immediate neighbour', 'FontWeight', 'normal');
% C: deviation curves
ax = nexttile(tl); hold(ax, 'on');
for bm = [30 60]
    R = res.C.("bin" + bm);
    plot(ax, R.LagHours, R.GroupDeviation, '-o', 'MarkerSize', 3, 'DisplayName', sprintf('%d-min windows', bm));
end
yline(ax, 0, ':', 'HandleVisibility', 'off'); legend(ax, 'Location', 'southeast');
xlabel(ax, 'Lag (h)'); ylabel(ax, {'Sleep-amount autocorrelation', 'minus renewal null (group mean)'});
title(ax, {'C  per0 DD windows (DESCRIPTIVE: the test', 'over-rejects and is not homeostat-specific)'}, 'FontWeight', 'normal'); box(ax, 'on');
% D: rho_S grid
ax = nexttile(tl);
taus = [0.5 1 2 4 8 16];
imagesc(ax, res.D.GroupRhoS); colorbar(ax); colormap(ax, parula);
set(ax, 'XTick', 1:6, 'XTickLabel', taus, 'YTick', 1:6, 'YTickLabel', taus);
xlabel(ax, 'Drain time constant during sleep, T_s (h)'); ylabel(ax, 'Build-up time constant during wake, T_w (h)');
title(ax, {sprintf('D  per0 DD: pressure proxy vs next sleep, max %.3f at (0.5, 0.5) h', res.D.Group.MaxRhoS), ...
    'shortest memory wins -- as a shared slow state did in simulation'}, 'FontWeight', 'normal');
% D2: per fly rho_S vs r_pre
ax = nexttile(tl); hold(ax, 'on');
scatter(ax, res.A.per0DD.PerFly.r_pre, res.D.PerFly.RhoS_atGroupBest, 28, [0.3 0.45 0.7], 'filled');
xlabel(ax, 'r_{pre} (sleep vs the preceding wake bout)'); ylabel(ax, '\rho_S at (0.5, 0.5) h');
title(ax, 'D vs A per fly: the proxy is mostly the preceding wake bout', 'FontWeight', 'normal');
xline(ax, 0, ':'); yline(ax, 0, ':'); box(ax, 'on'); axis(ax, 'square');
title(tl, 'Homeostat analyses (PRESPEC 2026-11-06 A-D, amended): per0 DD, Canton-S LD and DD', 'FontSize', 11);
exportgraphics(fh, fullfile('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts', 'homeostat_tests_2026-11-06.png'), Resolution=140);
