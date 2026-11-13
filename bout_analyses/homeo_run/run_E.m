cd('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode');
assert(isequal(fileparts(which('SleepVarianceTime')), pwd) && isequal(fileparts(which('RenewalSurrogate')), pwd));
load('C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\homeo_run\homeo_results.mat', 'flies');
F = flies.per0DD;
R = SleepVarianceTime(F, NumSurrogates=1000, Seed=1, FlyNames=[F.Name]);
disp(R.Group)
fprintf('mean logR at 1 2 4 8 16 24 h: %s\n', mat2str(round(R.Group.MeanLogR, 3)));
fprintf('per-fly slope: median %.3f, IQR %.3f to %.3f; flies with p_rev < 0.05: %d/%d; slopes < -0.10: %d\n', ...
    median(R.Slope), prctile(R.Slope, 25), prctile(R.Slope, 75), sum(R.PerFly.Slope_P_Reversion < 0.05), height(R.PerFly), sum(R.Slope < -0.10));
save('C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\homeo_run\E_results.mat', 'R');
fh = figure('Color', 'w', 'Position', [50 50 1100 450]); try theme(fh, 'light'); catch, end
tl = tiledlayout(fh, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, R.ScalesHours, R.LogR.', '-', 'Color', [0.75 0.75 0.75]);
plot(ax, R.ScalesHours, R.Group.MeanLogR, 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
plot(ax, R.ScalesHours, [-0.002 -0.009 -0.021 -0.052 -0.103 -0.16], 'b--', 'LineWidth', 1.2);
plot(ax, R.ScalesHours, [-0.1 -0.338 -0.668 -0.897 -1.064 -1.164], 'r--', 'LineWidth', 1.2);
plot(ax, R.ScalesHours, [-0.051 -0.09 -0.119 -0.142 -0.199 -0.236], 'g--', 'LineWidth', 1.2);
set(ax, 'XScale', 'log', 'XTick', R.ScalesHours); yline(ax, 0, ':');
xlabel(ax, 'Timescale (h)'); ylabel(ax, {'log( variance of the sleep surplus /', 'same bouts in random order )'});
legend(ax, {'flies'}, 'Location', 'southwest'); lg = ax.Legend; lg.String = lg.String(1);
title(ax, {'per0 DD, 48 flies: does the sleep surplus wander less than chance?', ...
  'black: mean; dashed: simulated groups -- blue independent, red homeostat, green shared state'}, 'FontWeight', 'normal');
box(ax, 'on');
ax = nexttile(tl); histogram(ax, R.Slope, 15, 'FaceColor', [0.3 0.45 0.7]); xline(ax, 0, ':'); xline(ax, -0.10, 'r--');
xlabel(ax, 'Per-fly slope of the curve (negative = mean reversion)'); ylabel(ax, 'Flies');
title(ax, sprintf('group mean slope %.3f, p(reversion) %.3f', R.Group.MeanSlope, R.Group.P_Reversion), 'FontWeight', 'normal'); box(ax, 'on');
exportgraphics(fh, 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\sleep_mean_reversion_2026-11-06.png', Resolution=140);
