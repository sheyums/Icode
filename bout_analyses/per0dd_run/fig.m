out = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\per0dd_run';
load(fullfile(out, 'per0dd_results.mat'), 'res');
fh = figure('Color', 'w', 'Position', [50 50 1350 760]);
try theme(fh, 'light'); catch, end
tl = tiledlayout(fh, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
stats = ["UpFraction","Bowley","Rho"];
ylabs = ["Fraction of changes that are increases", "Skewness of the changes (Bowley)", "Memory: rank correlation"];
nullOf = ["iaaft","iaaft","shuffle"]; centre = [0.5 0 0];
types = ["sleep","wake"];
lags = [1 2 4 8];
rs = RandStream('mt19937ar', 'Seed', 1);
for ti = 1:2
  R = res.(types(ti)); T = R.PerFly;
  for s = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    b = R.GroupNull.(nullOf(s)).(stats(s) + "Band");
    fill(ax, [lags fliplr(lags)], [b(:,1).' fliplr(b(:,2).')], [0.55 0.7 0.95], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
    for li = 1:4
      v = T.(stats(s))(T.Lag == lags(li));
      x = lags(li) * 2.^(0.12 * (rs.rand(size(v)) - 0.5));
      scatter(ax, x, v, 10, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.5);
    end
    plot(ax, lags, R.Group.(stats(s)), 'k-o', 'LineWidth', 1.6, 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    yline(ax, centre(s), ':k');
    set(ax, 'XScale', 'log', 'XTick', lags, 'XLim', [0.8 10]);
    if s == 3, pcol = R.Group.Rho_P_Shuffle; else, pcol = R.Group.(stats(s) + "_P_IAAFT"); end
    lo = min([b(:); R.Group.(stats(s))]); hi = max([b(:); R.Group.(stats(s))]); pad = 0.6*(hi - lo);
    ylim(ax, [lo - pad, hi + pad]);
    for li = 1:4
      text(ax, lags(li), hi + 0.45*pad, sprintf('p=%.3f', pcol(li)), 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    xlabel(ax, 'Bouts apart (lag k)'); ylabel(ax, ylabs(s));
    title(ax, sprintf('%s bouts', upper(extractBefore(types(ti), 2)) + extractAfter(types(ti), 1)), 'FontWeight', 'normal');
    box(ax, 'on');
  end
end
title(tl, {'per0, constant darkness, 48 flies: is there a time arrow in successive bouts?', ...
  'black = mean over flies; grey = individual flies (y-range clipped to the band); blue band = 95% of the mean under the null (IAAFT for the first two columns, shuffle for memory)'}, 'FontSize', 10);
exportgraphics(fh, fullfile(out, 'per0dd_time_arrow.png'), Resolution=150);
