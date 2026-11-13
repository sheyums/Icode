out = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\csld_run';
load(fullfile(out, 'csld_results.mat'), 'res');
fh = figure('Color', 'w', 'Position', [50 50 1300 760]);
try theme(fh, 'light'); catch, end
tl = tiledlayout(fh, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
types = ["sleep","wake"]; lags = [1 2 4 8]; rs = RandStream('mt19937ar', 'Seed', 2);
for ti = 1:2
  R = res.(types(ti)); T = R.PerFly; L = R.Settings.Lags;
  for panel = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    if panel == 1, st = 'UpFraction'; ylab = {'Fraction of increases','(data minus own null)'};
    else, st = 'Rho'; ylab = {'Memory: rank correlation','(data minus own null)'}; end
    for li = 1:4
      r = T(T.Lag == lags(li), :);
      nm = squeeze(mean(R.PerFlyNull.strata.(st)(:, L == lags(li), :), 3, 'omitnan'));
      d = r.(st) - nm; sig = r.([st '_P_Strata']) < 0.05;
      x = lags(li) * 2.^(0.15 * (rs.rand(size(d)) - 0.5));
      scatter(ax, x(~sig), d(~sig), 14, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.6);
      scatter(ax, x(sig), d(sig), 30, [0.85 0.3 0.1], 'filled');
      plot(ax, lags(li) + [-0.12 0.12]*lags(li), median(d)*[1 1], 'k-', 'LineWidth', 2.5);
      text(ax, lags(li), 0, '', 'Visible', 'off');
      ns(li) = sum(sig); %#ok<SAGROW>
    end
    yline(ax, 0, ':k');
    yl = ylim(ax);
    for li = 1:4
      text(ax, lags(li), yl(2) - 0.06*diff(yl), sprintf('%d/%d', ns(li), height(T)/4), 'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    set(ax, 'XScale', 'log', 'XTick', lags, 'XLim', [0.75 11]);
    xlabel(ax, 'Bouts apart (lag k)'); ylabel(ax, ylab);
    title(ax, sprintf('%s bouts', upper(extractBefore(types(ti), 2)) + extractAfter(types(ti), 1)), 'FontWeight', 'normal');
    box(ax, 'on');
  end
end
title(tl, {'Canton-S, 12:12 LD, 40 flies: order beyond each fly''s day/night profile', ...
  'dots = flies (orange: p < 0.05 vs its own profile-matched null); black bar = median; n/40 = flies significant (5% expected by chance: 2/40)'}, 'FontSize', 10);
exportgraphics(fh, fullfile(out, 'csld_time_arrow.png'), Resolution=150);
