D = 'C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data\per0 ld';
f = dir(fullfile(D, '*.txt')); names = string({f.name});
exps = extractBefore(names, regexpPattern('_(multibeam|4b)'));
ue = unique(exps, 'stable');
fh = figure('Color', 'w', 'Position', [50 50 1300 700]); try theme(fh, 'light'); catch, end
tl = tiledlayout(fh, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for e = ue
  idx = find(exps == e); P = [];
  for i = idx
    x = readmatrix(fullfile(D, names(i)), 'FileType', 'text'); x = x(:, 1); x = x(~isnan(x));
    nd = floor(numel(x) / 86400);
    y = reshape(x(1:nd*86400), 900, 96, nd); P(:, :, end+1) = squeeze(sum(y, 1)); %#ok<AGROW>
  end
  P = P(:, :, 2:end);
  ax = nexttile(tl); hold(ax, 'on');
  dayMean = squeeze(mean(P, 3));                 % 96 x nd
  imagesc(ax, (0:95)/4 + 0.125, 1:size(dayMean, 2), dayMean.'); axis(ax, 'ij');
  colormap(ax, flipud(gray)); xlim(ax, [0 24]); ylim(ax, [0.5 size(dayMean, 2) + 0.5]);
  [~, o] = sort(mean(dayMean, 2), 'descend');
  xline(ax, (o(1:2) - 1)/4, 'r-', 'LineWidth', 1);
  xlabel(ax, 'Hours from record start (each row = one day, fly mean)'); ylabel(ax, 'Day');
  title(ax, sprintf('%s (%d flies): largest bins at %.2f and %.2f h', e, numel(idx), (o(1:2) - 1)/4), ...
      'FontWeight', 'normal', 'Interpreter', 'none');
  set(ax, 'XTick', 0:3:24); box(ax, 'on');
end
title(tl, {'per0 LD: activity per 15 min by hour since the record started (darker = more beam crossings)', ...
  'Lights-on and lights-off startles appear as vertical stripes 12 h apart; a ZT0 start puts them at 0 and 12 h'}, 'FontSize', 10);
exportgraphics(fh, fullfile('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts', 'per0ld_phase_check_2026-11-06.png'), Resolution=130);
