R = load('N_results.mat'); o = R.out;
fig = figure('Color','w','Position',[50 50 1300 800]); tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
c = [.25 .45 .75; .8 .2 .1]; x = [1 1.5; 2 2.5; 3 3.5; 4 4.5];
for k = 1:2
    for m = 1:2
        ax = nexttile; hold on
        if m == 1, V = exp(o(k).L) / 60; yl = 'Typical sleep-bout length (min, geometric mean)'; else, V = o(k).A; yl = 'Sleep in that 12-h half (min)'; end
        mu = squeeze(mean(V, 1, 'omitnan')); se = squeeze(std(V, 0, 1, 'omitnan')) ./ sqrt(squeeze(sum(isfinite(V), 1)));
        for hh = 1:2
            errorbar(x(:, hh), mu(:, hh), se(:, hh), 'o', 'Color', c(hh, :), 'MarkerFaceColor', c(hh, :), 'LineWidth', 1.5, 'CapSize', 0);
        end
        plot(reshape(x.', 1, []), reshape(mu.', 1, []), '-', 'Color', [.7 .7 .7]);
        set(ax, 'XTick', 1.25:1:4.25, 'XTickLabel', {'day 1','day 2','day 3','day 4'}, 'FontSize', 11); xlim([.6 4.9])
        ylabel(yl); title(o(k).name, 'FontWeight', 'normal')
        if k == 1 && m == 1, legend({'first half (subjective day, 0-12 h)', 'second half (subjective night, 12-24 h)'}, 'Location', 'northeast'), end
    end
end
title(tl, {'Entry N: is the first day different? (mean \pm s.e. over flies, 2026-11-10)', ...
    'tim0: no half-day (LD-like) difference on day 1; bouts shorten steadily over the 4 days.   per0: day-1 night bouts longer.'}, 'FontSize', 13)
exportgraphics(fig, 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\tim0dd_first_day_2026-11-10.png', 'Resolution', 130)
