S = load('Q_results.mat'); R = S.R; ds = S.ds(:); c2 = load('Qcal2.mat'); c3 = load('Qcal3.mat');
rs = RandStream('mrg32k3a', 'Seed', 11); nB = 2000; Fb = NaN(nB, 3);
for b = 1:nB
    ix = []; for k = 1:3, m = find(ds == k); ix = [ix; m(randi(rs, numel(m), numel(m), 1))]; end %#ok<AGROW>
    Fb(b, :) = mean(R.fast(ix, 1:3), 1, 'omitnan');
end
fig = figure('Color','w','Position',[50 50 1650 560]); tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
% A: raw / slow / fast by lag
ax = nexttile; hold on
raw = mean(R.obs(:, 1:3), 'omitnan'); slow = mean(R.nullmean(:, 1:3), 'omitnan'); fast = mean(R.fast(:, 1:3), 'omitnan');
b = bar([raw; slow; fast].', 'grouped', 'EdgeColor', 'none'); b(1).FaceColor = [.35 .35 .35]; b(2).FaceColor = [.55 .7 .9]; b(3).FaceColor = [.8 .2 .1];
q = quantile(Fb, [.025 .975]); errorbar(b(3).XEndPoints, fast, fast - q(1, :), q(2, :) - fast, 'k.', 'CapSize', 0, 'LineWidth', 1.2)
yline(0); set(ax, 'XTick', 1:3, 'XTickLabel', {'1 step: sleep↔wake neighbours', '2 steps: sleep→sleep, wake→wake', '3 steps'}, 'FontSize', 10); xtickangle(10)
ylabel('Rank correlation (mean over 113 flies)'); legend(b, {'observed', 'slow part (kept by 2-h shuffle)', 'fast part, 95% range'}, 'Location', 'northeast')
title({'The fast link reaches the next bout only', 'Two bouts away, what remains is the slow layer'}, 'FontWeight', 'normal')
% B: kappa vs calibration
ax = nexttile; hold on
edges = -1.5:0.1:2.5;
histogram(c2.K, edges, 'Normalization', 'probability', 'FaceColor', [.55 .7 .9], 'EdgeColor', 'none');
histogram(c3.K, edges, 'Normalization', 'probability', 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none', 'FaceAlpha', .6);
ylim([0 0.45]); yl = ylim; patch([S.ci(1) S.ci(2) S.ci(2) S.ci(1)], [0 0 yl(2) yl(2)], [.8 .2 .1], 'FaceAlpha', .15, 'EdgeColor', 'none');
xline(S.kap, '-', sprintf('data %.2f', S.kap), 'Color', [.8 .2 .1], 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle'); xline(0.5, '--', 'decision line 0.5', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom')
xlim([-1.5 2.5]); set(ax, 'FontSize', 10)
xlabel('\kappa = share of the 1-step fast link left after 2 steps'); ylabel('Fraction of simulated groups')
legend({'simulated: fast state lasting several bouts', 'simulated: one-step link', 'data, 95% range'}, 'Location', 'northwest')
title({'Data sit with the one-step simulation', 'never called "one-step" when a multi-bout state was built in'}, 'FontWeight', 'normal')
% C: per dataset
ax = nexttile; hold on; dn = {'per0 DD','CS LD','CS DD'}; M = NaN(3, 2);
for k = 1:3, M(k, :) = mean(R.fast(ds == k, 1:2), 'omitnan'); end
b = bar(M, 'grouped', 'EdgeColor', 'none'); b(1).FaceColor = [.8 .2 .1]; b(2).FaceColor = [.95 .7 .6]; yline(0); ylim([-0.025 0.08])
set(ax, 'XTick', 1:3, 'XTickLabel', dn, 'FontSize', 11); ylabel('Fast part of rank correlation')
legend(b, {'1 step (sleep↔wake)', '2 steps (sleep→sleep, wake→wake)'}, 'Location', 'northeast')
title({'Per dataset (descriptive)', 'one-step link in all three; 2-step part small, sign varies'}, 'FontWeight', 'normal')
title(tl, 'Entry Q: sleep bouts are linked to the wake bouts beside them, not to the next sleep bout (2026-11-11; per0 DD + CS DD from 24 h, CS LD)', 'FontSize', 13)
exportgraphics(fig, 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\one_step_link_2026-11-11.png', 'Resolution', 120)
