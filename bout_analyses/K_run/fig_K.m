load K_results.mat
f = figure('Color','w','Position',[100 100 1250 470]); t = tiledlayout(1,2,'TileSpacing','compact');
nexttile; hold on
M = raw.M; x = (0:size(M,2)-1) * 0.5 + 0.25;
for i = 1:size(M,1), k = isfinite(M(i,:)); plot(x(k), exp(M(i,k)), '-', 'Color', [.6 .6 .6 .35]); end
g = exp(mean(M,1,'omitnan')); n = sum(isfinite(M),1); k = n >= 10;
plot(x(k), g(k), 'k-o', 'LineWidth', 2.5, 'MarkerFaceColor','k'); yline(1,'k:'); set(gca,'YScale','log')
xlabel('Days since the record start (12-h stretches)'); ylabel({'Typical sleep-bout length in the stretch','relative to the fly''s own average (log scale)'});
title({'per0 DD, 48 flies','grey = each fly; black = group (stretches with >= 10 flies)'},'FontWeight','normal')
set(gca,'Box','off','TickDir','out','FontSize',10)
nexttile; b = bar(rej / 100, 'grouped'); b(1).FaceColor = [.85 .4 .3]; b(2).FaceColor = [.4 .55 .8];
yline(0.05,'k:'); yline(0.095,'k--'); ylim([0 1.05])
set(gca,'XTickLabel',{'No change','Shared trend +10%/day','Each fly its own trend','First day 1.5x longer','A real 12-h state'},'Box','off','TickDir','out','FontSize',10)
xtickangle(20); ylabel('How often the test says "yes" (100 simulated groups)')
legend({'"Stretches vary beyond trends" test','"Shared trend" test'},'Location','northwest','Box','off')
title({'Check on simulated flies (no state in the first 4)','should say "yes" ~5%; dashed line = 9.5% limit'},'FontWeight','normal')
title(t,{'per0 DD: the half-day variation in sleep-bout length is mostly the first day of the record, not a recurring state', ...
    'The variation test is fooled by a first-day effect (26% false "yes"), so it can only be read descriptively'},'FontSize',11)
exportgraphics(f,'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\per0dd_trend_vs_state_2026-11-09.png','Resolution',150)
