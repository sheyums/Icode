sc = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad';
Lp = load(fullfile(sc,'L_run','L_results.mat'),'R24'); Lt = load(fullfile(sc,'O_run','OL_results.mat'),'R24');
Mp = load(fullfile(sc,'M_run','M_results.mat'),'out'); Mt = load(fullfile(sc,'O_run','OM_results.mat'),'out');
LL = {Lp.R24, Lt.R24}; MM = {Mp.out.per0DD, Mt.out.tim0DD};
gn = ["per0 DD, from 24 h (48 flies)", "tim0 DD, from 24 h (25 flies)"];
nlab = {'whole record','same 12-h stretch','same 2-h window'};
fig = figure('Color','w','Position',[50 50 1750 850]); tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for k = 1:2
    % 1: long-bout clustering (entry L)
    ax = nexttile; hold on; R = LL{k};
    for j = 1:3, g = R.Gn(:,1,j); q = quantile(g,[.025 .975]); plot([j j], q, '-', 'Color', [.6 .6 .6], 'LineWidth', 6); plot(j, mean(g), 'o', 'Color', [.4 .4 .4], 'MarkerFaceColor', [.4 .4 .4]); end
    ylim([0.7 1.6]); yline(R.G(1), '-', 'observed', 'Color', [.8 .2 .1], 'LineWidth', 2, 'LabelHorizontalAlignment','left'); yline(1, ':')
    set(ax,'XTick',1:3,'XTickLabel',nlab,'FontSize',10); xlim([.5 3.5]); xtickangle(20)
    ylabel({gn(k), 'Long sleep follows long sleep (\times chance)'}); title('Long sleep bouts cluster?', 'FontWeight','normal')
    % 2: coupling vs nulls (entry M)
    ax = nexttile; hold on; R = MM{k};
    for j = 1:3, g = R.Gn(:,1,j); q = quantile(g,[.025 .975]); plot([j j], q, '-', 'Color', [.6 .6 .6], 'LineWidth', 6); plot(j, mean(g), 'o', 'Color', [.4 .4 .4], 'MarkerFaceColor', [.4 .4 .4]); end
    ylim([-0.03 0.15]); yline(R.G(1), '-', 'observed', 'Color', [.8 .2 .1], 'LineWidth', 2, 'LabelHorizontalAlignment','left'); yline(0, ':')
    set(ax,'XTick',1:3,'XTickLabel',nlab,'FontSize',10); xlim([.5 3.5]); xtickangle(20)
    ylabel('Rank correlation, neighbouring sleep and wake'); title('Neighbouring sleep and wake bouts coupled?', 'FontWeight','normal')
    % 3: direction
    ax = nexttile; hold on; v = R.obs(:,2:3); ok = all(isfinite(v),2);
    plot(repmat([1;2],1,sum(ok)), v(ok,:).', '-', 'Color', [.75 .75 .85]); plot([1 2], mean(v(ok,:)), 'ko-', 'MarkerFaceColor','k', 'LineWidth', 2); yline(0,':')
    set(ax,'XTick',[1 2],'XTickLabel',{'wake → next sleep','sleep → next wake'},'FontSize',10); xlim([.6 2.4]); ylim([-0.3 0.5])
    ylabel('Rank correlation (one line per fly)'); title(sprintf('Same in both directions?  %.3f vs %.3f', mean(v(ok,1)), mean(v(ok,2))), 'FontWeight','normal')
    % 4: long-wake tail
    ax = nexttile; hold on; c = R.obs(:,6:7); ok = all(isfinite(c),2);
    plot(repmat([1;2],1,sum(ok)), c(ok,:).', '-', 'Color', [.85 .8 .75]); plot([1 2], mean(c(ok,:)), 'o-', 'Color', [.8 .2 .1], 'MarkerFaceColor', [.8 .2 .1], 'LineWidth', 2); yline(1,':','chance')
    set(ax,'XTick',[1 2],'XTickLabel',{'sleep AFTER long wake','sleep BEFORE long wake'},'FontSize',10); xlim([.6 2.4]); ylim([0 6])
    ylabel('How often that sleep is long (\times chance)'); title(sprintf('After a long wake, next sleep long?  diff %.2f', R.G(4)), 'FontWeight','normal')
end
title(tl, {'Entry O: does tim0 DD repeat the per0 DD results? (2026-11-10; bouts from 24 h; gray = shuffles that keep timing, 95% range)', ...
    'Replicates: neighbouring bouts coupled beyond 2-h shuffles.   Does not: long-bout clustering, long wake \rightarrow long sleep.   tim0 coupling is stronger sleep \rightarrow wake.'}, 'FontSize', 13)
exportgraphics(fig, 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\tim0_vs_per0_2026-11-10.png', 'Resolution', 120)
