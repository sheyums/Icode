S = load('M_results.mat'); o = S.out;
nm = ["per0DD","CSLD"]; tt = ["per0, constant darkness (from 24 h)","Wild type (CS), light-dark"];
fig = figure('Color','w','Position',[50 50 1500 900]); tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nlab = containers.Map({'P','S12','W2','Z'},{'whole record','same 12-h stretch','same 2-h window','same clock hour'});
for k = 1:2
    R = o.(nm(k)); nl = R.nulls;
    % A: coupling vs nulls
    ax = nexttile; hold on
    for j = 1:numel(nl)
        g = R.Gn(:,1,j); q = quantile(g,[.025 .975]);
        plot([j j], q, '-', 'Color', [.6 .6 .6], 'LineWidth', 6); plot(j, mean(g), 'o', 'Color', [.4 .4 .4], 'MarkerFaceColor', [.4 .4 .4]);
    end
    yl = ylim; ylim([yl(1) max(yl(2), R.G(1) * 1.2)]); yline(R.G(1), '-', 'observed', 'Color', [.8 .2 .1], 'LineWidth', 2, 'LabelHorizontalAlignment','left');
    set(ax,'XTick',1:numel(nl),'XTickLabel',arrayfun(@(s) nlab(char(s)), nl,'UniformOutput',false),'FontSize',11); xlim([.5 numel(nl)+.5]); xtickangle(20)
    ylabel('Rank correlation, neighbouring sleep and wake bouts'); title({tt(k), 'Neighbouring bouts vs shuffles that keep timing'}, 'FontWeight','normal')
    % B: direction
    ax = nexttile; hold on
    v = R.obs(:,2:3); ok = all(isfinite(v),2);
    plot(repmat([1;2],1,sum(ok)), v(ok,:).', '-', 'Color', [.75 .75 .85]);
    plot([1 2], mean(v(ok,:)), 'ko-', 'MarkerFaceColor','k', 'LineWidth', 2); yline(0,':')
    set(ax,'XTick',[1 2],'XTickLabel',{'wake → next sleep','sleep → next wake'},'FontSize',11); xlim([.6 2.4])
    ylabel('Rank correlation (one line per fly)'); title({'Same coupling in both directions?', sprintf('mean %.3f vs %.3f', mean(v(ok,1)), mean(v(ok,2)))}, 'FontWeight','normal')
    % C: long wake tail
    ax = nexttile; hold on
    c = R.obs(:,6:7); ok = all(isfinite(c),2);
    plot(repmat([1;2],1,sum(ok)), c(ok,:).', '-', 'Color', [.85 .8 .75]);
    plot([1 2], mean(c(ok,:)), 'o-', 'Color', [.8 .2 .1], 'MarkerFaceColor', [.8 .2 .1], 'LineWidth', 2); yline(1,':','chance')
    set(ax,'XTick',[1 2],'XTickLabel',{'sleep AFTER a long wake','sleep BEFORE a long wake'},'FontSize',11); xlim([.6 2.4]); ylim([0 6])
    ylabel('How often that sleep is long (× chance)'); title({'After a long wake, is the next sleep long?', sprintf('difference %.2f', R.G(4))}, 'FontWeight','normal')
end
title(tl, 'Entry M: neighbouring sleep and wake bouts (2026-11-09). "Long" = a fly''s top 10% of that bout type', 'FontSize', 14)
exportgraphics(fig, 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\sleep_wake_coupling_2026-11-09.png', 'Resolution', 130)
