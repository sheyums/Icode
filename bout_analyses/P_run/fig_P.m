sc = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad';
Mp = load(fullfile(sc,'M_run','M_results.mat'),'out'); Pp = load(fullfile(sc,'P_run','P_results.mat'),'out');
R = Pp.out.CSDD;
fig = figure('Color','w','Position',[50 50 1600 950]); tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nlab = {'whole record','same 12-h stretch','same 2-h window','same hour since start'};
% A: CS DD coupling vs nulls
ax = nexttile; hold on
for j = 1:4, g = R.Gn(:,1,j); q = quantile(g,[.025 .975]); plot([j j], q, '-', 'Color', [.6 .6 .6], 'LineWidth', 6); plot(j, mean(g), 'o', 'Color', [.4 .4 .4], 'MarkerFaceColor', [.4 .4 .4]); end
ylim([-0.08 0.06]); yline(R.G(1), '-', 'observed', 'Color', [.8 .2 .1], 'LineWidth', 2, 'LabelHorizontalAlignment','left'); yline(0, ':')
set(ax,'XTick',1:4,'XTickLabel',nlab,'FontSize',10); xlim([.5 4.5]); xtickangle(20)
ylabel('Rank correlation, neighbouring sleep and wake'); title({'CS DD from 24 h (25 flies)', 'Neighbouring bouts vs shuffles that keep timing'}, 'FontWeight','normal')
% B: direction
ax = nexttile; hold on; v = R.obs(:,2:3); ok = all(isfinite(v),2);
plot(repmat([1;2],1,sum(ok)), v(ok,:).', '-', 'Color', [.75 .75 .85]); plot([1 2], mean(v(ok,:)), 'ko-', 'MarkerFaceColor','k', 'LineWidth', 2); yline(0,':')
set(ax,'XTick',[1 2],'XTickLabel',{'wake → next sleep','sleep → next wake'},'FontSize',10); xlim([.6 2.4])
ylabel('Rank correlation (one line per fly)'); title({'CS DD: same in both directions?', sprintf('%.3f vs %.3f (p 0.16)', mean(v(ok,1)), mean(v(ok,2)))}, 'FontWeight','normal')
% C: long-wake tail
ax = nexttile; hold on; c = R.obs(:,6:7); ok = all(isfinite(c),2);
plot(repmat([1;2],1,sum(ok)), c(ok,:).', '-', 'Color', [.85 .8 .75]); plot([1 2], mean(c(ok,:)), 'o-', 'Color', [.8 .2 .1], 'MarkerFaceColor', [.8 .2 .1], 'LineWidth', 2); yline(1,':','chance')
set(ax,'XTick',[1 2],'XTickLabel',{'sleep AFTER long wake','sleep BEFORE long wake'},'FontSize',10); xlim([.6 2.4]); ylim([0 6])
ylabel('How often that sleep is long (\times chance)'); title({'CS DD: after a long wake, next sleep long?', sprintf('diff %.2f (descriptive: clock confound)', R.G(4))}, 'FontWeight','normal')
% D: slow vs fast split across datasets
ax = nexttile([1 2]); hold on
D = {Mp.out.per0DD, Mp.out.CSLD, R}; dn = {'per0 DD (no clock), from 24 h','CS LD (clock + light)','CS DD (clock, no light), from 24 h'};
B = zeros(3, 3); E = zeros(3, 2);
for k = 1:3, g = D{k}.Gn(:,1,3); B(k,:) = [D{k}.G(1), mean(g), D{k}.G(1) - mean(g)]; E(k,:) = quantile(D{k}.G(1) - g, [.025 .975]); end
b = bar(B, 'grouped', 'EdgeColor', 'none'); b(1).FaceColor = [.35 .35 .35]; b(2).FaceColor = [.55 .7 .9]; b(3).FaceColor = [.8 .2 .1];
errorbar(b(3).XEndPoints, B(:,3), B(:,3) - E(:,1), E(:,2) - B(:,3), 'k.', 'CapSize', 0, 'LineWidth', 1.2)
for j = 1:3, text(b(j).XEndPoints, B(:,j) + sign(B(:,j) + eps) * .006, compose('%.3f', B(:,j)), 'HorizontalAlignment','center', 'FontSize', 9), end
yline(0, '-'); set(ax, 'XTick', 1:3, 'XTickLabel', dn, 'FontSize', 11); ylim([-0.04 0.13])
ylabel('Rank correlation, neighbouring sleep and wake')
legend(b, {'observed coupling', 'slow part (kept by 2-h shuffle)', 'fast part (observed minus slow), 95% range'}, 'Location', 'northeast')
title('Split of the coupling: the fast part is similar everywhere; the slow part is positive without a clock, zero or negative with one', 'FontWeight', 'normal')
% E: text
ax = nexttile; axis off
text(0, .95, {'Reading', '', '\bullet Fast part (red): sleep and wake bouts next to each', '   other share something on < 2 h. Present in all three', '   (0.04-0.06), so it does not need a clock or light.', '', ...
    '\bullet Slow part (blue): in per0 a slow "tempo" lengthens', '   sleep AND wake together (+0.047). With a clock the', '   slow part is negative; untested guess: rest phase =', '   long sleep + short wake, active phase the reverse.', '', ...
    '\bullet Long wake \rightarrow long sleep: per0 0.54 (tested);', '   CS DD 0.56 (descriptive, clock not removed);', '   CS LD 0.26 (descriptive).'}, 'VerticalAlignment', 'top', 'FontSize', 11)
title(tl, 'Entry P: neighbouring sleep and wake bouts in CS DD, with per0 DD and CS LD for comparison (2026-11-11)', 'FontSize', 14)
exportgraphics(fig, 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\csdd_sleep_wake_coupling_2026-11-11.png', 'Resolution', 120)
