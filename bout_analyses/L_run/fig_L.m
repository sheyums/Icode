load L_results.mat; R = R24;
f = figure('Color','w','Position',[100 100 1200 450]); t = tiledlayout(1,2,'TileSpacing','compact');
nexttile; hold on; cols = [.7 .7 .7; .95 .6 .2; .3 .7 .4];
for j = 1:3, histogram(R.Gn(:,1,j),'BinWidth',.025,'FaceColor',cols(j,:),'EdgeColor','none','FaceAlpha',.7); end
xline(R.G(1),'r','LineWidth',2.5); xline(1,'k:'); xlim([.6 1.6])
legend({'Shuffled (no structure)','Each 12-h stretch''s bouts kept together','Each 2-h window''s bouts kept together','Observed flies'},'Location','northwest','Box','off')
xlabel({'How much more often a long sleep bout (top 10%) follows','a long one than chance, averaged over flies'}); ylabel('Reshuffled groups (of 1000)')
title({'per0 DD, 48 flies, first 24 h excluded','Long bouts still cluster; keeping 12-h stretches explains it'},'FontWeight','normal')
set(gca,'Box','off','TickDir','out','FontSize',10)
nexttile; b = bar(rej / 100); for j = 1:3, b(j).FaceColor = cols(j,:); end
yline(0.095,'k--'); ylim([0 1.05])
set(gca,'XTickLabel',{'No change','Shared decline 2%/day','Acclimation fading over ~1 day','Each fly its own trend','A real 3-h state'},'Box','off','TickDir','out','FontSize',10); xtickangle(20)
ylabel({'Share of 100 simulated groups','where the test says "beyond chance"'})
title({'Check on simulated flies: no false "yes" from slow changes (first 4);','a real 3-h state is caught 60% of the time by the 12-h test'},'FontWeight','normal')
title(t,'After the first day, per0 sleep-bout clustering is slow differences between 12-h stretches, not a state lasting hours','FontSize',11)
exportgraphics(f,'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\per0dd_after_day1_2026-11-09.png','Resolution',150)
