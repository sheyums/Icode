load I_results.mat
f = figure('Color','w','Position',[100 100 1250 430]); t = tiledlayout(1,3,'TileSpacing','compact');
pn = ["CSLD_day","CSLD_night","per0DD"]; tt = ["CS LD, day (ZT 0-12)","CS LD, night (ZT 12-24)","per0 DD (no clock, no light cycle)"];
cols = [.7 .7 .7; .35 .55 .85; .95 .6 .2]; js = [3 1 2];
nm = ["Shuffled (no structure)","Time of day kept (1-h bins)","Each 12-h segment's level kept"];
for i = 1:3
    nexttile; R = res.(pn(i)); hold on
    for q = 1:3
        g = squeeze(mean(R.nul(:,:,1,js(q)),1,'omitnan'));
        histogram(g, 'BinWidth', 0.03, 'FaceColor', cols(q,:), 'EdgeColor','none', 'FaceAlpha', .7);
    end
    xline(mean(R.obs(:,1),'omitnan'),'r','LineWidth',2.5); xline(1,'k:');
    title(tt(i),'FontWeight','normal'); xlim([0.55 1.9]);
    xlabel({'How much more often a long sleep bout (top 10%)','follows a long one than chance, averaged over flies'});
    if i==1, ylabel('Reshuffled groups (of 1000)'); legend([nm "Observed flies"],'Location','northeast','Box','off'); end
    set(gca,'Box','off','TickDir','out','FontSize',10)
end
title(t,{'Clustering of long sleep bouts: night = time of day; day = beyond time of day; per0 DD = differences between 12-h stretches', ...
    '(our flies: CS LD 40, per0 DD 48; each histogram: the same bouts reshuffled 1000 times, keeping the stated structure)'},'FontSize',11)
exportgraphics(f,'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\long_bout_clustering_source_2026-11-09.png','Resolution',150)
