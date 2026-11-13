load J_results.mat; I = load('..\I_run\I_results.mat');
f = figure('Color','w','Position',[100 100 1250 760]); t = tiledlayout(2,3,'TileSpacing','compact');
pn = ["CSLD_day","CSLD_night","per0DD"]; tt = ["CS LD, day","CS LD, night","per0 DD"];
for i = 1:3
    nexttile(i); R = res.(pn(i)); hold on
    g1 = squeeze(mean(I.res.(pn(i)).nul(:,:,1,1),1,'omitnan')); g4 = squeeze(mean(R.n4(:,:,1,1),1,'omitnan'));
    histogram(g1,'BinWidth',.03,'FaceColor',[.35 .55 .85],'EdgeColor','none','FaceAlpha',.7);
    histogram(g4,'BinWidth',.03,'FaceColor',[.3 .7 .4],'EdgeColor','none','FaceAlpha',.7);
    xline(mean(R.obs(:,1,1),'omitnan'),'r','LineWidth',2.5); xline(1,'k:'); xlim([.6 2])
    title(tt(i),'FontWeight','normal'); xlabel({'How much more often a long bout follows a long one','than chance (1 = no clustering)'});
    if i==1, ylabel('Reshuffled groups (of 1000)'); legend({'Time of day kept (entry I)','Same day AND same 2-h window kept','Observed'},'Location','northwest','Box','off'); end
    set(gca,'Box','off','TickDir','out','FontSize',10)
    nexttile(i+3); hold on
    n1 = squeeze(mean(R.n1,1,'omitnan')); m1 = mean(n1,1); ex = squeeze(mean(R.obs(:,:,1),1,'omitnan')) - m1;
    lo = prctile(n1 - m1, 2.5, 1); hi = prctile(n1 - m1, 97.5, 1);
    fill([1:4 4:-1:1],[lo fliplr(hi)],[.35 .55 .85],'FaceAlpha',.25,'EdgeColor','none');
    plot(1:4, ex, 'r-o','LineWidth',2,'MarkerFaceColor','r'); yline(0,'k:'); xticks(1:4); xlim([.7 4.3])
    xlabel('Bouts apart (1 = the next sleep bout)'); 
    if i==1, ylabel({'Clustering beyond time of day','(observed minus time-of-day null)'}); legend({'Time-of-day null, 95% range','Observed'},'Location','northeast','Box','off'); end
    set(gca,'Box','off','TickDir','out','FontSize',10)
end
title(t,{'Day and per0 DD: long sleep bouts cluster because bout length drifts over hours, not because one bout sets up the next (night: the two nulls disagree, see text)', ...
    'Top: keeping each day''s 2-h windows reproduces daytime and per0 clustering.  Bottom: the excess does not fall off after one bout, as memory would'},'FontSize',11)
exportgraphics(f,'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\long_bout_state_vs_memory_2026-11-09.png','Resolution',150)
