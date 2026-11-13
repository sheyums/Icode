load H_results.mat
f = figure('Color','w','Position',[100 100 1100 720]); t = tiledlayout(2,2,'TileSpacing','compact');
ph = ["day","night"]; tl = ["Long bouts (top 10%)","Short bouts (bottom 10%)"];
for i = 1:2, for k = 1:2
    nexttile; o = out.(ph(i));
    g = squeeze(mean(o.fitR(:,:,k),1,'omitnan')); g0 = squeeze(mean(o.p0R(:,:,k),1,'omitnan'));
    histogram(g0,20,'FaceColor',[.75 .75 .75],'EdgeColor','none'); hold on
    histogram(g,20,'FaceColor',[.35 .55 .85],'EdgeColor','none');
    xo = mean(o.obs(:,k),'omitnan'); xline(xo,'r','LineWidth',2.5); xline(1,'k:');
    title(sprintf('%s, %s', upper(extractBefore(ph(i),2))+extractAfter(ph(i),1), tl(k)),'FontWeight','normal');
    xlabel({'How much more often an extreme bout follows an extreme bout','than chance (1 = no clustering), averaged over flies'});
    ylabel('Simulated groups');
    if i==1 && k==1, legend({'His model, no memory (\phi=0)','His fitted model (Gaussian memory)','Observed flies'},'Location','northeast','Box','off'); end
    set(gca,'Box','off','TickDir','out','FontSize',10)
end, end
title(t,{'CS LD sleep bouts: long bouts follow long bouts more often than the Gaussian latent model allows; short bouts do not','(each histogram: 200 simulated groups with the same flies, bout counts and bout distributions)'},'FontSize',11)
exportgraphics(f,'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\extreme_bout_clustering_2026-11-09.png','Resolution',150)
