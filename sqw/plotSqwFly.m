function fh = plotSqwFly(F)
%PLOTSQWFLY  Symbol raster and per-phase numbers of one fly (SQWOBSERVABLEFLY output).
%
%   fh = plotSqwFly(F)
%
%   Top: one row per day, one column per minute from ZT0; Q light green,
%   W dark green, S blue, missing grey, unknown wake orange. Bottom: the
%   phase's count J, count EPR, c_obs and the chain's J and EPR.

arguments
    F (1,1) struct
end
fh = figure('Name', "sqw fly " + F.file, 'Color', 'w');
try, theme(fh, 'light'); catch, end                 % pin before any axes (Icode CLAUDE.md)
nd = F.T / 1440;
img = reshape(F.y, 1440, nd)';
cmap = [0.75 0.75 0.75; 0.55 0.85 0.45; 0.10 0.45 0.15; 0.20 0.35 0.85; 0.95 0.60 0.10];
ax = subplot(2, 1, 1, 'Parent', fh);
image(ax, (0:1439) / 60, 1:nd, img + 1);
colormap(ax, cmap); clim(ax, [0.5 5.5]);
xlabel(ax, 'ZT (h)'); ylabel(ax, 'day'); xlim(ax, [0 24]); xticks(ax, 0:6:24);
[~, nm] = fileparts(F.file);
title(ax, sprintf('%s  (Q light green, W dark green, S blue, missing grey)', nm), 'Interpreter', 'none');
for p = 1:F.nP
    if F.nP == 2 && p == 2, xline(ax, 12, 'k-'); end
end
ax2 = subplot(2, 1, 2, 'Parent', fh);
v = [[F.Q.J]; [F.Q.EPR]; [F.Q.c_obs]; [F.Q.J_chain]; [F.Q.EPR_chain]]';
bar(ax2, categorical(cellstr(F.phaseNames)), v);
legend(ax2, {'count J (cyc/h)', 'count EPR (nats/h)', 'c_{obs}', 'chain J', 'chain EPR'}, 'Location', 'eastoutside');
ylabel(ax2, 'value'); grid(ax2, 'on');
title(ax2, sprintf('reference OK: %s', mat2str(F.refOK')));
end
