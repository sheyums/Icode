function fh = plotSqwOnsager(O)
%PLOTSQWONSAGER  EPR against J with the fitted EPR = J^2/L per genotype (SQWONSAGERGROUP output).
%
%   fh = plotSqwOnsager(O)
%
%   One panel per genotype: the fly x phase points and the curve J^2/L
%   through the origin; the title carries L and its fly-bootstrap interval.

arguments
    O table
end
P = O.Properties.UserData.points;
if O.Properties.UserData.pointType == "count", kJ = "J"; kE = "EPR"; else, kJ = "J_chain"; kE = "EPR_chain"; end
fh = figure('Name', "Onsager (" + O.Properties.UserData.pointType + " points)", 'Color', 'w', ...
    'Position', [100 100 1100 650]);
try, theme(fh, 'light'); catch, end
n = height(O); nc = ceil(n / 2);
for g = 1:n
    ax = subplot(2, nc, g, 'Parent', fh); hold(ax, 'on');
    R = P(P.genotype == O.genotype(g), :);
    scatter(ax, R.(kJ), R.(kE), 16, 'filled', 'MarkerFaceAlpha', 0.6);
    jj = linspace(min([R.(kJ); 0]), max([R.(kJ); 0]), 200);
    plot(ax, jj, jj .^ 2 / O.L(g), 'k-', 'LineWidth', 1.5);
    xlabel(ax, 'J (cycles/h)'); ylabel(ax, 'EPR (nats/h)'); grid(ax, 'on'); box(ax, 'on');
    title(ax, sprintf('%s  L = %.3f [%.3f, %.3f]', O.genotype(g), O.L(g), O.L_lo(g), O.L_hi(g)));
end
end
