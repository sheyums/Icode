function X = exportFigA2(D, outDir, options)
%EXPORTFIGA2  Fig A2, EPR vs J and the Onsager conductance L (FIGURE_RECIPES K2; fig_planA.py): data and figure.
%
%   X = exportFigA2(D, outDir)
%   X = exportFigA2(D, outDir, NBoot=2000, RandomSeed=20, Plot=true)
%
%   D       PLANADATA output
%   outDir  folder for figA2_points.csv, figA2_curves.csv, figA2_L.csv, figA2.mat
%           and the figure (figA2_epr_j_onsager_observable.png / .pdf)
%
%   (a) one panel per genotype: counted EPR against counted J per fly x phase
%       (fly-phases with a Q/W reference and a chain with a cycle, finite values;
%       both phases pooled), with EPR = J^2 / L, L from OLS through the origin
%       (ONSAGERFIT). Light open, dark / DD filled.
%   (b) L per genotype with its fly-bootstrap 95% CI (flies resampled with all their
%       points; SQWONSAGERGROUP), dashed line at CS.
%   X  struct: points, curves (J grid and fitted EPR per genotype), L (the
%      SQWONSAGERGROUP table, incl. the bootstrap difference against CS), fig
%   Check (recipe): L = 0.272 (CS), 0.262 (CS DD), 0.286 (per0 LD), 0.352 (per0 DD),
%   0.243 (fumin), 0.145 (sss).
%
%   See also SQWONSAGERGROUP, ONSAGERFIT, EXPORTFIGA1.

arguments
    D (1,1) struct
    outDir (1,1) string
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
    options.Plot (1,1) logical = true
end
if ~isfolder(outDir), mkdir(outDir); end
S = planAStyle();
assert(isequal(D.order, S.order), 'exportFigA2:order', 'D must hold the genotypes in the order %s.', strjoin(S.order));
O = sqwOnsagerGroup(D.G, Points="count", Reference="CS", WaldSet=["per0LD" "per0DD" "fumin" "sss"], ...
    NBoot=options.NBoot, RandomSeed=options.RandomSeed);
P = O.Properties.UserData.points;
P.filled = P.phaseName ~= "light";
C = table();
for k = 1:6
    q = P(P.genotype == S.order(k), :);
    xx = linspace(min(-0.05, min(q.J)), max(0.4, max(q.J)), 100)';
    C = [C; table(repmat(S.order(k), 100, 1), xx, O.s(k) * xx .^ 2, 'VariableNames', {'genotype', 'J', 'EPR_fit'})]; %#ok<AGROW>
end
X.points = P; X.curves = C; X.L = O; X.fig = [];
X.boot = O.Properties.UserData.boot;

if options.Plot
    fh = planAFigure("Fig A2", 7.0, 5.0); X.fig = fh;
    tl = tiledlayout(fh, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:6
        ax = nexttile(tl, (floor((k - 1) / 3)) * 4 + mod(k - 1, 3) + 1); hold(ax, 'on');
        q = P(P.genotype == S.order(k), :); c = S.col(k, :);
        plot(ax, q.J(~q.filled), q.EPR(~q.filled), 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c, 'LineWidth', 0.7);
        plot(ax, q.J(q.filled), q.EPR(q.filled), 'o', 'MarkerSize', 3, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
        cc = C(C.genotype == S.order(k), :); plot(ax, cc.J, cc.EPR_fit, 'k-', 'LineWidth', 0.9);
        title(ax, S.name(k), 'Color', c, 'FontWeight', 'normal', 'FontSize', 8);
        xlabel(ax, 'J (cycles h^{-1})'); if mod(k - 1, 3) == 0, ylabel(ax, 'EPR (nats h^{-1})'); end
        text(ax, 0.04, 0.92, sprintf('L = %.2f', O.L(k)), 'Units', 'normalized', 'FontSize', 7);
    end
    ax = nexttile(tl, 4, [2 1]); hold(ax, 'on');
    for k = 1:6
        y = 6 - k;
        plot(ax, [O.L_lo(k) O.L_hi(k)], [y y], '-', 'Color', S.col(k, :), 'LineWidth', 1.6);
        plot(ax, O.L(k), y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', S.col(k, :), 'MarkerEdgeColor', S.col(k, :));
    end
    xline(ax, O.L(1), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    set(ax, 'YTick', 0:5, 'YTickLabel', fliplr(S.name), 'YLim', [-0.6 5.6]);
    xlabel(ax, {'Onsager conductance L', '(cycles^2 nats^{-1} h^{-1}; EPR = J^2/L)'});
    text(ax, -0.25, 1.02, 'b', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
    text(nexttile(tl, 1), -0.35, 1.12, 'a', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
    exportgraphics(fh, fullfile(outDir, "figA2_epr_j_onsager_observable.png"), 'Resolution', 300);
    exportgraphics(fh, fullfile(outDir, "figA2_epr_j_onsager_observable.pdf"), 'ContentType', 'vector');
end
writetable(P, fullfile(outDir, "figA2_points.csv"));
writetable(C, fullfile(outDir, "figA2_curves.csv"));
writetable(O, fullfile(outDir, "figA2_L.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "figA2.mat"), '-struct', 'Xs');
end
