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
%   (c) (2026-11-02) effective affinity per turn A_eff = EPR / |J| per fly x phase (the
%       points of (a) with J ~= 0), clipped at 3 for display, median with fly-bootstrap 95% CI;
%       top axis e^A, the odds of a forward over a backward turn. The sss panel of (a) notes
%       that its light (open) points are reversed.
%   X  struct: points (with A_eff), curves (J grid and fitted EPR per genotype), L (the
%      SQWONSAGERGROUP table, incl. the bootstrap difference against CS), Aeff, fig
%   Check (recipe): L = 0.272 (CS), 0.262 (CS DD), 0.286 (per0 LD), 0.352 (per0 DD),
%   0.243 (fumin), 0.145 (sss); A_eff medians 0.57, 0.51, 0.43, 0.64, 1.11, 0.70.
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
% (c) effective affinity per turn, A_eff = EPR / |J| (nats per cycle), the same points with J ~= 0
P.A_eff = P.EPR ./ abs(P.J); P.A_eff(P.J == 0) = NaN;
rsA = RandStream('twister', 'Seed', options.RandomSeed + 1); AE = table();
for k = 1:6
    a = P.A_eff(P.genotype == S.order(k) & isfinite(P.A_eff));
    [m, lo, hi] = medianCI(a, NBoot=options.NBoot, Stream=rsA);
    AE = [AE; table(S.order(k), numel(a), m, lo, hi, 'VariableNames', {'genotype', 'n', 'A_eff_med', 'lo', 'hi'})]; %#ok<AGROW>
end
X.points = P; X.curves = C; X.L = O; X.Aeff = AE; X.fig = [];
X.boot = O.Properties.UserData.boot;

if options.Plot
    fh = planAFigure("Fig A2", 7.2, 5.0); X.fig = fh;
    tl = tiledlayout(fh, 2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:6
        ax = nexttile(tl, (floor((k - 1) / 3)) * 5 + mod(k - 1, 3) + 1); hold(ax, 'on');
        q = P(P.genotype == S.order(k), :); c = S.col(k, :);
        plot(ax, q.J(~q.filled), q.EPR(~q.filled), 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c, 'LineWidth', 0.7);
        plot(ax, q.J(q.filled), q.EPR(q.filled), 'o', 'MarkerSize', 3, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
        cc = C(C.genotype == S.order(k), :); plot(ax, cc.J, cc.EPR_fit, 'k-', 'LineWidth', 0.9);
        if S.order(k) == "sss"   % 2026-11-02: sleepless runs the cycle backwards by day
            title(ax, {S.name(k), '(open = light: reversed)'}, 'Color', c, 'FontWeight', 'normal', 'FontSize', 7);
        else
            title(ax, S.name(k), 'Color', c, 'FontWeight', 'normal', 'FontSize', 8);
        end
        xlabel(ax, 'J (cycles h^{-1})'); if mod(k - 1, 3) == 0, ylabel(ax, 'EPR (nats h^{-1})'); end
        text(ax, 0.04, 0.92, sprintf('L = %.2f', O.L(k)), 'Units', 'normalized', 'FontSize', 7);
        ax.XAxis.TickValues = unique(round(linspace(ax.XLim(1), ax.XLim(2), 3), 1));
    end
    ax = nexttile(tl, 4, [2 1]); hold(ax, 'on');
    for k = 1:6
        y = 6 - k;
        plot(ax, [O.L_lo(k) O.L_hi(k)], [y y], '-', 'Color', S.col(k, :), 'LineWidth', 1.6);
        plot(ax, O.L(k), y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', S.col(k, :), 'MarkerEdgeColor', S.col(k, :));
    end
    xline(ax, O.L(1), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    set(ax, 'YTick', 0:5, 'YTickLabel', fliplr(S.name), 'YLim', [-0.6 6.75]);   % same rows as (c), which carries the e^A scale on top
    xlabel(ax, {'conductance L', '(cycles^2 nats^{-1} h^{-1})'});
    text(ax, -0.25, 1.02, 'b', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
    ax = nexttile(tl, 5, [2 1]); hold(ax, 'on');
    for k = 1:6
        y = 6 - k; a = P.A_eff(P.genotype == S.order(k) & isfinite(P.A_eff));
        plot(ax, min(max(a, 0), 3), y + rand(rsA, numel(a), 1) * 0.36 - 0.18, 'o', 'MarkerSize', 1.8, ...
            'MarkerFaceColor', 'none', 'MarkerEdgeColor', S.col(k, :), 'LineWidth', 0.4);
        plot(ax, [AE.lo(k) AE.hi(k)], [y y], 'k-', 'LineWidth', 1.4);
        plot(ax, AE.A_eff_med(k), y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', S.col(k, :), 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end
    set(ax, 'YTick', 0:5, 'YTickLabel', [], 'XLim', [0 3], 'YLim', [-0.6 6.75], 'Box', 'off');
    xlabel(ax, {'A_{eff} = EPR/|J|', '(nats per cycle)'});
    % top axis: odds forward : backward, e^A
    odds = [1 2 5 10 20]; yt = 5.75;   % the e^A scale drawn inside the axes, above the last row
    plot(ax, [0 3], [yt yt], 'k-', 'LineWidth', 0.6, 'Clipping', 'off');
    for o = odds
        plot(ax, log(o) * [1 1], yt - [0 0.12], 'k-', 'LineWidth', 0.6, 'Clipping', 'off');
        text(ax, log(o), yt + 0.08, num2str(o), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    text(ax, 1.5, yt + 0.45, 'odds forward : backward, e^{A}', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 6.5);
    text(ax, -0.1, 1.02, 'c', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
    text(nexttile(tl, 1), -0.35, 1.12, 'a', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
    exportgraphics(fh, fullfile(outDir, "figA2_epr_j_onsager_observable.png"), 'Resolution', 300);
    exportgraphics(fh, fullfile(outDir, "figA2_epr_j_onsager_observable.pdf"), 'ContentType', 'vector');
end
writetable(P, fullfile(outDir, "figA2_points.csv"));
writetable(C, fullfile(outDir, "figA2_curves.csv"));
writetable(O, fullfile(outDir, "figA2_L.csv"));
writetable(AE, fullfile(outDir, "figA2_Aeff.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "figA2.mat"), '-struct', 'Xs');
end
