function X = exportFigA3(D, outDir, options)
%EXPORTFIGA3  Fig A3, energy landscape of the observed chain (FIGURE_RECIPES K3; fig_planA.py): data and figure.
%
%   X = exportFigA3(D, outDir)
%   X = exportFigA3(D, outDir, RandomSeed=20, Plot=true)
%
%   D       PLANADATA output
%   outDir  folder for figA3_points.csv, figA3_summary.csv, figA3.mat and the figure
%           (figA3_landscape_observable.png / .pdf)
%
%   Fly-phases with a Q/W reference and a chain with a cycle (ENERGYLANDSCAPE of the
%   fly's observed first-order chain):
%     (a) E_W - E_Q = ln(pi_Q / pi_W), (b) E_S - E_Q = ln(pi_Q / pi_S), as strips with
%         median and fly-bootstrap 95% CI
%     (c) landscape prediction F_SW (circles) and F_WS (squares) against counted
%         routing Xn, Yn on the observed symbols (XYCOUNTS: >= 3 routes per
%         denominator, else not drawn), identity line
%   X  struct: points, summary, fig
%   Check (recipe): E_S - E_Q medians -1.78 (CS dark), +1.00 (fumin dark).
%
%   See also ENERGYLANDSCAPE, ROUTECOUNTSOBSERVED, EXPORTFIGA1.

arguments
    D (1,1) struct
    outDir (1,1) string
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
    options.Plot (1,1) logical = true
end
if ~isfolder(outDir), mkdir(outDir); end
S = planAStyle(); B = S.blocks; nB = height(B);
assert(isequal(D.order, S.order), 'exportFigA3:order', 'D must hold the genotypes in the order %s.', strjoin(S.order));
P = table(); e21 = cell(nB, 1); e31 = cell(nB, 1);
for b = 1:nB
    T = D.G(B.g(b)).perFly; T = T(T.phase == B.phase(b) & T.refOK & T.hasCycle, :);
    e21{b} = T.E21; e31{b} = T.E31;
    P = [P; table(repmat(B.label(b), height(T), 1), T.genotype, T.file, T.fly, T.phaseName, ...
        repmat(B.filled(b), height(T), 1), T.E21, T.E31, T.F_SW, T.Xn_cnt, T.F_WS, T.Yn_cnt, ...
        'VariableNames', {'block', 'genotype', 'file', 'fly', 'phase', 'filled', 'E21', 'E31', ...
        'F_SW', 'Xn_cnt', 'F_WS', 'Yn_cnt'})]; %#ok<AGROW>
end
P.block = erase(erase(erase(P.block, "\it"), "\rm"), "^");
rs = RandStream('twister', 'Seed', options.RandomSeed);
X.fig = [];
if options.Plot
    fh = planAFigure("Fig A3", 7.0, 2.8); X.fig = fh;
    tl = tiledlayout(fh, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    a1 = nexttile(tl); a2 = nexttile(tl); a3 = nexttile(tl);
else
    fh = figure('Visible', 'off'); a1 = axes(fh); a2 = axes(fh); a3 = [];
end
T1 = planAStrip(a1, S, e21, Stream=rs, Gap=0.4, LabelStyle="short");
T2 = planAStrip(a2, S, e31, Stream=rs, Gap=0.4, LabelStyle="short");
X.summary = [[table(repmat("E21", height(T1), 1), 'VariableNames', {'quantity'}), T1]; ...
             [table(repmat("E31", height(T2), 1), 'VariableNames', {'quantity'}), T2]];
X.summary.label = erase(erase(erase(X.summary.label, "\it"), "\rm"), "^");
X.points = P;
if options.Plot
    for a = [a1 a2], a.XTickLabelRotation = 90; a.XAxis.FontSize = 6; end
    ylabel(a1, 'E_W - E_Q = ln(\pi_Q/\pi_W) (nats)'); ylabel(a2, 'E_S - E_Q = ln(\pi_Q/\pi_S) (nats)');
    hold(a3, 'on');
    mk = ["o" "s"]; pk = ["F_SW" "F_WS"]; ck = ["Xn_cnt" "Yn_cnt"];
    for b = 1:nB
        q = P(P.block == erase(erase(erase(B.label(b), "\it"), "\rm"), "^"), :); c = S.col(B.g(b), :);
        if B.filled(b), mf = c; else, mf = 'w'; end
        for j = 1:2
            ok = isfinite(q.(pk(j))) & isfinite(q.(ck(j)));
            plot(a3, q.(pk(j))(ok), q.(ck(j))(ok), mk(j), 'MarkerSize', 2.2, 'MarkerFaceColor', mf, ...
                'MarkerEdgeColor', c, 'LineWidth', 0.5);
        end
    end
    plot(a3, [-1 1], [-1 1], 'k-', 'LineWidth', 0.7);
    set(a3, 'XLim', [-1 1], 'YLim', [-1 1], 'DataAspectRatio', [1 1 1]);
    xlabel(a3, 'landscape prediction'); ylabel(a3, 'counted routing (observed)');
    text(a3, -0.95, -0.8, {'o  S\rightarrowW (X_n)', [char(9633) '  W\rightarrowS (Y_n)']}, 'FontSize', 6);
    lab = 'abc'; ax = [a1 a2 a3];
    for k = 1:3, text(ax(k), -0.12, 1.04, lab(k), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold'); end
    exportgraphics(fh, fullfile(outDir, "figA3_landscape_observable.png"), 'Resolution', 300);
    exportgraphics(fh, fullfile(outDir, "figA3_landscape_observable.pdf"), 'ContentType', 'vector');
else
    close(fh);
end
writetable(X.points, fullfile(outDir, "figA3_points.csv"));
writetable(X.summary, fullfile(outDir, "figA3_summary.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "figA3.mat"), '-struct', 'Xs');
end
