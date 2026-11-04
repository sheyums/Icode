function X = exportSleepVsCycleFigure(Tobs, Thid, outDir)
%EXPORTSLEEPVSCYCLEFIGURE  C14 / Act 3: cycle flux and direction against sleep, observed (top) and hidden (bottom) (sqw22h.py figure).
%
%   X = exportSleepVsCycleFigure(Tobs, Thid, outDir)
%
%   Tobs    SLEEPVSCYCLEROWS(D) (observed: sleep, c_obs, J, cycle)
%   Thid    SLEEPVSCYCLEROWS(D, Level="hidden", ...) (hidden: model sleep 60 pi_S, d, J_hmm, cycle)
%   outDir  folder for fig22_J_vs_sleep_obs_hidden.png / .pdf and fig22_points.csv
%
%   The seven PRESPEC 22 primary blocks (CS L, D; per0 DD; fumin L, D; sss L, D). (a) observed cycle
%   flux J (cycle fly-phases; no-cycle ones as crosses at 0), (b) observed direction c_obs, (c) hidden
%   J, (d) hidden net forward turns per sleep onset d; x = sleep (min per hour of the phase; hidden:
%   from the model). Open symbols light, filled dark / DD. Descriptive; the tests are SLEEPVSCYCLETESTS.
%   X  struct: points (table), fig
%
%   Example (from the Icode\sqw folder):
%       jsonl = fullfile("..", "hmm_audit", "sqw20b_fit2_theta0.5.jsonl");
%       D = buildPerFlyData("D:\flies", jsonl, Cache="perfly_cache.mat");
%       To = sleepVsCycleRows(D); Th = sleepVsCycleRows(D, Level="hidden", DataRoot="D:\flies", JsonlFile=jsonl);
%       X = exportSleepVsCycleFigure(To, Th, "figures");
%
%   See also SLEEPVSCYCLEROWS, SLEEPVSCYCLETESTS, NEWPAPERFIGURE.

arguments
    Tobs table
    Thid table
    outDir (1,1) string
end
if ~isfolder(outDir), mkdir(outDir); end
S = paperFigureStyle();
prim = {"CS", 1; "CS", 2; "per0DD", 1; "fumin", 1; "fumin", 2; "sss", 1; "sss", 2};
col = containers.Map(["CS" "per0DD" "fumin" "sss"], {[42 120 214] / 255, [237 161 0] / 255, [212 61 81] / 255, [106 61 154] / 255});
nm = containers.Map(["CS" "per0DD" "fumin" "sss"], ["CS" "\itper^0\rm DD" "\itfumin\rm" "\itsss\rm"]);
pan = {Tobs, "J", 'observed cycle flux J (cycles h^{-1})'; Tobs, "c_obs", 'observed direction per transition c_{obs}'; ...
       Thid, "J", 'hidden cycle flux J (cycles h^{-1})'; Thid, "d", {'hidden net forward turns', 'per sleep onset'}};
fh = newPaperFigure("C14 sleep vs cycle", 7.0, 5.8); X.fig = fh;
tl = tiledlayout(fh, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); pts = table();
for c = 1:4
    T = pan{c, 1}; v = pan{c, 2}; ax = nexttile(tl); hold(ax, 'on');
    for k = 1:size(prim, 1)
        g = prim{k, 1}; p = prim{k, 2}; U = T(T.genotype == g & T.phase == p, :); cc = col(g);
        if p == 2 || g == "per0DD", mf = cc; else, mf = 'w'; end
        if v == "J", V = U(U.cycle, :); else, V = U(isfinite(U.(v)), :); end
        plot(ax, V.sleep, V.(v), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', mf, 'MarkerEdgeColor', cc, 'LineWidth', 0.6);
        if v == "J", Z = U(~U.cycle, :); plot(ax, Z.sleep, zeros(height(Z), 1), 'x', 'MarkerSize', 4, 'Color', cc, 'LineWidth', 0.6); end
        pts = [pts; table(repmat(c, height(V), 1), repmat(g, height(V), 1), repmat(p, height(V), 1), V.fly, V.sleep, V.(v), ...
            'VariableNames', {'panel', 'genotype', 'phase', 'fly', 'sleep', 'value'})]; %#ok<AGROW>
    end
    yline(ax, 0, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5); ylabel(ax, pan{c, 3});
    if c > 2, xlabel(ax, 'sleep (min per hour of the phase)'); end
    text(ax, -0.16, 1.03, char('a' + c - 1), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
    if c == 1, title(ax, 'observed chain', 'FontSize', 8.5, 'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.01]); end
    if c == 3, title(ax, 'hidden model (post hoc; x = model sleep)', 'FontSize', 8.5, 'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.01]); end
end
ax = nexttile(tl, 1); h = gobjects(1, 7); ks = ["CS" "per0DD" "fumin" "sss"];
for k = 1:4, h(k) = plot(ax, NaN, NaN, 's', 'MarkerFaceColor', col(ks(k)), 'MarkerEdgeColor', col(ks(k))); end
h(5) = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k'); h(6) = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
h(7) = plot(ax, NaN, NaN, 'x', 'Color', 'k');
lg = legend(ax, h, [values(nm, cellstr(ks)) {'light', 'dark / DD', 'no cycle'}], 'Orientation', 'horizontal', 'NumColumns', 7, 'Box', 'off', 'FontSize', 7);
lg.Layout.Tile = 'north';
X.points = pts;
exportgraphics(fh, fullfile(outDir, "fig22_J_vs_sleep_obs_hidden.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "fig22_J_vs_sleep_obs_hidden.pdf"), 'ContentType', 'vector');
writetable(pts, fullfile(outDir, "fig22_points.csv"));
end
