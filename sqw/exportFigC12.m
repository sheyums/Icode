function X = exportFigC12(jsonFile, outDir, options)
%EXPORTFIGC12  Fig C12, circuit silencing with Kir2.1 (FIGURE_RECIPES section L; fig_kir.py; PRESPEC 20y).
%
%   X = exportFigC12(jsonFile, outDir)
%   X = exportFigC12(jsonFile, outDir, RandomSeed=1)
%
%   jsonFile  hmm_audit/y_kir.json (one record per fly: g = empty / dfb / vnc, dead, and per phase
%             (light, dark) in 'main': sleep, P (QS, SW, ...), E31, J, ok (Q/W reference), cycle)
%   outDir    folder for figC12_points.csv, figC12_summary.csv, figC12.mat and the figure
%             (figC12_kir_silencing.png / .pdf)
%
%   Five panels, each light | dark, groups control (empty>Kir, #8a8a8a), dFB (#2a78d6), VNC (#eb6834):
%     sleep (min/h); P(Q->S) per min; P(S->W) per min; depth of the sleep well -(E_S - E_Q) = -E31
%     (nats); cycle flux J (count, cycles/h)
%   Dead flies dropped; fly x phases without a Q/W reference excluded (sleep excepted; none here);
%   E31 and J only for fly x phases with a cycle. Per-fly points (open light, filled dark), median with
%   fly-bootstrap 95% CI (2000; MEDIANCI). '*' = Holm < 0.05 in the 20y ledger entry, flags copied
%   (as fig_kir.py), not recomputed. Depth panel ylim -0.6..3.3; values above are counted and named
%   on the panel rather than hard-coded.
%
%   See also PLANASTYLE, PLANAFIGURE, MEDIANCI.

arguments
    jsonFile (1,1) string {mustBeFile}
    outDir (1,1) string
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 1
end
if ~isfolder(outDir), mkdir(outDir); end
Y = jsondecode(fileread(jsonFile)); if iscell(Y), Y = [Y{:}]'; end
Y = Y(~[Y.dead]);
G = ["empty" "dfb" "vnc"]; col = [138 138 138; 42 120 214; 235 104 52] / 255;
lab = ["control (empty>Kir)" "dFB silenced" "VNC silenced"]; ph = ["light" "dark"];
% key, y label, title, ledger Holm < 0.05 as [group phase] (group 2 = dfb, 3 = vnc; phase 1 light, 2 dark)
P = {"sleep", 'sleep (min per hour)', 'Sleep amount', [3 1; 2 2; 3 2];
     "QS", 'P(quiet wake \rightarrow sleep) per min', 'Falling asleep', [3 1; 2 2; 3 2];
     "SW", 'P(sleep \rightarrow active wake) per min', 'Waking up', [2 1; 3 1; 2 2; 3 2];
     "depth", {'depth of the sleep well,', 'E_Q - E_S (nats)'}, {'Sleep vs quiet wake', '(occupancy)'}, [2 1; 3 1; 2 2; 3 2];
     "J", 'cycle flux J (cycles h^{-1})', {'Drive around', ['W' char(8594) 'Q' char(8594) 'S' char(8594) 'W']}, zeros(0, 2)};

pts = table();
for i = 1:numel(Y)
    m = Y(i).main; if iscell(m), m = [m{:}]'; end
    for p = 1:2
        q = m(p);
        pts = [pts; table(string(Y(i).g), string(Y(i).f), ph(p), logical(q.ok), logical(q.cycle), q.sleep, q.P.QS, q.P.SW, ...
            -num(q.E31), num(q.J), 'VariableNames', {'group', 'file', 'phase', 'ok', 'cycle', 'sleep', 'QS', 'SW', 'depth', 'J'})]; %#ok<AGROW>
    end
end
rs = RandStream('twister', 'Seed', options.RandomSeed);
fh = planAFigure("Fig C12", 7.2, 2.9); X.fig = fh;
tl = tiledlayout(fh, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
summ = table(); S = planAStyle(); hl = gobjects(3, 1);
for a = 1:5
    key = P{a, 1}; ax = nexttile(tl); hold(ax, 'on'); x0 = 0; ticks = zeros(1, 2);
    for p = 1:2
        for g = 1:3
            x = x0 + (g - 1) * 0.9;
            r = pts.group == G(g) & pts.phase == ph(p) & (key == "sleep" | pts.ok);
            if ismember(key, ["J" "depth"]), r = r & pts.cycle; end
            v = pts.(key)(r);
            if p == 2, mf = col(g, :); else, mf = 'w'; end
            hh = plot(ax, x + rand(rs, numel(v), 1) * 0.4 - 0.2, v, 'o', 'MarkerSize', 2.2, 'MarkerFaceColor', mf, ...
                'MarkerEdgeColor', col(g, :), 'LineWidth', 0.5);
            if a == 1 && p == 2, hl(g) = hh; end
            [md, lo, hi] = medianCI(v, NBoot=S.NBoot, Stream=rs);
            plot(ax, [x - 0.3 x + 0.3], [md md], 'k-', 'LineWidth', 1.2); plot(ax, [x x], [lo hi], 'k-', 'LineWidth', 0.8);
            star = any(all(P{a, 4} == [g p], 2));
            if star, text(ax, x, 1, '*', 'Units', 'data', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                    'FontSize', 11, 'Tag', 'star'); end
            summ = [summ; table(key, G(g), ph(p), numel(v), md, lo, hi, star, ...
                'VariableNames', {'quantity', 'group', 'phase', 'n', 'median', 'lo', 'hi', 'ledgerHolmStar'})]; %#ok<AGROW>
        end
        ticks(p) = x0 + 0.9; x0 = x0 + 3.3;
    end
    set(ax, 'XTick', ticks, 'XTickLabel', ph, 'XLim', [-0.6 x0 - 0.9]);
    ylabel(ax, P{a, 2}); title(ax, P{a, 3}, 'FontWeight', 'normal', 'FontSize', 8.5);
    if key == "J", yline(ax, 0, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5); end
    if key == "depth"
        ylim(ax, [-0.6 3.3]);
        out = pts(pts.ok & pts.cycle & pts.depth > 3.3, :);
        if ~isempty(out)
            txt = strjoin(compose("%s %s at %.1f", out.group, out.phase, out.depth), ", ");
            text(ax, 0.02, 0.01, {'not shown (above 3.3):', char(replace(txt, "empty", "control"))}, 'Units', 'normalized', ...
                'VerticalAlignment', 'bottom', 'FontSize', 5.5);
        end
    end
    yl = ylim(ax);   % stars at the top of each panel
    set(findobj(ax, 'Tag', 'star'), 'Units', 'data');
    for s = findobj(ax, 'Tag', 'star')', s.Position(2) = yl(2); end
end
lg = legend(hl, lab, 'NumColumns', 3, 'Box', 'off', 'FontSize', 7.5); lg.Layout.Tile = 'north';
xlabel(tl, 'Each point is one fly (open: light phase, filled: dark phase); bar: median with 95% CI; * differs from control.', ...
    'FontSize', 6.5, 'FontName', S.font);
X.points = pts; X.summary = summ;
exportgraphics(fh, fullfile(outDir, "figC12_kir_silencing.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "figC12_kir_silencing.pdf"), 'ContentType', 'vector');
writetable(pts, fullfile(outDir, "figC12_points.csv"));
writetable(summ, fullfile(outDir, "figC12_summary.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "figC12.mat"), '-struct', 'Xs');
end

function v = num(v)
% JSON null (a fly-phase without a cycle has no landscape) decodes to []
if isempty(v), v = NaN; end
end
