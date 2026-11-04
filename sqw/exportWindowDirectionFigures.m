function X = exportWindowDirectionFigures(R, dataRoot, outDir)
%EXPORTWINDOWDIRECTIONFIGURES  PRESPEC 25 SI figures: window rasters, counts against the null, and the cumulative walk (fig25.py, fig25_walk.py).
%
%   X = exportWindowDirectionFigures(R, dataRoot, outDir)
%
%   R         WINDOWDIRECTIONTESTS output
%   dataRoot  folder with the count-file subfolders CS, CSDD, per0LD, per0DD, sss (for the walks)
%   outDir    folder for the figures and fig25_walk_stats.csv
%
%   fig25_raster_6h / _4h   one row per fly (sorted by net direction: forward minus backward windows),
%       one column per window in time order: blue = measurably forward, orange = measurably backward,
%       light grey = not measurable, white = excluded (missing minutes); bar above: dark (filled) and
%       light (open) windows; columns grouped by day. Under each panel: forward and backward counts and
%       the clock-matched null mean of backward windows.
%   fig25_counts   (a) measurably backward and (b) forward windows, (c) backward share of measurable
%       windows, per genotype and window size: dot = observed (filled 6 h, open 4 h), grey bar = 95% of
%       the clock-matched null, tick = null mean (median for the share).
%   fig25_cumulative_walk   (a) each fly's CYCLEWALK (flies with >= 60% scored minutes), mean in black,
%       dark phase shaded (LD); (b) histogram of the net forward count (N_SW - N_WS) per 6-h window
%       (excluded windows left out), clipped at +-6, with the share of windows above and below 0.
%   X   struct: walkStats (table: genotype, n, hours, end_mean, frac_end_pos, frac_end_neg), figs
%
%   Example (from the Icode\sqw folder):
%       R = windowDirectionTests("D:\flies", Cache="windows.mat");
%       X = exportWindowDirectionFigures(R, "D:\flies", "figures");
%
%   See also WINDOWDIRECTIONTESTS, CYCLEWALK.

arguments
    R (1,1) struct
    dataRoot (1,1) string
    outDir (1,1) string
end
if ~isfolder(outDir), mkdir(outDir); end
S = paperFigureStyle(); B = R.blocks;
G = ["CS" "CSDD" "per0LD" "per0DD" "sss"]; ld = ["CS" "per0LD" "sss"];
NAME = ["CS, light–dark" "CS, constant dark" "\itper^0\rm, light–dark" "\itper^0\rm, constant dark" "\itsleepless\rm, light–dark"];
SH = ["CS LD" "CS DD" "\itper^0\rm LD" "\itper^0\rm DD" "\itsss\rm LD"];
fwd = [42 120 214] / 255; bwd = [235 104 52] / 255; none = [236 236 236] / 255;
X.figs = gobjects(1, 4);
% ---- rasters
for L = [360 240]
    sz = L / 60 + "h"; perDay = 1440 / L;
    nWin = arrayfun(@(g) max(cellfun(@numel, R.flies{B.key == g + "|" + sz}.cls)), G);
    fh = newPaperFigure("fig25 raster " + sz, 7.2, 3.6); tl = tiledlayout(fh, 1, sum(nWin), 'TileSpacing', 'compact', 'Padding', 'compact');
    col0 = 1;
    for k = 1:numel(G)
        g = G(k); b = find(B.key == g + "|" + sz); fl = R.flies{b}; nw = nWin(k);
        M = 9 * ones(height(fl), nw); for i = 1:height(fl), c = fl.cls{i}; M(i, 1:numel(c)) = c; end
        net = sum(M == 1, 2) - sum(M == -1, 2); nb = sum(M == -1, 2);
        [~, o] = sortrows([net -nb], [-1 -2]); M = M(o, :);
        img = zeros([size(M) 3]); cmap = {bwd, none, fwd}; vals = [-1 0 1];
        for v = 1:3, for ch = 1:3, layer = img(:, :, ch); layer(M == vals(v)) = cmap{v}(ch); img(:, :, ch) = layer; end, end
        for ch = 1:3, layer = img(:, :, ch); layer(M == 9) = 1; img(:, :, ch) = layer; end
        ax = nexttile(tl, col0, [1 nw]); col0 = col0 + nw; hold(ax, 'on');
        image(ax, [0.5 nw - 0.5], [0.5 size(M, 1) - 0.5], img);
        for d = 1:floor(nw / perDay) - 1, xline(ax, d * perDay, 'Color', 'w', 'LineWidth', 1.2); end
        for w = 0:nw - 1
            dark = ~ismember(g, ld) || mod(w, perDay) >= perDay / 2;
            if dark, fc = [0.25 0.25 0.25]; else, fc = 'w'; end
            rectangle(ax, 'Position', [w -1.6 1 1], 'FaceColor', fc, 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.4, 'Clipping', 'off');
        end
        set(ax, 'YDir', 'reverse', 'YLim', [-1.8 size(M, 1)], 'XLim', [0 nw], 'YTick', [], 'XTick', (0:floor(nw / perDay) - 1) * perDay + perDay / 2, ...
            'XTickLabel', "d" + (1:floor(nw / perDay)), 'TickLength', [0 0], 'XColor', 'k', 'Box', 'off'); ax.XAxis.FontSize = 6.5; ax.YAxis.Visible = 'off';
        title(ax, NAME(k), 'FontSize', 7.5, 'FontWeight', 'normal');
        xlabel(ax, {sprintf('fwd %d | back %d', B.F(b), B.B(b)), sprintf('null back \\approx %.0f', B.B_null_mean_clock(b))}, 'FontSize', 6.5);
        if k == 1, ax.YAxis.Visible = 'on'; ax.YColor = 'none'; ylabel(ax, 'flies (sorted by net direction)', 'Color', 'k'); end
    end
    ax = nexttile(tl, 1); h = [patch(ax, NaN, NaN, fwd) patch(ax, NaN, NaN, bwd) patch(ax, NaN, NaN, none)];
    lg = legend(ax, h, {'measurably forward (W\rightarrowQ\rightarrowS\rightarrowW)', 'measurably backward', 'no measurable cycle'}, ...
        'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 7); lg.Layout.Tile = 'north';
    xlabel(tl, sprintf(['Each cell is one %s window on one day. Bar above each panel: dark windows (filled) and light windows ' ...
        '(open). Columns are grouped by day (d1, d2, ...).'], sz), 'FontSize', 6.5, 'FontName', S.font);
    save2(fh, outDir, "fig25_raster_" + sz); X.figs(1 + (L == 240)) = fh;
end
% ---- counts
fh = newPaperFigure("fig25 counts", 7.2, 2.9); tl = tiledlayout(fh, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = [repelem(G', 2) repmat(["6h"; "4h"], numel(G), 1)]; ys = (size(rows, 1):-1:1)';
spec = {"B", 'measurably backward windows'; "F", 'measurably forward windows'; "frac", 'backward share of measurable windows'};
for c = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    for r = 1:size(rows, 1)
        q = B(B.key == rows(r, 1) + "|" + rows(r, 2), :); y = ys(r);
        if spec{c, 1} == "frac", lo = q.fracB_null_lo_clock; hi = q.fracB_null_hi_clock; mid = q.fracB_null_med_clock; obs = q.frac_B;
        else, w = spec{c, 1}; lo = q.(w + "_null_lo_clock"); hi = q.(w + "_null_hi_clock"); mid = q.(w + "_null_mean_clock"); obs = q.(w); end
        if isfinite(lo), plot(ax, [lo hi], [y y], '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 6); end
        if isfinite(mid), plot(ax, [mid mid], y + [-0.3 0.3], '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1); end
        if ~isfinite(obs)
            x0 = 0; if c == 3, x0 = 0.5; end
            text(ax, x0, y, '  none measurable', 'FontSize', 6, 'Color', [0.4 0.4 0.4]); continue
        end
        if rows(r, 2) == "6h", mf = 'k'; else, mf = 'w'; end
        plot(ax, obs, y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', mf, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
    xlabel(ax, spec{c, 2}); set(ax, 'YLim', [0.4 numel(ys) + 0.6], 'YTick', flipud(ys), 'TickLength', [0 0]);
    if c == 1, set(ax, 'YTickLabel', flipud(SH(arrayfun(@(g) find(G == g), rows(:, 1)))' + ", " + rows(:, 2)), 'TickLabelInterpreter', 'tex');
    else, set(ax, 'YTickLabel', []); end
    if c == 3, xline(ax, 0.5, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5); xlim(ax, [0 1]); end
    text(ax, 0, 1.06, char('a' + c - 1), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
end
xlabel(tl, {['Dots: observed (filled 6-h windows, open 4-h). Grey bar: 95% range from 500 simulated records per fly of a ' ...
    'constant process with each fly''s own'], 'clock-time-matched transition rates; tick: its mean. Dotted line in (c): no preferred direction.'}, ...
    'FontSize', 6.5, 'FontName', S.font);
save2(fh, outDir, "fig25_counts"); X.figs(3) = fh;
% ---- cumulative walk + per-window net counts
fh = newPaperFigure("fig25 walk", 7.2, 5.4); tl = tiledlayout(fh, 5, numel(G), 'TileSpacing', 'compact', 'Padding', 'compact');
st = table();
for k = 1:numel(G)
    g = G(k); files = sqwListFiles(fullfile(dataRoot, g)); cc = S.col(S.order == g, :); W = {};
    for i = 1:numel(files)
        y = sqwFlyInput(files(i), ismember(g, ld)); if mean(y >= 1 & y <= 3) < 0.6, continue, end
        W{end + 1} = cycleWalk(y); %#ok<AGROW>
    end
    T = min(cellfun(@numel, W)); M = cell2mat(cellfun(@(w) w(1:T)', W', 'UniformOutput', false)); t = (0:T - 1) / 60;
    ax = nexttile(tl, k, [3 1]); hold(ax, 'on'); yl = [min(cellfun(@min, W)) max(cellfun(@max, W))];
    if ismember(g, ld)
        for d = 0:floor(max(cellfun(@numel, W)) / 1440), patch(ax, 24 * d + [12 24 24 12], yl([1 1 2 2]), [0.9 0.9 0.9], 'EdgeColor', 'none'); end
    end
    for i = 1:numel(W), plot(ax, (0:numel(W{i}) - 1) / 60, W{i}, '-', 'Color', [cc 0.35], 'LineWidth', 0.5); end
    plot(ax, t, mean(M, 1), 'k-', 'LineWidth', 1.6); yline(ax, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    set(ax, 'XLim', [0 t(end)], 'XTick', 0:24:t(end)); xlabel(ax, 'hours'); title(ax, NAME(k), 'FontSize', 7.5, 'FontWeight', 'normal');
    if k == 1, ylabel(ax, {'cumulative net forward cycles', '(one line per fly; black: mean)'}); end
    e = M(:, end); st = [st; table(g, size(M, 1), t(end), mean(e), mean(e > 0), mean(e < 0), 'VariableNames', ...
        {'genotype', 'n', 'hours', 'end_mean', 'frac_end_pos', 'frac_end_neg'})]; %#ok<AGROW>
    fl = R.flies{B.key == g + "|6h"}; net = [];
    for i = 1:height(fl), c = fl.cls{i}; d = fl.nSW{i} - fl.nWS{i}; net = [net; d(c ~= 9)']; end %#ok<AGROW>
    net = min(max(net, -6), 6); xv = -6:6; hh = arrayfun(@(v) mean(net == v), xv);
    ax2 = nexttile(tl, 3 * numel(G) + k, [2 1]); hold(ax2, 'on');
    for v = 1:numel(xv)
        if xv(v) < 0, fc = bwd; elseif xv(v) > 0, fc = fwd; else, fc = [0.75 0.75 0.75]; end
        bar(ax2, xv(v), hh(v), 0.85, 'FaceColor', fc, 'EdgeColor', 'none');
    end
    set(ax2, 'XTick', [-6 -3 0 3 6], 'XTickLabel', {[char(8804) char(8722) '6'], [char(8722) '3'], '0', '3', [char(8805) '6']}, 'XLim', [-7 7], ...
        'TickLabelInterpreter', 'none', 'XTickLabelRotation', 0);
    xlabel(ax2, {'net forward cycles', 'in one 6-h window'}); if k == 1, ylabel(ax2, 'fraction of windows'); end
    text(ax2, 0.97, 0.95, {sprintf('fwd %.0f%%', 100 * mean(net > 0)), sprintf('back %.0f%%', 100 * mean(net < 0))}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 6.5);
    if k == 1, text(ax, -0.5, 1.05, 'a', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold'); text(ax2, -0.5, 1.05, 'b', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold'); end
end
xlabel(tl, ['Grey shading: dark phase. Each fly takes small forward and backward steps; the forward bias appears only as the ' ...
    'steps accumulate over days.'], 'FontSize', 6.5, 'FontName', S.font);
save2(fh, outDir, "fig25_cumulative_walk"); X.figs(4) = fh;
X.walkStats = st; writetable(st, fullfile(outDir, "fig25_walk_stats.csv"));
end

function save2(fh, outDir, name)
exportgraphics(fh, fullfile(outDir, name + ".png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, name + ".pdf"), 'ContentType', 'vector');
end
