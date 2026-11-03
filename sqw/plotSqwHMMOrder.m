function fh = plotSqwHMMOrder(R, names)
%PLOTSQWHMMORDER  HMM order figure for one or more SQWBATCHHMM batches (fig_hmm_order.py, fig_emission.py).
%
%   fh = plotSqwHMMOrder(R)                       % one batch
%   fh = plotSqwHMMOrder({R1, R2, ...}, ["CS" "\itsss\rm" ...])
%
%   R      one SQWBATCHHMM result (struct) or a cell array of them, one per group
%   names  (optional) one label per batch, TeX allowed; default each batch's R.label (its folder name)
%   fh     the figure handle (save it with EXPORTGRAPHICS)
%
%   Top: per batch, the median over flies with an identified H1 (deg3 false) of the labelled
%   K = 3 emission matrix P(symbol | state). Bottom: fraction of flies whose BIC is lowest at
%   each K, with n = minutes (primary) and n = symbol changes (sensitivity), as fig_hmm_order.py;
%   and, for K = 4, the fraction of flies with 1, 2 or 3 hidden states emitting mostly S
%   (two would be the "deep sleep" reading of the extra state). Theme pinned light.
%
%   Example (from the Icode\sqw folder):
%       R1 = sqwBatchHMM("D:\flies\CS_symbols"); R2 = sqwBatchHMM("D:\flies\per0DD_symbols");
%       fh = plotSqwHMMOrder({R1, R2}, ["CS" "\itper^0\rm DD"]);
%       exportgraphics(fh, "hmm_order.png", Resolution=300);
%
%   See also SQWBATCHHMM.

if isstruct(R), R = {R}; end
if nargin < 2, names = string(cellfun(@(r) r.label, R, 'UniformOutput', false)); end
names = string(names); G = numel(R); S = paperFigureStyle();
KC = [201 214 232; 31 78 121; 143 182 224] / 255; KS = [2 3 4];
SC = [0.93 0.80 0.86; 0.80 0.37 0.55; 0.45 0.10 0.27];
one = G == 1;
if one, nc = 4; fh = newPaperFigure("HMM order", 7.0, 1.9); tl = tiledlayout(fh, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
else, nc = 3 * (G > 3) + 3; fh = newPaperFigure("HMM order", 7.0, 4.2); tl = tiledlayout(fh, 2, nc, 'TileSpacing', 'compact', 'Padding', 'compact'); end
for g = 1:G
    T = R{g}.perFly; ax = nexttile(tl, g);
    if ~ismember('deg3', T.Properties.VariableNames), axis(ax, 'off'); continue, end
    U = T(~T.deg3, :); M = zeros(3); sym = ["Q" "W" "S"];
    for h = 1:3, for v = 1:3, M(h, v) = median(U.("B_H" + h + "_" + sym(v))); end, end
    image(ax, repmat(1 - M, 1, 1, 3)); axis(ax, 'image'); hold(ax, 'on');
    for h = 1:3
        for v = 1:3
            if M(h, v) > 0.55, tc = 'w'; else, tc = 'k'; end
            text(ax, v, h, sprintf('%.2f', M(h, v)), 'HorizontalAlignment', 'center', 'FontSize', 6.5, 'Color', tc);
        end
    end
    set(ax, 'XTick', 1:3, 'XTickLabel', {'Q', 'W', 'S'}, 'YTick', 1:3, 'TickLength', [0 0]);
    if g == 1, set(ax, 'YTickLabel', {'H_1', 'H_2', 'H_3'}); ylabel(ax, 'K = 3 state'); else, set(ax, 'YTickLabel', []); end
    title(ax, sprintf('%s (n %d)', names(g), height(U)), 'FontWeight', 'normal', 'FontSize', 7.5);
end
if one, w = 1; first = 2; else, w = nc / 3; first = nc + 1; end
lab = names + " (" + string(cellfun(@(r) height(r.perFly), R)) + ")";
cols = {"bestK_T", "n = minutes (primary)"; "bestK_chg", "n = symbol changes"};
for c = 1:3
    ax = nexttile(tl, first + (c - 1) * w, [1 w]); hold(ax, 'on');
    for g = 1:G
        T = R{g}.perFly;
        if c < 3, v = T.(cols{c, 1}); levels = KS; col = KC;
        elseif ismember('nS_K4', T.Properties.VariableNames), v = T.nS_K4; levels = 1:3; col = SC;
        else, continue
        end
        left = 0;
        for k = 1:numel(levels)
            f = mean(v == levels(k));
            patch(ax, left + [0 f f 0], g + [-0.4 -0.4 0.4 0.4], col(k, :), 'EdgeColor', 'w', 'LineWidth', 1.2);
            if f >= 0.08
                if (c < 3 && levels(k) == 3) || (c == 3 && k == 3), tc = 'w'; else, tc = 'k'; end
                text(ax, left + f / 2, g, sprintf('%.0f%%', 100 * f), 'HorizontalAlignment', 'center', 'FontSize', 6.5, 'Color', tc);
            end
            left = left + f;
        end
    end
    set(ax, 'YDir', 'reverse', 'YTick', 1:G, 'XLim', [0 1], 'YLim', [0.4 G + 0.6], 'TickLabelInterpreter', 'tex');
    if c == 1, set(ax, 'YTickLabel', lab); else, set(ax, 'YTickLabel', []); end
    if c < 3
        title(ax, cols{c, 2}, 'FontWeight', 'normal', 'FontSize', 7.5); xlabel(ax, 'fraction of flies, lowest BIC');
    else
        title(ax, 'K = 4: states emitting mostly S', 'FontWeight', 'normal', 'FontSize', 7.5); xlabel(ax, 'fraction of flies');
    end
end
ax = nexttile(tl, first); h = gobjects(1, 6);
for k = 1:3, h(k) = patch(ax, NaN, NaN, KC(k, :), 'EdgeColor', 'none'); h(k + 3) = patch(ax, NaN, NaN, SC(k, :), 'EdgeColor', 'none'); end
lg = legend(ax, h, {'K = 2', 'K = 3', 'K = 4', '1 S-state', '2 S-states', '3 S-states'}, 'Box', 'off', 'FontSize', 6.5, ...
    'Orientation', 'horizontal', 'NumColumns', 6);
lg.Layout.Tile = 'south';
end
