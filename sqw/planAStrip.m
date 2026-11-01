function T = planAStrip(ax, S, v, options)
%PLANASTRIP  One strip panel of a Plan A figure: per-fly points, median and fly-bootstrap 95% CI per block.
%
%   T = planAStrip(ax, S, v)
%   T = planAStrip(ax, S, v, Hmm=h, Zero=true, Stream=rs, Gap=0.5, LabelStyle="two-line")
%
%   S     PLANASTYLE; one strip per row of S.blocks, light before dark
%   v     (nBlock x 1) cell, the per-fly values of each block
%   Hmm   (nBlock x 1) cell of HMM values per block (non-degenerate flies);
%         their median is drawn as a grey triangle beside the block, labelled
%         "HMM (model-based estimate), median" (not a bound; 2026-11-01)
%   Zero  draw a grey line at 0
%   Stream  RandStream for jitter and bootstrap (default twister seeded 20)
%   Gap   extra space after a dark / DD block
%
%   Points: light open, dark / DD filled, in the genotype colour. Median: black
%   bar; CI: black line (MEDIANCI: non-finite values dropped, >= 3 needed).
%   T  table per block: label, x, n (finite), median, lo, hi, nHmm, hmmMedian.
%
%   See also EXPORTFIGA1, EXPORTFIGA3.

arguments
    ax
    S (1,1) struct
    v cell
    options.Hmm cell = {}
    options.Zero (1,1) logical = true
    options.Stream = RandStream('twister', 'Seed', 20)
    options.Gap (1,1) double = 0.5
    options.LabelStyle (1,1) string {mustBeMember(options.LabelStyle, ["two-line" "short"])} = "two-line"
end
rs = options.Stream; B = S.blocks; nB = height(B);
x = 0; T = table();
hold(ax, 'on');
for b = 1:nB
    c = S.col(B.g(b), :); vals = v{b}(:);
    if B.filled(b), mf = c; else, mf = 'w'; end
    plot(ax, x + (rand(rs, numel(vals), 1) * 0.36 - 0.18), vals, 'o', 'MarkerSize', 2.6, ...
        'MarkerFaceColor', mf, 'MarkerEdgeColor', c, 'LineWidth', 0.6);
    [m, lo, hi] = medianCI(vals, NBoot=S.NBoot, Level=S.Level, Stream=rs);
    plot(ax, [x - 0.3 x + 0.3], [m m], 'k-', 'LineWidth', 1.2);
    plot(ax, [x x], [lo hi], 'k-', 'LineWidth', 0.8);
    r = table(B.label(b), x, nnz(isfinite(vals)), m, lo, hi, 0, NaN, ...
        'VariableNames', {'label', 'x', 'n', 'median', 'lo', 'hi', 'nHmm', 'hmmMedian'});
    if ~isempty(options.Hmm)
        h = options.Hmm{b}(:); h = h(isfinite(h));
        if ~isempty(h)
            r.nHmm = numel(h); r.hmmMedian = median(h);
            plot(ax, x + 0.36, r.hmmMedian, '<', 'MarkerSize', 4, 'MarkerFaceColor', [0.35 0.35 0.35], ...
                'MarkerEdgeColor', [0.35 0.35 0.35]);
        end
    end
    T = [T; r]; %#ok<AGROW>
    x = x + 1 + options.Gap * B.filled(b);
end
if options.Zero, hl = yline(ax, 0, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5); uistack(hl, 'bottom'); end
lab = B.label;
if options.LabelStyle == "two-line"
    lab = replace(lab, [" light" " dark"], ["\newline light" "\newline dark"]);   % a newline character would split one label into two ticks
else
    lab = replace(lab, [" light" " dark"], [" L" " D"]);
end
set(ax, 'XTick', T.x, 'XTickLabel', lab, 'XLim', [T.x(1) - 0.6, T.x(end) + 0.6], 'TickLabelInterpreter', 'tex');
ax.XAxis.FontSize = 6.5;
end
