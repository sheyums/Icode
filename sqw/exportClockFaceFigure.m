function X = exportClockFaceFigure(jsonFile, outDir)
%EXPORTCLOCKFACEFIGURE  Fig C10, the clock drives the cycle (FIGURE_RECIPES section L; fig_clockface.py; PRESPEC 20w-B).
%
%   X = exportClockFaceFigure(jsonFile, outDir)
%
%   jsonFile  hmm_audit/planB_20wB.json (field 'real': per genotype J (12 two-hour clock bins, pooled
%             over flies), AJ, phJ, pJ, nullJ95)
%   outDir    folder for figC10_bins.csv, figC10_amplitude.csv, figC10.mat and the figure
%             (figC10_clock_face.png / .pdf)
%
%   (a) 2 x 2: CS LD, CS DD, per0 LD, per0 DD. Points: pooled J per 2-h bin at bin centres 1, 3, ..,
%       23 h; line: least-squares 24-h cosine on the 12 bins (design [1 cos(2 pi b/12) sin(2 pi b/12)],
%       b = 0..11, drawn against t - 1); grey band: cosine mean +/- nullJ95 (the 95th-percentile
%       amplitude a random clock phase per fly gives); shading: dark phase (LD) or the whole panel
%       (DD); ylim -0.02..0.30; the verdict word of 20w-B per panel.
%   (b) per genotype: grey bar 0..nullJ95, dot = AJ, text p.
%   The verdict words are 20w-B's (Holm over CS, CS DD, per0 LD: 0.015, 0.03, 0.059; per0 DD is the
%   control), copied as fig_clockface.py does, not recomputed. X.check compares the port's cosine
%   amplitude with the stored AJ.
%
%   Example (from the Icode\sqw folder):
%       X = exportClockFaceFigure(fullfile("..", "hmm_audit", "planB_20wB.json"), "figures");
%
%   See also PAPERFIGURESTYLE, NEWPAPERFIGURE.

arguments
    jsonFile (1,1) string {mustBeFile}
    outDir (1,1) string
end
if ~isfolder(outDir), mkdir(outDir); end
R = jsondecode(fileread(jsonFile)); R = R.real;
S = paperFigureStyle();
G = ["CS" "CSDD" "per0LD" "per0DD"]; col = S.col(arrayfun(@(g) find(S.order == g), G), :);
nameP = ["CS, light" + char(8211) + "dark", "CS, constant dark", "\itper^0\rm, light" + char(8211) + "dark", ...
    "\itper^0\rm, constant dark"];
nameS = ["CS LD" "CS DD" "\itper^0\rm LD" "\itper^0\rm DD"];
verd = ["clock-locked rhythm" "clock-locked rhythm" "not detected" "none (control)"];
ld = [true false true false];
b = (0:11)'; t = 2 * b + 1; Xd = [ones(12, 1) cos(2 * pi * b / 12) sin(2 * pi * b / 12)];
bins = table(); amp = table();
for k = 1:4
    J = R.(G(k)).J(:); beta = Xd \ J;
    A = hypot(beta(2), beta(3));
    bins = [bins; table(repmat(G(k), 12, 1), t, J, Xd * beta, 'VariableNames', {'genotype', 'binCentre_h', 'J', 'cosineFit'})]; %#ok<AGROW>
    amp = [amp; table(G(k), R.(G(k)).n, beta(1), A, R.(G(k)).AJ, R.(G(k)).phJ, R.(G(k)).pJ, R.(G(k)).nullJ95, verd(k), ...
        'VariableNames', {'genotype', 'n', 'mean', 'A_port', 'AJ', 'phJ', 'pJ', 'nullJ95', 'verdict'})]; %#ok<AGROW>
end
X.bins = bins; X.amplitude = amp; X.check = max(abs(amp.A_port - amp.AJ));

fh = newPaperFigure("Fig C10", 7.0, 4.6); X.fig = fh;
tl = tiledlayout(fh, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
tt = linspace(0, 24, 200)';
for k = 1:4
    ax = nexttile(tl, (floor((k - 1) / 2)) * 3 + mod(k - 1, 2) + 1); hold(ax, 'on');
    r = amp(k, :); B = bins(bins.genotype == G(k), :);
    beta = Xd \ B.J;
    if ld(k), patch(ax, [12 24 24 12], [-0.02 -0.02 0.30 0.30], [0.88 0.88 0.88], 'EdgeColor', 'none');
    else, patch(ax, [0 24 24 0], [-0.02 -0.02 0.30 0.30], [0.94 0.94 0.94], 'EdgeColor', 'none'); end
    patch(ax, [0 24 24 0], r.mean + r.nullJ95 * [-1 -1 1 1], [0.75 0.75 0.75], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
    plot(ax, tt, beta(1) + beta(2) * cos(2 * pi * (tt - 1) / 24) + beta(3) * sin(2 * pi * (tt - 1) / 24), '-', ...
        'Color', col(k, :), 'LineWidth', 2);
    plot(ax, B.binCentre_h, B.J, 'o', 'MarkerSize', 4, 'MarkerFaceColor', col(k, :), 'MarkerEdgeColor', 'w', 'LineWidth', 0.6);
    set(ax, 'XLim', [0 24], 'XTick', 0:6:24, 'YLim', [-0.02 0.30], 'Layer', 'top');
    if ld(k), xlabel(ax, 'clock time (h) (ZT)'); else, xlabel(ax, 'clock time (h) from start of each day'); end
    if mod(k, 2) == 1, ylabel(ax, {'cycle flux J (cycles h^{-1})', 'W\rightarrowQ\rightarrowS\rightarrowW, all flies pooled'}); end
    title(ax, nameP(k), 'FontWeight', 'normal', 'FontSize', 8.5, 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02]);
    top = k == 1 || k == 3;
    if top, text(ax, 0.03, 0.97, verd(k), 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontAngle', 'italic', 'FontSize', 7);
    else, text(ax, 0.98, 0.04, verd(k), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'FontAngle', 'italic', 'FontSize', 7); end
end
ax = nexttile(tl, 3, [2 1]); hold(ax, 'on');
for k = 1:4
    y = 4 - k; r = amp(k, :);
    plot(ax, [0 r.nullJ95], [y y], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 6);
    plot(ax, r.AJ, y, 'o', 'MarkerSize', 7, 'MarkerFaceColor', col(k, :), 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
    text(ax, max(r.AJ, r.nullJ95) + 0.008, y, sprintf('p = %.3f', r.pJ), 'FontSize', 7);
end
set(ax, 'YTick', 0:3, 'YTickLabel', fliplr(nameS), 'XLim', [0 0.12], 'YLim', [-0.6 3.6], 'TickLength', [0 0]);
xlabel(ax, {'size of the daily rhythm in J', '(cosine amplitude, cycles h^{-1})'});
title(ax, {'grey bar: largest rhythm (95%) a random', 'clock phase per fly gives from the same data'}, ...
    'FontWeight', 'normal', 'FontSize', 6.5);
text(nexttile(tl, 1), -0.3, 1.12, 'a', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
text(ax, -0.3, 1.06, 'b', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
xlabel(tl, ['Shading: dark phase (light' char(8211) 'dark) or constant darkness. Grey band in (a): mean ' char(177) ' the null amplitude.'], ...
    'FontSize', 6.5, 'FontName', S.font);
exportgraphics(fh, fullfile(outDir, "figC10_clock_face.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "figC10_clock_face.pdf"), 'ContentType', 'vector');
writetable(bins, fullfile(outDir, "figC10_bins.csv"));
writetable(amp, fullfile(outDir, "figC10_amplitude.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "figC10.mat"), '-struct', 'Xs');
end
