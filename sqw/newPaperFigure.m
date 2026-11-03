function fh = newPaperFigure(name, widthIn, heightIn)
%NEWPAPERFIGURE  New figure for a Plan A panel set: theme pinned light, size in inches, Plan A font.
%
%   fh = newPaperFigure(name, widthIn, heightIn)
%
%   The theme is pinned BEFORE any axes exist so children inherit it (root
%   CLAUDE.md: a new R2025a+ figure takes an unpredictable theme). Paper size
%   matches the figure, so PDF export keeps the stated size.

S = paperFigureStyle();
fh = figure('Name', name, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 widthIn heightIn], ...
    'PaperUnits', 'inches', 'PaperSize', [widthIn heightIn], 'PaperPosition', [0 0 widthIn heightIn]);
try theme(fh, 'light'); catch, end
set(fh, 'DefaultAxesFontName', S.font, 'DefaultTextFontName', S.font, 'DefaultAxesFontSize', S.size, ...
    'DefaultTextFontSize', S.size, 'DefaultAxesTickDir', 'in', 'DefaultAxesLineWidth', 0.6, ...
    'DefaultAxesBox', 'off', 'DefaultAxesColor', 'w', 'DefaultAxesXColor', 'k', 'DefaultAxesYColor', 'k');
end
