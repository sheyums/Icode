out = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\ar1_run';
sets = ["per0dd_sleep_all","per0dd_wake_all","csld_sleep_light","csld_sleep_dark","csld_wake_light","csld_wake_dark"];
lab = ["per0 DD sleep","per0 DD wake","CS LD sleep, light","CS LD sleep, dark","CS LD wake, light","CS LD wake, dark"];
fams = ["ShiftedExp","AsymLaplace","ShiftWeibull","ShiftLogNormal","GH","GenNorm","Laplace","Normal"];
famLab = ["Shifted exponential","Asymmetric Laplace","Shifted Weibull","Shifted lognormal","GH family","Generalized normal","Laplace","Normal"];
cols = [0.30 0.45 0.70; 0.85 0.45 0.20; 0.40 0.65 0.35; 0.60 0.40 0.70; 0.80 0.25 0.35; 0.55 0.55 0.55; 0.75 0.70 0.30; 0.20 0.20 0.20];
fh = figure('Color', 'w', 'Position', [40 40 1500 900]);
try theme(fh, 'light'); catch, end
tl = tiledlayout(fh, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
rs = RandStream('mt19937ar', 'Seed', 3);
T = cell(1, 6);
for i = 1:6, T{i} = readtable(fullfile(out, sets(i) + "_perfly.csv"), 'TextType', 'string'); T{i} = T{i}(isfinite(T{i}.B), :); end
panels = {'B', 'Memory B: slope of next bout on this bout', 'B_P'; ...
          'Hetero', {'Spread dependence: rank corr. of |\epsilon| with this bout'}, 'Hetero_P'; ...
          'GenBelowXmin', {'Regenerated bouts below the floor','(300 s sleep, 2 s wake), fraction'}, ''};
for p = 1:3
    ax = nexttile(tl, p); hold(ax, 'on');
    for i = 1:6
        v = T{i}.(panels{p, 1});
        x = i + 0.25 * (rs.rand(size(v)) - 0.5);
        if isempty(panels{p, 3}), sig = false(size(v)); else, sig = T{i}.(panels{p, 3}) < 0.05; end
        scatter(ax, x(~sig), v(~sig), 12, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.6);
        scatter(ax, x(sig), v(sig), 24, [0.85 0.3 0.1], 'filled');
        plot(ax, i + [-0.3 0.3], median(v) * [1 1], 'k-', 'LineWidth', 2.5);
        if ~isempty(panels{p, 3})
            text(ax, i, 0, '', 'Visible', 'off');
        end
    end
    if p < 3, yline(ax, 0, ':k'); end
    set(ax, 'XTick', 1:6, 'XTickLabel', lab, 'XLim', [0.4 6.6]); xtickangle(ax, 30);
    ylabel(ax, panels{p, 2}); box(ax, 'on');
end
title(nexttile(tl, 1), 'A  Memory (orange: fly p < 0.05)', 'FontWeight', 'normal');
title(nexttile(tl, 2), 'B  Does the noise spread with the previous bout?', 'FontWeight', 'normal');
title(nexttile(tl, 3), 'C  Can the fitted model regenerate the bouts?', 'FontWeight', 'normal');
% D: BIC winners
ax = nexttile(tl, 4); hold(ax, 'on');
M = zeros(6, numel(fams));
for i = 1:6, for k = 1:numel(fams), M(i, k) = sum(T{i}.BICWinner == fams(k)); end, end
keepF = any(M > 0, 1);
bh = bar(ax, 1:6, M(:, keepF), 'stacked', 'EdgeColor', 'none');
c = cols(keepF, :); for k = 1:numel(bh), bh(k).FaceColor = c(k, :); end
legend(ax, famLab(keepF), 'Location', 'eastoutside', 'FontSize', 8);
set(ax, 'XTick', 1:6, 'XTickLabel', lab); xtickangle(ax, 30); ylabel(ax, 'Flies (BIC winner for \epsilon)');
title(ax, 'D  Which family describes \epsilon best?', 'FontWeight', 'normal'); box(ax, 'on');
% E: residual skew vs no-dynamics skew
ax = nexttile(tl, 5); hold(ax, 'on');
mk = 'osd^vp';
for i = 1:6
    scatter(ax, T{i}.ResidSkew_Null, T{i}.ResidSkew, 18, 'filled', 'Marker', mk(i), 'DisplayName', lab(i));
end
lim = [0 max(cellfun(@(t) max([t.ResidSkew; t.ResidSkew_Null]), T))];
plot(ax, lim, lim, 'k-', 'HandleVisibility', 'off');
xlabel(ax, 'Skewness of \epsilon with the bouts in random order'); ylabel(ax, 'Skewness of \epsilon, real order');
legend(ax, 'Location', 'northwest', 'FontSize', 8); box(ax, 'on');
title(ax, 'E  \epsilon has the shape the bouts give with no dynamics', 'FontWeight', 'normal');
% F: SD vs no-dynamics SD
ax = nexttile(tl, 6); hold(ax, 'on');
for i = 1:6
    scatter(ax, T{i}.ResidSD_Null / 60, T{i}.ResidSD / 60, 18, 'filled', 'Marker', mk(i), 'DisplayName', lab(i));
end
lim = [0 max(cellfun(@(t) max([t.ResidSD; t.ResidSD_Null]), T)) / 60];
plot(ax, lim, lim, 'k-', 'HandleVisibility', 'off');
xlabel(ax, 'SD of \epsilon, random order (min)'); ylabel(ax, 'SD of \epsilon, real order (min)');
box(ax, 'on'); title(ax, 'F  ... and almost the same width', 'FontWeight', 'normal');
title(tl, {'x(t+1) = A + B x(t) + \epsilon, fitted per fly in bout space: per0 constant darkness (48 flies), Canton-S 12:12 LD by phase (40 flies)', ...
    'per0: shuffle null; CS: null permuting bouts within 1-h ZT bins. Black bars = medians.'}, 'FontSize', 10);
exportgraphics(fh, fullfile(out, 'ar1_noise.png'), Resolution=150);
