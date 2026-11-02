function X = exportFigDwellDirection(files, outDir, options)
%EXPORTFIGDWELLDIRECTION  Direction-dependent dwell times in CS at night (fig_dwell_direction.py; descriptive, as PRESPEC 20v).
%
%   X = exportFigDwellDirection(files, outDir)
%   X = exportFigDwellDirection(files, outDir, NBoot=1000, RandomSeed=3)
%
%   files   CS count files (a folder: SQWLISTFILES)
%   1-min symbols (SQWFLYINPUT: theta 0.5, light/dark references, missing = 0). A run of Q (or W)
%   counts if it is not the record's first or last run, its entry and exit neighbours are
%   classified symbols (not missing), and the run with both neighbours lies in the dark phase.
%   Q runs: W->Q->S ("into sleep") vs S->Q->W ("out of sleep"); W runs: Q->W->S vs S->W->Q.
%   (a, b) survival P(dwell >= t), t = 1..30 min, flies pooled, fly-bootstrap 95% band; (c) per
%   fly mean dwell (flies with >= 2 runs of a kind), the two directions joined per fly.
%   Writes fig_dwell_direction_runs.csv, _survival.csv, _perfly.csv, .mat and
%   fig_dwell_direction_CS_night.png / .pdf.
%   Check (ledger 2026-11-01): Q into / out n 649 / 433, per-fly mean 2.68 vs 2.38 min, into longer
%   in 64%; W out / into n 585 / 418, 4.43 vs 3.47 min, out longer in 80% (25 flies).
%
%   See also SQWFLYINPUT, MINUTEPHASES.

arguments
    files
    outDir (1,1) string
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 1000
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 3
end
if ~isfolder(outDir), mkdir(outDir); end
files = sqwListFiles(files); nF = numel(files); runs = cell(nF, 1);
for i = 1:nF
    [y, ph] = sqwFlyInput(files(i), true);
    ch = [1; find(diff(y) ~= 0) + 1; numel(y) + 1]; r = zeros(0, 4);
    for k = 2:numel(ch) - 2
        a = ch(k); b = ch(k + 1); s = y(a);              % run = minutes a .. b-1
        if ~ismember(s, [1 2]), continue, end
        e = y(a - 1); x = y(b);
        if ~ismember(e, 1:3) || ~ismember(x, 1:3) || ~all(ph(a - 1:b) == 2), continue, end
        r(end + 1, :) = [s e x b - a]; %#ok<AGROW>
    end
    runs{i} = r;
end
pairs = {1, [2 3], 'into sleep  W\rightarrowQ\rightarrowS'; 1, [3 2], 'out of sleep  S\rightarrowQ\rightarrowW'; ...
         2, [1 3], 'into sleep  Q\rightarrowW\rightarrowS'; 2, [3 1], 'out of sleep  S\rightarrowW\rightarrowQ'};
T = (1:30)'; rs = RandStream('twister', 'Seed', options.RandomSeed);
surv = @(d) arrayfun(@(t) mean(d >= t), T);
SV = table(); PF = table(); per = cell(nF, 4); n = zeros(4, 1);
for j = 1:4
    for i = 1:nF
        r = runs{i}; per{i, j} = r(r(:, 1) == pairs{j, 1} & r(:, 2) == pairs{j, 2}(1) & r(:, 3) == pairs{j, 2}(2), 4);
    end
    pool = vertcat(per{:, j}); n(j) = numel(pool); sv = surv(pool); bs = zeros(numel(T), options.NBoot);
    for b = 1:options.NBoot
        pick = randi(rs, nF, nF, 1); d = vertcat(per{pick, j});
        if isempty(d), bs(:, b) = NaN; else, bs(:, b) = surv(d); end
    end
    lo = NaN(numel(T), 1); hi = lo;
    for t = 1:numel(T), q = quantile7(bs(t, isfinite(bs(t, :))), [0.025 0.975]); lo(t) = q(1); hi(t) = q(2); end
    SV = [SV; table(repmat(j, numel(T), 1), T, sv, lo, hi, 'VariableNames', {'kind', 't_min', 'S', 'lo', 'hi'})]; %#ok<AGROW>
    m = cellfun(@(d) mean(d), per(:, j)); m(cellfun(@numel, per(:, j)) < 2) = NaN;
    PF.("mean_" + j) = m;
end
X.runs = runs; X.n = n; X.survival = SV; X.perFly = PF; X.pairs = pairs(:, 3);
st = table();
for s = 1:2
    a = PF.("mean_" + (2 * s - 1)); b = PF.("mean_" + (2 * s)); m = isfinite(a) & isfinite(b);
    st = [st; table(extract("QW", s), nnz(m), n(2 * s - 1), n(2 * s), median(a(m)), median(b(m)), mean(a(m) > b(m)), ...
        'VariableNames', {'state', 'flies', 'n_into', 'n_out', 'median_into_min', 'median_out_min', 'frac_into_longer'})]; %#ok<AGROW>
end
X.summary = st;

S = planAStyle(); cols = [42 120 214; 235 104 52] / 255;
fh = planAFigure("dwell direction", 7.0, 2.6); X.fig = fh;
tl = tiledlayout(fh, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:2
    ax = nexttile(tl); hold(ax, 'on'); h = gobjects(2, 1);
    for jj = 1:2
        j = 2 * (s - 1) + jj; q = SV(SV.kind == j, :);
        patch(ax, [T; flipud(T)], [q.lo; flipud(q.hi)], cols(jj, :), 'FaceAlpha', 0.18, 'EdgeColor', 'none');
        h(jj) = plot(ax, T, q.S, '-', 'Color', cols(jj, :), 'LineWidth', 2);
    end
    set(ax, 'YScale', 'log', 'XLim', [1 30], 'YLim', [5e-3 1.05]);
    xlabel(ax, 'time in state (min)'); ylabel(ax, 'fraction of stays lasting \geq t');
    nmS = ["quiet wake (Q)" "active wake (W)"]; title(ax, sprintf('%s, CS at night', nmS(s)), 'FontSize', 8.5, 'FontWeight', 'normal');
    legend(ax, h, {sprintf('%s (n %d)', pairs{2 * s - 1, 3}, n(2 * s - 1)), sprintf('%s (n %d)', pairs{2 * s, 3}, n(2 * s))}, ...
        'Box', 'off', 'FontSize', 6.3, 'Location', 'northeast');
end
ax = nexttile(tl); hold(ax, 'on');
for s = 1:2
    a = PF.("mean_" + (2 * s - 1)); b = PF.("mean_" + (2 * s)); m = isfinite(a) & isfinite(b); x0 = 2 * (s - 1);
    plot(ax, [x0; x0 + 0.8] * ones(1, nnz(m)), [a(m)'; b(m)'], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
    plot(ax, x0 * ones(nnz(m), 1), a(m), 'o', 'MarkerSize', 3, 'MarkerFaceColor', cols(1, :), 'MarkerEdgeColor', cols(1, :));
    plot(ax, (x0 + 0.8) * ones(nnz(m), 1), b(m), 'o', 'MarkerSize', 3, 'MarkerFaceColor', cols(2, :), 'MarkerEdgeColor', cols(2, :));
end
set(ax, 'XTick', [0.4 2.4], 'XTickLabel', {'Q', 'W'}, 'XLim', [-0.5 3.3]);
ylabel(ax, 'mean stay per fly (min)'); title(ax, 'per fly: into vs out of sleep', 'FontSize', 8.5, 'FontWeight', 'normal');
xlabel(tl, ['DESCRIPTIVE. Blue: stays on the way into sleep; orange: stays on the way out of sleep. ' ...
    'Bands: 95% bootstrap over flies. Grey lines join one fly''s two means.'], 'FontSize', 6, 'FontName', S.font);
exportgraphics(fh, fullfile(outDir, "fig_dwell_direction_CS_night.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "fig_dwell_direction_CS_night.pdf"), 'ContentType', 'vector');
writetable(SV, fullfile(outDir, "fig_dwell_direction_survival.csv"));
writetable(PF, fullfile(outDir, "fig_dwell_direction_perfly.csv"));
writetable(st, fullfile(outDir, "fig_dwell_direction_summary.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "fig_dwell_direction.mat"), '-struct', 'Xs');
end
