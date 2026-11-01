function X = figC8_wake_hazard(outDir, nullReps)
% FIGC8_WAKE_HAZARD  Fig C8: wake-bout hazard of per0 DD and CS DD at the four fixed bin counts
% (PLANS C8; ledger entries "20r replication: per0 DD wake hazard rise at the four fixed bin counts"
% and "20r-A: CS DD wake hazard rise (blind)"). Run BY NAME from the Icode folder.
%
%   X = figC8_wake_hazard(outDir)
%   X = figC8_wake_hazard(outDir, nullReps)     % default 999, as the ledgered runs
%
% Recomputes, with the settings of RISE_TEST_20R unchanged (null = FitHyperexponentialMLE
% MaxComponents 2, discrete, SamplingInterval 1, RandomSeed 1, its K = 2 SurvivalHandle;
% BoutHazard Spacing "log", MinAtRisk 10, Alpha 0.05, Bootstrap 0, NullReplicates 999,
% RandomSeed 1), from the saved reference bouts:
%   per0 DD  per0DD_bouts_20r_replication.mat (SHA-256 b14ba340...), n = 4017
%   CS DD    CSDD_bouts_20r.mat (SHA-256 80402183...), n = 3309
% With the same code and seeds the RiseNullP values must equal the ledger's; X.check compares.
%
% Figure (theme pinned light): rows per0 DD, CS DD; columns 12 / 16 / 20 / 24 bins. Points =
% life-table hazard per bin (per minute), bars = Wilson 95% band (pointwise), dashed = the
% two-exponential null's hazard on the same bins, open circles = trough and the later
% maximum the rise statistic uses; text = RiseRatio and RiseNullP. Each row states the verdict
% by the 2026-09-15 rule. Writes figC8_wake_hazard.png / .pdf, figC8_hazard_bins.csv (every bin)
% and figC8_summary.csv (one row per data set x bin count).

assert(isequal(fileparts(which('BoutHazard')), pwd), ...
    'Run from the Icode folder: BoutHazard resolves to %s', which('BoutHazard'));
assert(isequal(fileparts(which('FitHyperexponentialMLE')), pwd), ...
    'FitHyperexponentialMLE resolves to %s', which('FitHyperexponentialMLE'));
if nargin < 2 || isempty(nullReps), nullReps = 999; end
if ~isfolder(outDir), mkdir(outDir); end
here = fullfile(pwd, 'analysis_20r');
src = {'per0 DD', fullfile(here, 'per0DD_bouts_20r_replication.mat'), [237 161 0] / 255, ...
       [0.001 0.009 0.059 0.086];
       'CS DD', fullfile(here, 'CSDD_bouts_20r.mat'), [98 80 214] / 255, ...
       [1.000 1.000 0.994 0.140]};                        % ledger RiseNullP, 12 / 16 / 20 / 24 bins
nb = [12 16 20 24];
bins = table(); summ = table(); H = cell(2, 4);
for d = 1:2
    B = load(src{d, 2});
    x = B.Wake(isfinite(B.Wake) & B.Wake >= B.xminWake); x = x(:);
    F = FitHyperexponentialMLE(x, B.xminWake, MaxComponents=2, DistributionType="discrete", ...
        SamplingInterval=1, RandomSeed=1, Verbose=false);
    f2 = F.AllFits(2);
    for i = 1:4
        h = BoutHazard(x, B.xminWake, SamplingInterval=1, NumBins=nb(i), Spacing="log", MinAtRisk=10, ...
            Alpha=0.05, Bootstrap=0, NullSurvivalHandle=f2.SurvivalHandle, NullReplicates=nullReps, RandomSeed=1);
        h.ModelHazard = nullHazard(h, f2.SurvivalHandle);
        H{d, i} = h;
        k = numel(h.BinCenters);
        bins = [bins; table(repmat(string(src{d, 1}), k, 1), repmat(nb(i), k, 1), h.BinCenters(:), ...
            h.Hazard(:) * 60, h.Lower(:) * 60, h.Upper(:) * 60, h.ModelHazard(:) * 60, h.AtRisk(:), h.Events(:), ...
            'VariableNames', {'dataset', 'NumBins', 'BinCenter_s', 'Hazard_per_min', 'Lower', 'Upper', ...
            'NullHazard_per_min', 'AtRisk', 'Events'})]; %#ok<AGROW>
        summ = [summ; table(string(src{d, 1}), numel(x), nb(i), h.BinCenters(h.TroughIndex), ...
            h.BinCenters(h.PeakIndex), h.RiseRatio, h.RiseNullQuantiles(1), h.RiseNullQuantiles(2), ...
            h.RiseNullP, src{d, 4}(i), h.NullDropped, 'VariableNames', {'dataset', 'n', 'NumBins', ...
            'trough_s', 'peak_s', 'RiseRatio', 'nullMedian', 'null95', 'RiseNullP', 'ledgerP', 'NullDropped'})]; %#ok<AGROW>
    end
    X.null(d) = struct('dataset', string(src{d, 1}), 'tau', f2.Tau(:)', 'q', f2.WeightsObserved(:)', 'logL', f2.LogLik);
end
summ.matchesLedger = abs(round(summ.RiseNullP, 3) - summ.ledgerP) < 1e-9;
X.bins = bins; X.summary = summ; X.check = all(summ.matchesLedger);
for d = 1:2
    p = summ.RiseNullP(summ.dataset == src{d, 1});
    if any(summ.NullDropped(summ.dataset == src{d, 1}) > 0), v = "not computable";
    elseif all(p < 0.05), v = "supported";
    elseif any(p < 0.05), v = "supported but bin-sensitive";
    else, v = "not supported"; end
    X.verdict.(matlab.lang.makeValidName(src{d, 1})) = v;
end

fh = figure('Name', 'Fig C8', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 7.0 3.9], ...
    'PaperUnits', 'inches', 'PaperSize', [7.0 3.9], 'PaperPosition', [0 0 7.0 3.9]);
try theme(fh, 'light'); catch, end
set(fh, 'DefaultAxesFontName', 'Times New Roman', 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 7, 'DefaultTextFontSize', 7, 'DefaultAxesTickDir', 'in', 'DefaultAxesLineWidth', 0.6);
tl = tiledlayout(fh, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for d = 1:2
    c = src{d, 3};
    for i = 1:4
        h = H{d, i}; ax = nexttile(tl); hold(ax, 'on');
        ok = isfinite(h.Hazard) & h.Hazard > 0;
        ok = ok(:); t = h.BinCenters(:) / 60; hz = h.Hazard(:) * 60; lo = h.Lower(:) * 60; hi = h.Upper(:) * 60;
        errorbar(ax, t(ok), hz(ok), hz(ok) - lo(ok), hi(ok) - hz(ok), 'o', 'Color', c, ...
            'MarkerFaceColor', c, 'MarkerSize', 3, 'CapSize', 0, 'LineWidth', 0.7);
        plot(ax, t, h.ModelHazard(:) * 60, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
        plot(ax, t([h.TroughIndex h.PeakIndex]), hz([h.TroughIndex h.PeakIndex]), 'ko', 'MarkerSize', 6, 'LineWidth', 0.8);
        set(ax, 'XScale', 'log', 'YScale', 'log', 'Box', 'off', 'XLim', [0.02 500], 'YLim', [0.005 1.5], ...
            'XTick', [0.1 1 10 100], 'XTickLabel', {'0.1', '1', '10', '100'}, ...
            'YTick', [0.01 0.1 1], 'YTickLabel', {'0.01', '0.1', '1'}, 'XMinorTick', 'on', 'YMinorTick', 'on');
        r = summ(summ.dataset == src{d, 1} & summ.NumBins == nb(i), :);
        text(ax, 0.97, 0.97, sprintf('RiseRatio %.2f\nRiseNullP %.3f', r.RiseRatio, r.RiseNullP), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 6.5);
        if d == 1, title(ax, sprintf('%d bins', nb(i)), 'FontWeight', 'normal', 'FontSize', 8); end
        if i == 1
            ylabel(ax, {sprintf('%s (n = %d)', src{d, 1}, r.n), 'wake hazard (min^{-1})'}, 'Color', c);
            text(ax, -0.42, 1.06, char('a' + d - 1), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
        end
        if d == 2, xlabel(ax, 'time awake (min)'); end
    end
end
title(tl, {'Wake-bout hazard (points, Wilson 95% band); dashed: the two-exponential null; circles: trough and the later maximum', ...
    sprintf('Rise test (2026-09-15 rule): per0 DD %s; CS DD %s', X.verdict.per0DD, X.verdict.CSDD)}, ...
    'FontSize', 7, 'FontName', 'Times New Roman');
X.fig = fh;
exportgraphics(fh, fullfile(outDir, 'figC8_wake_hazard.png'), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, 'figC8_wake_hazard.pdf'), 'ContentType', 'vector');
writetable(bins, fullfile(outDir, 'figC8_hazard_bins.csv'));
writetable(summ, fullfile(outDir, 'figC8_summary.csv'));
disp(summ); fprintf('ledger RiseNullP reproduced at every count: %d\n', X.check);
end

function mh = nullHazard(h, S)
% The null's hazard on the kept bins, by BoutHazard's own formula (interval probability
% (S(left) - S(right)) / S(left), then -log1p(-p) / width). Computed here, not through
% BoutHazard's SurvivalHandle option: FitHyperexponentialMLE's handle returns a COLUMN for a
% row of edges, so there (sl - sr) ./ sl ./ width broadcasts to an n x n matrix (reported).
L = h.BinEdges(1:end - 1); R = h.BinEdges(2:end);
C = sqrt(max(L, h.SamplingInterval / 2) .* R);
keep = ismembertol(C, h.BinCenters, 1e-12);
L = L(keep); R = R(keep); assert(numel(L) == numel(h.BinCenters), 'figC8: bins not recovered');
sl = reshape(S(L), [], 1); sr = reshape(S(R), [], 1);
p = (sl - sr) ./ sl; p(~(sl > 0)) = NaN;
mh = (-log1p(-min(p, 1)) ./ (R(:) - L(:)))';
end
