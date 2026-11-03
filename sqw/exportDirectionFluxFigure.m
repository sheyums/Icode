function X = exportDirectionFluxFigure(D, outDir, options)
%EXPORTDIRECTIONFLUXFIGURE  Fig A1, direction and flux at the observable level (FIGURE_RECIPES K1; fig_planA.py): data and figure.
%
%   X = exportDirectionFluxFigure(D, outDir)
%   X = exportDirectionFluxFigure(D, outDir, RandomSeed=20, Plot=true, HmmMarkers=false)
%
%   D       BUILDPERFLYDATA output (six genotypes in the paper's order)
%   outDir  folder for figA1_points.csv, figA1_summary.csv, figA1.mat and the
%           figure (figA1_direction_flux_observable.png / .pdf; with HmmMarkers,
%           figA1_direction_flux_SI_hmm.png / .pdf)
%
%   Per fly x phase with a Q/W reference (the others excluded):
%     (a) c_obs, (b) counted J (cycles/h), (c) counted EPR (nats/h, rule of three)
%   The main text is observed-level only (Sheyum, 2026-11-01; FIGURE_RECIPES L note), so
%   the figure has no HMM markers by default. HmmMarkers=true draws the SI version: in
%   (b) and (c) a grey triangle beside each block is the median over NON-DEGENERATE flies
%   of the HMM's value, averaged over 200 FFBS paths per fly (HMMPOSTERIORFLUXEPR): "HMM
%   (model-based estimate), median" -- not a bound. The summary carries those medians
%   (nHmm, hmmMedian) either way.
%   Each block: points (light open, dark / DD filled), median, fly-bootstrap 95% CI.
%
%   X  struct: points (table, one row per fly x phase plotted), summary (one row per
%      panel x block: n, median, CI, nHmm, hmmMedian), fig (handle or empty)
%   Check (recipe): CS dark c_obs +0.147, J +0.167, EPR 0.116 (20s-A ledger).
%
%   Example (from the Icode\sqw folder):
%       dataRoot = "D:\flies";          % one subfolder of 1-s count files per genotype
%       jsonl = fullfile("..", "hmm_audit", "sqw20b_fit2_theta0.5.jsonl");
%       D = buildPerFlyData(dataRoot, jsonl, Cache="perfly_cache.mat");
%       X = exportDirectionFluxFigure(D, "figures");
%       Xsi = exportDirectionFluxFigure(D, "figures", HmmMarkers=true);   % SI version
%
%   See also BUILDPERFLYDATA, PLOTSTRIPPANEL, EXPORTEPRONSAGERFIGURE.

arguments
    D (1,1) struct
    outDir (1,1) string
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
    options.Plot (1,1) logical = true
    options.HmmMarkers (1,1) logical = false
end
if ~isfolder(outDir), mkdir(outDir); end
S = paperFigureStyle(); B = S.blocks; nB = height(B);
assert(isequal(D.order, S.order), 'exportDirectionFluxFigure:order', 'D must hold the genotypes in the order %s.', strjoin(S.order));

P = table(); v = struct('c_obs', {cell(nB, 1)}, 'J', {cell(nB, 1)}, 'EPR', {cell(nB, 1)});
h = struct('J', {cell(nB, 1)}, 'EPR', {cell(nB, 1)});
for b = 1:nB
    T = D.G(B.g(b)).perFly; T = T(T.phase == B.phase(b) & T.refOK, :);
    for k = ["c_obs" "J" "EPR"], v.(k){b} = T.(k); end
    nd = ~T.hmm_deg;
    h.J{b} = T.J_hmm(nd); h.EPR{b} = T.EPR_hmm(nd);
    P = [P; table(repmat(B.label(b), height(T), 1), T.genotype, T.file, T.fly, T.phaseName, ...
        repmat(B.filled(b), height(T), 1), T.c_obs, T.J, T.EPR, T.rule_of_three, T.hmm_deg, T.J_hmm, T.EPR_hmm, ...
        'VariableNames', {'block', 'genotype', 'file', 'fly', 'phase', 'filled', 'c_obs', 'J', 'EPR', ...
        'rule_of_three', 'hmm_deg', 'J_hmm', 'EPR_hmm'})]; %#ok<AGROW>
end
P.block = tex2plain(P.block);

rs = RandStream('twister', 'Seed', options.RandomSeed);
X.fig = [];
if options.Plot
    fh = newPaperFigure("Fig A1", 7.0, 6.6); X.fig = fh;
    tl = tiledlayout(fh, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); ax = gobjects(3, 1);
    for k = 1:3, ax(k) = nexttile(tl); end
else
    ax = gobjects(3, 1);
    fh = figure('Visible', 'off'); for k = 1:3, ax(k) = axes(fh); end   % summaries only
end
T1 = plotStripPanel(ax(1), S, v.c_obs, Stream=rs);
T2 = plotStripPanel(ax(2), S, v.J, Hmm=h.J, Stream=rs);
T3 = plotStripPanel(ax(3), S, v.EPR, Hmm=h.EPR, Stream=rs, Zero=false);
X.summary = [addPanel(T1, "c_obs"); addPanel(T2, "J"); addPanel(T3, "EPR")];
% sss in (a, b): the median is 0 (many flies have no qualifying transitions, value exactly 0), so the
% per-fly MEAN with a fly-bootstrap 95% CI is drawn too, left of the median (2026-11-02 recipe note)
X.summary.mean = NaN(height(X.summary), 1); X.summary.mean_lo = X.summary.mean; X.summary.mean_hi = X.summary.mean;
mcol = [176 48 106] / 255; Tq = {T1, T2}; kq = ["c_obs" "J"];
for q = 1:2
    for b = reshape(find(S.order(S.blocks.g) == "sss"), 1, [])
        vv = v.(kq(q)){b}; vv = vv(isfinite(vv));
        bm = zeros(S.NBoot, 1);
        for i = 1:S.NBoot, bm(i) = mean(vv(randi(rs, numel(vv), numel(vv), 1))); end
        ci = quantile7(bm, [0.025 0.975]); xq = Tq{q}.x(b) - 0.38;
        r = X.summary.quantity == kq(q) & X.summary.label == S.blocks.label(b);
        X.summary.mean(r) = mean(vv); X.summary.mean_lo(r) = ci(1); X.summary.mean_hi(r) = ci(2);
        if options.Plot
            plot(ax(q), [xq xq], ci, '-', 'Color', mcol, 'LineWidth', 0.8);
            plot(ax(q), xq, mean(vv), 'd', 'MarkerSize', 3.6, 'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
        end
    end
end
X.summary.label = tex2plain(X.summary.label);
X.points = P;
if options.Plot
    ylabel(ax(1), {'observable current c_{obs}', '(P(Q before sleep) - P(Q after waking))'});
    ylabel(ax(2), {'cycle flux J (cycles h^{-1})', '(forward W\rightarrowQ\rightarrowS\rightarrowW > 0)'});
    ylabel(ax(3), {'entropy production (nats h^{-1})', 'count estimate, rule of three'});
    for k = 1:3, text(ax(k), -0.105, 1.04, char('a' + k - 1), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold'); end
    hold(ax(1), 'on');
    l1 = plot(ax(1), NaN, NaN, 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
    l2 = plot(ax(1), NaN, NaN, 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    l3 = plot(ax(1), NaN, NaN, 'k-', 'LineWidth', 1.2);
    l5 = plot(ax(1), NaN, NaN, 'd', 'MarkerSize', 3.6, 'MarkerFaceColor', mcol, 'MarkerEdgeColor', 'w');
    hl = [l1 l2 l3 l5]; stem = "figA1_direction_flux_observable";
    lt = {'light (one point per fly)', 'dark / DD', 'median, 95% CI', ...
        '\itsss\rm: mean, 95% CI (median is 0: many flies have no qualifying transitions)'};
    if options.HmmMarkers
        hl(5) = plot(ax(1), NaN, NaN, '<', 'MarkerSize', 4, 'MarkerFaceColor', [0.35 0.35 0.35], 'MarkerEdgeColor', [0.35 0.35 0.35]);
        lt{5} = 'HMM (model-based estimate), median'; stem = "figA1_direction_flux_SI_hmm";
    else
        delete(findobj(fh, 'Tag', 'hmm'));        % main text: observed level only
    end
    lg = legend(ax(1), hl, lt, 'NumColumns', 3, 'Box', 'off', 'FontSize', 6.3);
    lg.Layout.Tile = 'north';
    exportgraphics(fh, fullfile(outDir, stem + ".png"), 'Resolution', 300);
    exportgraphics(fh, fullfile(outDir, stem + ".pdf"), 'ContentType', 'vector');
else
    close(fh);
end
writetable(X.points, fullfile(outDir, "figA1_points.csv"));
writetable(X.summary, fullfile(outDir, "figA1_summary.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "figA1.mat"), '-struct', 'Xs');
end

function T = addPanel(T, name)
T = [table(repmat(name, height(T), 1), 'VariableNames', {'quantity'}), T];
end

function s = tex2plain(s)
s = erase(erase(erase(string(s), "\it"), "\rm"), "^");
end
