function X = exportFigA1(D, outDir, options)
%EXPORTFIGA1  Fig A1, direction and flux at the observable level (FIGURE_RECIPES K1; fig_planA.py): data and figure.
%
%   X = exportFigA1(D, outDir)
%   X = exportFigA1(D, outDir, RandomSeed=20, Plot=true)
%
%   D       PLANADATA output (six genotypes in the paper's order)
%   outDir  folder for figA1_points.csv, figA1_summary.csv, figA1.mat and the
%           figure (figA1_direction_flux_observable.png / .pdf)
%
%   Per fly x phase with a Q/W reference (the others excluded):
%     (a) c_obs, (b) counted J (cycles/h), (c) counted EPR (nats/h, rule of three)
%   In (b) and (c) a grey triangle beside each block is the median over
%   NON-DEGENERATE flies of the HMM's value, averaged over 200 FFBS paths per fly
%   (PLANAFLYHMM): "HMM (model-based estimate), median" -- not a bound.
%   Each block: points (light open, dark / DD filled), median, fly-bootstrap 95% CI.
%
%   X  struct: points (table, one row per fly x phase plotted), summary (one row per
%      panel x block: n, median, CI, nHmm, hmmMedian), fig (handle or empty)
%   Check (recipe): CS dark c_obs +0.147, J +0.167, EPR 0.116 (20s-A ledger).
%
%   See also PLANADATA, PLANASTRIP, EXPORTFIGA2.

arguments
    D (1,1) struct
    outDir (1,1) string
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
    options.Plot (1,1) logical = true
end
if ~isfolder(outDir), mkdir(outDir); end
S = planAStyle(); B = S.blocks; nB = height(B);
assert(isequal(D.order, S.order), 'exportFigA1:order', 'D must hold the genotypes in the order %s.', strjoin(S.order));

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
    fh = planAFigure("Fig A1", 7.0, 6.6); X.fig = fh;
    tl = tiledlayout(fh, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); ax = gobjects(3, 1);
    for k = 1:3, ax(k) = nexttile(tl); end
else
    ax = gobjects(3, 1);
    fh = figure('Visible', 'off'); for k = 1:3, ax(k) = axes(fh); end   % summaries only
end
T1 = planAStrip(ax(1), S, v.c_obs, Stream=rs);
T2 = planAStrip(ax(2), S, v.J, Hmm=h.J, Stream=rs);
T3 = planAStrip(ax(3), S, v.EPR, Hmm=h.EPR, Stream=rs, Zero=false);
X.summary = [addPanel(T1, "c_obs"); addPanel(T2, "J"); addPanel(T3, "EPR")];
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
    l4 = plot(ax(1), NaN, NaN, '<', 'MarkerSize', 4, 'MarkerFaceColor', [0.35 0.35 0.35], 'MarkerEdgeColor', [0.35 0.35 0.35]);
    lg = legend(ax(1), [l1 l2 l3 l4], {'light (one point per fly)', 'dark / DD', 'median, 95% CI', ...
        'HMM (model-based estimate), median'}, 'NumColumns', 4, 'Box', 'off', 'FontSize', 6.5);
    lg.Layout.Tile = 'north';
    exportgraphics(fh, fullfile(outDir, "figA1_direction_flux_observable.png"), 'Resolution', 300);
    exportgraphics(fh, fullfile(outDir, "figA1_direction_flux_observable.pdf"), 'ContentType', 'vector');
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
