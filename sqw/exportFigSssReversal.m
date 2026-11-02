function X = exportFigSssReversal(D, outDir, options)
%EXPORTFIGSSSREVERSAL  Sleepless runs the cycle backwards by day (fig_sss_reversal.py; FIGURE_RECIPES note 2026-11-02).
%
%   X = exportFigSssReversal(D, outDir)
%   X = exportFigSssReversal(D, outDir, RandomSeed=2)
%
%   D   PLANADATA output. Per fly x phase (CS and sss, light and dark) with a Q/W reference and a
%       non-degenerate 20b fit: the observed chain's exact cycle flux J_chain (circle) joined to the
%       HMM's flux J_hmm (triangle, mean over FFBS paths; the check). Bar = mean over flies; text =
%       share of flies with J < 0 (observed | HMM) and the two-sided Wilcoxon signed-rank p of J
%       against 0 (WILCOXONP, as scipy's default). Descriptive; the pre-registered sign test is 20d.
%   Writes fig_sss_reversal_points.csv, fig_sss_reversal_summary.csv, fig_sss_reversal.mat and
%   fig_sss_reversal.png / .pdf.
%   Check (PLANS C2): sss light n 50, chain J mean -0.038, J < 0 in 64%, p 0.004; HMM J -0.062,
%   p 7e-5; sss dark chain p 0.59.
%   CAUTION: a fly whose observed chain has no cycle has J_chain ~ +-1e-15 (round-off), so its
%   sign, the "% J < 0" and (through tiny ranks) p depend on floating-point detail and differ
%   between implementations (sss light: 13 of 50 such flies). The figure text mirrors the Python;
%   the summary adds n_noCycle, neg_obs_cycleOnly (|J| < 1e-12 counted as 0) and p_obs_zeroed
%   (those set to exactly 0, dropped by the Wilcoxon as scipy's zero_method 'wilcox').
%
%   See also PLANADATA, WILCOXONP, EXPORTFIGA1.

arguments
    D (1,1) struct
    outDir (1,1) string
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 2
end
if ~isfolder(outDir), mkdir(outDir); end
S = planAStyle(); rs = RandStream('twister', 'Seed', options.RandomSeed);
B = table(["CS"; "CS"; "sss"; "sss"], [1; 2; 1; 2], ["CS\newlinelight"; "CS\newlinedark"; ...
    "\itsss\rm\newlinelight"; "\itsss\rm\newlinedark"], 'VariableNames', {'genotype', 'phase', 'label'});
pts = table(); summ = table();
for k = 1:4
    T = D.G(S.order == B.genotype(k)).perFly;
    T = T(T.phase == B.phase(k) & T.refOK & ~T.hmm_deg & isfinite(T.J_hmm), :);
    jo = T.J_chain; jh = T.J_hmm;
    pts = [pts; table(repmat(B.genotype(k), height(T), 1), repmat(B.phase(k), height(T), 1), T.file, jo, jh, ...
        'VariableNames', {'genotype', 'phase', 'file', 'J_chain', 'J_hmm'})]; %#ok<AGROW>
    jz = jo; jz(abs(jz) < 1e-12) = 0;    % no-cycle chains (|J| < 1e-12, as hasCycle): J is round-off
    summ = [summ; table(B.genotype(k), B.phase(k), numel(jo), mean(jo), mean(jh), mean(jo < 0), mean(jh < 0), ...
        wilcoxonP(jo), wilcoxonP(jh), nnz(jz == 0), mean(jz < 0), wilcoxonP(jz), 'VariableNames', {'genotype', 'phase', 'n', ...
        'mean_obs', 'mean_hmm', 'neg_obs', 'neg_hmm', 'p_obs', 'p_hmm', 'n_noCycle', 'neg_obs_cycleOnly', 'p_obs_zeroed'})]; %#ok<AGROW>
end
X.points = pts; X.summary = summ;

fh = planAFigure("sss reversal", 5.2, 3.0); X.fig = fh; ax = axes(fh); hold(ax, 'on');
for k = 1:4
    q = pts(pts.genotype == B.genotype(k) & pts.phase == B.phase(k), :); c = S.col(S.order == B.genotype(k), :);
    if B.phase(k) == 2, mf = c; else, mf = 'w'; end
    x0 = 3 * (k - 1);
    plot(ax, [x0; x0 + 1] * ones(1, height(q)), [q.J_chain'; q.J_hmm'], '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.5);
    scatter(ax, x0 + rand(rs, height(q), 1) * 0.24 - 0.12, q.J_chain, 9, 'o', 'MarkerFaceColor', mf, 'MarkerEdgeColor', c, 'LineWidth', 0.7);
    scatter(ax, x0 + 1 + rand(rs, height(q), 1) * 0.24 - 0.12, q.J_hmm, 9, '^', 'MarkerFaceColor', mf, 'MarkerEdgeColor', c, 'LineWidth', 0.7);
    plot(ax, x0 + [-0.3 0.3], mean(q.J_chain) * [1 1], 'k-', 'LineWidth', 1.4);
    plot(ax, x0 + 1 + [-0.3 0.3], mean(q.J_hmm) * [1 1], 'k-', 'LineWidth', 1.4);
    r = summ(k, :);
    text(ax, x0 + 0.5, 0.60, {sprintf('J < 0 in %.0f%% | %.0f%%', 100 * r.neg_obs, 100 * r.neg_hmm), ...
        sprintf('mean p %s | %s', g1(r.p_obs), g1(r.p_hmm))}, 'HorizontalAlignment', 'center', 'FontSize', 6);
end
yline(ax, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
set(ax, 'YLim', [-0.5 0.75], 'XTick', 3 * (0:3) + 0.5, 'XTickLabel', B.label, 'TickLabelInterpreter', 'tex', 'XLim', [-0.7 10.7]);
ylabel(ax, {'cycle flux J (cycles h^{-1})', '> 0: W\rightarrowQ\rightarrowS\rightarrowW,  < 0: reversed'});
h1 = scatter(ax, NaN, NaN, 9, 'o', 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', [0.4 0.4 0.4]);
h2 = scatter(ax, NaN, NaN, 9, '^', 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', [0.4 0.4 0.4]);
legend(ax, [h1 h2], {'observed S/Q/W chain', 'hidden Markov model (check)'}, 'Box', 'off', 'FontSize', 6.5, 'Location', 'southwest');
title(ax, 'Sleepless flies run the cycle backwards on average by day', 'FontSize', 9, 'FontWeight', 'normal');
xlabel(ax, {'Per fly: circle = observed chain, triangle = hidden Markov model, joined by a line; bar = mean.', ...
    'Percentages: flies with J < 0 (observed | HMM); p: Wilcoxon test of J against 0.'}, 'FontSize', 5.8);
exportgraphics(fh, fullfile(outDir, "fig_sss_reversal.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "fig_sss_reversal.pdf"), 'ContentType', 'vector');
writetable(pts, fullfile(outDir, "fig_sss_reversal_points.csv"));
writetable(summ, fullfile(outDir, "fig_sss_reversal_summary.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "fig_sss_reversal.mat"), '-struct', 'Xs');
end

function s = g1(p)
% one significant figure, as Python's {:.1g}
s = sprintf('%.1g', p);
end
