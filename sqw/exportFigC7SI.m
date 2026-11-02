function X = exportFigC7SI(c7Mat, H, outDir, options)
%EXPORTFIGC7SI  SI figure: sleep-length memory seen at the observed AND hidden level (fig_R1_main.py; figS_C7_memory_with_hmm_check).
%
%   X = exportFigC7SI(c7Mat, H, outDir)
%   X = exportFigC7SI(c7Mat, H, outDir, NBoot=2000, SigFrom="ledger")
%
%   c7Mat   figC7.mat written by EXPORTFIGC7 (observed per-fly D, blocks with the null band, bins)
%   H       PLANAROUTEHIDDEN table
%   (a) observed per-fly D (as Fig C7a). (b) hidden per-fly D: per fly the mean over 200 FFBS paths;
%       grey = 2.5-97.5% of the mean over flies under the null (each fly's records simulated from its
%       own fit, NSim times, decoded with 20 paths; sqw10_report.py); p = 2 min((1 + #>=)/(1 + N),
%       (1 + #<=)/(1 + N)), Holm over the 9 blocks. (c) CS dark: P(quiet exit) per sleep-length bin,
%       hidden (H1 exit, records of 20 paths pooled per fly) and observed (Q exit), flies weighted
%       equally (>= 10 flies), fly-bootstrap 95% CIs; dashed: the landscape prediction, mean over CS
%       flies of (F_SW + 1)/2 of the fit's dark matrix.
%   Stars and fading come from the 20k LEDGER by default (SigFrom="ledger"): with 400 simulations the
%   smallest p is 2/401 = 0.005, so CS dark's hidden Holm (9 x 0.005 = 0.045) sits at that floor.
%   X.hiddenBlocks carries the port's own values beside the ledger's.
%   Writes figS_C7_hidden_blocks.csv, figS_C7_hidden_perfly.csv, figS_C7_bins.csv, figS_C7.mat and
%   figS_C7_memory_with_hmm_check.png / .pdf.
%
%   See also EXPORTFIGC7, PLANAROUTEHIDDEN, ROUTEDURATIONHIDDEN.

arguments
    c7Mat (1,1) string {mustBeFile}
    H table
    outDir (1,1) string
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.SigFrom (1,1) string {mustBeMember(options.SigFrom, ["ledger" "port"])} = "ledger"
end
if ~isfolder(outDir), mkdir(outDir); end
C7 = load(c7Mat); ob = C7.blocks; S = planAStyle(); rs = RandStream('twister', 'Seed', 5);
% ---- hidden blocks, 20k order (Holm family)
hb = ob(:, {'genotype', 'phase', 'phaseName'});
led = [24 0.041 -0.031 0.195 1; 25 -0.275 0.017 0.00499 0.0449; 29 -0.053 -0.010 0.559 1; 18 -0.144 -0.013 0.354 1; ...
       27 -0.024 -0.008 0.813 1; 25 -0.039 0.009 0.409 1; 23 -0.027 -0.000 0.579 1; 18 0.021 0.000 0.868 1; 23 0.004 0.013 0.883 1];
hp = table();
for b = 1:height(hb)
    F = H(H.genotype == hb.genotype(b), :); p = hb.phase(b);
    D = cellfun(@(v) v(p), F.hidD); ok = isfinite(D);
    Sm = cell2mat(cellfun(@(s) s(:, p)', F.sims, 'UniformOutput', false));     % fly x sim
    sm = mean(Sm, 1, 'omitnan')'; real = mean(D(ok)); N = numel(sm);
    q = quantile7(sm, [0.025 0.975]);
    hb.n(b) = nnz(ok); hb.D(b) = real; hb.nullMean(b) = mean(sm); hb.null_lo(b) = q(1); hb.null_hi(b) = q(2);
    hb.p(b) = min(1, 2 * min((1 + sum(sm >= real)) / (1 + N), (1 + sum(sm <= real)) / (1 + N)));
    hp = [hp; table(repmat(hb.genotype(b), nnz(ok), 1), repmat(hb.phaseName(b), nnz(ok), 1), F.file(ok), D(ok), ...
        'VariableNames', {'genotype', 'phase', 'file', 'D'})]; %#ok<AGROW>
end
hb.holm = holmAdjust(hb.p);
hb.ledger_n = led(:, 1); hb.ledger_D = led(:, 2); hb.ledger_nullMean = led(:, 3); hb.ledger_p = led(:, 4); hb.ledger_holm = led(:, 5);
if options.SigFrom == "ledger", hb.sig = hb.ledger_holm < 0.05; else, hb.sig = hb.holm < 0.05; end
% ---- panel (c): CS dark bins, hidden from the path records
edges = [5 10 20 40 80 Inf]; F = H(H.genotype == "CS", :); bins = table();
for a = 1:numel(edges) - 1
    props = [];
    for j = 1:height(F)
        r = F.records{j};
        if isempty(r), continue, end
        q = [r.q_out]'; L = [r.sleep_min]'; m = [r.phase]' == 2 & ~isnan(q) & L >= edges(a) & L < edges(a + 1);
        if any(m), props(end + 1, 1) = mean(q(m)); end %#ok<AGROW>
    end
    if numel(props) < 10, continue, end
    bs = zeros(options.NBoot, 1);
    for i = 1:options.NBoot, bs(i) = mean(props(randi(rs, numel(props), numel(props), 1))); end
    ci = quantile7(bs, [0.025 0.975]);
    bins = [bins; table(edges(a), edges(a + 1), sqrt(edges(a) * min(edges(a + 1), 2 * edges(a))), mean(props), ci(1), ci(2), ...
        numel(props), 'VariableNames', {'from_min', 'to_min', 'x_min', 'Pquiet', 'lo', 'hi', 'nFlies'})]; %#ok<AGROW>
end
pred = mean(F.Fdark, 'omitnan');
X.hiddenBlocks = hb; X.hiddenPoints = hp; X.hiddenBins = bins; X.observedBins = C7.bins; X.landscapePred = pred;

% ---- figure
fh = planAFigure("Fig S C7", 7.0, 2.9); X.fig = fh;
tl = tiledlayout(fh, 1, 14, 'TileSpacing', 'loose', 'Padding', 'compact');
order = [1 2 5 6 7 3 4 8 9];
nameS = containers.Map(["CS" "per0LD" "per0DD" "fumin" "sss"], ["CS" "\itper^0\rm LD" "\itper^0\rm DD" "\itfmn\rm" "\itsss\rm"]);
lev = {ob, C7.points, '(a) observed behaviour'; hb, hp, '(b) hidden states'};
for c = 1:2
    ax = nexttile(tl, [1 5]); hold(ax, 'on'); B = lev{c, 1}; P = lev{c, 2}; labs = strings(1, 9);
    for x = 1:9
        b = order(x); g = B.genotype(b); col = S.col(S.order == g, :); xx = x - 1;
        v = P.D(P.genotype == g & P.phase == B.phaseName(b));
        patch(ax, xx + [-0.33 0.33 0.33 -0.33], [B.null_lo(b) B.null_lo(b) B.null_hi(b) B.null_hi(b)], [0.86 0.86 0.86], 'EdgeColor', 'none');
        if B.phaseName(b) == "light", mf = 'w'; else, mf = col; end
        sc = scatter(ax, xx + rand(rs, numel(v), 1) * 0.36 - 0.18, min(max(v, -1.5), 1.5), 7, 'MarkerFaceColor', mf, 'MarkerEdgeColor', col, 'LineWidth', 0.6);
        if ~B.sig(b), sc.MarkerFaceAlpha = 0.55; sc.MarkerEdgeAlpha = 0.55; end
        plot(ax, xx + [-0.28 0.28], mean(v) * [1 1], '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.6);
        if B.sig(b), text(ax, xx, 1.38, '*', 'HorizontalAlignment', 'center', 'FontSize', 11); end
        if B.phaseName(b) == "whole", labs(x) = nameS(g); else, labs(x) = nameS(g) + " " + extractBefore(upper(B.phaseName(b)), 2); end
    end
    yline(ax, 0, 'Color', [0.53 0.53 0.53], 'LineWidth', 0.5);
    set(ax, 'XTick', 0:8, 'XTickLabel', labs, 'XTickLabelRotation', 60, 'YLim', [-1.55 1.55], 'XLim', [-0.6 8.6], 'TickLabelInterpreter', 'tex');
    ax.XAxis.FontSize = 6.5;
    if c == 1, ylabel(ax, 'D = \langleln t\rangle_{quiet exit} - \langleln t\rangle_{active exit}'); end
    title(ax, lev{c, 3}, 'FontSize', 8, 'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02]);
end
ax = nexttile(tl, [1 4]); hold(ax, 'on'); col = S.col(1, :); ob2 = C7.bins;
h1 = errorbar(ax, bins.x_min, bins.Pquiet, bins.Pquiet - bins.lo, bins.hi - bins.Pquiet, 'o-', 'Color', col, 'MarkerFaceColor', col, ...
    'MarkerSize', 4, 'LineWidth', 1.0, 'CapSize', 3);
h2 = errorbar(ax, ob2.x_min, ob2.Pquiet, ob2.Pquiet - ob2.lo, ob2.hi - ob2.Pquiet, 's--', 'Color', col, 'MarkerFaceColor', 'w', ...
    'MarkerSize', 3.5, 'LineWidth', 1.0, 'CapSize', 3);
h3 = yline(ax, pred, '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 0.9);
set(ax, 'XScale', 'log', 'XLim', [4.5 160], 'XTick', [5 10 20 40 80], 'XTickLabel', {'5', '10', '20', '40', '80'}, 'YLim', [0 0.8], 'XMinorTick', 'off');
xlabel(ax, 'sleep-bout length (min)'); ylabel(ax, 'P(quiet exit)');
title(ax, '(c) CS, dark', 'FontSize', 8, 'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02]);
legend(ax, [h1 h2 h3], {'hidden (H1 exit)', 'observed (Q exit)', 'landscape (Markov) prediction'}, 'Box', 'off', 'FontSize', 6, 'Location', 'northeast');
xlabel(tl, {['(a, b) Points: flies (open light, filled dark/DD; faded where not significant); bar: mean; grey: 95% range of the ' ...
    'mean without dependence; *: Holm p < 0.05; |D| clipped at 1.5.'], ['(c) Mean over flies of each fly''s proportion per bin ' ...
    '(>= 10 flies), bootstrap 95% intervals; dashed: prediction from the fitted transition matrix, constant in bout length.']}, ...
    'FontSize', 5.6, 'Color', [0.2 0.2 0.2], 'FontName', S.font);
exportgraphics(fh, fullfile(outDir, "figS_C7_memory_with_hmm_check.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "figS_C7_memory_with_hmm_check.pdf"), 'ContentType', 'vector');
writetable(hb, fullfile(outDir, "figS_C7_hidden_blocks.csv"));
writetable(hp, fullfile(outDir, "figS_C7_hidden_perfly.csv"));
writetable(bins, fullfile(outDir, "figS_C7_bins.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "figS_C7.mat"), '-struct', 'Xs');
end
