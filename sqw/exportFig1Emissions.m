function X = exportFig1Emissions(jsonlFile, outDir)
%EXPORTFIG1EMISSIONS  Fig 1 panel: the HMM's hidden states coincide with the observed symbols (fig_emission.py; note 2026-11-02).
%
%   X = exportFig1Emissions(jsonlFile, outDir)
%
%   jsonlFile  hmm_audit/sqw20b_fit2_theta0.5.jsonl (READSTORED20BFITS)
%   Per genotype (CS, per0 LD, per0 DD, fumin, sss): the element-wise median over flies with a
%   NON-DEGENERATE fit of the labelled K = 3 emission matrix B3 (rows H1, H2, H3 = quiet / active /
%   sleep-like hidden state; columns observed Q, W, S), greyscale 0-1 with values printed.
%   Writes fig1_hmm_emissions.csv (one row per genotype x state), fig1_hmm_emissions.mat and
%   fig1_hmm_emissions.png / .pdf.
%   Check (recipe): diagonals Q 0.80-0.94, W 0.91-0.96, S 1.00.
%
%   See also READSTORED20BFITS.

arguments
    jsonlFile (1,1) string {mustBeFile}
    outDir (1,1) string
end
if ~isfolder(outDir), mkdir(outDir); end
F = readStored20bFits(jsonlFile); S = planAStyle();
G = ["CS" "per0LD" "per0DD" "fumin" "sss"]; nm = ["CS" "\itper^0\rm LD" "\itper^0\rm DD" "\itfmn\rm" "\itsss\rm"];
X.M = zeros(3, 3, 5); X.n = zeros(1, 5); T = table();
for k = 1:5
    f = F([F.g] == G(k) & ~[F.deg]);
    X.n(k) = numel(f); X.M(:, :, k) = median(cat(3, f.B), 3);
    T = [T; table(repmat(G(k), 3, 1), repmat(X.n(k), 3, 1), ["H1"; "H2"; "H3"], X.M(:, 1, k), X.M(:, 2, k), X.M(:, 3, k), ...
        'VariableNames', {'genotype', 'n', 'state', 'Q', 'W', 'S'})]; %#ok<AGROW>
end
X.table = T; X.diag = [squeeze(X.M(1, 1, :)) squeeze(X.M(2, 2, :)) squeeze(X.M(3, 3, :))];
fh = planAFigure("Fig 1 emissions", 7.0, 1.9); X.fig = fh;
tl = tiledlayout(fh, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:5
    ax = nexttile(tl); M = X.M(:, :, k);
    image(ax, repmat(1 - M, 1, 1, 3)); axis(ax, 'image'); hold(ax, 'on');
    for i = 1:3
        for j = 1:3
            if M(i, j) > 0.55, tc = 'w'; else, tc = 'k'; end
            text(ax, j, i, sprintf('%.2f', M(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 6.5, 'Color', tc);
        end
    end
    set(ax, 'XTick', 1:3, 'XTickLabel', {'Q', 'W', 'S'}, 'YTick', 1:3, 'TickLength', [0 0]);
    if k == 1, set(ax, 'YTickLabel', {'H_1', 'H_2', 'H_3'}); ylabel(ax, 'hidden state'); else, set(ax, 'YTickLabel', []); end
    title(ax, sprintf('%s (n %d)', nm(k), X.n(k)), 'FontWeight', 'normal', 'FontSize', 7.5);
end
xlabel(tl, 'observed symbol (median emission probability over flies)', 'FontSize', 7, 'FontName', S.font);
exportgraphics(fh, fullfile(outDir, "fig1_hmm_emissions.png"), 'Resolution', 300);
exportgraphics(fh, fullfile(outDir, "fig1_hmm_emissions.pdf"), 'ContentType', 'vector');
writetable(T, fullfile(outDir, "fig1_hmm_emissions.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "fig1_hmm_emissions.mat"), '-struct', 'Xs');
end
