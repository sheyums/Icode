function X = exportFigA4(L, perFlyDir, outDir, options)
%EXPORTFIGA4  Fig A4, beneath the observed cycle (FIGURE_RECIPES K4; fig_beneath.py): data and figure.
%
%   X = exportFigA4(L, perFlyDir, outDir)
%   X = exportFigA4(L, perFlyDir, outDir, RandomSeed=1, Plot=true)
%
%   L          PLANALIKELIHOODS table (CS and fumin)
%   perFlyDir  folder holding the per-fly CSVs that SQWMEMORYGROUP and SQWTURGROUP
%              write: <g>_memory_perfly.csv, CSDD_memory_halves_perfly.csv,
%              <g>_tur_perfly.csv, for CS, CSDD, per0LD, per0DD, fumin, sss
%   outDir     folder for figA4_a_dBIC.csv, figA4_b_shares.csv, figA4_c_excess.csv,
%              figA4_d_tur.csv, figA4_rows.csv (one row per plotted group: median,
%              CI), figA4.mat and the figure (figA4_beneath_the_surface.png / .pdf)
%
%   (a) 20n: per fly dBIC = BIC(chain) - BIC(HMM), n = T (> 0: hidden layer preferred).
%   (b) 20o: flies with HMM gain > 3 nats: R_mem = (LL_MC2 - LL_MC1)/gain (filled) and
%       R_drift = 1 - (LL_HMM6 - LL_MC1_6)/gain (open); values > 5 drawn at 5 (and
%       < -1.5 at -1.5); medians and CIs from the unclipped values.
%   (c) 20p: excess irreversibility E = D3 - median D3 under null A per fly, per block,
%       median and CI; CS DD as its two 12-h halves (amendment A, null A2). Filled =
%       the block passes the pre-registered rule (20p ledger verdict, reproduced by
%       VALIDATE_20PQ: per0 DD, fumin dark; CS DD's halves do not).
%   (d) 20q: B_c / sigma_chain per block, tau 3 h, amendment-D set (use_180); filled =
%       primary DD blocks (CS DD, per0 DD); line at 1 (the TUR bound).
%   Medians with a fly-bootstrap 95% CI (MEDIANCI, 2000 resamples). Genotype order rule.
%   Sign check (recipe): dBIC medians -16.6 (CS), -30.4 (fumin); 28% and 23% above 0.
%
%   See also PLANALIKELIHOODS, SQWMEMORYGROUP, SQWTURGROUP, EXPORTFIGA1.

arguments
    L table
    perFlyDir (1,1) string {mustBeFolder}
    outDir (1,1) string
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 1
    options.Plot (1,1) logical = true
end
if ~isfolder(outDir), mkdir(outDir); end
S = planAStyle(); rs = RandStream('twister', 'Seed', options.RandomSeed);
ci = @(v) medianCIrow(v, S, rs);
gi = @(g) find(S.order == g);

% ---- (a), (b)
A = L(:, {'genotype', 'fly0', 'file', 'T', 'deg', 'MC1_12', 'HMM_12', 'dBIC_T'});
Bt = L(L.gain > 3, {'genotype', 'fly0', 'file', 'gain', 'MC1_12', 'HMM_12', 'MC2_12', 'MC1_6', 'HMM_6', 'R_mem', 'R_drift'});
rows = table();
for g = ["CS" "fumin"]
    v = A.dBIC_T(A.genotype == g);
    r = ci(v); r.panel = "a"; r.group = g; r.what = "dBIC_T"; r.fracAbove0 = mean(v > 0);
    r.gainPer1000min = median(1000 * (L.HMM_12(L.genotype == g) - L.MC1_12(L.genotype == g)) ./ L.T(L.genotype == g));
    rows = [rows; r]; %#ok<AGROW>
    for w = ["R_mem" "R_drift"]
        r = ci(Bt.(w)(Bt.genotype == g)); r.panel = "b"; r.group = g; r.what = w; r.fracAbove0 = NaN; r.gainPer1000min = NaN;
        rows = [rows; r]; %#ok<AGROW>
    end
end

% ---- (c) 20p
specC = {"CS light", "CS", "memory", "light", "E_A", false; "CS dark", "CS", "memory", "dark", "E_A", false; ...
    "CS DD, ZT 0-12", "CSDD", "memory_halves", "ZT0-12", "E_A2", false; ...
    "CS DD, ZT 12-24", "CSDD", "memory_halves", "ZT12-24", "E_A2", false; ...
    "per0 LD light", "per0LD", "memory", "light", "E_A", false; "per0 LD dark", "per0LD", "memory", "dark", "E_A", false; ...
    "per0 DD", "per0DD", "memory", "whole", "E_A", true; "fmn light", "fumin", "memory", "light", "E_A", false; ...
    "fmn dark", "fumin", "memory", "dark", "E_A", true; "sss light", "sss", "memory", "light", "E_A", false; ...
    "sss dark", "sss", "memory", "dark", "E_A", false};
Cp = table(); vc = cell(size(specC, 1), 1);
for j = 1:size(specC, 1)
    [lab, g, src, ph, col, filled] = specC{j, :};
    T = readtable(fullfile(perFlyDir, g + "_" + src + "_perfly.csv"), 'TextType', 'string');
    T = T(T.phase == ph & T.hasRef == 1, :); vc{j} = T.(col);
    Cp = [Cp; table(repmat(lab, height(T), 1), T.genotype, T.file, T.fly, T.phase, T.(col), ...
        'VariableNames', {'row', 'genotype', 'file', 'fly', 'phase', 'E'})]; %#ok<AGROW>
    r = ci(vc{j}); r.panel = "c"; r.group = lab; r.what = "E (nats/h)"; r.fracAbove0 = mean(vc{j} > 0);
    r.gainPer1000min = NaN; r.filled = filled; rows = [rows; r]; %#ok<AGROW>
end

% ---- (d) 20q
specD = {"CS light", "CS", "light", false; "CS dark", "CS", "dark", false; "CS DD", "CSDD", "whole", true; ...
    "per0 LD light", "per0LD", "light", false; "per0 LD dark", "per0LD", "dark", false; "per0 DD", "per0DD", "whole", true; ...
    "fmn light", "fumin", "light", false; "fmn dark", "fumin", "dark", false; "sss light", "sss", "light", false; ...
    "sss dark", "sss", "dark", false};
Dp = table(); vd = cell(size(specD, 1), 1);
for j = 1:size(specD, 1)
    [lab, g, ph, filled] = specD{j, :};
    T = readtable(fullfile(perFlyDir, g + "_tur_perfly.csv"), 'TextType', 'string');
    T = T(T.phase == ph & T.use_180 == 1, :); vd{j} = T.B_c_180 ./ T.sig_mc;
    Dp = [Dp; table(repmat(lab, height(T), 1), T.genotype, T.file, T.fly, T.phase, T.B_c_180, T.sig_mc, vd{j}, ...
        'VariableNames', {'row', 'genotype', 'file', 'fly', 'phase', 'B_c_180', 'sig_mc', 'ratio'})]; %#ok<AGROW>
    r = ci(vd{j}); r.panel = "d"; r.group = lab; r.what = "B_c/sigma_chain"; r.fracAbove0 = mean(vd{j} > 1);
    r.gainPer1000min = NaN; r.filled = filled; rows = [rows; r]; %#ok<AGROW>
end
rows = movevars(rows, {'panel', 'group', 'what'}, 'Before', 1);
rows = renamevars(rows, 'fracAbove0', 'fracAbove0or1');
X.a = A; X.b = Bt; X.c = Cp; X.d = Dp; X.rows = rows; X.fig = [];

if options.Plot
    fh = planAFigure("Fig A4", 7.0, 5.6); X.fig = fh;
    tl = tiledlayout(fh, 2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
    % (a)
    ax = nexttile(tl); hold(ax, 'on');
    for k = 1:2
        g = ["CS" "fumin"]; v = A.dBIC_T(A.genotype == g(k)); c = S.col(gi(g(k)), :);
        plot(ax, k - 1 + (rand(rs, numel(v), 1) * 0.3 - 0.15), v, 'o', 'MarkerSize', 3, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
        r = rows(rows.panel == "a" & rows.group == g(k), :);
        plot(ax, k - 1 + [-0.28 0.28], [r.median r.median], 'k-', 'LineWidth', 1.2); plot(ax, [k k] - 1, [r.lo r.hi], 'k-', 'LineWidth', 0.8);
    end
    yline(ax, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
    set(ax, 'XTick', [0 1], 'XTickLabel', S.name([1 5]), 'XLim', [-0.6 1.6]);
    ylabel(ax, {'BIC(plain chain) - BIC(HMM)', '(> 0: hidden layer preferred)'});
    ra = rows(rows.panel == "a", :);
    xlabel(ax, sprintf('hidden layer preferred by BIC in %.0f%% (CS), %.0f%% (fmn) of flies;\nlikelihood gain ~%.1f / %.1f nats per 1000 min', ...
        100 * ra.fracAbove0or1(1), 100 * ra.fracAbove0or1(2), ra.gainPer1000min(1), ra.gainPer1000min(2)), 'FontSize', 6.5);
    title(ax, 'a   Hidden layer vs plain S/Q/W chain (20n)', 'FontSize', 8.5, 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [-0.3 1.02]);
    % (b)
    ax = nexttile(tl); hold(ax, 'on');
    for k = 1:2
        g = ["CS" "fumin"]; c = S.col(gi(g(k)), :);
        for w = ["R_mem" "R_drift"]
            off = -0.18 + 0.36 * (w == "R_drift"); x = k - 1 + off; v = Bt.(w)(Bt.genotype == g(k));
            if w == "R_mem", mf = c; else, mf = 'w'; end
            plot(ax, x + (rand(rs, numel(v), 1) * 0.12 - 0.06), min(max(v, -1.5), 5), 'o', 'MarkerSize', 2.6, ...
                'MarkerFaceColor', mf, 'MarkerEdgeColor', c, 'LineWidth', 0.6);
            r = rows(rows.panel == "b" & rows.group == g(k) & rows.what == w, :);
            plot(ax, x + [-0.12 0.12], [r.median r.median], 'k-', 'LineWidth', 1.2);
            plot(ax, [x x], [max(r.lo, -1.5) min(r.hi, 5)], 'k-', 'LineWidth', 0.8);
        end
    end
    yline(ax, 1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6); yline(ax, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
    set(ax, 'XTick', [0 1], 'XTickLabel', S.name([1 5]), 'XLim', [-0.6 1.6], 'YLim', [-1.6 5.2]);
    ylabel(ax, {'share of the HMM''s likelihood gain', 'matched by ...'});
    l1 = plot(ax, NaN, NaN, 'o', 'MarkerSize', 3, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0.3 0.3 0.3]);
    l2 = plot(ax, NaN, NaN, 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.3 0.3 0.3]);
    legend(ax, [l1 l2], {'one more minute of memory (2nd-order chain)', '6-h drift of the rates'}, 'Box', 'off', ...
        'FontSize', 6.3, 'Location', 'southoutside');
    text(ax, 0.98, 0.98, 'values > 5 shown at 5', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 6);
    title(ax, 'b   What the HMM''s gain is made of (20o)', 'FontSize', 8.5, 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [-0.3 1.02]);
    % (c), (d)
    forest(nexttile(tl), rows(rows.panel == "c", :), specC(:, 2), S, gi, 0, ...
        'extra irreversibility from 3-min memory, E (nats h^{-1})', {'filled: passes the', 'pre-registered rule'}, ...
        'c   Does the memory carry dissipation? (20p)', []);
    forest(nexttile(tl), rows(rows.panel == "d", :), specD(:, 2), S, gi, 1, ...
        'cycle precision / its observed dissipation, B_c/\sigma_{chain}', {'filled: primary (DD);', 'open: LD, descriptive'}, ...
        'd   Is the cycle more precise than it pays for? (20q)', [0.2 1.5]);
    exportgraphics(fh, fullfile(outDir, "figA4_beneath_the_surface.png"), 'Resolution', 300);
    exportgraphics(fh, fullfile(outDir, "figA4_beneath_the_surface.pdf"), 'ContentType', 'vector');
end
writetable(X.a, fullfile(outDir, "figA4_a_dBIC.csv")); writetable(X.b, fullfile(outDir, "figA4_b_shares.csv"));
writetable(X.c, fullfile(outDir, "figA4_c_excess.csv")); writetable(X.d, fullfile(outDir, "figA4_d_tur.csv"));
writetable(X.rows, fullfile(outDir, "figA4_rows.csv"));
Xs = rmfield(X, 'fig'); save(fullfile(outDir, "figA4.mat"), '-struct', 'Xs');
end

function r = medianCIrow(v, S, rs)
[m, lo, hi] = medianCI(v, NBoot=S.NBoot, Level=S.Level, Stream=rs);
r = table(nnz(isfinite(v)), m, lo, hi, 'VariableNames', {'n', 'median', 'lo', 'hi'});
r.filled = false;
end

function forest(ax, R, geno, S, gi, ref, xlab, note, ttl, xl)
hold(ax, 'on'); n = height(R);
for j = 1:n
    c = S.col(gi(geno{j}), :); y = n - j;
    plot(ax, [R.lo(j) R.hi(j)], [y y], '-', 'Color', c, 'LineWidth', 1.4);
    if R.filled(j), mf = c; else, mf = 'w'; end
    plot(ax, R.median(j), y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', mf, 'MarkerEdgeColor', c, 'LineWidth', 1.2);
end
if ref == 1
    xline(ax, 1, 'k-', 'LineWidth', 0.9); text(ax, 1.02, n - 0.6, 'TUR bound', 'FontSize', 6.5);
else
    xline(ax, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
end
lab = replace(R.group, ["per0" "fmn" "sss"], ["\itper^0\rm" "\itfmn\rm" "\itsss\rm"]);
set(ax, 'YTick', 0:n - 1, 'YTickLabel', flipud(lab), 'YLim', [-0.7 n - 0.3], 'TickLabelInterpreter', 'tex');
ax.YAxis.FontSize = 7;
if ~isempty(xl), xlim(ax, xl); end
xlabel(ax, xlab);
text(ax, 0.98, 0.02, note, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 6.3);
title(ax, ttl, 'FontSize', 8.5, 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [-0.3 1.02]);
end
