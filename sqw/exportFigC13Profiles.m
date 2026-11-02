function X = exportFigC13Profiles(dataRoot, refMat, outDir, options)
%EXPORTFIGC13PROFILES  C13 / PRESPEC 21 descriptive figure: 1-s activity at the two ends of wake bouts (fig_planB_21.py).
%
%   X = exportFigC13Profiles(dataRoot, refMat, outDir)
%   X = exportFigC13Profiles(dataRoot, refMat, outDir, NBoot=2000, RandomSeed=2026, Reference="")
%
%   dataRoot  count files per genotype: per0DD, CSDD, CS, per0LD, fumin, sss (all 57)
%   refMat    analysis_20r/per0DD_bouts_20r_replication.mat: per0 DD traces are cut to its
%             traceLen (matched by file-name prefix), as planB_21.traces
%   Every trace is cut to whole days. Wake bouts between consecutive bridged sleep bouts
%   (SLEEPBOUTS) of at least 2L s; per fly, the active fraction in b-s bins over the first L s after
%   waking (m_start) and over the last L s before the next onset, time-reversed (m_end); flies with
%   fewer than 5 such bouts dropped. Settings L300 (L 300 s, 10-s bins; primary) and L120 (120, 5).
%   Profiles: mean over flies of each fly's mean; d = m_start - m_end, median over flies;
%   fly-bootstrap 95% bands (NBoot, percentiles type 7).
%   Figures: fig_planB21_profiles_primary (per0 DD, CS DD) and _supp (CS, per0 LD, fumin, sss).
%   Reference: hmm_audit/fig_planB21_profiles.json; if given, X.check compares n, m_start, m_end
%   and d (deterministic) with the Python's.
%
%   See also SLEEPBOUTS, READCOUNTS.

arguments
    dataRoot (1,1) string {mustBeFolder}
    refMat (1,1) string {mustBeFile}
    outDir (1,1) string
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 2026
    options.Reference (1,1) string = ""
end
if ~isfolder(outDir), mkdir(outDir); end
G = ["per0DD" "CSDD" "CS" "per0LD" "fumin" "sss"];
names = containers.Map(cellstr(G), {'per0 DD', 'CS DD', 'CS LD', 'per0 LD', 'fumin LD', 'sss LD'});
cfg = struct('name', {'L300', 'L120'}, 'L', {300, 120}, 'b', {10, 5});
R = load(refMat, 'files', 'traceLen'); rs = RandStream('twister', 'Seed', options.RandomSeed);
X = struct();
for gg = G
    files = sqwListFiles(fullfile(dataRoot, gg)); [~, fn, fe] = fileparts(files); fn = fn + fe;
    if gg == "per0DD"   % the reference order and trace lengths of the 20r replication set
        sel = strings(numel(R.files), 1); tl = double(R.traceLen(:));
        for i = 1:numel(R.files), sel(i) = files(startsWith(fn, string(R.files{i}))); end
        files = sel;
    else
        tl = Inf(numel(files), 1);
    end
    P = cell(numel(files), 2);
    parfor i = 1:numel(files)
        x = readCounts(files(i)); x = x(1:min(numel(x), tl(i))); x = x(1:floor(numel(x) / 86400) * 86400);
        P(i, :) = {flyProfiles(x, 300, 10), flyProfiles(x, 120, 5)};
    end
    for s = 1:2
        ok = ~cellfun(@isempty, P(:, s)); Q = vertcat(P{ok, s});
        Sm = vertcat(Q.ms); Em = vertcat(Q.me); Dm = Sm - Em; n = size(Sm, 1);
        bm = zeros(options.NBoot, size(Sm, 2)); be = bm; bd = bm;
        for k = 1:options.NBoot
            i = randi(rs, n, n, 1); bm(k, :) = mean(Sm(i, :), 1); be(k, :) = mean(Em(i, :), 1); bd(k, :) = median(Dm(i, :), 1);
        end
        X.(cfg(s).name).(gg) = struct('n', n, 'nBouts', [Q.n], 'ms', mean(Sm, 1), 'me', mean(Em, 1), 'd', median(Dm, 1), ...
            'ms_b', q7(bm), 'me_b', q7(be), 'd_b', q7(bd), 'excluded', nnz(~ok));
    end
end
X.check = table();
if strlength(options.Reference) > 0
    J = jsondecode(fileread(options.Reference));
    for s = 1:2
        for gg = G
            a = X.(cfg(s).name).(gg); b = J.(cfg(s).name).(gg);
            X.check = [X.check; table(string(cfg(s).name), gg, a.n, b.n, max(abs(a.ms - b.ms(:)')), max(abs(a.me - b.me(:)')), ...
                max(abs(a.d - b.d(:)')), 'VariableNames', {'setting', 'genotype', 'n', 'n_py', 'max_dms', 'max_dme', 'max_dd'})]; %#ok<AGROW>
        end
    end
    writetable(X.check, fullfile(outDir, "fig_planB21_check.csv"));
end
makeFig(X, cfg, ["per0DD" "CSDD"], names, fullfile(outDir, "fig_planB21_profiles_primary"), ...
    'DESCRIPTIVE (PRESPEC 21 / C13): 1-s activity at the two ends of wake bouts, DD. Bands: fly bootstrap 95%.');
makeFig(X, cfg, ["CS" "per0LD" "fumin" "sss"], names, fullfile(outDir, "fig_planB21_profiles_supp"), ...
    'DESCRIPTIVE, supplementary (PRESPEC 21 / C13): LD genotypes. Bands: fly bootstrap 95%.');
rows = table();
for s = 1:2
    for gg = G
        a = X.(cfg(s).name).(gg); t = ((1:numel(a.ms))' - 0.5) * cfg(s).b;
        rows = [rows; table(repmat(string(cfg(s).name), numel(t), 1), repmat(gg, numel(t), 1), repmat(a.n, numel(t), 1), t, ...
            a.ms', a.ms_b(1, :)', a.ms_b(2, :)', a.me', a.me_b(1, :)', a.me_b(2, :)', a.d', a.d_b(1, :)', a.d_b(2, :)', ...
            'VariableNames', {'setting', 'genotype', 'nFlies', 'tau_s', 'm_start', 'ms_lo', 'ms_hi', 'm_end', 'me_lo', 'me_hi', ...
            'd', 'd_lo', 'd_hi'})]; %#ok<AGROW>
    end
end
writetable(rows, fullfile(outDir, "fig_planB21_profiles.csv"));
save(fullfile(outDir, "fig_planB21_profiles.mat"), '-struct', 'X');
end

function p = flyProfiles(x, L, b)
% wake bouts between consecutive sleep bouts lasting >= 2L s; >= 5 bouts required
a = double(x > 0); sb = sleepBouts(x); st = []; en = [];
for k = 1:size(sb, 1) - 1
    w0 = sb(k, 1) - 1 + sb(k, 2); w1 = sb(k + 1, 1) - 1;      % 0-based first wake second, first sleep second
    if w1 - w0 < 2 * L, continue, end
    st(end + 1, :) = mean(reshape(a(w0 + 1:w0 + L), b, []), 1); %#ok<AGROW>
    en(end + 1, :) = mean(reshape(flipud(a(w1 - L + 1:w1)), b, []), 1); %#ok<AGROW>
end
if size(st, 1) >= 5, p = struct('ms', mean(st, 1), 'me', mean(en, 1), 'n', size(st, 1)); else, p = []; end
end

function q = q7(B)
q = zeros(2, size(B, 2));
for j = 1:size(B, 2), q(:, j) = quantile7(B(:, j), [0.025 0.975])'; end
end

function makeFig(X, cfg, groups, names, stem, ttl)
blue = [42 120 214] / 255; orange = [235 104 52] / 255;
fh = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 15 2.9 * numel(groups) + 0.8]);
try theme(fh, 'light'); catch, end
tl = tiledlayout(fh, numel(groups), 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:numel(groups)
    g = groups(r);
    for s = 1:2
        a = X.(cfg(s).name).(g); t = ((1:numel(a.ms)) - 0.5) * cfg(s).b;
        ax = nexttile(tl); hold(ax, 'on');
        patch(ax, [t fliplr(t)], [a.ms_b(1, :) fliplr(a.ms_b(2, :))], blue, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
        patch(ax, [t fliplr(t)], [a.me_b(1, :) fliplr(a.me_b(2, :))], orange, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
        h1 = plot(ax, t, a.ms, '-', 'Color', blue, 'LineWidth', 2); h2 = plot(ax, t, a.me, '-', 'Color', orange, 'LineWidth', 2);
        set(ax, 'XLim', [0 cfg(s).L * 1.18]); grid(ax, 'on');
        ylabel(ax, {sprintf('%s (n %d)', names(char(g)), a.n), 'active fraction'});
        if r == 1, title(ax, sprintf('profiles, L %d s, %d-s bins', cfg(s).L, cfg(s).b), 'FontWeight', 'normal'); end
        if r == 1 && s == 1, legend(ax, [h1 h2], {'after waking', 'before sleep onset (reversed)'}, 'Box', 'off', 'Location', 'northeast'); end
        if r == numel(groups), xlabel(ax, '\tau, seconds from the sleep boundary'); end
        ax = nexttile(tl); hold(ax, 'on');
        patch(ax, [t fliplr(t)], [a.d_b(1, :) fliplr(a.d_b(2, :))], blue, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
        plot(ax, t, a.d, '-', 'Color', blue, 'LineWidth', 2); yline(ax, 0, 'Color', [0.32 0.32 0.31]);
        set(ax, 'XLim', [0 cfg(s).L]); grid(ax, 'on'); ylabel(ax, 'd = after - before');
        if r == 1, title(ax, sprintf('median d(\\tau), L %d s', cfg(s).L), 'FontWeight', 'normal'); end
        if r == numel(groups), xlabel(ax, '\tau, seconds from the sleep boundary'); end
    end
end
title(tl, ttl, 'FontSize', 10);
exportgraphics(fh, stem + ".png", 'Resolution', 200);
exportgraphics(fh, stem + ".pdf", 'ContentType', 'vector');
close(fh);
end
