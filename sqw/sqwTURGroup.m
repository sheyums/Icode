function G = sqwTURGroup(files, genotype, options)
%SQWTURGROUP  Thermodynamic uncertainty relation (PRESPEC 20q) for a group of flies: per fly-phase table, phase summary, CSV, plot.
%
%   G = sqwTURGroup(files, genotype)
%   G = sqwTURGroup(files, genotype, LightDark=true, HMMFits=[], NNullMC=200, NNullHMM=100, ...
%                   Taus=[180 120 360], RuleOfThree=true, RandomSeed=16000, OutDir="", Plot=true)
%
%   files, genotype, LightDark, HMMFits as SQWMEMORYGROUP. Taus(1) is the primary
%   block length (3 h); the others are sensitivity.
%
%   G  struct:
%     perFly   one row per fly x phase: genotype, file, fly, phase, hasRef, deg,
%              SQWTURFLY's exact columns (sig_mc, B_mc, B_mc_ct, kappa, sig_hmm,
%              B_hmm, kind, r3, ...), J (counted, /h), and per tau (suffix _<tau>):
%              D, n (blocks), B_obs, medB (median B of the chain records), nSim
%              (finite chain records), p1 = (1 + #{B* >= B_obs}) / (1 + nSim),
%              B_c = B_obs B_mc / medB / kappa (bias-corrected, deflated by the
%              discrete-time inflation), inH (B_obs inside the HMM records' central 95%;
%              NaN if degenerate or < 50 finite HMM records), use (in the summary)
%     summary  one row per phase over flies with a Q/W reference, excluding fly-phases
%              without a usable cycle (amendment D: kind ~= "cycle" after the rule of
%              three) and those with no estimate (B_obs not finite or < 100 finite chain
%              records), each counted (nNoCycle_* by kind, nNoEstimate):
%              S1 median B_obs/medB, Wilcoxon p of its log, fraction p1 < 0.05 and the
%              binomial p; S2 median B_c/sig_mc, Wilcoxon p of its log, fraction > 1;
%              HMM_enough (sig_hmm >= B_c, non-degenerate; Inf counts as >=), HMM_inside;
%              Q_mc = B_mc/sig_mc, Q_obs = B_c/sig_mc, Q_hmm = B_c/sig_hmm (0 where sig_hmm
%              is infinite or 0) medians; kappa max; medians of J, D, blocks, B_obs,
%              sig_mc; S1 and S2 medians at the other taus. Holm (primary DD blocks
%              only, 20q) is the caller's (HOLMADJUST).
%     flies    the SQWTURFLY structs
%   Writes <genotype>_tur_perfly.csv and _tur_summary.csv to OutDir.
%
%   See also SQWTURFLY, SQWMEMORYGROUP, HOLMADJUST.

arguments
    files string
    genotype (1,1) string
    options.LightDark (1,1) logical = true
    options.HMMFits = []
    options.NNullMC (1,1) double {mustBeNonnegative, mustBeInteger} = 200
    options.NNullHMM (1,1) double {mustBeNonnegative, mustBeInteger} = 100
    options.Taus (1,:) double {mustBePositive, mustBeInteger} = [180 120 360]
    options.RuleOfThree (1,1) logical = true
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 16000
    options.OutDir (1,1) string = ""
    options.Plot (1,1) logical = true
end
files = sqwListFiles(files);
n = numel(files); out = cell(n, 1); fits = options.HMMFits; o = options;
if ~isempty(fits) && numel(fits) ~= n, error('sqw:fits', 'HMMFits must have one entry per file (%d).', n); end
parfor k = 1:n
    [y, ph, nP, refs, names] = sqwFlyInput(files(k), o.LightDark);
    if isempty(fits), h = []; deg = NaN; else, h = fits(k); deg = h.deg; end
    F = sqwTURFly(y, ph, nP, HMM=h, NNullMC=o.NNullMC, NNullHMM=o.NNullHMM, Taus=o.Taus, ...
        RandomSeed=o.RandomSeed + k, RuleOfThree=o.RuleOfThree);
    F.rows = flyRows(F, genotype, fileName(files(k)), k, names, refs, deg);
    out{k} = F;
end
G.flies = vertcat(out{:});
G.perFly = vertcat(G.flies.rows);
G.summary = summarise(G.perFly, options.Taus);
if strlength(options.OutDir) > 0
    if ~isfolder(options.OutDir), mkdir(options.OutDir); end
    writetable(G.perFly, fullfile(options.OutDir, genotype + "_tur_perfly.csv"));
    writetable(G.summary, fullfile(options.OutDir, genotype + "_tur_summary.csv"));
end
if options.Plot, plotTUR(G, genotype, options.Taus(1)); end
end

function T = flyRows(F, genotype, file, k, names, refs, deg)
nP = numel(names);
T = [table(repmat(genotype, nP, 1), repmat(file, nP, 1), repmat(k, nP, 1), names, ...
    arrayfun(@(nm) isfinite(refs.(nm)), names), repmat(deg, nP, 1), [F.obs.J]', ...
    'VariableNames', {'genotype', 'file', 'fly', 'phase', 'hasRef', 'deg', 'J'}), F.ex];
for q = 1:numel(F.Taus)
    s = "_" + F.Taus(q);
    for p = 1:nP
        Bo = F.obs(p).B(q); Bs = F.simMC(:, p, q); Bs = Bs(isfinite(Bs));
        Bh = F.simH(:, p, q); Bh = Bh(isfinite(Bh)); e = F.ex(p, :);
        T.("D" + s)(p) = F.obs(p).D(q); T.("n" + s)(p) = F.obs(p).n(q); T.("B_obs" + s)(p) = Bo;
        T.("nSim" + s)(p) = numel(Bs);
        if isempty(Bs), md = NaN; else, md = median(Bs); end
        T.("medB" + s)(p) = md;
        T.("p1" + s)(p) = (1 + sum(Bs >= Bo)) / (1 + numel(Bs));
        T.("B_c" + s)(p) = Bo * e.B_mc / md / e.kappa;
        if numel(Bh) >= 50 && ~(deg == 1)
            qh = quantile7(Bh, [0.025 0.975]); T.("inH" + s)(p) = double(qh(1) <= Bo && Bo <= qh(2));
        else
            T.("inH" + s)(p) = NaN;
        end
        cyc = isfinite(e.sig_mc) && e.sig_mc > 1e-9 && isfinite(e.B_mc_ct) && e.B_mc_ct > 1e-12;   % amendment D, after E
        T.cycle(p) = cyc;
        T.("use" + s)(p) = T.hasRef(p) && cyc && isfinite(Bo) && numel(Bs) >= 100;
    end
end
end

function S = summarise(T, taus)
S = table(); t0 = "_" + taus(1);
for ph = unique(T.phase, 'stable')'
    A = T(T.phase == ph & T.hasRef, :); r = table(); r.phase = ph;
    r.nNoCycle_i = sum(A.kind == "i" & ~A.cycle); r.nNoCycle_ii = sum(A.kind == "ii"); r.nNoCycle_iii = sum(A.kind == "iii");
    r.nRuleOfThree = sum(A.r3);
    r.nNoEstimate = sum(A.cycle & ~A.("use" + t0));
    D = A(A.("use" + t0), :); r.n = height(D);
    a = D.("B_obs" + t0) ./ D.("medB" + t0); b = D.("B_c" + t0) ./ D.sig_mc;
    r.S1_med = median(a); r.S1_p = wilcoxonP(log(a));
    r.S1_frac = mean(D.("p1" + t0) < 0.05); r.S1_bin = binomGreater(sum(D.("p1" + t0) < 0.05), height(D), 0.05);
    r.S2_med = median(b); r.S2_p = wilcoxonP(log(b)); r.S2_frac = mean(log(b) > 0);
    nd = D(D.deg ~= 1, :); bc = nd.("B_c" + t0); sh = nd.sig_hmm;
    r.HMM_enough = mean(sh >= bc); r.HMM_inside = mean(nd.("inH" + t0), 'omitnan');
    r.Q_mc = median(D.B_mc ./ D.sig_mc); r.Q_obs = median(b);
    qh = zeros(height(nd), 1); okh = isfinite(sh) & sh > 0; qh(okh) = bc(okh) ./ sh(okh); r.Q_hmm = median(qh);
    r.kappa_max = max(D.kappa); r.J = median(D.J); r.D = median(D.("D" + t0)); r.blocks = fix(median(D.("n" + t0)));
    r.B_obs = median(D.("B_obs" + t0)); r.sig_mc = median(D.sig_mc);
    for q = 2:numel(taus)
        s = "_" + taus(q); Dq = A(A.("use" + s), :);
        r.("S1" + s) = median(Dq.("B_obs" + s) ./ Dq.("medB" + s));
        r.("S2" + s) = median(Dq.("B_c" + s) ./ Dq.sig_mc);
    end
    S = [S; r]; %#ok<AGROW>
end
end

function plotTUR(G, genotype, tau)
s = "_" + tau; T = G.perFly(G.perFly.("use" + s), :); phs = unique(T.phase, 'stable');
fh = figure('Name', genotype + " TUR (20q)", 'Color', 'w');
try theme(fh, 'light'); catch, end
tl = tiledlayout(fh, 1, numel(phs), 'TileSpacing', 'compact');
for k = 1:numel(phs)
    D = T(T.phase == phs(k), :); ax = nexttile(tl); hold(ax, 'on');
    x = D.("B_obs" + s) ./ D.("medB" + s); yv = D.("B_c" + s) ./ D.sig_mc;
    scatter(ax, x, yv, 30, [0.2 0.4 0.7], 'filled');
    set(ax, 'XScale', 'log', 'YScale', 'log');
    xline(ax, 1, 'k-'); yline(ax, 1, 'k-');
    xlabel(ax, 'S1: B_{obs} / median B (chain records)'); ylabel(ax, 'S2: B_c / \sigma_{MC}');
    title(ax, sprintf('%s %s (n %d): S1 %.2f, S2 %.2f', genotype, phs(k), height(D), median(x), median(yv)));
    box(ax, 'on');
end
title(tl, sprintf('Each point a fly, tau = %d min. Above 1 on S2: the TUR would demand dissipation the observed chain lacks', tau));
end

function s = fileName(f)
[~, a, b] = fileparts(f); s = a + b;
end
