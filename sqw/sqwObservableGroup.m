function G = sqwObservableGroup(files, genotype, options)
%SQWOBSERVABLEGROUP  Plan A observable table for one genotype (PRESPEC 20s-A, A1 and A3; sqw17 + sqw17_report).
%
%   G = sqwObservableGroup(files, genotype)
%   G = sqwObservableGroup(files, genotype, LightDark=true, Theta=0.5, SixHour=false, ...
%                          NBoot=2000, RandomSeed=20, OutDir="", Plot=false)
%
%   files      string array of count files, one per fly (or a folder: every
%              *.txt in it or below, sorted by full path, as sqw.files)
%   genotype   label written into every row ("CS", "per0DD", ...)
%   LightDark  true for LD (light / dark phases), false for DD (one phase "whole")
%   SixHour    also the four 6-h ZT blocks (C3; LD only)
%   OutDir     if given, CSV files for Origin are written there:
%                <genotype>_observable_perfly.csv   one row per fly x phase
%                <genotype>_observable_summary.csv  one row per phase
%                <genotype>_observable_sixhour_*.csv (SixHour)
%   Plot       draw PLOTSQWGROUP
%   RuleOfThree  passed to SQWOBSERVABLEFLY (default true, convention of 2026-11-02)
%
%   G   struct:
%     genotype, files, flies (SQWOBSERVABLEFLY outputs)
%     perFly    table, one row per fly x phase, every quantity plus refOK, hasCycle
%     summary   table, one row per phase (rules below)
%     sixPerFly, sixSummary   (SixHour) the same for the 6-h blocks
%
%   Summary rules, those of sqw17_report.py (PRESPEC 20s-A):
%     n          fly-phases WITH a Q/W reference (the others are excluded)
%     nCycle     of those, chains with a cycle (|J_chain| >= 1e-12)
%     c_obs, J, EPR (count), E21, E31: median and fly-bootstrap 95% CI
%                (MEDIANCI: non-finite values dropped, >= 3 needed)
%     *_fracPos  share of the n flies with value > 0 (NaN counts as not > 0)
%     chain J / A / EPR medians over the nCycle flies, non-finite dropped;
%     EPR_chain_inf = share with infinite chain EPR, and EPR_chain_med_withInf
%     keeps those flies (Inf ranks highest) -- the ledger's median drops them
%     dXn, dYn   median (counted - predicted) routing, rhoXn, rhoYn Spearman,
%                over cycle flies with both values finite
%   One RandStream (RandomSeed) runs through the whole summary.
%
%   See also SQWOBSERVABLEFLY, SQWONSAGERGROUP, PLOTSQWGROUP.

arguments
    files
    genotype (1,1) string
    options.LightDark (1,1) logical = true
    options.Theta (1,1) double {mustBePositive} = 0.5
    options.UnknownWake (1,1) string {mustBeMember(options.UnknownWake, ["missing" "unknown"])} = "missing"
    options.SixHour (1,1) logical = false
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
    options.OutDir (1,1) string = ""
    options.Plot (1,1) logical = false
    options.RuleOfThree (1,1) logical = true
end
files = string(files);
if isscalar(files) && isfolder(files)
    d = dir(fullfile(files, '**', '*.txt'));        % recursive, as sqw.files
    files = sort(string(fullfile({d.folder}, {d.name})))';
end
files = files(:);
if isempty(files), error('sqwObservableGroup:noFiles', 'No count files given.'); end

G.genotype = genotype; G.files = files; G.options = options;
nF = numel(files);
flies = cell(nF, 1);
for i = 1:nF
    flies{i} = sqwObservableFly(files(i), LightDark=options.LightDark, Theta=options.Theta, ...
        UnknownWake=options.UnknownWake, SixHour=options.SixHour, RuleOfThree=options.RuleOfThree);
end
G.flies = vertcat(flies{:});

rs = RandStream('twister', 'Seed', options.RandomSeed);
G.perFly = flyTable(G.flies, genotype, 'Q', 'phaseNames');
G.summary = summarise(G.perFly, options.NBoot, rs);
G.sixPerFly = table(); G.sixSummary = table();        % same fields in every group output
if options.SixHour
    G.sixPerFly = flyTable(G.flies, genotype, 'six', '');
    G.sixSummary = summarise(G.sixPerFly, options.NBoot, rs);
end
if strlength(options.OutDir) > 0
    if ~isfolder(options.OutDir), mkdir(options.OutDir); end
    writetable(dropMatrices(G.perFly), fullfile(options.OutDir, genotype + "_observable_perfly.csv"));
    writetable(G.summary, fullfile(options.OutDir, genotype + "_observable_summary.csv"));
    if options.SixHour
        writetable(dropMatrices(G.sixPerFly), fullfile(options.OutDir, genotype + "_observable_sixhour_perfly.csv"));
        writetable(G.sixSummary, fullfile(options.OutDir, genotype + "_observable_sixhour_summary.csv"));
    end
end
if options.Plot
    plotSqwGroup(G);
end
end

% ---------------------------------------------------------------------------
function T = flyTable(flies, genotype, which, ~)
rows = cell(numel(flies), 1);
for i = 1:numel(flies)
    f = flies(i);
    if which == "six", Q = f.six.Q; names = f.six.phaseNames; else, Q = f.Q; names = f.phaseNames; end
    [~, nm, ext] = fileparts(f.file);
    t = struct2table(rmfield(Q, {'N', 'routes'}), 'AsArray', true);
    t.N = reshape(permute(cat(3, Q.N), [3 1 2]), numel(Q), 9);       % N(i,j) in column 3(j-1)+i
    t.routes = vertcat(Q.routes);
    t = [table(repmat(genotype, numel(Q), 1), repmat(i, numel(Q), 1), repmat(string(nm) + ext, numel(Q), 1), ...
        (1:numel(Q))', names(:), 'VariableNames', {'genotype', 'fly', 'file', 'phase', 'phaseName'}), t]; %#ok<AGROW>
    rows{i} = t;
end
T = vertcat(rows{:});
end

function T = dropMatrices(T)
% CSV for Origin: matrices split into named columns.
N = T.N; R = T.routes; T = removevars(T, {'N', 'routes'});
lab = 'QWS';
for j = 1:3
    for i = 1:3
        if i ~= j, T.("N_" + lab(i) + lab(j)) = N(:, 3 * (j - 1) + i); end
    end
end
nm = ["N312" "N32" "N213" "N23"];
for k = 1:4, T.(nm(k)) = R(:, k); end
end

function S = summarise(T, nBoot, rs)
phases = unique(T.phase)';
S = table();
for p = phases
    R = T(T.phase == p & T.refOK, :);
    C = R(R.hasCycle, :);
    s = table();
    s.phase = p; s.phaseName = string(T.phaseName(find(T.phase == p, 1)));
    s.n = height(R); s.nCycle = height(C);
    for k = ["c_obs" "J" "EPR"]
        [s.(k + "_med"), s.(k + "_lo"), s.(k + "_hi")] = medianCI(R.(k), NBoot=nBoot, Stream=rs);
        s.(k + "_fracPos") = mean(R.(k) > 0);
    end
    s.rule3_frac = mean(R.rule_of_three);
    s.J_chain_med = median(C.J_chain(isfinite(C.J_chain)));
    s.A_chain_med = median(C.A_chain(isfinite(C.A_chain)));
    s.EPR_chain_med = median(C.EPR_chain(isfinite(C.EPR_chain)));
    s.EPR_chain_inf = mean(~isfinite(C.EPR_chain));
    s.EPR_chain_med_withInf = median(C.EPR_chain);
    for k = ["E21" "E31"]
        [s.(k + "_med"), s.(k + "_lo"), s.(k + "_hi")] = medianCI(C.(k), NBoot=nBoot, Stream=rs);
    end
    okX = isfinite(C.Xn_cnt) & isfinite(C.F_SW); okY = isfinite(C.Yn_cnt) & isfinite(C.F_WS);
    s.dXn_med = median(C.Xn_cnt(okX) - C.F_SW(okX)); s.rhoXn = rho5(C.F_SW(okX), C.Xn_cnt(okX)); s.nXn = nnz(okX);
    s.dYn_med = median(C.Yn_cnt(okY) - C.F_WS(okY)); s.rhoYn = rho5(C.F_WS(okY), C.Yn_cnt(okY)); s.nYn = nnz(okY);
    S = [S; s]; %#ok<AGROW>
end
end

function r = rho5(a, b)
% The report computes Spearman only with more than 4 pairs.
if numel(a) > 4, r = spearmanRho(a, b); else, r = NaN; end
end
