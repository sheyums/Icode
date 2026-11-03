function R = sqwBatchHMM(input, options)
%SQWBATCHHMM  Per-fly phase-switched HMMs at K = 2, 3, 4 on the S/Q/W symbols, with BIC order selection (PRESPEC 20b; fig_hmm_order.py).
%
%   R = sqwBatchHMM(folder)
%   R = sqwBatchHMM(folder, LightDark=true, Ks=[2 3 4], Starts=20, RandomSeed=20, OutDir=..., Name="hmm", Plot=true)
%
%   input      a folder of symbol files written by SQWBATCHSYMBOLS (*_symbols.csv: lab coding
%              S = 0, Q = 1, W = 2, NaN = missing, with a phase column), or a folder / list of 1-s
%              count files (scored here by SQWSYMBOLS, theta 0.5, as SQWBATCHSYMBOLS)
%   LightDark  count-file input only: true = LD (light / dark matrices), false = DD (one)
%   Ks         numbers of hidden states to fit (default 2, 3, 4)
%   Starts, RandomSeed, MaxIter, Tol   passed to HMMFITPHASED (20 random starts, seed 20, as 20b)
%   OutDir     default: the input folder for symbol files, "<folder>_hmm" beside a count folder
%   Resume     true (default): each fly's fits are saved in OutDir/<Name>_fits/ as they finish and
%              reused on the next call, so an interrupted batch continues where it stopped
%
%   Per fly and K: the maximum-likelihood fit, its log-likelihood, and
%     BIC = -2 ll + k ln n,  k = (K-1) + nP K (K-1) + 2K   (HMMBIC)
%   with n = minutes (primary) and n = number of symbol changes (sensitivity). bestK_T / bestK_chg
%   = the K with the lowest BIC; flag_T / flag_chg = some K ~= 3 beats K = 3 by more than 10 (the
%   PRESPEC 20 flag rule). The K = 3 fit is labelled H1, H2, H3 (LABELSTATES: H3 = most S, H2 = more
%   W, H1 = rest) and deg3 = H1 not identified (ISDEGENERATEH1). For K = 2 and 4, dom_K2 / dom_K4
%   give the symbol each hidden state emits most, e.g. "S+S+Q+W" (two sleep-like states), and
%   nS_K4 / nQ_K4 / nW_K4 count them. nNear_K = starts ending within 0.5 nats of the best (a
%   value of 1 means the optimum was found once; it is not evidence of identification).
%
%   Writes <Name>_perfly.csv (one row per fly, names not codes: B_H1_Q = P(Q | H1),
%   P_light_H1_H2 = P(H1 -> H2) in the light phase, ...), <Name>.mat (R, including every raw fit)
%   and <Name>_order.png / .pdf (K = 3 emissions, BIC selection, dominant-symbol patterns).
%   R  struct: label (the input folder's name, used in figures), perFly (table), fits (nFly x numel(Ks) cell of HMMFITPHASED structs), K3 (labelled
%      fits), Ks, options, outDir, fig.
%   The fits use MATLAB's random streams, so a fly can land on a different (higher or lower)
%   optimum than the Python's 20b; compare likelihoods, not parameters.
%
%   Example (from the Icode\sqw folder):
%       X = sqwBatchSymbols("D:\flies\CS");     % symbol files first (or pass the count folder directly)
%       R = sqwBatchHMM(X.outDir);              % writes hmm_perfly.csv, hmm.mat, hmm_order.png there
%       mean(R.perFly.bestK_T == 3)             % share of flies where BIC picks K = 3
%
%   See also SQWBATCHSYMBOLS, HMMFITPHASED, HMMBIC, LABELSTATES, PLOTSQWHMMORDER.

arguments
    input
    options.LightDark (1,1) logical = true
    options.Ks (1,:) double {mustBePositive, mustBeInteger} = [2 3 4]
    options.Starts (1,1) double {mustBePositive, mustBeInteger} = 20
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
    options.MaxIter (1,1) double {mustBePositive, mustBeInteger} = 1000
    options.Tol (1,1) double {mustBePositive} = 1e-7
    options.OutDir (1,1) string = ""
    options.Name (1,1) string = "hmm"
    options.Resume (1,1) logical = true
    options.Plot (1,1) logical = true
    options.MaxFlies (1,1) double = Inf
end
[src, files, outDir] = resolveInput(input, options.OutDir);
files = files(1:min(numel(files), options.MaxFlies));
if ~isfolder(outDir), mkdir(outDir); end
fitDir = fullfile(outDir, options.Name + "_fits"); if ~isfolder(fitDir), mkdir(fitDir); end
o = options; Ks = o.Ks; nF = numel(files); fits = cell(nF, numel(Ks)); meta = cell(nF, 1); t0 = tic;
parfor i = 1:nF
    [y, ph, nP, name] = readFly(files(i), src, o.LightDark);
    cache = fullfile(fitDir, name + ".mat"); F = [];
    if o.Resume && isfile(cache)
        C = load(cache); if isequal(C.Ks, Ks) && C.Starts == o.Starts && C.RandomSeed == o.RandomSeed && C.T == numel(y), F = C.F; end
    end
    if isempty(F)
        F = cell(1, numel(Ks));
        for k = 1:numel(Ks)
            F{k} = hmmFitPhased(y, ph, nP, Ks(k), Starts=o.Starts, RandomSeed=o.RandomSeed, MaxIter=o.MaxIter, Tol=o.Tol);
        end
        parsave(cache, F, Ks, o.Starts, o.RandomSeed, numel(y));
    end
    fits(i, :) = F;
    meta{i} = struct('name', name, 'y', y, 'nP', nP);
end
fprintf('sqwBatchHMM: %d flies x K = %s, %.0f s\n', nF, mat2str(Ks), toc(t0));
[R.perFly, R.K3] = tabulate(fits, meta, Ks);
R.label = batchLabel(input, options.Name); R.fits = fits; R.Ks = Ks; R.options = options; R.outDir = outDir; R.files = files;
writetable(R.perFly, fullfile(outDir, options.Name + "_perfly.csv"));
if options.Plot
    R.fig = plotSqwHMMOrder(R, options.Name);
    exportgraphics(R.fig, fullfile(outDir, options.Name + "_order.png"), 'Resolution', 300);
    exportgraphics(R.fig, fullfile(outDir, options.Name + "_order.pdf"), 'ContentType', 'vector');
end
Rs = R; if isfield(Rs, 'fig'), Rs = rmfield(Rs, 'fig'); end
save(fullfile(outDir, options.Name + ".mat"), '-struct', 'Rs');
end

function s = batchLabel(input, name)
% the input folder's name without "_symbols" (the batch name if a file list was given)
s = string(input);
if isscalar(s) && isfolder(s), [~, s] = fileparts(strip(s, "right", filesep)); s = erase(s, "_symbols"); else, s = name; end
end

function [src, files, outDir] = resolveInput(input, outDir)
s = string(input);
if isscalar(s) && isfolder(s)
    d = dir(fullfile(s, '*_symbols.csv'));
    if ~isempty(d)
        src = "symbols"; files = sort(string(fullfile({d.folder}, {d.name}))');
        if strlength(outDir) == 0, outDir = s; end
        return
    end
end
src = "counts"; files = sqwListFiles(input);
if strlength(outDir) == 0
    if ~(isscalar(s) && isfolder(s)), error('sqw:outDir', 'Give OutDir when the input is a file list.'); end
    in = char(s); if any(in(end) == '/\'), in = in(1:end - 1); end
    outDir = string(in) + "_hmm";
end
end

function [y, ph, nP, name] = readFly(file, src, lightDark)
[~, nm, ex] = fileparts(file);
if src == "symbols"
    t = readtable(file);
    if ~all(ismember({'symbol', 'phase'}, t.Properties.VariableNames)), error('sqw:symbolFile', '%s has no symbol / phase column.', file); end
    y = labCodeToSqw(t.symbol); ph = t.phase; name = erase(string(nm), "_symbols");
    if ~all(ismember(ph, [1 2])), error('sqw:symbolFile', '%s: phase must be 1 or 2.', file); end
    nP = max(ph);
else
    [y, ph, nP] = sqwFlyInput(file, lightDark); name = erase(string(nm) + ex, ".txt");
end
end

function parsave(f, F, Ks, Starts, RandomSeed, T) %#ok<INUSD>
save(f, 'F', 'Ks', 'Starts', 'RandomSeed', 'T');
end

function [T, K3] = tabulate(fits, meta, Ks)
sym = ["Q" "W" "S"]; nF = size(fits, 1); rows = cell(nF, 1); K3 = cell(nF, 1);
for i = 1:nF
    m = meta{i}; n = numel(m.y); nP = m.nP; r = table(m.name, n, nP, 'VariableNames', {'file', 'T', 'nP'});
    bT = NaN(1, numel(Ks)); bC = bT;
    for k = 1:numel(Ks)
        f = fits{i, k}; K = Ks(k);
        [bT(k), ~, nE] = hmmBic(f.ll, K, nP, n, m.y); bC(k) = hmmBic(f.ll, K, nP, [], m.y);
        r.("ll_K" + K) = f.ll; r.("BIC_T_K" + K) = bT(k); r.("BIC_chg_K" + K) = bC(k);
        r.("converged_K" + K) = f.converged; r.("nNear_K" + K) = nnz(f.llStarts >= f.ll - 0.5);
    end
    r.nchg = nE;
    [~, a] = min(bT); r.bestK_T = Ks(a); [~, a] = min(bC); r.bestK_chg = Ks(a);
    k3 = find(Ks == 3);
    if isempty(k3)
        r.flag_T = false; r.flag_chg = false;
    else
        r.flag_T = any(bT([1:k3 - 1, k3 + 1:end]) < bT(k3) - 10); r.flag_chg = any(bC([1:k3 - 1, k3 + 1:end]) < bC(k3) - 10);
        for K = [2 4]
            if any(Ks == K), r.("dBIC_T_K" + K + "m3") = bT(Ks == K) - bT(k3); end
        end
        f = fits{i, k3}; [P, B, p0] = labelStates(f.P, f.B, f.pi);
        K3{i} = struct('pi', p0, 'P', P, 'B', B); r.deg3 = isDegenerateH1(B);
        for h = 1:3, for v = 1:3, r.("B_H" + h + "_" + sym(v)) = B(h, v); end, end
        if nP == 2, pn = ["light" "dark"]; else, pn = "whole"; end
        for p = 1:nP, for a2 = 1:3, for b2 = 1:3, r.("P_" + pn(p) + "_H" + a2 + "_H" + b2) = P(a2, b2, p); end, end, end
    end
    for K = [2 4]
        if ~any(Ks == K), continue, end
        [~, d] = max(fits{i, Ks == K}.B, [], 2); d = sort(mod(d, 3));     % S (3 -> 0) first, then Q, W
        SQW = sym([3 1 2]);
        r.("dom_K" + K) = strjoin(SQW(d' + 1), "+");
        if K == 4, r.nS_K4 = nnz(d == 0); r.nQ_K4 = nnz(d == 1); r.nW_K4 = nnz(d == 2); end
    end
    rows{i} = r;
end
T = vertcat(rows{:});
end
