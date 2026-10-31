function Bres = run_20rB(stageDir, boutsMat, label, workDir, nRep)
% RUN_20RB  20r-B renewal control (PRESPEC 2026-10-31 20r-B, with the extractor of
% "20r amendment: stored per0 DD bouts lost"). Run BY NAME from the Icode folder.
%
%   stageDir : folder holding ONLY this genotype's cut count files (the 20r input)
%   boutsMat : the observed bouts saved by extract_20r for those files
%   label    : 'per0 DD' or 'CS DD'
%   workDir  : empty scratch folder for the synthetic count files
%   nRep     : number of synthetic pooled datasets (199 as pre-specified)
%
% Per fly, from the retained trace (Teliminate 0, so the first traceLen samples
% of the cut file): maximal runs of count > 0 (active) and count == 0
% (inactive, all lengths). Synthetic trace: start in the state of the fly's
% first second; alternate active / inactive runs drawn iid with replacement
% from that fly's own runs; stop at the record length, truncating the last run.
% The 25 synthetic traces are written as count files (1 for active seconds, 0
% for inactive: only count > 0 versus 0 enters sleeptrace) and scored by
% extract_20r -- the identical pipeline. The pooled synthetic wake bouts are
% subsampled without replacement to the observed n if larger. Statistic:
% BoutHazard RiseRatio at NumBins {12,16,20,24}, Spacing "log", MinAtRisk 10.
% Replicate b uses rng(2000 + b, 'twister') for every draw it makes.
% p_ren = (1 + #{RiseRatio* >= observed}) / (nRep + 1).

assert(isequal(fileparts(which('BoutHazard')), pwd), ...
    'Run from the Icode folder: BoutHazard resolves to %s', which('BoutHazard'));
here = fileparts(mfilename('fullpath'));
nb = [12 16 20 24];

O = load(boutsMat);
xObs = O.Wake(isfinite(O.Wake) & O.Wake >= O.xminWake);
nObs = numel(xObs);
rrObs = arrayfun(@(k) BoutHazard(xObs, O.xminWake, SamplingInterval=1, NumBins=k, ...
    Spacing="log", MinAtRisk=10).RiseRatio, nb);

% ---------------------------------------------------- per-fly run lengths
nFly = numel(O.files);
runs = struct('first', {}, 'act', {}, 'ina', {}, 'len', {});
for i = 1:nFly
    c = readmatrix(fullfile(stageDir, O.files{i})); c = c(:);
    assert(numel(c) >= O.traceLen(i), 'fly %d: file shorter than its retained trace', i);
    c = c(1:O.traceLen(i));                  % the retained trace (Teliminate 0)
    a = c > 0;
    d = diff([~a(1); a]);  st = find(d ~= 0);
    len = diff([st; numel(a) + 1]);
    state = a(st);
    runs(i).first = a(1);
    runs(i).act = len(state);
    runs(i).ina = len(~state);
    runs(i).len = numel(a);
end

% ------------------------------------------------ pre-run histogram check
rng(2000, 'twister');
syn = synth(runs(1));
[sa, si] = runLengths(syn);
fprintf('%s pre-run check, fly 1 (%s):\n', label, O.files{1});
q = [0.1 0.25 0.5 0.75 0.9 0.99];
fprintf('  active runs   real n %5d, quantiles %s; synthetic n %5d, quantiles %s\n', ...
    numel(runs(1).act), mat2str(quantile(runs(1).act, q)), numel(sa), mat2str(quantile(sa, q)));
fprintf('  inactive runs real n %5d, quantiles %s; synthetic n %5d, quantiles %s\n', ...
    numel(runs(1).ina), mat2str(quantile(runs(1).ina, q)), numel(si), mat2str(quantile(si, q)));
if nRep == 0, Bres = runs; return, end      % check-only mode: read this before the run

% ---------------------------------------------------------- replicates
rrRen = nan(nRep, 4); nSyn = nan(nRep, 1);
if ~isfolder(workDir), mkdir(workDir); end
for b = 1:nRep
    rng(2000 + b, 'twister');
    delete(fullfile(workDir, '*.txt'));
    for i = 1:nFly
        s = synth(runs(i));
        fid = fopen(fullfile(workDir, sprintf('syn_%02d.txt', i)), 'w');
        fprintf(fid, '%d\n', s); fclose(fid);
    end
    Bs = evalc_extract(workDir, fullfile(workDir, 'syn.mat'));
    xs = Bs.Wake(isfinite(Bs.Wake) & Bs.Wake >= Bs.xminWake);
    nSyn(b) = numel(xs);
    if nSyn(b) > nObs
        xs = xs(randperm(nSyn(b), nObs));
    end
    rrRen(b, :) = arrayfun(@(k) BoutHazard(xs, Bs.xminWake, SamplingInterval=1, ...
        NumBins=k, Spacing="log", MinAtRisk=10).RiseRatio, nb);
    if mod(b, 20) == 0, fprintf('  replicate %d / %d\n', b, nRep); end
end

pRen = (1 + sum(rrRen >= rrObs, 1)) / (nRep + 1);
nBeyond = sum(pRen < 0.05);
if nBeyond >= 3
    verdict = 'rise BEYOND renewal';
elseif sum(pRen >= 0.05) >= 3
    verdict = 'EXPLAINED by renewal and scoring';
else
    verdict = 'inconclusive (bin-sensitive)';
end

fprintf('%s 20r-B: observed n %d; synthetic pooled n median %d (min %d, max %d); %d of %d subsampled\n', ...
    label, nObs, median(nSyn), min(nSyn), max(nSyn), sum(nSyn > nObs), nRep);
for k = 1:4
    fprintf(['  bins %2d: observed RiseRatio %.3f; renewal median %.3f, 95th %.3f; ' ...
        'p_ren %.3f\n'], nb(k), rrObs(k), median(rrRen(:, k)), prctile(rrRen(:, k), 95), pRen(k));
end
fprintf('  verdict: %s\n', verdict);

Bres = struct('label', label, 'NumBins', nb, 'nObs', nObs, 'rrObs', rrObs, 'rrRen', rrRen, ...
    'nSyn', nSyn, 'pRen', pRen, 'verdict', verdict);
save(fullfile(here, sprintf('run_20rB_%s.mat', strrep(label, ' ', ''))), 'Bres', '-v7');
end

function s = synth(r)
% Alternating renewal from the fly's own runs, starting in its first state.
s = zeros(r.len, 1); pos = 0; st = r.first;
while pos < r.len
    if st
        L = r.act(randi(numel(r.act)));
    else
        L = r.ina(randi(numel(r.ina)));
    end
    L = min(L, r.len - pos);
    s(pos + 1:pos + L) = st;
    pos = pos + L; st = ~st;
end
end

function [act, ina] = runLengths(s)
a = s > 0;
d = diff([~a(1); a]); st = find(d ~= 0);
len = diff([st; numel(a) + 1]);
act = len(a(st)); ina = len(~a(st));
end

function B = evalc_extract(dirIn, outMat)
% extract_20r prints a report per call; keep the replicate loop's log readable.
B = []; %#ok<NASGU>
evalc('B = extract_20r(dirIn, outMat);');
end
