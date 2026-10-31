function C = run_20rC(boutsMat, eFile, genotype)
% RUN_20RC  20r-C, EXPLORATORY (PRESPEC 2026-10-31 20r-C, with step (4) of
% "20r amendment: stored per0 DD bouts lost"). Run BY NAME from the Icode folder.
%
%   boutsMat : bouts saved by extract_20r (the 20r reference set for the genotype)
%   eFile    : hmm_audit/e_perfly_20p.csv
%   genotype : the csv's genotype label, e.g. 'per0DD'
%
% Pooled fit: FitHyperErlangMLE with Shapes [1 1 5] (the FINDINGS shapes),
% discrete, SamplingInterval 1, RandomSeed 1. For a wake bout of k bins the
% responsibility of branch j is w_j [S_j(k-1) - S_j(k)] / sum_i w_i [S_i(k-1) -
% S_i(k)] -- the engine's convention that a bout of k bins lies in (k-1, k] --
% with Theta's untruncated weights w and Erlang survivals S_j; the truncation
% normalizer cancels. Per fly: the mean responsibility of the m = 5 branch over
% its wake bouts. Spearman correlation with E_nats_per_h and with
% E_wakeonly_nats_per_h, matched by file name, flies with qw_reference False
% excluded. No decision rule: report rho, p, n.

assert(isequal(fileparts(which('FitHyperErlangMLE')), pwd), ...
    'Run from the Icode folder: FitHyperErlangMLE resolves to %s', which('FitHyperErlangMLE'));

B = load(boutsMat);
keep = isfinite(B.Wake) & B.Wake >= B.xminWake;
x = B.Wake(keep);
H = FitHyperErlangMLE(x, B.xminWake, SamplingInterval=1, Shapes=[1 1 5], RandomSeed=1);
mm = H.Shapes; J = numel(mm);
th = H.Theta(:)'; rates = th(1:J); w = th(J+1:2*J);
series = find(mm > 1);
assert(isscalar(series), 'expected exactly one series branch');

Sj = @(t, j) erlangSF(t, rates(j), mm(j));
nFly = size(B.Wake, 1);
resp = nan(nFly, 1); nb = zeros(nFly, 1);
for i = 1:nFly
    k = B.Wake(i, keep(i, :)); k = k(:);
    p = zeros(numel(k), J);
    for j = 1:J
        p(:, j) = w(j) * (Sj(k - 1, j) - Sj(k, j));
    end
    r = p(:, series) ./ sum(p, 2);
    resp(i) = mean(r); nb(i) = numel(k);
end

T = readtable(eFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
T.Properties.VariableNames{1} = 'genotype';   % the header carries a byte-order mark
T = T(T.genotype == string(genotype), :);
files = string(B.files(:));
E = nan(nFly, 1); Ew = nan(nFly, 1); ref = false(nFly, 1);
for i = 1:nFly
    m = find(T.file == files(i));
    if isscalar(m)
        ref(i) = strcmpi(string(T.qw_reference(m)), "True");
        E(i) = T.E_nats_per_h(m); Ew(i) = T.E_wakeonly_nats_per_h(m);
    end
end
use = ref & isfinite(E);
[rho1, p1] = corr(resp(use), E(use), 'Type', 'Spearman');
[rho2, p2] = corr(resp(use), Ew(use), 'Type', 'Spearman');

fprintf('20r-C (EXPLORATORY), %s: hyper_erlang Shapes %s, logL %.4f\n', genotype, mat2str(mm), H.LogLik);
fprintf('  rates %s, untruncated w %s\n', mat2str(rates, 4), mat2str(w, 4));
fprintf('  reported: %s\n', strjoin(compose('%s=%.4g', string(H.ParamNames(:)), H.Params(:)), ', '));
if isfield(H.Diagnostics, 'GuardOK'), fprintf('  GuardOK %d\n', H.Diagnostics.GuardOK); end
fprintf('  flies matched with a Q/W reference: %d of %d\n', nnz(use), nFly);
fprintf('  per-fly mean series-branch responsibility: median %.3f, range %.3f-%.3f\n', ...
    median(resp), min(resp), max(resp));
fprintf('  Spearman with E_nats_per_h:          rho %.3f, p %.3f, n %d\n', rho1, p1, nnz(use));
fprintf('  Spearman with E_wakeonly_nats_per_h: rho %.3f, p %.3f, n %d\n', rho2, p2, nnz(use));

C = struct('genotype', genotype, 'Shapes', mm, 'Theta', th, 'logL', H.LogLik, ...
    'resp', resp, 'nBouts', nb, 'E', E, 'Ewake', Ew, 'used', use, ...
    'rhoE', rho1, 'pE', p1, 'rhoEwake', rho2, 'pEwake', p2, 'n', nnz(use), 'files', {B.files});
save(fullfile(fileparts(mfilename('fullpath')), sprintf('run_20rC_%s.mat', genotype)), 'C', '-v7');
end

function S = erlangSF(t, lam, m)
% S(t) = exp(-lam t) sum_{i<m} (lam t)^i / i!, exact for integer m; S(t) = 1 for t <= 0.
x = lam * max(t, 0);
term = ones(size(x)); s = term;
for i = 1:m-1
    term = term .* x / i; s = s + term;
end
S = exp(-x) .* s;
end
