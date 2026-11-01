function fit = hmmFitPhased(y, ph, nP, K, options)
%HMMFITPHASED  Maximum-likelihood phase-switched HMM by batched Baum-Welch from random starts (sqw_ref.em_fit; PRESPEC 20b).
%
%   fit = hmmFitPhased(y, ph, nP, K)
%   fit = hmmFitPhased(y, ph, nP, K, Starts=20, MaxIter=1000, Tol=1e-7, RandomSeed=20, Start=[])
%
%   y    symbols 1 = Q, 2 = W, 3 = S, 0 = missing, 4 = unknown wake
%   ph   phase of each minute 1..nP (the step INTO minute t uses P(:,:,ph(t)))
%   K    number of hidden states
%   Start  optional struct with fields pi (K x 1), P (K x K x nP), B (K x 3):
%          ONE fixed start instead of random ones (e.g. a stored fit)
%
%   fit  struct: ll (best final log-likelihood), pi, P (K x K x nP), B (K x 3),
%        iterations, llStarts (final ll of every start), converged
%
%   All starts are iterated together and the loop stops when EVERY start's
%   relative gain (ll - ll_old) < Tol |ll|, or at MaxIter, exactly as sqw_ref
%   (so a start does not stop early on its own). Starts: pi, every row of every
%   P and every row of B ~ Dirichlet(1); the first floor(Starts/2) starts use
%   P <- 0.5 P + 0.5 I. numpy's draws cannot be reproduced: verify against
%   HMMEMSTEP (exact) and by reaching the same optimum; a higher likelihood
%   than the Python's is legitimate. Transition counts go to the phase of the
%   arrival minute; emissions are pooled over phases.
%
%   See also HMMEMSTEP, LABELSTATES, HMMBIC.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    K (1,1) double {mustBePositive, mustBeInteger}
    options.Starts (1,1) double {mustBePositive, mustBeInteger} = 20
    options.MaxIter (1,1) double {mustBePositive, mustBeInteger} = 1000
    options.Tol (1,1) double {mustBePositive} = 1e-7
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
    options.Start = []
end
if numel(y) ~= numel(ph), error('sqw:length', 'y and ph must have the same length.'); end
if ~all(ismember(y, 0:4)), error('sqw:symbols', 'Symbols must be 0..4.'); end
T = numel(y);

% ---------------------------------------------------------------- starts
if isempty(options.Start)
    S = options.Starts;
    rs = RandStream('twister', 'Seed', options.RandomSeed);
    dir1 = @(varargin) dirichlet1(rs, varargin{:});
    pi0 = dir1(K, S);                                   % K x S
    P = dir1(K, K * nP * S);                           % each column one Dirichlet row
    P = permute(reshape(P, K, K, nP, S), [2 1 3 4]);    % P(i,:,p,s) is a Dirichlet row
    half = floor(S / 2);
    P(:, :, :, 1:half) = 0.5 * P(:, :, :, 1:half) + 0.5 * eye(K);
    B = permute(reshape(dir1(3, K * S), 3, K, S), [2 1 3]);   % K x 3 x S
else
    S = 1; st = options.Start;
    pi0 = st.pi(:); P = reshape(st.P, K, K, nP, 1); B = reshape(st.B, K, 3, 1);
end

% ---------------------------------------------------------------- fixed pieces
onehot = double([y == 1, y == 2, y == 3]);             % T x 3
u = y == 4; anyU = any(u);
idx = arrayfun(@(p) find(ph(2:end) == p) + 1, 1:nP, 'UniformOutput', false);
llOld = -Inf(1, S);
for it = 1:options.MaxIter
    E = emissionBatch(B, y);                            % K x S x T
    [A, c] = fwdBatch(pi0, P, E, ph);
    ll = sum(log(c), 1);                                % 1 x S
    Bk = ones(K, S, T);
    for t = T - 1:-1:1
        Pt = squeeze4(P(:, :, ph(t + 1), :));           % K x K x S
        v = reshape(E(:, :, t + 1) .* Bk(:, :, t + 1), 1, K, S);
        Bk(:, :, t) = reshape(sum(Pt .* v, 2), K, S) ./ c(t + 1, :);
    end
    G = A .* Bk; G = G ./ sum(G, 1);
    R = E .* Bk ./ reshape(c', 1, S, T);
    Ap = permute(A, [1 3 2]); Rp = permute(R, [3 1 2]);   % K x T x S, T x K x S
    for p = 1:nP
        Xi = pagemtimes(Ap(:, idx{p} - 1, :), Rp(idx{p}, :, :)) .* squeeze4(P(:, :, p, :));
        P(:, :, p, :) = reshape(Xi ./ max(sum(Xi, 2), 1e-300), K, K, 1, S);
    end
    Gp = permute(G, [1 3 2]);                           % K x T x S
    Bn = pagemtimes(Gp, onehot);                        % K x 3 x S
    if anyU
        g = reshape(sum(Gp(:, u, :), 2), K, 1, S);
        qw = B(:, 1, :) + B(:, 2, :); qw(qw <= 0) = 1;
        Bn(:, 1, :) = Bn(:, 1, :) + g .* B(:, 1, :) ./ qw;
        Bn(:, 2, :) = Bn(:, 2, :) + g .* B(:, 2, :) ./ qw;
    end
    B = Bn ./ max(sum(Bn, 2), 1e-300);
    pi0 = G(:, :, 1);
    if all(ll - llOld < options.Tol * abs(ll))
        break
    end
    llOld = ll;
end
E = emissionBatch(B, y);
[~, c] = fwdBatch(pi0, P, E, ph);
ll = sum(log(c), 1); ll(~isfinite(ll)) = -Inf;
[best, s] = max(ll);
fit = struct('ll', best, 'pi', pi0(:, s), 'P', P(:, :, :, s), 'B', B(:, :, s), ...
    'iterations', it, 'llStarts', ll, 'converged', it < options.MaxIter, 'K', K, 'nP', nP);
end

% ---------------------------------------------------------------------------
function X = dirichlet1(rs, k, n)
% n draws from Dirichlet(ones(k,1)), as columns (k x n).
g = -log(rand(rs, k, n));
X = g ./ sum(g, 1);
end

function M = squeeze4(M)
% K x K x 1 x S  ->  K x K x S
M = reshape(M, size(M, 1), size(M, 2), []);
end

function E = emissionBatch(B, y)
% K x S x T emission likelihoods (HMMEMISSION for every start).
[K, ~, S] = size(B);
Bx = cat(2, ones(K, 1, S), B, B(:, 1, :) + B(:, 2, :));   % columns: missing, Q, W, S, unknown
E = permute(Bx(:, y + 1, :), [1 3 2]);
end

function [A, c] = fwdBatch(pi0, P, E, ph)
[K, S, T] = size(E);
A = zeros(K, S, T); c = zeros(T, S);
a = pi0 .* E(:, :, 1); c(1, :) = sum(a, 1); A(:, :, 1) = a ./ c(1, :);
for t = 2:T
    Pt = squeeze4(P(:, :, ph(t), :));
    a = reshape(sum(reshape(A(:, :, t - 1), K, 1, S) .* Pt, 1), K, S) .* E(:, :, t);
    c(t, :) = sum(a, 1); A(:, :, t) = a ./ c(t, :);
end
end
