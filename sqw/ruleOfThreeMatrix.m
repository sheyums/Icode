function [Q, info] = ruleOfThreeMatrix(P, N, options)
%RULEOFTHREEMATRIX  Rule-of-three fill of one-way edges in a transition matrix (sqw_ref.rule_of_three_matrix; section 16, 2026-11-02).
%
%   [Q, info] = ruleOfThreeMatrix(P, N)
%   [Q, info] = ruleOfThreeMatrix(P, N, Zero=1e-10)
%
%   P   (K x K) row-stochastic transition matrix of ONE phase (observed chain: CHAINFITPHASED;
%       HMM: the K = 3 fit)
%   N   (K x K) transition counts of the same phase (observed chain: observed counts; HMM: mean
%       posterior-path counts, need not be integers). Only off-diagonal entries are used.
%
%   As sqw_ref, edge by edge in the order a = 1..K, b = 1..K: an edge a -> b is ONE-WAY if
%   P(a,b) < Zero and P(b,a) >= Zero in the ORIGINAL P. If k = N(b,a) > 0,
%       n_a = sum_{j ~= a} N(a,j) / max(1 - Q(a,a), 1e-12)     (or max(k, 1) if that sum is 0)
%       row = n_a Q(a,:); row(b) = min(3, k); Q(a,:) = row / sum(row)
%   where Q is the matrix AS UPDATED SO FAR: a row with two one-way edges is renormalised twice,
%   and the second fill sees the first (so the result can depend slightly on the edge order; this
%   mirrors the Python and is noted in PORT_VALIDATION). If k = 0 the edge is left one-way.
%   Pairs unseen in BOTH directions stay at zero.
%
%   info  struct: changed (true if any edge was filled, sqw_ref's second output), filled (K x K
%         logical), added (K x K, the counts set), stillOneWay (K x K logical: one-way edges with
%         N(b,a) = 0, left one-way; MODELENERGETICS then gives EPR = Inf there).
%
%   See also MODELENERGETICSR3, MODELENERGETICS, EPR_COUNT_BASED.

arguments
    P (:,:) double {mustBeNonnegative}
    N (:,:) double {mustBeNonnegative}
    options.Zero (1,1) double {mustBePositive} = 1e-10
end
K = size(P, 1);
if size(P, 2) ~= K || ~isequal(size(N), [K K]), error('sqw:size', 'P and N must both be K x K.'); end
Z = P < options.Zero; Q = P;
info.filled = false(K); info.added = zeros(K); info.stillOneWay = false(K);
for a = 1:K
    for b = 1:K
        if a == b || ~(Z(a, b) && ~Z(b, a)), continue, end
        k = N(b, a);
        if k <= 0, info.stillOneWay(a, b) = true; continue, end
        off = sum(N(a, [1:a - 1, a + 1:K]));
        if off > 0, na = off / max(1 - Q(a, a), 1e-12); else, na = max(k, 1); end
        row = Q(a, :) * na; row(b) = min(3, k); Q(a, :) = row / sum(row);
        info.filled(a, b) = true; info.added(a, b) = min(3, k);
    end
end
info.changed = any(info.filled(:));
end
