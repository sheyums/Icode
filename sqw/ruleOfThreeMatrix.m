function [Q, info] = ruleOfThreeMatrix(P, N, options)
%RULEOFTHREEMATRIX  Rule-of-three fill of one-way edges in a transition matrix (sqw_ref.rule_of_three_matrix; section 16, revision 663a0c2).
%
%   [Q, info] = ruleOfThreeMatrix(P, N)
%   [Q, info] = ruleOfThreeMatrix(P, N, Zero=1e-10)
%
%   P   (K x K) row-stochastic transition matrix of ONE phase (observed chain: CHAINFITPHASED;
%       HMM: the K = 3 fit)
%   N   (K x K) transition counts of the same phase (observed chain: observed counts; HMM: mean
%       posterior-path counts, need not be integers). Only off-diagonal entries are used.
%
%   An edge a -> b is ONE-WAY if P(a,b) < Zero and P(b,a) >= Zero. Per row a, from the ORIGINAL P:
%       n_a = sum_{j ~= a} N(a,j) / (1 - P(a,a))                    (row a's total steps, recovered once)
%       row = n_a P(a,:); row(b) = min(3, N(b,a)) for EVERY one-way edge b of the row with N(b,a) > 0
%       Q(a,:) = row / sum(row)                                      (renormalised once)
%   so the result does not depend on the order of the edges. A row with no observed exit
%   (sum_{j ~= a} N(a,j) = 0) or with P(a,a) = 1 is left UNCHANGED: the state is never left, there
%   is no cycle through it (MODELENERGETICSR3 then reports no cycle). Pairs unseen in BOTH
%   directions stay at zero. (The first version, 2026-11-02, renormalised after each fill and used a
%   fallback count for rows without an exit; both were removed in sqw_ref 663a0c2.)
%
%   info  struct: changed (true if any edge was filled, sqw_ref's second output), filled (K x K
%         logical), added (K x K, the counts set), stillOneWay (K x K logical: one-way edges left
%         one-way, because N(b,a) = 0 or the row has no observed exit), noExit (K x 1 logical).
%
%   See also MODELENERGETICSR3, MODELENERGETICS, EPR_COUNT_BASED.

arguments
    P (:,:) double {mustBeNonnegative}
    N (:,:) double {mustBeNonnegative}
    options.Zero (1,1) double {mustBePositive} = 1e-10
end
K = size(P, 1);
if size(P, 2) ~= K || ~isequal(size(N), [K K]), error('sqw:size', 'P and N must both be K x K.'); end
Z = P < options.Zero; Q = P; off = ~eye(K);
oneWay = Z & ~Z' & off;
info.filled = false(K); info.added = zeros(K); info.stillOneWay = false(K); info.noExit = false(K, 1);
for a = 1:K
    b = find(oneWay(a, :));
    if isempty(b), continue, end
    exits = sum(N(a, off(a, :)));
    if exits <= 0 || P(a, a) >= 1
        info.noExit(a) = true; info.stillOneWay(a, b) = true; continue
    end
    k = N(b, a)'; f = b(k > 0); info.stillOneWay(a, b(k <= 0)) = true;
    if isempty(f), continue, end
    row = P(a, :) * (exits / (1 - P(a, a)));
    row(f) = min(3, N(f, a)');
    Q(a, :) = row / sum(row);
    info.filled(a, f) = true; info.added(a, f) = row(f);
end
info.changed = any(info.filled(:));
end
