function [Q, info] = ruleOfThreeMatrix(P, N, options)
%RULEOFTHREEMATRIX  Rule-of-three fill of one-way edges in a transition matrix (sqw_ref.rule_of_three_matrix; convention 2026-11-02).
%
%   [Q, info] = ruleOfThreeMatrix(P, N)
%   [Q, info] = ruleOfThreeMatrix(P, N, Zero=1e-10)
%
%   P   (K x K) row-stochastic transition matrix of ONE phase (observed chain: CHAINFITPHASED;
%       HMM: the K = 3 fit), states in the matrix's own order
%   N   (K x K) transition counts of the same phase (observed chain: observed counts; HMM: mean
%       posterior-path counts, need not be integers). Only off-diagonal entries are used.
%
%   For every edge a -> b with P(a,b) < Zero while P(b,a) >= Zero (a ONE-WAY edge), the unseen
%   count is set to min(3, N(b,a)) and row a is renormalised on its count scale:
%       n_a = sum_{j ~= a} N(a,j) / (1 - P(a,a))      (row a's total steps)
%       Q(a,:) = (n_a P(a,:) + sum over a's one-way edges b of min(3, N(b,a)) e_b) / (n_a + that sum)
%   All of a row's one-way edges are filled together, from the ORIGINAL P, so the result does not
%   depend on the order of the edges. Pairs unseen in BOTH directions stay at zero; rows without a
%   one-way edge are returned unchanged.
%
%   info  struct: filled (K x K logical, the edges given a count), added (K x K, the counts added),
%         nSteps (K x 1, n_a; NaN where not needed), stillOneWay (K x K logical: one-way edges
%         left one-way because N(b,a) = 0 or row a has no steps to scale, n_a not finite and > 0).
%         A stillOneWay edge keeps MODELENERGETICS' infinite EPR; it is reported, not hidden.
%
%   See also MODELENERGETICSR3, MODELENERGETICS, EPR_COUNT_BASED.

arguments
    P (:,:) double {mustBeNonnegative}
    N (:,:) double {mustBeNonnegative}
    options.Zero (1,1) double {mustBePositive} = 1e-10
end
K = size(P, 1);
if size(P, 2) ~= K || ~isequal(size(N), [K K]), error('sqw:size', 'P and N must both be K x K.'); end
z = options.Zero; Z = P < z; off = ~eye(K);
oneWay = Z & ~Z' & off;                                 % a -> b unseen, b -> a seen
Q = P; info.filled = false(K); info.added = zeros(K); info.nSteps = NaN(K, 1); info.stillOneWay = false(K);
for a = 1:K
    b = find(oneWay(a, :));
    if isempty(b), continue, end
    na = sum(N(a, off(a, :))) / (1 - P(a, a));
    c = min(3, N(b, a)');
    if ~(isfinite(na) && na > 0)
        info.stillOneWay(a, b) = true; continue
    end
    info.nSteps(a) = na;
    m = na * P(a, :); m(b) = m(b) + c;                  % P(a,b) < Zero: count na*P(a,b) is ~0 and kept
    Q(a, :) = m / sum(m);
    info.filled(a, b) = c > 0; info.added(a, b) = c; info.stillOneWay(a, b) = c == 0;
end
end
