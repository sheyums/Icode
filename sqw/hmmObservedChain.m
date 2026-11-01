function [T, d] = hmmObservedChain(P, B)
%HMMOBSERVEDCHAIN  Markov chain on (hidden state, symbol) pairs for the HMM's OBSERVED current (sqw_ref.hmm_observed_chain; PRESPEC 20q).
%
%   [T, d] = hmmObservedChain(P, B)
%
%   P  (K x K) one phase's hidden transition matrix
%   B  (K x 3) emissions, columns Q, W, S
%
%   T  (3K x 3K): state (z, y) has index 3 (z - 1) + y, and
%      T(3(z-1)+y, 3(z'-1)+y') = P(z, z') B(z', y')
%   d  (3K x 3K): increments of the counted cycle current on the SYMBOLS,
%      d(3(z-1)+y, 3(z'-1)+y') = CYCLEINCREMENTS(y, y')
%
%   TUREXACT(T, d) then gives the exact J, D and B = J^2/D of the current an
%   observer of this HMM's symbols would count. sqw_ref fixes K = 3; any K works here.
%
%   See also TUREXACT, CYCLEINCREMENTS.

arguments
    P (:,:) double {mustBeNonnegative}
    B (:,3) double {mustBeNonnegative}
end
K = size(P, 1);
if size(P, 2) ~= K || size(B, 1) ~= K, error('hmmObservedChain:size', 'P must be K x K and B K x 3.'); end
T = zeros(3 * K);
for z = 1:K
    row = reshape((P(z, :)' .* B)', 1, []);      % (z', y') in order 3 (z'-1) + y'
    T(3 * (z - 1) + (1:3), :) = repmat(row, 3, 1);
end
d = repmat(cycleIncrements(), K, K);
end
