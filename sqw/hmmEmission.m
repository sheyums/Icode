function E = hmmEmission(B, y)
%HMMEMISSION  Emission likelihood of each minute under each state (sqw_ref.emission_lik).
%
%   E = hmmEmission(B, y)
%
%   B   (K x 3) emission probabilities, columns Q, W, S
%   y   (T x 1) symbols 1 = Q, 2 = W, 3 = S, 0 = missing, 4 = wake of unknown type
%   E   (T x K): E(t,k) = B(k, y(t)); a missing minute contributes 1 for every
%       state; an unknown-wake minute contributes B(k,Q) + B(k,W).

arguments
    B (:,3) double {mustBeNonnegative}
    y (:,1) double
end
if ~all(ismember(y, 0:4))
    error('sqw:symbols', 'Symbols must be 1 = Q, 2 = W, 3 = S, 0 = missing or 4 = unknown wake.');
end
K = size(B, 1);
Bx = [ones(1, K); B'; (B(:, 1) + B(:, 2))'];      % rows: missing, Q, W, S, unknown
E = Bx(y + 1, :);
end
