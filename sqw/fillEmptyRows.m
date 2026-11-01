function P = fillEmptyRows(P, Allow)
%FILLEMPTYROWS  Make every row of a fitted transition matrix simulable (sqw_ref.sim_rows; 20p amendment B).
%
%   P = fillEmptyRows(P)
%   P = fillEmptyRows(P, Allow)
%
%   P      (K x K x nP) transition matrices
%   Allow  (K x K) logical, allowed transitions (default all); for the pair
%          chain pass PAIRCHAINSTRUCTURE's Allow
%
%   FOR SIMULATION ONLY. A row whose sum is < 0.5 (zero expected count: a state
%   never left within that phase) is replaced by its count-start row, uniform
%   over the allowed transitions, so a simulated record cannot get stuck in a
%   state the data never leave. Such rows have zero weight in the likelihood,
%   so no fitted value changes. Other rows are returned unchanged.
%
%   See also CHAINFITPHASED, CHAIN2FITPHASED, HMMSIMULATEPHASED.

arguments
    P (:,:,:) double
    Allow logical = true(size(P, 1), size(P, 2))
end
K = size(P, 1);
if size(P, 2) ~= K || ~isequal(size(Allow), [K K])
    error('fillEmptyRows:size', 'P must be K x K x nP and Allow K x K.');
end
if any(sum(Allow, 2) == 0), error('fillEmptyRows:allow', 'Every row of Allow needs an allowed transition.'); end
U = double(Allow) ./ sum(Allow, 2);
for p = 1:size(P, 3)
    for k = 1:K
        if sum(P(k, :, p)) < 0.5, P(k, :, p) = U(k, :); end
    end
end
end
