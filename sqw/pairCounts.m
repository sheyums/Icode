function N = pairCounts(z, ph, nP)
%PAIRCOUNTS  Single-step change counts on a hidden path (sqw_ref.pair_counts).
%
%   N = pairCounts(z, ph, nP)
%
%   z   (T x 1) complete path in 1..3;  N (3 x 3 x nP), N(i,j,p) = number of
%       i -> j changes on the run-compressed path, by phase of the ARRIVAL
%       minute (global, no window seams). Same counting as
%       OBSERVABLEPAIRCOUNTS, which on a complete sequence is identical.

arguments
    z (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
end
if any(z < 1 | z > 3), error('pairCounts:states', 'z must hold states 1..3 only.'); end
N = observablePairCounts(z, ph, nP);
end
