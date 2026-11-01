function N = observablePairCounts(y, ph, nP)
%OBSERVABLEPAIRCOUNTS  Symbol-change counts on the observed S/Q/W sequence (sqw_ref.observable_pair_counts; PRESPEC 20f).
%
%   N = observablePairCounts(y, ph, nP)
%
%   y    symbols 1 = Q, 2 = W, 3 = S (0 missing, 4 wake of unknown type)
%   ph   phase of each minute (MINUTEPHASES), 1..nP
%   N    (3 x 3 x nP) N(i,j,p) = number of i -> j changes whose ARRIVAL minute
%        is in phase p
%
%   The sequence is run-compressed: only changes count, never stays. A
%   change counts only between two ADJACENT minutes that both carry Q, W or
%   S, so a change across a missing or unknown-wake minute is not counted.
%   The indices are those of the hidden states (Q = H1, W = H2, S = H3), so
%   JCYCLECOUNT and EPR_COUNT_BASED apply unchanged.
%
%   See also JCYCLECOUNT, EPR_COUNT_BASED, COBSMINUTES.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
end
checkSeq(y, ph, nP);

t = (2:numel(y))';
a = y(t - 1); b = y(t);
use = a >= 1 & a <= 3 & b >= 1 & b <= 3 & a ~= b;
N = accumarray([a(use), b(use), ph(t(use))], 1, [3 3 nP]);
end

function checkSeq(y, ph, nP)
if numel(y) ~= numel(ph)
    error('sqw:length', 'y and ph must have the same length.');
end
if ~all(ismember(y, 0:4))
    error('sqw:symbols', 'Symbols must be 1 = Q, 2 = W, 3 = S, 0 = missing or 4 = unknown wake.');
end
if any(ph > nP)
    error('sqw:phase', 'ph holds a phase above nP.');
end
end
