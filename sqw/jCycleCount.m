function J = jCycleCount(N, hours)
%JCYCLECOUNT  Counted cycle current, cycles per hour (sqw_ref.j_cycle_count).
%
%   J = jCycleCount(N, hours)
%
%   N      (3 x 3 x nP) directed change counts (OBSERVABLEPAIRCOUNTS, or a
%          hidden path's pair counts); index 1 = Q/H1, 2 = W/H2, 3 = S/H3
%   hours  (nP x 1) hours in each phase
%   J(p) = (N(3,2,p) - N(2,3,p)) / hours(p), forward-positive (S -> W is the
%          forward leg of W -> Q -> S -> W). Linear in the counts, so the mean
%          over posterior paths equals J of the mean counts.
%
%   Same numerator as EPR_COUNT_BASED's J_cycle.

arguments
    N (3,3,:) double
    hours (:,1) double {mustBePositive}
end
nP = size(N, 3);
if numel(hours) ~= nP
    error('jCycleCount:hours', 'hours needs one entry per phase (%d).', nP);
end
J = (squeeze(N(3, 2, :)) - squeeze(N(2, 3, :))) ./ hours;
J = J(:);
end
