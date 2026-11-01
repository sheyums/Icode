function s = roundTripStats(counts, hours)
%ROUNDTRIPSTATS  Summary of round trips for one fly and phase (sqw_ref.round_trip_stats; PRESPEC 20g).
%
%   s = roundTripStats(counts, hours)
%
%   counts  (1 x 4) for one record (observed level), or (nPaths x 4) for
%           posterior paths (hidden level: every statistic is computed per
%           path, then averaged; paths with no round trip are skipped)
%           columns [forward reverse direct quiet_both]
%   s       struct: n (mean count), frac (1 x 4), rate (1 x 4, per hour),
%           p_q_entry, p_q_exit, dependence (1 x 4) = fraction - product of the
%           marginals, frac_paths_used. Empty ([]) if no path has a round trip.

arguments
    counts (:,4) double
    hours (1,1) double {mustBePositive}
end
n = sum(counts, 2); ok = n > 0;
if ~any(ok), s = []; return, end
Fr = counts(ok, :) ./ n(ok);
qi = Fr(:, 1) + Fr(:, 4); qo = Fr(:, 2) + Fr(:, 4);
prod = [qi .* (1 - qo), (1 - qi) .* qo, (1 - qi) .* (1 - qo), qi .* qo];
s = struct('n', mean(n), 'frac', mean(Fr, 1), 'rate', mean(counts, 1) / hours, ...
    'p_q_entry', mean(qi), 'p_q_exit', mean(qo), 'dependence', mean(Fr - prod, 1), ...
    'frac_paths_used', mean(ok));
end
