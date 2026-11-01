function C = routeCounts(z, ph, nP)
%ROUTECOUNTS  Route counts per phase of arrival on one complete path (sqw_ref.route_counts).
%
%   C = routeCounts(z, ph, nP)
%
%   C   (nP x 4) counts in ROUTES order [N312 N32 N213 N23] (ROUTEEVENTS)
%
%   See also ROUTECOUNTSOBSERVED, XYCOUNTS.

arguments
    z (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
end
if numel(z) ~= numel(ph), error('sqw:length', 'z and ph must have the same length.'); end
[r, t] = routeEvents(z);
C = accumarray([ph(t), r], 1, [nP 4]);
end
