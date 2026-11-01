function C = routeCountsObserved(y, ph, nP)
%ROUTECOUNTSOBSERVED  Route counts on the observed symbols (sqw_ref.routes_observed; PRESPEC 20s-A).
%
%   C = routeCountsObserved(y, ph, nP)
%
%   y   symbols 1 = Q, 2 = W, 3 = S; 0 (missing) and 4 (unknown wake) split
%       the record, and each run of observed minutes is scanned separately,
%       so no route is assembled across a gap
%   C   (nP x 4) counts in ROUTES order [N312 N32 N213 N23], by phase of the
%       arrival minute. Equals ROUTECOUNTS on a complete record.
%
%   See also ROUTEEVENTS, XYCOUNTS.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
end
if numel(y) ~= numel(ph), error('sqw:length', 'y and ph must have the same length.'); end
ok = y >= 1 & y <= 3;
d = diff([false; ok; false]);
s = find(d == 1); e = find(d == -1) - 1;
C = zeros(nP, 4);
for k = 1:numel(s)
    [r, t] = routeEvents(y(s(k):e(k)));
    if isempty(r), continue, end
    C = C + accumarray([ph(s(k) + t - 1), r], 1, [nP 4]);
end
end
