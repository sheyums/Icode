function [route, t] = routeEvents(z)
%ROUTEEVENTS  Sleep-wake routes on a run-compressed state or symbol path (sqw_ref.route_events).
%
%   [route, t] = routeEvents(z)
%
%   z      (T x 1) path in 1..3 (1 = H1/Q, 2 = H2/W, 3 = H3/S); no missing minutes
%   route  (nRoutes x 1) index in ROUTES order: 1 = N312 (S -> Q -> W),
%          2 = N32 (S -> W), 3 = N213 (W -> Q -> S), 4 = N23 (W -> S)
%   t      (nRoutes x 1) minute of ARRIVAL in the destination (1-based)
%
%   Excursions S -> Q -> S and W -> Q -> W are not routes. For a sequence
%   with missing minutes use ROUTECOUNTSOBSERVED, which treats each run of
%   observed minutes separately.

arguments
    z (:,1) double
end
if any(z < 1 | z > 3)
    error('routeEvents:states', 'z must hold states 1..3 only.');
end
ch = [1; find(diff(z) ~= 0) + 1];
st = z(ch);
route = zeros(0, 1); t = zeros(0, 1);
for k = 2:numel(st)
    a = st(k - 1); b = st(k);
    if a == 3 && b == 2
        r = 2;
    elseif a == 2 && b == 3
        r = 4;
    elseif k >= 3 && a == 1 && st(k - 2) == 3 && b == 2
        r = 1;
    elseif k >= 3 && a == 1 && st(k - 2) == 2 && b == 3
        r = 3;
    else
        continue
    end
    route(end + 1, 1) = r; %#ok<AGROW>
    t(end + 1, 1) = ch(k); %#ok<AGROW>
end
end
