function ph = minutePhases(T, options)
%MINUTEPHASES  Phase index of each minute (sqw_ref.phases).
%
%   ph = minutePhases(T)                    % "12h"
%   ph = minutePhases(T, Scheme="6h")
%
%   T        number of minutes (the record starts at ZT0)
%   Scheme   "12h": 1 = ZT 0-12, 2 = ZT 12-24
%            "6h" : 1..4 = ZT 0-6, 6-12, 12-18, 18-24
%            "none": all 1 (DD, one block)
%
%   ph(t) selects the transition matrix for the step INTO minute t, so ph(1)
%   is never used for a transition, and every transition, route or pair
%   count is attributed to the phase of its ARRIVAL minute.

arguments
    T (1,1) double {mustBeNonnegative, mustBeInteger}
    options.Scheme (1,1) string {mustBeMember(options.Scheme, ["12h" "6h" "none"])} = "12h"
end

t = mod((0:T - 1)', 1440);
switch options.Scheme
    case "12h", ph = 1 + double(t >= 720);
    case "6h",  ph = 1 + floor(t / 360);
    case "none", ph = ones(T, 1);
end
end
