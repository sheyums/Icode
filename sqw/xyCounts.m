function [Xn, Yn] = xyCounts(C, options)
%XYCOUNTS  Counted routing from route counts (sqw_ref.xy_counts).
%
%   [Xn, Yn] = xyCounts(C)
%   [Xn, Yn] = xyCounts(C, MinRoutes=3)
%
%   C    (n x 4) route counts in ROUTES order [N312 N32 N213 N23]
%   Xn = (N312 - N32)/(N312 + N32)    sleep -> wake, +1 = all via quiet wake
%   Yn = (N213 - N23)/(N213 + N23)    wake -> sleep
%   NaN where the denominator is below MinRoutes.

arguments
    C (:,4) double {mustBeNonnegative}
    options.MinRoutes (1,1) double {mustBeNonnegative} = 3
end
dX = C(:, 1) + C(:, 2); dY = C(:, 3) + C(:, 4);
Xn = (C(:, 1) - C(:, 2)) ./ dX;
Yn = (C(:, 3) - C(:, 4)) ./ dY;
Xn(dX < options.MinRoutes) = NaN;
Yn(dY < options.MinRoutes) = NaN;
end
