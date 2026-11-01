function [Xn, Yn, c] = routingPredicted(P)
%ROUTINGPREDICTED  Routing predicted by a labelled 3x3 matrix (sqw_ref.xn_pred, yn_pred, c_pred).
%
%   [Xn, Yn, c] = routingPredicted(P)
%
%   P   (3 x 3) per-step matrix; 1 = H1/Q (quiet wake), 2 = H2/W (active wake),
%       3 = H3/S (sleep)
%
%   Branching out of H1: b12 = P12/(P12 + P13) to wake, b13 = 1 - b12 to sleep.
%   Sleep -> wake:  r = P31 b12 / P32,   Xn = (r - 1)/(r + 1)
%                   (+1 all via quiet wake, -1 all direct)
%   Wake -> sleep:  r = P21 b13 / P23,   Yn = (r - 1)/(r + 1)
%   c = (Yn - Xn)/2, the signed cycle current: sign(c) = sign(A) = sign(J).
%   Equal to ENERGYLANDSCAPE's F_SW and F_WS. 0/0 gives NaN, as in sqw_ref.

arguments
    P (3,3) double {mustBeNonnegative}
end
b12 = P(1, 2) / (P(1, 2) + P(1, 3));
b13 = P(1, 3) / (P(1, 2) + P(1, 3));
a = P(3, 1) * b12; Xn = (a - P(3, 2)) / (a + P(3, 2));
a = P(2, 1) * b13; Yn = (a - P(2, 3)) / (a + P(2, 3));
c = (Yn - Xn) / 2;
end
