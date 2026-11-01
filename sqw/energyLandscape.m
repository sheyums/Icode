function L = energyLandscape(P)
%ENERGYLANDSCAPE  Energy landscape of a labelled 3x3 matrix (sqw_ref.landscape; Methods Eqs. 5-16, PRESPEC 20l).
%
%   L = energyLandscape(P)
%
%   P   (3 x 3) labelled matrix, 1 = H1/Q, 2 = H2/W, 3 = H3/S
%
%   L   struct, energies in nats (k_B T = 1), relative to H1 = quiet wake:
%       E2_minus_E1 = ln(pi_1/pi_2), E3_minus_E1 = ln(pi_1/pi_3)
%       b12 = P12/(P12 + P13), b13 = 1 - b12            branching out of H1
%       E31c = -ln(P31 b12), E21c = -ln(P21 b13)         branching-corrected barriers
%       E32 = -ln P32, E23 = -ln P23                     direct barriers
%       F_SW = tanh((E32 - E31c)/2), F_WS = tanh((E23 - E21c)/2)
%       pi   stationary distribution (STATIONARY_DIST)
%   F_SW and F_WS equal ROUTINGPREDICTED's Xn and Yn. A zero probability gives
%   an infinite barrier, which is reported as such.

arguments
    P (3,3) double {mustBeNonnegative}
end
p = stationary_dist(P);
b12 = P(1, 2) / (P(1, 2) + P(1, 3)); b13 = 1 - b12;
L.E2_minus_E1 = log(p(1) / p(2));
L.E3_minus_E1 = log(p(1) / p(3));
L.b12 = b12; L.b13 = b13;
L.E31c = -log(P(3, 1) * b12);
L.E21c = -log(P(2, 1) * b13);
L.E32 = -log(P(3, 2));
L.E23 = -log(P(2, 3));
L.F_SW = tanh((L.E32 - L.E31c) / 2);
L.F_WS = tanh((L.E23 - L.E21c) / 2);
L.pi = p;
end
