function E = modelEnergetics(P, options)
%MODELENERGETICS  Exact cycle current, affinity and entropy production of a 3-state chain (sqw_ref.model_energetics).
%
%   E = modelEnergetics(P)
%   E = modelEnergetics(P, Zero=1e-10, StepsPerHour=60)
%
%   P   (3 x 3) labelled per-minute transition matrix, rows sum to 1 (an all-zero
%       row = a state never visited in this phase: its pi is 0, and with S
%       absent there is no cycle);
%       index 1 = H1/Q, 2 = H2/W, 3 = H3/S
%
%   E   struct, at P's stationary distribution pi, rates per hour:
%       J    = pi_2 P(2,1) - pi_1 P(1,2)        forward-positive (W -> Q -> S -> W),
%              equal on every edge of the 3-cycle
%       EPR  = sum_{i<j} (pi_i P_ij - pi_j P_ji) ln(pi_i P_ij / (pi_j P_ji))
%       A    = ln[P(2,1) P(1,3) P(3,2) / (P(1,2) P(3,1) P(2,3))]   nats per cycle
%       pi   (3 x 1) stationary distribution (STATIONARY_DIST, QR solve)
%       oneWay  true if some edge has P < Zero in exactly one direction
%
%   An edge with P < Zero in ONE direction makes EPR (and A) infinite: an
%   irreversible step dissipates without bound per unit flux. The infinite
%   value is reported, never floored. An edge below Zero in BOTH directions
%   contributes nothing. For a 3-state stationary chain sign(J) = sign(A);
%   |J| < 1e-12 means no cycle (report the block as "no cycle").
%
%   Zero (1e-10) is sqw_ref's numerical zero for FITTED matrices. It is not
%   the 1e-7 observability floor (P_min_bio) of COMPUTE_HMM_ENERGETICS, which
%   masked barriers in the older per-window analysis.
%
%   See also STATIONARY_DIST, ENERGYLANDSCAPE, ROUTINGPREDICTED.

arguments
    P (3,3) double {mustBeNonnegative}
    options.Zero (1,1) double {mustBeNonnegative} = 1e-10
    options.StepsPerHour (1,1) double {mustBePositive} = 60
end
checkStochastic(P);

h = options.StepsPerHour;
Z = P < options.Zero;
p = stationary_dist(P);
E.J = h * (p(2) * P(2, 1) - p(1) * P(1, 2));
E.A = log(P(2, 1) * P(1, 3) * P(3, 2) / (P(1, 2) * P(3, 1) * P(2, 3)));
E.pi = p;
E.oneWay = any(Z(:) ~= reshape(Z', [], 1));
if E.oneWay
    E.EPR = Inf;
    return
end
s = 0;
for i = 1:3
    for j = i + 1:3
        if Z(i, j), continue, end
        f = p(i) * P(i, j); r = p(j) * P(j, i);
        s = s + (f - r) * log(f / r);
    end
end
E.EPR = h * s;
end

function checkStochastic(P)
% A row of zeros is a state the chain never left in this phase (CHAINFITPHASED returns it
% that way, e.g. a fly with no sleep minute in its light half); it is allowed. Any other
% row must sum to 1.
s = sum(P, 2);
if any(abs(s - 1) > 1e-8 & s > 0)
    error('sqw:notStochastic', 'Rows of P must sum to 1 or be empty (max deviation %.3g).', max(abs(s - 1)));
end
end
