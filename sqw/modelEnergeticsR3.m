function E = modelEnergeticsR3(P, N, options)
%MODELENERGETICSR3  MODELENERGETICS after the rule-of-three fill of one-way edges (sqw_ref.model_energetics_r3; revision 663a0c2).
%
%   E = modelEnergeticsR3(P, N)
%   E = modelEnergeticsR3(P, N, Zero=1e-10, StepsPerHour=60)
%
%   P, N  one phase's 3 x 3 transition matrix and its transition counts (RULEOFTHREEMATRIX)
%   E     MODELENERGETICS of RULEOFTHREEMATRIX(P, N): J, A, EPR (per hour), pi; plus
%         r3 (true if any edge was filled), r3info (RULEOFTHREEMATRIX's info), Pr3 (the filled
%         matrix, which ENERGYLANDSCAPE takes in place of P), and noCycle (see below).
%   No-cycle convention: if |J| < 1e-12 on the filled matrix AND the EPR or A would be non-finite
%   (e.g. a state that is never left, so the stationary mass sits in it), EPR = 0 and A = NaN,
%   noCycle = true. Finite round-off values are left as computed. The count-based EPR
%   (EPR_COUNT_BASED) is a different quantity and is unchanged by this convention.
%
%   See also RULEOFTHREEMATRIX, MODELENERGETICS, ENERGYLANDSCAPE.

arguments
    P (3,3) double {mustBeNonnegative}
    N (3,3) double {mustBeNonnegative}
    options.Zero (1,1) double {mustBePositive} = 1e-10
    options.StepsPerHour (1,1) double {mustBePositive} = 60
end
[Q, info] = ruleOfThreeMatrix(P, N, Zero=options.Zero);
E = modelEnergetics(Q, Zero=options.Zero, StepsPerHour=options.StepsPerHour);
E.r3 = info.changed; E.r3info = info; E.Pr3 = Q;
E.noCycle = abs(E.J) < 1e-12 && ~(isfinite(E.EPR) && isfinite(E.A));
if E.noCycle, E.EPR = 0; E.A = NaN; end
end
