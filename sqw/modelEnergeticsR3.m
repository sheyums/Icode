function E = modelEnergeticsR3(P, N, options)
%MODELENERGETICSR3  MODELENERGETICS after the rule-of-three fill of one-way edges (sqw_ref.model_energetics_r3; convention 2026-11-02).
%
%   E = modelEnergeticsR3(P, N)
%   E = modelEnergeticsR3(P, N, Zero=1e-10, StepsPerHour=60)
%
%   P, N  one phase's 3 x 3 transition matrix and its transition counts (RULEOFTHREEMATRIX)
%   E     MODELENERGETICS of RULEOFTHREEMATRIX(P, N): J, A, EPR (per hour), pi; plus
%         r3 (true if any edge was filled), r3info (RULEOFTHREEMATRIX's info) and Pr3 (the
%         filled matrix, which ENERGYLANDSCAPE takes in place of P).
%   No one-way edge gives infinity any more, unless RULEOFTHREEMATRIX could not fill it
%   (r3info.stillOneWay). The count-based EPR (EPR_COUNT_BASED) is a different quantity and is
%   unchanged by this convention.
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
end
