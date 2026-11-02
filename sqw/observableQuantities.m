function Q = observableQuantities(y, ph, nP, Pchain, options)
%OBSERVABLEQUANTITIES  Plan A per-phase quantities at the observed S/Q/W level (sqw_ref.observable_quantities; PRESPEC 20s-A).
%
%   Q = observableQuantities(y, ph, nP, Pchain)
%   Q = observableQuantities(y, ph, nP, Pchain, RuleOfThree=false)     % the pre-2026-11-02 convention
%
%   y       symbols (SQWSYMBOLS), ph phases 1..nP (MINUTEPHASES)
%   Pchain  (3 x 3 x nP) the fly's observed first-order chain (CHAINFITPHASED)
%
%   Q   (nP x 1) struct, per phase p, with hours = minutes in phase p / 60:
%       c_obs          COBSMINUTES
%       J              count J, cycles/h (OBSERVABLEPAIRCOUNTS -> JCYCLECOUNT)
%       EPR            count EPR, nats/h, rule-of-three lower bound (EPR_COUNT_BASED 'bound')
%       rule_of_three  true if some pair was seen in one direction only
%       J_chain        MODELENERGETICS of Pchain(:,:,p) as fitted (defines the cycle, |J_chain| >= 1e-12)
%       A_chain, EPR_chain   MODELENERGETICSR3 of Pchain(:,:,p) with the phase's counts N: a one-way
%                      edge gets the rule-of-three count (convention of 2026-11-02, sqw_ref section 16)
%                      instead of an infinite A or EPR; r3_chain = true where that happened
%       E21, E31       ENERGYLANDSCAPE E2-E1, E3-E1 (nats; E1 = quiet wake) of the same filled matrix
%       F_SW, F_WS     landscape-predicted routing, same matrix
%   RuleOfThree=false gives MODELENERGETICS and ENERGYLANDSCAPE of the fitted matrix (Inf kept), as
%   every number before 2026-11-02 and the golden files.
%       Xn_cnt, Yn_cnt counted routing on the observed symbols (ROUTECOUNTSOBSERVED -> XYCOUNTS)
%       hours, N (3 x 3 counts), routes (1 x 4 counts)
%   Sign: forward W -> Q -> S -> W positive throughout.
%
%   Exclusions are the caller's (per fly x phase, counted and reported):
%   no Q/W reference, and for chain quantities a chain without a cycle
%   (|J_chain| < 1e-12).

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    Pchain (3,3,:) double
    options.RuleOfThree (1,1) logical = true
end
if size(Pchain, 3) ~= nP, error('observableQuantities:P', 'Pchain needs one matrix per phase.'); end
N = observablePairCounts(y, ph, nP);
c = cObsMinutes(y, ph, nP);
R = routeCountsObserved(y, ph, nP);
[X, Y] = xyCounts(R);
hours = accumarray(ph, 1, [nP 1]) / 60;
J = jCycleCount(N, hours);

Q = repmat(struct(), nP, 1);
for p = 1:nP
    Np = N(:, :, p);
    E0 = modelEnergetics(Pchain(:, :, p));
    if options.RuleOfThree, E = modelEnergeticsR3(Pchain(:, :, p), Np); Pl = E.Pr3; else, E = E0; E.r3 = false; Pl = Pchain(:, :, p); end
    L = energyLandscape(Pl);
    Q(p).c_obs = c(p);
    Q(p).J = J(p);
    Q(p).EPR = epr_count_based(Np, hours(p) * 3600, 'bound');
    Q(p).rule_of_three = oneWay(Np);
    Q(p).J_chain = E0.J;
    Q(p).A_chain = E.A;
    Q(p).EPR_chain = E.EPR;
    Q(p).r3_chain = E.r3;
    Q(p).E21 = L.E2_minus_E1;
    Q(p).E31 = L.E3_minus_E1;
    Q(p).F_SW = L.F_SW;
    Q(p).F_WS = L.F_WS;
    Q(p).Xn_cnt = X(p);
    Q(p).Yn_cnt = Y(p);
    Q(p).hours = hours(p);
    Q(p).N = Np;
    Q(p).routes = R(p, :);
end
end

function tf = oneWay(N)
% sqw_ref.epr_count_bound's "substituted": some pair seen in exactly one direction.
u = [N(1, 2) N(1, 3) N(2, 3)]; v = [N(2, 1) N(3, 1) N(3, 2)];
tf = any(xor(u == 0, v == 0));
end
