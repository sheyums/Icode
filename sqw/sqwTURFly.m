function F = sqwTURFly(y, ph, nP, options)
%SQWTURFLY  Thermodynamic uncertainty relation for one fly: data, exact chain and HMM, matched nulls (PRESPEC 20q, per fly).
%
%   F = sqwTURFly(y, ph, nP)
%   F = sqwTURFly(y, ph, nP, HMM=fit, NNullMC=200, NNullHMM=100, Taus=[180 120 360], ...
%                 RandomSeed=0, RuleOfThree=true)
%
%   y, ph, nP  as SQWMEMORYFLY
%   HMM        the fly's LABELLED K = 3 fit (pi, P, B); [] skips HMM quantities
%   RuleOfThree  amendment E: a phase whose fitted chain has a ONE-WAY edge
%              (kind i) gets b* = min(3, k) counts on the unseen direction,
%              k = the observed reverse count; its exact quantities and its
%              null records are recomputed from that matrix
%
%   F  struct:
%     obs    (nP x 1) TURESTIMATE of the data (J, D, B, n per tau)
%     ex     (nP x 1) table, exact per phase, per hour: sig_mc (MODELENERGETICS EPR of
%            the fitted first-order chain), J_mc, D_mc, B_mc (TUREXACT), B_mc_ct
%            (TUREXACTCT), kappa = B_mc / B_mc_ct; sig_hmm, J_hmm, D_hmm, B_hmm,
%            B_hmm_ct (the HMM's OBSERVED current, HMMOBSERVEDCHAIN);
%            kind: "cycle" (usable), "i" one-way edge, "ii" exactly symmetric counts,
%            "iii" no cycle (amendment E); r3 true where the rule of three was applied
%     simMC  (NNullMC x nP x nTau) B of records simulated from the chain
%     simH   (NNullHMM x nP x nTau) B of records simulated from the HMM
%     Taus   block lengths (minutes)
%
%   Every simulated record gets the fly's missing-minute mask and goes through
%   TURESTIMATE as the data. Chain rows never left are filled for simulation
%   (FILLEMPTYROWS) before anything is computed from the chain, as 20q did.
%
%   See also TURESTIMATE, TUREXACT, TUREXACTCT, HMMOBSERVEDCHAIN, SQWTURGROUP.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    options.HMM = []
    options.NNullMC (1,1) double {mustBeNonnegative, mustBeInteger} = 200
    options.NNullHMM (1,1) double {mustBeNonnegative, mustBeInteger} = 100
    options.Taus (1,:) double {mustBePositive, mustBeInteger} = [180 120 360]
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    options.RuleOfThree (1,1) logical = true
end
T = numel(y);
if numel(ph) ~= T, error('sqw:length', 'y and ph must have the same length.'); end
rs = RandStream('twister', 'Seed', options.RandomSeed);
miss = y == 0; taus = options.Taus; nT = numel(taus); I3 = eye(3);
[~, pm, Pm] = chainFitPhased(y, ph, nP); Pm = fillEmptyRows(Pm);
haveH = ~isempty(options.HMM);
if haveH, options.HMM.P = fillEmptyRows(options.HMM.P); end   % see SQWMEMORYFLY: HMM rows never left

F.Taus = taus;
F.obs = turEstimate(y, ph, nP, Taus=taus);
F.simMC = sims(pm, Pm, I3, ph, nP, miss, options.NNullMC, taus, rs);
if haveH, F.simH = sims(options.HMM.pi(:), options.HMM.P, options.HMM.B, ph, nP, miss, options.NNullHMM, taus, rs);
else, F.simH = NaN(0, nP, nT); end

ex = table();
for p = 1:nP
    r = chainExact(Pm(:, :, p));
    if haveH
        [Th, dh] = hmmObservedChain(options.HMM.P(:, :, p), options.HMM.B);
        [r.J_hmm, r.D_hmm, r.B_hmm] = turExact(Th, dh); [~, ~, r.B_hmm_ct] = turExactCT(Th, dh);
        r.sig_hmm = modelEnergetics(options.HMM.P(:, :, p)).EPR;
    else
        r.J_hmm = NaN; r.D_hmm = NaN; r.B_hmm = NaN; r.B_hmm_ct = NaN; r.sig_hmm = NaN;
    end
    C = phaseCounts(y, ph, p);
    if hasCycle(r), r.kind = "cycle"; else, r.kind = cycleKind(C, y, ph, p); end
    r.r3 = false;
    if options.RuleOfThree && r.kind == "i"
        Pr = ruleOfThree(C); Pm3 = Pm; Pm3(:, :, p) = Pr;
        h = r; r = chainExact(Pr);
        for f = ["J_hmm" "D_hmm" "B_hmm" "B_hmm_ct" "sig_hmm" "kind"], r.(f) = h.(f); end
        r.r3 = true;
        S = sims(pm, Pm3, I3, ph, nP, miss, options.NNullMC, taus, rs);
        F.simMC(:, p, :) = S(:, p, :);
    end
    ex = [ex; struct2table(r)]; %#ok<AGROW>
end
F.ex = ex;
end

function r = chainExact(P)
r.sig_mc = modelEnergetics(P).EPR;
[r.J_mc, r.D_mc, r.B_mc] = turExact(P); [~, ~, r.B_mc_ct] = turExactCT(P);
r.kappa = r.B_mc / r.B_mc_ct;
end

function tf = hasCycle(r)
% amendment D: sigma and the continuous-time bound must be finite and nonzero
tf = isfinite(r.sig_mc) && r.sig_mc > 1e-9 && isfinite(r.B_mc_ct) && r.B_mc_ct > 1e-12;
end

function S = sims(pi0, P, B, ph, nP, miss, n, taus, rs)
S = NaN(n, nP, numel(taus));
for k = 1:n
    [~, ys] = hmmSimulatePhased(pi0, P, B, ph, rs);
    ys(miss) = 0;
    R = turEstimate(ys, ph, nP, Taus=taus);
    for p = 1:nP, S(k, p, :) = R(p).B; end
end
end

function C = phaseCounts(y, ph, p)
t = find(ph(2:end) == p & y(2:end) >= 1 & y(2:end) <= 3 & y(1:end - 1) >= 1 & y(1:end - 1) <= 3) + 1;
C = accumarray([y(t - 1), y(t)], 1, [3 3]);
end

function k = cycleKind(C, y, ph, p)
seen = unique(y(ph == p & y >= 1 & y <= 3));
off = [1 2; 1 3; 2 3];
if numel(seen) < 3 || any(C(sub2ind([3 3], off(:, 1), off(:, 2))) == 0 & C(sub2ind([3 3], off(:, 2), off(:, 1))) == 0)
    k = "iii";
elseif any((C(sub2ind([3 3], off(:, 1), off(:, 2))) > 0) ~= (C(sub2ind([3 3], off(:, 2), off(:, 1))) > 0))
    k = "i";
else
    k = "ii";
end
end

function P = ruleOfThree(C)
for a = 1:3
    for b = a + 1:3
        if C(a, b) == 0 && C(b, a) > 0, C(a, b) = min(3, C(b, a));
        elseif C(b, a) == 0 && C(a, b) > 0, C(b, a) = min(3, C(a, b)); end
    end
end
P = C ./ sum(C, 2);
end
