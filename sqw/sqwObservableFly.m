function F = sqwObservableFly(source, options)
%SQWOBSERVABLEFLY  Plan A observable analysis of one fly (sqw17.py's per-fly step, from sqw_ref functions).
%
%   F = sqwObservableFly(file)
%   F = sqwObservableFly(counts, LightDark=false, Theta=0.5, UnknownWake="missing", ...
%                        SixHour=false, Plot=false)
%
%   source       a count file (READCOUNTS) or a vector of 1-s counts from ZT0
%   LightDark    true for LD records (r_f and transition matrices per light /
%                dark half-day); false for DD (one r_f, one matrix: phase "whole")
%   Theta        Q/W threshold (SQWSYMBOLS), default 0.5
%   UnknownWake  "missing" (ledger) or "unknown" (SQWSYMBOLS)
%   SixHour      also fit and summarise the four 6-h ZT blocks (LD only; C3)
%   Plot         draw the fly's symbol raster and its per-phase numbers
%
%   F   struct:
%     file, T (minutes), y (symbols), ph, nP, phaseNames
%     refs        SQWSYMBOLS references; refOK (nP x 1) true if the phase has one
%     chain       struct ll, pi, P (3 x 3 x nP), iterations, converged
%     Q           (nP x 1) OBSERVABLEQUANTITIES; Q(p).refOK and Q(p).hasCycle added
%     six         (if SixHour) the same for phases "6h": ph, chain, Q, refOK
%   A fly-phase without a reference is NOT removed here: refOK says so, and
%   the group step excludes it. hasCycle = |J_chain| >= 1e-12.
%
%   See also SQWOBSERVABLEGROUP, PLOTSQWFLY, OBSERVABLEQUANTITIES.

arguments
    source
    options.LightDark (1,1) logical = true
    options.Theta (1,1) double {mustBePositive} = 0.5
    options.UnknownWake (1,1) string {mustBeMember(options.UnknownWake, ["missing" "unknown"])} = "missing"
    options.SixHour (1,1) logical = false
    options.Plot (1,1) logical = false
end
if ischar(source) || isstring(source)
    F.file = string(source);
    x = readCounts(source);
else
    F.file = "";
    x = double(source(:));
end
[y, refs, info] = sqwSymbols(x, Theta=options.Theta, LightDark=options.LightDark, ...
    UnknownWake=options.UnknownWake);
F.T = numel(y); F.y = y; F.refs = refs; F.nBouts = size(info.bouts, 1); F.nWakeEvents = numel(info.events);
if options.LightDark
    F.ph = minutePhases(F.T); F.nP = 2; F.phaseNames = ["light"; "dark"];
    F.refOK = [isfinite(refs.light); isfinite(refs.dark)];
else
    F.ph = minutePhases(F.T, Scheme="none"); F.nP = 1; F.phaseNames = "whole";
    F.refOK = isfinite(refs.whole);
end
if F.T == 0
    error('sqwObservableFly:short', 'The record holds less than one whole day.');
end
F.chain = fitChain(y, F.ph, F.nP);
F.Q = annotate(observableQuantities(y, F.ph, F.nP, F.chain.P), F.refOK);

if options.SixHour
    if ~options.LightDark
        error('sqwObservableFly:sixDD', 'SixHour needs an LD record (ZT blocks).');
    end
    s.ph = minutePhases(F.T, Scheme="6h");
    s.refOK = F.refOK([1 1 2 2]);
    s.chain = fitChain(y, s.ph, 4);
    s.Q = annotate(observableQuantities(y, s.ph, 4, s.chain.P), s.refOK);
    s.phaseNames = ["ZT0-6"; "ZT6-12"; "ZT12-18"; "ZT18-24"];
    F.six = s;
end
if options.Plot
    plotSqwFly(F);
end
end

function c = fitChain(y, ph, nP)
[c.ll, c.pi, c.P, info] = chainFitPhased(y, ph, nP);
c.iterations = info.iterations; c.converged = info.converged;
end

function Q = annotate(Q, refOK)
for p = 1:numel(Q)
    Q(p).refOK = refOK(p);
    Q(p).hasCycle = abs(Q(p).J_chain) >= 1e-12;
end
end
