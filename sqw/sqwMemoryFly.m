function F = sqwMemoryFly(y, ph, nP, options)
%SQWMEMORYFLY  Block irreversibility of one fly against matched first-order-chain nulls (PRESPEC 20p, per fly).
%
%   F = sqwMemoryFly(y, ph, nP)
%   F = sqwMemoryFly(y, ph, nP, HMM=fit, NNullA=200, NNullB=200, NModel=50, RandomSeed=0, Halves=false)
%
%   y    (T x 1) symbols (SQWSYMBOLS; missing = 0), whole days
%   ph   (T x 1) phases 1..nP: 12-h light/dark for LD records, all 1 for DD
%   HMM  struct with pi (3 x 1), P (3 x 3 x nP), B (3 x 3): the fly's LABELLED K = 3
%        fit (20b), for the HMM's predicted Delta3; [] skips the HMM prediction
%   Halves  (DD records, amendment A) also run the estimator in the two 12-h
%        clock halves, with null A2 (two 12-h chain matrices) and null B
%
%   Every simulated record gets the fly's missing-minute mask and goes through
%   BLOCKKLSTATS exactly as the data (sample-size-matched null).
%     null A   NNullA records from the first-order chain fitted in ph (CHAINFITPHASED)
%     null B   NNullB records from the chain with four 6-h matrices (drift guard, 20o)
%     HMM, MC2 NModel records each from the HMM and the second-order chain
%              (predictions, not nulls)
%   Empty chain rows are filled for simulation only (FILLEMPTYROWS).
%
%   F  struct:
%     obs     (nP x 1) BLOCKKLSTATS of the data
%     phase   (nP x 1) table: d2 d3 d4 D3 D4 D3QW (observed, nats/h);
%             E_A = D3 - median D3 under null A, p_A = (1 + #{D3* >= D3}) / (1 + NNullA)
%             (one-sided: excess irreversibility); E_B, p_B likewise under null B;
%             E_D4, E_D3QW (each minus its null-A median);
%             pred_HMM, in_HMM (median simulated D3; observed inside the central 95%),
%             pred_MC2, in_MC2
%     sims    struct of (nSim x nP) arrays: A_D3, A_D4, A_D3QW, B_D3, HMM_D3, MC2_D3
%     halves  (Halves only) table for the two 12-h halves: d2 d3 d4 D3, E_A2, p_A2, E_B, p_B
%     chain   struct: P (first-order, 12-h/whole), P6 (6-h), P2 (second order), each as fitted
%
%   Percentiles are type 7 (QUANTILE7), as numpy.
%
%   See also BLOCKKLSTATS, SQWMEMORYGROUP, CHAINFITPHASED, CHAIN2FITPHASED.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    options.HMM = []
    options.NNullA (1,1) double {mustBeNonnegative, mustBeInteger} = 200
    options.NNullB (1,1) double {mustBeNonnegative, mustBeInteger} = 200
    options.NModel (1,1) double {mustBeNonnegative, mustBeInteger} = 50
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    options.Halves (1,1) logical = false
end
T = numel(y);
if numel(ph) ~= T, error('sqw:length', 'y and ph must have the same length.'); end
rs = RandStream('twister', 'Seed', options.RandomSeed);
miss = y == 0; I3 = eye(3);
ph6 = minutePhases(T, Scheme="6h");
% A hidden state never left within a phase (e.g. H3 in the light of a fly with no S
% there) has an all-zero row. Filled uniform for simulation, as 20p amendment B does
% for the chain. sqw_ref.simulate instead jumps to the LAST state, which for H3 keeps a
% simulated fly asleep through the whole phase (PORT_VALIDATION, stage 3).
if ~isempty(options.HMM), options.HMM.P = fillEmptyRows(options.HMM.P); end

[~, pa, Pa] = chainFitPhased(y, ph, nP);
[~, pb, Pb] = chainFitPhased(y, ph6, 4);
[~, p2, P2, i2] = chain2FitPhased(y, ph, nP);
F.chain = struct('P', Pa, 'P6', Pb, 'P2', P2);
Pa = fillEmptyRows(Pa); Pb = fillEmptyRows(Pb); P2 = fillEmptyRows(P2, i2.Allow);

F.obs = blockKLStats(y, ph, nP);
[A3, A4, Aq] = simStats(pa, Pa, I3, ph, ph, nP, miss, options.NNullA, rs);
B3 = simStats(pb, Pb, I3, ph6, ph, nP, miss, options.NNullB, rs);
if isempty(options.HMM), H3 = NaN(0, nP);
else, H3 = simStats(options.HMM.pi(:), options.HMM.P, options.HMM.B, ph, ph, nP, miss, options.NModel, rs); end
M3 = simStats(p2, P2, i2.B, ph, ph, nP, miss, options.NModel, rs);
F.sims = struct('A_D3', A3, 'A_D4', A4, 'A_D3QW', Aq, 'B_D3', B3, 'HMM_D3', H3, 'MC2_D3', M3);

o = struct2table(F.obs);
for p = 1:nP
    o.E_A(p) = o.D3(p) - median(A3(:, p)); o.p_A(p) = pval(A3(:, p), o.D3(p));
    o.E_B(p) = o.D3(p) - median(B3(:, p)); o.p_B(p) = pval(B3(:, p), o.D3(p));
    o.E_D4(p) = o.D4(p) - median(A4(:, p)); o.E_D3QW(p) = o.D3QW(p) - median(Aq(:, p));
    [o.pred_HMM(p), o.in_HMM(p)] = predict(H3(:, p), o.D3(p));
    [o.pred_MC2(p), o.in_MC2(p)] = predict(M3(:, p), o.D3(p));
end
F.phase = o;

if options.Halves
    ph2 = minutePhases(T, Scheme="12h");
    [~, pa2, Pa2] = chainFitPhased(y, ph2, 2); Pa2 = fillEmptyRows(Pa2);
    h = struct2table(blockKLStats(y, ph2, 2)); h = h(:, {'d2', 'd3', 'd4', 'D3'});
    S2 = simStats(pa2, Pa2, I3, ph2, ph2, 2, miss, options.NNullA, rs);
    SB = simStats(pb, Pb, I3, ph6, ph2, 2, miss, options.NNullB, rs);
    for p = 1:2
        h.E_A2(p) = h.D3(p) - median(S2(:, p)); h.p_A2(p) = pval(S2(:, p), h.D3(p));
        h.E_B(p) = h.D3(p) - median(SB(:, p)); h.p_B(p) = pval(SB(:, p), h.D3(p));
    end
    F.halves = h;
end
end

function [S3, S4, Sq] = simStats(pi0, P, B, simPh, statPh, nP, miss, n, rs)
S3 = NaN(n, nP); S4 = S3; Sq = S3;
for r = 1:n
    [~, ys] = hmmSimulatePhased(pi0, P, B, simPh, rs);
    ys(miss) = 0;
    s = blockKLStats(ys, statPh, nP);
    S3(r, :) = [s.D3]; S4(r, :) = [s.D4]; Sq(r, :) = [s.D3QW];
end
end

function p = pval(s, o)
p = (1 + sum(s >= o)) / (1 + numel(s));
end

function [m, in] = predict(s, o)
if isempty(s), m = NaN; in = NaN; return, end
m = median(s); q = quantile7(s, [0.025 0.975]);
in = double(q(1) <= o && o <= q(2));
end
