function test_sqw_r3()
%TEST_SQW_R3  Rule-of-three fill of one-way edges (RULEOFTHREEMATRIX, MODELENERGETICSR3): stress tests.
%
%   Run by name from the sqw folder:   cd <Icode>/sqw; test_sqw_r3
%
%    1. no one-way edge: Q is P bit for bit, energetics equal MODELENERGETICS
%    2. hand example: one unseen edge, n_a = 20 / 0.1 = 200, N(b,a) = 7 -> Q(a,:) = [180 3 20] / 203;
%       N(b,a) = 2 -> the added count is 2, not 3
%    3. pairs unseen both ways stay 0; rows without a one-way edge unchanged; every row sums to 1
%    4. two one-way edges in one row are filled TOGETHER against the original row (sqw_ref 663a0c2;
%       hand: [34 2 3 2 4] / 45); relabelling the states permutes the answer (no order dependence)
%    5. when P is the empirical matrix of integer counts N, n_a is each row's actual number of steps,
%       so every row with an observed exit = (N + added counts) normalised, however many fills
%    6. N(b,a) = 0: the edge is left one-way and reported, EPR stays Inf (J ~= 0); a row with no
%       observed exit, or P(a,a) = 1, is left UNCHANGED (no fallback count), noExit reported
%    7. modelEnergeticsR3 on a one-way chain: EPR and A finite, EPR > 0, Pr3 = Q, r3 true, J sign
%       as the forward cycle
%    8. size checks: non-square P, N of another size refused
%    9. observableQuantities wiring: on a record with a one-way chain edge, A_chain / EPR_chain / E21 / E31
%       = modelEnergeticsR3 and energyLandscape of the filled matrix, r3_chain true, J_chain from the
%       FITTED matrix; RuleOfThree=false restores modelEnergetics (Inf) exactly
%   10. hmmPosteriorFluxEPR wiring with B = I (paths = the record): EPR_hmm_model = modelEnergeticsR3(P, counts),
%       r3_hmm true; RuleOfThree=false gives the old Inf
%   11. no-cycle convention: an absorbing chain (|J| < 1e-12, EPR Inf) -> EPR = 0, A = NaN, noCycle;
%       a reversible chain with finite round-off J and EPR is left exactly as MODELENERGETICS gives it

here = fileparts(mfilename('fullpath'));
assert(isequal(fileparts(which('ruleOfThreeMatrix')), here), 'Run from the sqw folder.');
nP = 0; nF = 0;
fprintf('\nRunning rule-of-three tests...\n\n');

%% 1
P = [0.8 0.1 0.1; 0.2 0.7 0.1; 0.05 0.15 0.8]; N = 10 * ones(3);
[Q, info] = ruleOfThreeMatrix(P, N); E = modelEnergeticsR3(P, N); E0 = modelEnergetics(P);
ok = isequal(Q, P) && ~any(info.filled(:)) && E.EPR == E0.EPR && E.A == E0.A && E.J == E0.J && ~E.r3;
[nP, nF] = rep(ok, nP, nF, 1, 'no one-way edge: P unchanged, energetics as modelEnergetics', '');

%% 2
P = [0.9 0 0.1; 0.05 0.9 0.05; 0.1 0.1 0.8]; N = zeros(3); N(1, 3) = 20; N(2, 1) = 7; N(2, 3) = 4; N(3, 1) = 9; N(3, 2) = 9;
[Q, info] = ruleOfThreeMatrix(P, N); N2 = N; N2(2, 1) = 2; Q2 = ruleOfThreeMatrix(P, N2);
ok = max(abs(Q(1, :) - [180 3 20] / 203)) < 1e-15 && isequal(Q(2:3, :), P(2:3, :)) ...
    && info.added(1, 2) == 3 && max(abs(Q2(1, :) - [180 2 20] / 202)) < 1e-15;
[nP, nF] = rep(ok, nP, nF, 2, 'hand example: [180 3 20]/203, and min(3, N(b,a)) = 2 when only 2 seen', mat2str(Q(1, :), 6));

%% 3
P = [0.9 0 0.1; 0 0.95 0.05; 0.1 0.1 0.8]; N = [0 0 12; 0 0 6; 9 8 0];      % 1<->2 unseen both ways
[Q, info] = ruleOfThreeMatrix(P, N);
ok = Q(1, 2) == 0 && Q(2, 1) == 0 && isequal(Q, P) && ~any(info.filled(:)) && max(abs(sum(Q, 2) - 1)) < 1e-15;
[nP, nF] = rep(ok, nP, nF, 3, 'both-way-unseen pair stays 0; nothing else touched', mat2str(Q, 4));

%% 4
Pm = [0.9 0 0 0.1; 0.05 0.9 0.05 0; 0 0.1 0.8 0.1; 0.05 0 0.05 0.9];   % 1 -> 2 one-way; 1 <-> 3 and 2 <-> 4 unseen both ways
Nm = [0 0 0 8; 4 0 6 0; 0 5 0 7; 3 0 2 0];
[Qm, im] = ruleOfThreeMatrix(Pm, Nm);
two = [0.85 0 0 0.05 0.1; 0.05 0.9 0 0 0.05; 0.05 0.05 0.85 0.05 0; 0.04 0.05 0.05 0.86 0; 0.1 0 0 0.05 0.85];
Nt = round(40 * two); Nt(2, 1) = 2; Nt(3, 1) = 5;          % row 1 counts [34 0 0 2 4]: edges to 2 and 3 one-way
[Qt, it] = ruleOfThreeMatrix(two, Nt); pr = [3 5 1 4 2]; Qp = ruleOfThreeMatrix(two(pr, pr), Nt(pr, pr));
ok = it.filled(1, 2) && it.filled(1, 3) && it.added(1, 2) == 2 && it.added(1, 3) == 3 && max(abs(Qt(1, :) - [34 2 3 2 4] / 45)) < 1e-14 ...
    && max(abs(Qp - Qt(pr, pr)), [], 'all') < 1e-15 ...
    && im.filled(1, 2) && ~im.filled(1, 3) && Qm(1, 3) == 0 && max(abs(sum(Qt, 2) - 1)) < 1e-14;
[nP, nF] = rep(ok, nP, nF, 4, 'two one-way edges filled together against the original row; relabelling permutes the answer', mat2str(Qt(1, :), 4));

%% 5
rs = RandStream('twister', 'Seed', 3); ok = true;
for r = 1:50
    C = randi(rs, [0 30], 3); C(~eye(3) & rand(rs, 3) < 0.25) = 0; C(1:4:end) = C(1:4:end) + 5;
    P = C ./ sum(C, 2); [Q, info] = ruleOfThreeMatrix(P, C);
    D = C + info.added; one = sum(C .* ~eye(3), 2) > 0;   % every row with an observed exit
    ok = ok && max(abs(Q(one, :) - D(one, :) ./ sum(D(one, :), 2)), [], 'all') < 1e-12;
end
[nP, nF] = rep(ok, nP, nF, 5, 'empirical P: n_a = actual row steps, row = (N + added) normalised for every row with an exit (50 matrices)', '');

%% 6
P = [0.9 0 0.1; 0.05 0.9 0.05; 0.1 0.1 0.8]; N = [0 0 20; 0 0 4; 9 9 0];   % 2->1 seen in P, but N(2,1) = 0
[Q, info] = ruleOfThreeMatrix(P, N); E = modelEnergeticsR3(P, N);
Pa = [1 0 0; 0.1 0.8 0.1; 0.1 0.1 0.8]; Na = [0 0 0; 3 0 3; 3 3 0]; [Qa, ia] = ruleOfThreeMatrix(Pa, Na);
ok = info.stillOneWay(1, 2) && isequal(Q, P) && isinf(E.EPR) && ~E.noCycle && isequal(Qa, Pa) && ia.noExit(1) && ia.stillOneWay(1, 2) && ia.stillOneWay(1, 3) && ~ia.changed;
[nP, nF] = rep(ok, nP, nF, 6, 'N(b,a) = 0 left one-way (EPR Inf, J ~= 0); row without an exit left unchanged, reported', '');

%% 7
P = [0.90 0.00 0.10; 0.08 0.90 0.02; 0.00 0.06 0.94]; N = [0 0 15; 12 0 3; 0 9 0];   % forward W->Q->S->W only on 2 edges
E = modelEnergeticsR3(P, N); E0 = modelEnergetics(P);
ok = isinf(E0.EPR) && isfinite(E.EPR) && E.EPR > 0 && isfinite(E.A) && E.r3 && isequal(E.Pr3, ruleOfThreeMatrix(P, N)) && E.J > 0 && E.A > 0;
[nP, nF] = rep(ok, nP, nF, 7, 'modelEnergeticsR3: finite EPR > 0 and A where modelEnergetics gave Inf; J > 0 for a forward cycle', ...
    sprintf('EPR %g, A %g, J %g', E.EPR, E.A, E.J));

%% 8
cnt = 0;
try, ruleOfThreeMatrix(ones(3, 2) / 2, ones(3, 2)); catch, cnt = cnt + 1; end
try, ruleOfThreeMatrix(eye(3), ones(4)); catch, cnt = cnt + 1; end
[nP, nF] = rep(cnt == 2, nP, nF, 8, 'size checks', sprintf('%d of 2 refused', cnt));

%% 9-10. wiring
try
    rs = RandStream('twister', 'Seed', 11); T = 2880; ph = minutePhases(T);
    P = repmat([0.90 0 0.10; 0.08 0.90 0.02; 0.02 0.06 0.92], 1, 1, 2);       % Q -> W never: one-way edge
    [~, y] = hmmSimulatePhased([1 1 1]' / 3, P, eye(3), ph, rs);
    [~, ~, Pc] = chainFitPhased(y, ph, 2); N = observablePairCounts(y, ph, 2);
    Q = observableQuantities(y, ph, 2, Pc); Q0 = observableQuantities(y, ph, 2, Pc, RuleOfThree=false); ok = true;
    for p = 1:2
        e = modelEnergeticsR3(Pc(:, :, p), N(:, :, p)); L = energyLandscape(e.Pr3); e0 = modelEnergetics(Pc(:, :, p));
        ok = ok && Q(p).A_chain == e.A && Q(p).EPR_chain == e.EPR && Q(p).E21 == L.E2_minus_E1 && Q(p).E31 == L.E3_minus_E1 ...
            && Q(p).r3_chain && Q(p).J_chain == e0.J && isfinite(Q(p).EPR_chain) && isinf(Q0(p).EPR_chain) && ~Q0(p).r3_chain ...
            && Q0(p).A_chain == e0.A;
    end
    [nP, nF] = rep(ok, nP, nF, 9, 'observableQuantities: r3 for A, EPR, landscape; J_chain fitted; RuleOfThree=false = old', '');
    f = struct('pi', [1 1 1]' / 3, 'P', P, 'B', eye(3), 'deg', false);
    H = hmmPosteriorFluxEPR(y, ph, 2, f, NPaths=5, RandomSeed=1); H0 = hmmPosteriorFluxEPR(y, ph, 2, f, NPaths=5, RandomSeed=1, RuleOfThree=false);
    ok = true;
    for p = 1:2
        e = modelEnergeticsR3(P(:, :, p), N(:, :, p));
        ok = ok && abs(H(p).EPR_hmm_model - e.EPR) < 1e-12 && H(p).r3_hmm && isinf(H0(p).EPR_hmm_model) && H(p).J_hmm_model == H0(p).J_hmm_model;
    end
    [nP, nF] = rep(ok, nP, nF, 10, 'hmmPosteriorFluxEPR: EPR_hmm_model = r3 on mean path counts; RuleOfThree=false = old Inf', '');
catch err, [nP, nF] = rep(false, nP, nF, 9, 'wiring', err.message); end

%% 11. no-cycle convention
Pa = [1 0 0; 0.1 0.8 0.1; 0.1 0.1 0.8]; Na = [0 0 0; 3 0 3; 3 3 0]; Ea = modelEnergeticsR3(Pa, Na); E0a = modelEnergetics(Pa);
Ps = [0.8 0.1 0.1; 0.1 0.8 0.1; 0.1 0.1 0.8]; Es = modelEnergeticsR3(Ps, 5 * ones(3)); E0s = modelEnergetics(Ps);
ok = isinf(E0a.EPR) && abs(Ea.J) < 1e-12 && Ea.EPR == 0 && isnan(Ea.A) && Ea.noCycle ...
    && ~Es.noCycle && Es.EPR == E0s.EPR && Es.A == E0s.A && Es.J == E0s.J && isfinite(Es.EPR);
[nP, nF] = rep(ok, nP, nF, 11, 'no cycle: absorbing chain -> EPR 0, A NaN; reversible round-off left as computed', ...
    sprintf('absorbing EPR %g A %g; reversible EPR %g J %g', Ea.EPR, Ea.A, Es.EPR, Es.J));

fprintf('\n%d passed, %d failed\n', nP, nF);
if nF > 0, error('test_sqw_r3:failed', '%d test(s) failed.', nF); end
end

function [nP, nF] = rep(ok, nP, nF, k, what, why)
if ok
    nP = nP + 1; fprintf('[PASS] %2d  %s\n', k, what);
else
    nF = nF + 1; fprintf('[FAIL] %2d  %s\n        %s\n', k, what, why);
end
end
