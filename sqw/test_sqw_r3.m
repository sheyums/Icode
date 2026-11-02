function test_sqw_r3()
%TEST_SQW_R3  Rule-of-three fill of one-way edges (RULEOFTHREEMATRIX, MODELENERGETICSR3): stress tests.
%
%   Run by name from the sqw folder:   cd <Icode>/sqw; test_sqw_r3
%
%    1. no one-way edge: Q is P bit for bit, energetics equal MODELENERGETICS
%    2. hand example: one unseen edge, n_a = 20 / 0.1 = 200, N(b,a) = 7 -> Q(a,:) = [180 3 20] / 203;
%       N(b,a) = 2 -> the added count is 2, not 3
%    3. pairs unseen both ways stay 0; rows without a one-way edge unchanged; every row sums to 1
%    4. two one-way edges in one row, both filled; relabelling the states permutes the answer
%       (no order dependence)
%    5. when P is the empirical matrix of integer counts N, n_a is row a's actual number of steps,
%       so Q = (N + added counts) normalised by row
%    6. edges that cannot be filled are reported (N(b,a) = 0; absorbing row P(a,a) = 1), the
%       EPR stays Inf there, nothing is silently floored
%    7. modelEnergeticsR3 on a one-way chain: EPR and A finite, EPR > 0, Pr3 = Q, r3 true, J sign
%       as the forward cycle
%    8. size checks: non-square P, N of another size refused

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
ok = max(abs(Q(1, :) - [180 3 20] / 203)) < 1e-15 && isequal(Q(2:3, :), P(2:3, :)) && abs(info.nSteps(1) - 200) < 1e-12 ...
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
two(1, 2:3) = 0; Nt = round(40 * two); Nt(2, 1) = 2; Nt(3, 1) = 5;          % row 1: edges to 2 and 3 one-way
[Qt, it] = ruleOfThreeMatrix(two, Nt); pr = [3 5 1 4 2]; Qp = ruleOfThreeMatrix(two(pr, pr), Nt(pr, pr));
ok = it.filled(1, 2) && it.filled(1, 3) && it.added(1, 2) == 2 && it.added(1, 3) == 3 && max(abs(Qp - Qt(pr, pr)), [], 'all') < 1e-15 ...
    && im.filled(1, 2) && ~im.filled(1, 3) && Qm(1, 3) == 0 && max(abs(sum(Qt, 2) - 1)) < 1e-14;
[nP, nF] = rep(ok, nP, nF, 4, 'two one-way edges in a row filled together; relabelling permutes the answer', mat2str(Qt(1, :), 4));

%% 5
rs = RandStream('twister', 'Seed', 3); ok = true;
for r = 1:50
    C = randi(rs, [0 30], 3); C(~eye(3) & rand(rs, 3) < 0.25) = 0; C(1:4:end) = C(1:4:end) + 5;
    P = C ./ sum(C, 2); [Q, info] = ruleOfThreeMatrix(P, C);
    D = C + info.added; ok = ok && max(abs(Q - D ./ sum(D, 2)), [], 'all') < 1e-12 ...
        && all(abs(info.nSteps(~isnan(info.nSteps)) - sum(C(~isnan(info.nSteps), :), 2)) < 1e-9);
end
[nP, nF] = rep(ok, nP, nF, 5, 'empirical P: n_a = actual row steps, Q = (N + added) by row (50 random count matrices)', '');

%% 6
P = [0.9 0 0.1; 0.05 0.9 0.05; 0.1 0.1 0.8]; N = [0 0 20; 0 0 4; 9 9 0];   % 2->1 seen in P, but N(2,1) = 0
[Q, info] = ruleOfThreeMatrix(P, N); E = modelEnergeticsR3(P, N);
Pa = [1 0 0; 0.1 0.8 0.1; 0.1 0.1 0.8]; Na = [0 0 0; 3 0 3; 3 3 0]; [Qa, ia] = ruleOfThreeMatrix(Pa, Na);
ok = info.stillOneWay(1, 2) && isequal(Q, P) && isinf(E.EPR) && ia.stillOneWay(1, 2) && ia.stillOneWay(1, 3) && isequal(Qa, Pa);
[nP, nF] = rep(ok, nP, nF, 6, 'unfillable edges reported (N(b,a) = 0; absorbing row), EPR stays Inf', '');

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
