function test_BoutIrreversibility(varargin)
%TEST_BOUTIRREVERSIBILITY  Suite for BoutIrreversibility.
%
%   test_BoutIrreversibility            all tests
%   test_BoutIrreversibility('fast')    skip the calibration tests (9-12)
%
% Ground truth that shares no code with the function: hand-computed
% statistics, the exact antisymmetry under time reversal, and simulated
% processes whose reversibility is known by construction --
%   reversible:   a Gaussian AR(1) mapped through a monotone transform to
%                 exponential bouts (any monotone map of a reversible
%                 process is reversible);
%   irreversible: a linear AR(1) with EXPONENTIAL innovations (Weiss 1975:
%                 a linear process with non-Gaussian innovations is not
%                 reversible) -- jumps up, decays down.
% Claude-written, like the code it tests: passing means self-consistency
% plus agreement with those constructions, not an independent review.

fast = nargin > 0 && strcmpi(varargin{1}, 'fast');
assert(isequal(fileparts(which('BoutIrreversibility')), pwd), ...
    'Run from the folder holding BoutIrreversibility.m (resolution check).');
t0 = tic;
res = {};
res{end+1} = run_test('1: hand example (lag 1)', @t1);
res{end+1} = run_test('2: time reversal negates skew, flips UpFraction, keeps Rho', @t2);
res{end+1} = run_test('3: pairs never span blocks', @t3);
res{end+1} = run_test('4: input forms agree; interior NaN errors', @t4);
res{end+1} = run_test('5: deterministic, global RNG untouched, nulls independent', @t5);
res{end+1} = run_test('6: IAAFT keeps values and memory; shuffle kills memory', @t6);
res{end+1} = run_test('7: one-fly group p equals the per-fly p', @t7);
res{end+1} = run_test('8: MinPairs drops a fly from the group', @t8);
if ~fast
    res{end+1} = run_test('9: calibration, no memory (iid exponential)', @t9);
    res{end+1} = run_test('10: calibration, reversible WITH memory', @t10);
    res{end+1} = run_test('11: power on an irreversible process', @t11);
    res{end+1} = run_test('12: group calibration, reversible with memory', @t12);
end
res{end+1} = run_test('13: plot is light-themed', @t13);
ok = cellfun(@(r) r, res);
fprintf('\n=== SUMMARY ===\n  Passed: %d\n  Failed: %d\n  Total : %d\n  %.0f s\n', ...
    sum(ok), sum(~ok), numel(ok), toc(t0));
end

% ---------------------------------------------------------------------
function ok = run_test(name, fn)
try
    fn();
    fprintf('  PASS  %s\n', name);
    ok = true;
catch err
    fprintf('  FAIL  %s\n        %s\n', name, err.message);
    ok = false;
end
end

function close_to(a, b, tol, what)
if ~(abs(a - b) <= tol)
    error('%s: got %.10g, expected %.10g', what, a, b);
end
end

% ---------------------------------------------------------------------
function t1
x = [1 3 2 6 4];
R = BoutIrreversibility(x, Lags=1, NumSurrogates=19, MinPairs=2);
T = R.PerFly;
assert(T.NPairs == 4);
% D = [2 -1 4 -2]; sorted [-2 -1 2 4], Fm = [.125 .375 .625 .875]
% Q1 = -1.5, Q2 = 0.5, Q3 = 3  ->  (3 - 1.5 - 1)/4.5
close_to(T.Bowley, 0.5/4.5, 1e-12, 'Bowley');
d = [2 -1 4 -2]; dc = d - mean(d);
close_to(T.Moment, mean(dc.^3)/mean(dc.^2)^1.5, 1e-12, 'Moment');
close_to(T.UpFraction, 0.5, 0, 'UpFraction');
% ranks of [1 3 2 6] and [3 2 6 4] are [1 3 2 4] and [2 1 4 3]: rho = 0
close_to(T.Rho, 0, 1e-12, 'Rho');
end

function t2
rs = RandStream('mt19937ar', 'Seed', 2);
blocks = {round(300 + 600*rs.rand(1, 80).^3), round(300 + 900*rs.rand(1, 60))};
A = BoutIrreversibility({blocks}, NumSurrogates=9, MinPairs=2);
rev = cellfun(@fliplr, blocks, 'UniformOutput', false);
B = BoutIrreversibility({rev}, NumSurrogates=9, MinPairs=2);
for l = 1:height(A.PerFly)
    close_to(B.PerFly.Bowley(l), -A.PerFly.Bowley(l), 1e-12, 'Bowley');
    close_to(B.PerFly.Moment(l), -A.PerFly.Moment(l), 1e-12, 'Moment');
    close_to(B.PerFly.UpFraction(l), 1 - A.PerFly.UpFraction(l), 1e-12, 'UpFraction');
    close_to(B.PerFly.Rho(l), A.PerFly.Rho(l), 1e-12, 'Rho');
end
end

function t3
R = BoutIrreversibility({{[1 2 3], [100 101 102]}}, Lags=[1 2], ...
    NumSurrogates=9, MinPairs=1);
assert(isequal(R.PerFly.NPairs, [4; 2]), 'pair counts');
assert(R.PerFly.UpFraction(1) == 1, 'UpFraction');
% all lag-1 differences are exactly 1 -- only possible if 3 -> 100 was not
% paired -- so both skewnesses are undefined
assert(isnan(R.PerFly.Moment(1)) && isnan(R.PerFly.Bowley(1)), 'constant D');
% as ONE block the 97 enters and D is no longer constant
R2 = BoutIrreversibility({[1 2 3 100 101 102]}, Lags=1, NumSurrogates=9, MinPairs=1);
assert(R2.PerFly.NPairs == 5 && isfinite(R2.PerFly.Moment), 'one block spans the gap');
end

function t4
rs = RandStream('mt19937ar', 'Seed', 4);
a = round(300 + 500*rs.rand(1, 50));
b = round(300 + 500*rs.rand(1, 40));
M = [a; [b NaN(1, 10)]];
o = {'Lags', [1 2], 'NumSurrogates', 19, 'MinPairs', 5};
R1 = BoutIrreversibility(M, o{:});
R2 = BoutIrreversibility({a, b}, o{:});
R3 = BoutIrreversibility({{a}, {b}}, o{:});
assert(isequaln(R1.PerFly{:, 2:end}, R2.PerFly{:, 2:end}), 'matrix vs cell');
assert(isequaln(R2.PerFly{:, 2:end}, R3.PerFly{:, 2:end}), 'cell vs cell of cells');
R4 = BoutIrreversibility(a(:), o{:});
assert(height(R4.PerFly) == 2, 'column vector is ONE fly');
M(1, 20) = NaN;
try
    BoutIrreversibility(M, o{:});
    error('no error on interior NaN');
catch err
    assert(strcmp(err.identifier, 'BoutIrreversibility:InteriorNaN'), err.message);
end
end

function t5
rs = RandStream('mt19937ar', 'Seed', 5);
x = round(300 + 600*rs.rand(1, 120));
rng(123); before = rand(1, 3); rng(123);
A = BoutIrreversibility(x, NumSurrogates=49, MinPairs=5);
after = rand(1, 3);
assert(isequal(before, after), 'global stream was used');
B = BoutIrreversibility(x, NumSurrogates=49, MinPairs=5);
assert(isequaln(A.PerFly, B.PerFly), 'not deterministic');
C = BoutIrreversibility(x, NumSurrogates=49, MinPairs=5, Null="iaaft");
assert(isequaln(A.PerFly.Bowley_P_IAAFT, C.PerFly.Bowley_P_IAAFT), ...
    'IAAFT p depends on whether the shuffle null ran');
assert(all(isnan(C.PerFly.Bowley_P_Shuffle)), 'skipped null not NaN');
end

function t6
x = reversible_ar(0.6, 200, 6);
R = BoutIrreversibility(x, Lags=1, NumSurrogates=200, MinPairs=5, ReturnNull=true);
rhoObs = R.PerFly.Rho;
rI = R.PerFlyNull.iaaft.Rho(:);
rS = R.PerFlyNull.shuffle.Rho(:);
fprintf('        observed Rho %.3f; IAAFT median %.3f; shuffle median %.3f\n', ...
    rhoObs, median(rI), median(rS));
assert(rhoObs > 0.35, 'simulated memory missing');
assert(abs(median(rI) - rhoObs) < 0.12, 'IAAFT lost the memory');
assert(abs(median(rS)) < 0.05, 'shuffle kept memory');
end

function t7
x = reversible_ar(0.3, 150, 7);
R = BoutIrreversibility(x, NumSurrogates=99, MinPairs=5);
assert(~any(startsWith(R.PerFly.Properties.VariableNames, 'Moment_P')), ...
    'Moment must carry no p-value');
for s = ["Bowley","UpFraction"]
    for n = ["IAAFT","Shuffle"]
        assert(isequaln(R.Group.(s + "_P_" + n), R.PerFly.(s + "_P_" + n)), s + n);
    end
end
assert(isequaln(R.Group.Rho_P_Shuffle, R.PerFly.Rho_P_Shuffle));
end

function t8
a = reversible_ar(0.3, 100, 8);
b = reversible_ar(0.3, 12, 9);
R = BoutIrreversibility({a, b}, Lags=[1 8], NumSurrogates=19, MinPairs=6);
% fly 2: lag 1 has 11 pairs (kept), lag 8 has 4 (dropped)
assert(isequal(R.Group.NFlies, [2; 1]), 'NFlies');
assert(isnan(R.PerFly.Bowley(4)) && isfinite(R.PerFly.Bowley(3)), 'NaN rows');
close_to(R.Group.Bowley(2), R.PerFly.Bowley(2), 1e-12, 'group mean over kept flies');
end

% ---------------------------------------------------------------------
% Calibration: rejection rates at 0.05, per null and statistic, lags 1-8.
function t9
% rates: lag x [Bowley UpFraction Rho] x [IAAFT shuffle]. 0.095 is about
% 0.05 + 3 binomial SDs at 200 reps.
rates = calib(@(i) iid_exp(150, 900 + i), 200, 199, 'iid exponential, n = 150');
assert(all(rates(:, 1:2, :) <= 0.095, 'all'), 'a skew null over-rejects with no memory');
assert(all(rates(:, 3, 2) <= 0.095), 'memory test over-rejects with no memory');
end

function t10
rates = calib(@(i) reversible_ar(0.5, 150, 1000 + i), 200, 199, ...
    'reversible, latent AR phi = 0.5, n = 150');
iaaft = rates(:, 1:2, 1);
assert(all(iaaft(:) <= 0.095), 'IAAFT over-rejects on a reversible process');
assert(rates(1, 3, 2) > 0.9, 'memory test misses rho ~0.48');
end

function t11
rates = calib(@(i) irreversible_ar(0.5, 150, 2000 + i), 200, 199, ...
    'IRREVERSIBLE, AR phi = 0.5 with exponential innovations, n = 150');
assert(rates(1, 2, 1) > 0.7, 'IAAFT UpFraction has little power at lag 1');
end

function t12
reps = 100; nf = 10; R0 = 199;
rej = zeros(1, 2);
for i = 1:reps
    flies = arrayfun(@(f) reversible_ar(0.5, 80, 5000 + 100*i + f), 1:nf, ...
        'UniformOutput', false);
    R = BoutIrreversibility(flies, Lags=1, NumSurrogates=R0, MinPairs=10, Null="iaaft");
    rej = rej + ([R.Group.Bowley_P_IAAFT, R.Group.UpFraction_P_IAAFT] < 0.05);
end
rate = rej / reps;
fprintf('        group (10 flies x 80 bouts, %d groups) IAAFT rejection: Bowley %.3f Up %.3f\n', ...
    reps, rate);
assert(all(rate <= 0.11), 'group test over-rejects');
assert(any(rej > 0), 'group test never rejects in 100 null groups: check it');
end

function t13
x = reversible_ar(0.3, 100, 13);
R = BoutIrreversibility({x, x(end:-1:1)}, NumSurrogates=19, MinPairs=5, Plot=true);
fh = R.Figure;
c = onCleanup(@() close(fh));
assert(isequal(fh.Color, [1 1 1]), 'figure not white');
ax = findobj(fh, 'Type', 'axes');
assert(numel(ax) == 3, 'three panels');
assert(all(arrayfun(@(a) all(a.Color > 0.9), ax)), 'axes not light');
end

% ---------------------------------------------------------------------
function rates = calib(gen, reps, R0, label)
lags = [1 2 4 8];
rej = zeros(4, 3, 2);   % lag x [Bowley UpFraction Rho] x [IAAFT shuffle]
for i = 1:reps
    R = BoutIrreversibility(gen(i), Lags=lags, NumSurrogates=R0, MinPairs=10);
    T = R.PerFly;
    P = cat(3, [T.Bowley_P_IAAFT, T.UpFraction_P_IAAFT, NaN(4, 1)], ...
               [T.Bowley_P_Shuffle, T.UpFraction_P_Shuffle, T.Rho_P_Shuffle]);
    rej = rej + (P < 0.05);
end
rates = rej / reps;
fprintf('        %s, %d flies, %d surrogates; rejection rate at 0.05\n', label, reps, R0);
fprintf('        lag | Bowley IAAFT  shuffle | Up IAAFT  shuffle | Rho shuffle\n');
for l = 1:4
    fprintf('        %3d |       %.3f    %.3f |   %.3f    %.3f |       %.3f\n', ...
        lags(l), rates(l, 1, 1), rates(l, 1, 2), rates(l, 2, 1), rates(l, 2, 2), ...
        rates(l, 3, 2));
end
end

function x = iid_exp(n, seed)
rs = RandStream('mt19937ar', 'Seed', seed);
x = 300 + round(-600 * log(rs.rand(n, 1)));
end

function x = reversible_ar(phi, n, seed)
% Gaussian AR(1), stationary start, mapped to exponential bouts (mean 600 s
% above a 300 s floor) through its own CDF: reversible by construction.
rs = RandStream('mt19937ar', 'Seed', seed);
z = zeros(n, 1);
z(1) = rs.randn;
for t = 2:n
    z(t) = phi * z(t-1) + sqrt(1 - phi^2) * rs.randn;
end
upper = 0.5 * erfc(z / sqrt(2));          % 1 - Phi(z), accurate in the tail
x = 300 + round(-600 * log(upper));
end

function x = irreversible_ar(phi, n, seed)
% Linear AR(1) with exponential innovations: big jumps up, geometric decay.
% Irreversible (Weiss 1975). Burn-in of 200.
rs = RandStream('mt19937ar', 'Seed', seed);
m = n + 200;
y = zeros(m, 1);
for t = 2:m
    y(t) = phi * y(t-1) - log(rs.rand);
end
y = y(201:end);
x = 300 + round(600 * (1 - phi) * y);
end
