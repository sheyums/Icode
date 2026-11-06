function test_BoutAR1Noise(varargin)
%TEST_BOUTAR1NOISE  Suite for BoutAR1Noise.
%
%   test_BoutAR1Noise            all tests
%   test_BoutAR1Noise('fast')    skip the calibration tests (6-9)
%
% Ground truth sharing no code with the function: polyfit for the line,
% kstest2 for the KS distance (Statistics Toolbox; skipped if absent), and
% simulated processes whose properties are known by construction:
%   iid exponential bouts   B = 0, |eps| independent of x(t);
%   EAR(1)                  linear AR(1) with iid noise: B > 0, |eps|
%                           independent of x(t), y(t+1) >= B*y(t) always;
%   Gaussian-copula AR(1)   AR(1) on a latent Gaussian mapped to exponential
%                           bouts (the transformed-AR model class): B > 0 and
%                           the spread of eps grows with x(t);
%   LD-like, no memory      bouts iid GIVEN time of day (profile only).
% Claude-written, like the code it tests.

fast = nargin > 0 && strcmpi(varargin{1}, 'fast');
assert(isequal(fileparts(which('BoutAR1Noise')), pwd), ...
    'Run from the folder holding BoutAR1Noise.m (resolution check).');
t0 = tic;
res = {};
res{end+1} = run_test('1: A, B and residuals equal polyfit', @t1);
res{end+1} = run_test('2: KS distance equals kstest2', @t2);
res{end+1} = run_test('3: EAR(1) has no pair below half of B; the copula AR(1) does', @t3);
res{end+1} = run_test('4: Xmin required; blocks respected; deterministic; global RNG untouched', @t4);
res{end+1} = run_test('5: strata input checked', @t5);
res{end+1} = run_test('10: regeneration counts bouts below the floor, never clips', @t10);
if ~fast
    res{end+1} = run_test('6: calibration, iid exponential bouts', @t6);
    res{end+1} = run_test('7: EAR(1), B = 0.3: B found, spread flat', @t7);
    res{end+1} = run_test('8: Gaussian-copula AR(1): spread grows with x(t)', @t8);
    res{end+1} = run_test('9: LD-like, no memory: strata null calibrated, shuffle not', @t9);
end
ok = cellfun(@(r) r, res);
fprintf('\n=== SUMMARY ===\n  Passed: %d\n  Failed: %d\n  Total : %d\n  %.0f s\n', ...
    sum(ok), sum(~ok), numel(ok), toc(t0));
end

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

% ---------------------------------------------------------------------
function t1
x = iid_exp(120, 1);
R = BoutAR1Noise(x, Xmin=300, NumSurrogates=19, ReturnResiduals=true);
pf = polyfit(x(1:end-1), x(2:end), 1);
T = R.PerFly;
assert(abs(T.B - pf(1)) < 1e-10 && abs(T.A - pf(2)) < 1e-7, 'line');
e = x(2:end) - polyval(pf, x(1:end-1));
assert(max(abs(R.Residuals{1}(:, 3) - e)) < 1e-7, 'residuals');
assert(T.NPairs == 119 && abs(T.MeanBout - mean(x)) < 1e-9);
end

function t2
if isempty(which('kstest2')), fprintf('        (kstest2 absent: skipped)\n'); return, end
src = fileread(which('BoutAR1Noise'));
d = fullfile(tempdir, 'ar1noise_ks'); if ~isfolder(d), mkdir(d); end
put(d, src, 'ks_dist'); put(d, src, 'count_le');
addpath(d); c = onCleanup(@() rmpath(d));
rs = RandStream('mt19937ar', 'Seed', 2);
for trial = 1:20
    x = round(rs.randn(50 + trial, 1) * 3);          % heavy ties
    y = round(rs.randn(300, 1) * 3 + 0.5);
    [~, ~, k] = kstest2(x, y);
    assert(abs(ks_dist(x, sort(y)) - k) < 1e-12, 'trial %d', trial);
end
end

function put(d, src, name)
% Copy one local function of the source into its own file.
starts = regexp(src, '(?m)^function ');
k = find(~cellfun(@isempty, regexp(arrayfun(@(s) src(s:min(end, s + 80)), starts, ...
    'UniformOutput', false), ['\<' name '\('], 'once')), 1);
stop = numel(src);
if k < numel(starts), stop = starts(k + 1) - 1; end
fid = fopen(fullfile(d, [name '.m']), 'w'); fwrite(fid, src(starts(k):stop)); fclose(fid);
end

function t3
y = ear1(0.4, 600, 3) - 300;
c = y(2:end); a = y(1:end-1);
assert(all(c >= 0.4 * a - 1), 'EAR(1) violates y(t+1) >= B*y(t) beyond rounding');
R = BoutAR1Noise(y + 300, Xmin=300, NumSurrogates=99);
fprintf('        EAR(1) B = 0.4, n = 600: fitted B %.3f, BelowHalfB %.3f (shuffled: %.3f)\n', ...
    R.PerFly.B, R.PerFly.BelowHalfB, R.PerFly.BelowHalfB_NullMean);
assert(R.PerFly.BelowHalfB == 0, 'EAR(1) data put pairs below half of B');
% The transformed-Gaussian AR(1), with a similar fitted B, does put pairs
% there: that is what makes the check discriminate. (Shuffled data are no
% contrast: their Bhat is ~0, so the threshold is ~0.)
Rc = BoutAR1Noise(copula_ar(0.5, 600, 3), Xmin=300, NumSurrogates=19);
fprintf('        Gaussian-copula AR(1) phi 0.5: fitted B %.3f, BelowHalfB %.3f\n', ...
    Rc.PerFly.B, Rc.PerFly.BelowHalfB);
assert(Rc.PerFly.BelowHalfB > 0.05, 'copula AR(1) should put pairs below half of B');
end

function t4
x = iid_exp(100, 4);
try
    BoutAR1Noise(x, NumSurrogates=9);
    error('no error without Xmin');
catch err
    assert(strcmp(err.identifier, 'BoutAR1Noise:Xmin'), err.message);
end
R = BoutAR1Noise({{x(1:50), x(51:100)}}, Xmin=300, NumSurrogates=9);
assert(R.PerFly.NPairs == 98, 'a pair spanned two blocks');
rng(7); before = rand(1, 2); rng(7);
A = BoutAR1Noise(x, Xmin=300, NumSurrogates=49);
assert(isequal(before, rand(1, 2)), 'global stream used');
B = BoutAR1Noise(x, Xmin=300, NumSurrogates=49);
assert(isequaln(A.PerFly, B.PerFly), 'not deterministic');
end

function t5
x = iid_exp(100, 5);
try
    BoutAR1Noise({x}, Xmin=300, Strata={ones(99, 1)}, NumSurrogates=9);
    error('no error');
catch err
    assert(strcmp(err.identifier, 'BoutAR1Noise:Strata'), err.message);
end
R = BoutAR1Noise({x}, Xmin=300, Strata={(1:100).'}, NumSurrogates=19);
% every bout alone in its stratum: surrogates equal the data
assert(R.PerFly.B_P == 1 && R.PerFly.Hetero_P == 1, 'immovable strata must give p = 1');
end

function t10
% Independent bouts: Bhat ~ 0, so x(t+1) ~ A + eps* = a resampled bout and
% almost nothing falls below the floor. EAR(1): its innovations are
% non-negative, y(t+1) - B*y(t) >= 0, so the resampled residuals keep the
% floor BY THEMSELVES -- measured 0 below 300 s (the first version of this
% test wrongly expected memory alone to break the floor). The copula AR(1)
% with a similar B does break it: that contrast is the check.
xi = iid_exp(400, 10);
Ri = BoutAR1Noise(xi, Xmin=300, NumSurrogates=19);
Re = BoutAR1Noise(ear1(0.4, 400, 11), Xmin=300, NumSurrogates=19);
Rc = BoutAR1Noise(copula_ar(0.5, 400, 12), Xmin=300, NumSurrogates=19);
fprintf(['        below 300 s -- iid: %.4f (B %.3f), EAR(1): %.4f (B %.3f), copula: %.4f (B %.3f); ' ...
    'GenKS/ref iid %.3f/%.3f\n'], Ri.PerFly.GenBelowXmin, Ri.PerFly.B, Re.PerFly.GenBelowXmin, ...
    Re.PerFly.B, Rc.PerFly.GenBelowXmin, Rc.PerFly.B, Ri.PerFly.GenKS, Ri.PerFly.GenKS_Ref);
assert(Ri.PerFly.GenBelowXmin < 0.02, 'iid bouts regenerate below the floor');
assert(Re.PerFly.GenBelowXmin < 0.002, 'EAR(1) residuals should keep the floor');
assert(Rc.PerFly.GenBelowXmin > 0.03, 'the copula AR(1) should break the floor');
assert(Rc.PerFly.GenMinBout < 300, 'nothing below 300 would mean clipping');
% the floor is Xmin, not 300: wake-like bouts with a 2 s floor
w = 2 + round(-60 * log(RandStream('mt19937ar', 'Seed', 13).rand(300, 1)));
Rw = BoutAR1Noise(w, Xmin=2, NumSurrogates=19);
Rw300 = BoutAR1Noise(w, Xmin=300, NumSurrogates=19);
fprintf('        wake-like (mean 60 s): below 2 s %.4f (B %.3f); same bouts against 300 s: %.4f\n', ...
    Rw.PerFly.GenBelowXmin, Rw.PerFly.B, Rw300.PerFly.GenBelowXmin);
assert(Rw300.PerFly.GenBelowXmin > 0.5 && Rw.PerFly.GenBelowXmin < Rw300.PerFly.GenBelowXmin, ...
    'the floor used is not Xmin');
end

% ---------------------------------------------------------------------
function t6
r = calib(@(i) iid_exp(150, 100 + i), 200, 199, 'iid exponential, n = 150', []);
assert(all(r <= 0.095), 'a test over-rejects on iid bouts');
end

function t7
r = calib(@(i) ear1(0.3, 150, 200 + i), 200, 199, 'EAR(1) B = 0.3, n = 150', []);
assert(r(1) > 0.5, 'B missed');
assert(r(2) <= 0.095, 'flat spread flagged as growing');
end

function t8
r = calib(@(i) copula_ar(0.5, 150, 300 + i), 200, 199, 'Gaussian-copula AR(1) phi = 0.5, n = 150', []);
assert(r(1) > 0.5, 'B missed');
assert(r(2) > 0.3, 'growing spread not detected');
end

function t9
[rS, ~] = calib(@(i) ld_fly(400 + i), 200, 199, 'LD-like, no memory, SHUFFLE null', [], true);
[rT, ~] = calib(@(i) ld_fly(400 + i), 200, 199, 'LD-like, no memory, STRATA null', 'strata', true);
assert(all(rT <= 0.095), 'strata null over-rejects with the profile alone');
fprintf('        (shuffle null on the same flies, B: %.3f -- the profile alone)\n', rS(1));
end

function [rates, B] = calib(gen, reps, R0, label, mode, ld)
if nargin < 6, ld = false; end
rej = zeros(1, 3);
B = NaN(reps, 1);
for i = 1:reps
    if ld
        [x, st] = gen(i);
        if isequal(mode, 'strata')
            R = BoutAR1Noise({x}, Xmin=300, Strata={st}, NumSurrogates=R0);
        else
            R = BoutAR1Noise({x}, Xmin=300, NumSurrogates=R0);
        end
    else
        R = BoutAR1Noise(gen(i), Xmin=300, NumSurrogates=R0);
    end
    T = R.PerFly;
    B(i) = T.B;
    rej = rej + ([T.B_P, T.Hetero_P, T.ShapeKS_P] < 0.05);
end
rates = rej / reps;
fprintf('        %s, %d flies, %d surrogates: rejection B %.3f, Hetero %.3f, ShapeKS %.3f; median B %.3f\n', ...
    label, reps, R0, rates, median(B));
end

% ---------------------------------------------------------------------
function x = iid_exp(n, seed)
rs = RandStream('mt19937ar', 'Seed', seed);
x = 300 + round(-600 * log(rs.rand(n, 1)));
end

function x = ear1(B, n, seed)
% EAR(1) (Gaver-Lewis): y(t) = B*y(t-1) + (0 w.p. B, else Exp(mean m)).
rs = RandStream('mt19937ar', 'Seed', seed);
m = 600;
y = zeros(n + 200, 1);
y(1) = -m * log(rs.rand);
for t = 2:numel(y)
    y(t) = B * y(t-1);
    if rs.rand >= B, y(t) = y(t) - m * log(rs.rand); end
end
x = 300 + round(y(201:end));
end

function x = copula_ar(phi, n, seed)
rs = RandStream('mt19937ar', 'Seed', seed);
z = zeros(n, 1);
z(1) = rs.randn;
for t = 2:n
    z(t) = phi * z(t-1) + sqrt(1 - phi^2) * rs.randn;
end
x = 300 + round(-600 * log(0.5 * erfc(z / sqrt(2))));
end

function [sleep, strat] = ld_fly(seed)
% Alternating wake/sleep on a real clock from ZT0, 5 days; sleep means 15
% min by day, rising 20 -> 80 min through the night; NO memory.
rs = RandStream('mt19937ar', 'Seed', seed);
t = 0; sleep = []; strat = [];
while t < 5 * 86400
    zt = mod(t / 3600, 24);
    mw = 25 * (zt < 12) + 8 * (zt >= 12);
    t = t + 2 + round(60 * mw * -log(rs.rand));
    zt = mod(t / 3600, 24);
    if zt < 12, ms = 15; else, ms = 20 + 60 * (zt - 12) / 12; end
    d = 300 + round(60 * ms * -log(rs.rand));
    sleep(end+1, 1) = d; %#ok<AGROW>
    strat(end+1, 1) = floor(zt); %#ok<AGROW>
    t = t + d;
end
end
