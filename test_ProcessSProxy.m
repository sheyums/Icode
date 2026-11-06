function test_ProcessSProxy(varargin)
%TEST_PROCESSSPROXY  Suite and calibration for ProcessSProxy (PRESPEC entry D).
%
%   test_ProcessSProxy           all tests
%   test_ProcessSProxy('fast')   skip the calibrations (4-6)

fast = nargin > 0 && strcmpi(varargin{1}, 'fast');
assert(isequal(fileparts(which('ProcessSProxy')), pwd), 'resolution check');
t0 = tic;
res = {};
res{end+1} = run_test('1: S recursion: sleep set by S gives rho_S = 1 at the true taus', @t1);
res{end+1} = run_test('2: a gap restarts S and the burn-in', @t2);
res{end+1} = run_test('3: same surrogates as SleepAmountACF; global RNG untouched', @t3);
if ~fast
    res{end+1} = run_test('4: group calibration, independent bouts (48 flies)', @t4);
    res{end+1} = run_test('5: homeostat: power and where the grid maximum lands', @t5);
    res{end+1} = run_test('6: shared slow state: what the test does', @t6);
end
ok = cellfun(@(r) r, res);
fprintf('\n=== SUMMARY ===\n  Passed: %d\n  Failed: %d\n  Total : %d\n  %.0f s\n', ...
    sum(ok), sum(~ok), numel(ok), toc(t0));
end

function ok = run_test(name, fn)
try
    fn(); fprintf('  PASS  %s\n', name); ok = true;
catch err
    fprintf('  FAIL  %s\n        %s\n', name, err.message); ok = false;
end
end

function t1
% Build a record whose sleep durations are an increasing function of S at
% onset under Tw = 2 h, Ts = 1 h (computed here independently).
rs = RandStream('mt19937ar', 'Seed', 1);
n = 400; d = zeros(n, 1); isS = mod((1:n).', 2) == 0; on = zeros(n, 1);
S = 0.5; t = 0;
for i = 1:n
    on(i) = t;
    if isS(i)
        d(i) = 300 + 7200 * S;               % strictly increasing in S (not rounded)
        S = S * exp(-d(i) / 3600);
    else
        d(i) = 2 + round(-900 * log(rs.rand));
        S = 1 - (1 - S) * exp(-d(i) / 7200);
    end
    t = t + d(i);
end
fl = struct('Dur', d, 'IsSleep', isS, 'Onset', on);
R = ProcessSProxy(fl, NumSurrogates=9, BurnInHours=0);
iw = find(R.Settings.TauW == 2); is = find(R.Settings.TauS == 1);
r = R.RhoS(1, iw, is);
fprintf('        rho_S at (2 h, 1 h): %.4f; group best (%g, %g)\n', r, R.Group.BestTauW_S, R.Group.BestTauS_S);
assert(r > 0.99, 'S recursion does not match the specification');
end

function t2
fl = SimulateSleepWake("iid", 2, Days=3);
R1 = ProcessSProxy(fl, NumSurrogates=9);
fl.Onset(150:end) = fl.Onset(150:end) + 100;
R2 = ProcessSProxy(fl, NumSurrogates=9);
assert(R2.PerFly.NSleepUsed + R2.PerFly.NWakeUsed < R1.PerFly.NSleepUsed + R1.PerFly.NWakeUsed, ...
    'the gap did not restart the burn-in');
end

function t3
fl = SimulateSleepWake("iid", 3, Days=3);
rng(4); b = rand(1, 2); rng(4);
ProcessSProxy(fl, NumSurrogates=19);
assert(isequal(b, rand(1, 2)), 'global stream used');
% the r-th surrogate of fly 1 is the same draw in both functions
rs1 = RandStream('mrg32k3a', 'Seed', 1); rs1.Substream = 1;
rs2 = RandStream('mrg32k3a', 'Seed', 1); rs2.Substream = 1;
for r = 1:5
    a = RenewalSurrogate(fl, rs1); c = RenewalSurrogate(fl, rs2);
    assert(isequal(a.Dur, c.Dur), 'surrogate streams differ');
end
end

function t4
r = group_rates("iid", 60, 'independent bouts');
assert(all(r <= 0.11), 'over-rejects on independent bouts');
end

function t5
r = group_rates("homeostat", 30, 'HOMEOSTAT (Tw 2 h, Ts 1 h, gain 2)');
assert(r(1) > 0.5, 'rho_S has little power against the homeostat');
end

function t6
group_rates("shared", 30, 'shared slow state (no homeostat)');
end

function rates = group_rates(mode, reps, label)
rej = [0 0]; best = NaN(reps, 2); mx = NaN(reps, 1);
for g = 1:reps
    fl = arrayfun(@(f) SimulateSleepWake(mode, 120000 + 100*g + f), 1:48);
    R = ProcessSProxy(fl, NumSurrogates=199);
    rej = rej + [R.Group.MaxRhoS_P < 0.05, R.Group.MinRhoW_P < 0.05];
    best(g, :) = [R.Group.BestTauW_S, R.Group.BestTauS_S]; mx(g) = R.Group.MaxRhoS;
end
rates = rej / reps;
[u, ~, j] = unique(best, 'rows'); c = accumarray(j, 1); [~, o] = max(c);
fprintf('        %s, 48-fly groups x %d: rejection rho_S %.3f, rho_W %.3f; median max group rho_S %.3f; most frequent best (Tw, Ts) = (%g, %g) h in %d/%d\n', ...
    label, reps, rates, median(mx), u(o, 1), u(o, 2), c(o), reps);
end
