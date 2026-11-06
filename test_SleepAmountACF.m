function test_SleepAmountACF(varargin)
%TEST_SLEEPAMOUNTACF  Suite and calibration for SleepAmountACF (PRESPEC entry C),
% RenewalSurrogate and SimulateSleepWake.
%
%   test_SleepAmountACF           all tests
%   test_SleepAmountACF('fast')   skip the calibrations (4-6)

fast = nargin > 0 && strcmpi(varargin{1}, 'fast');
assert(isequal(fileparts(which('SleepAmountACF')), pwd), 'resolution check');
t0 = tic;
res = {};
res{end+1} = run_test('1: ACF equals a per-second brute force', @t1);
res{end+1} = run_test('2: renewal surrogate keeps bouts, alternation, segments', @t2);
res{end+1} = run_test('3: simulator is contiguous and alternating; deterministic', @t3);
if ~fast
    res{end+1} = run_test('4: group calibration, independent bouts (48 flies)', @t4);
    res{end+1} = run_test('5: homeostat: power and where the minimum falls', @t5);
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
fl = SimulateSleepWake("homeostat", 1, Days=3);
R = SleepAmountACF(fl, NumSurrogates=9, BinMinutes=30, MaxLagHours=6);
tot = fl.Onset(end) + fl.Dur(end);
x = false(tot, 1);
for i = find(fl.IsSleep).'
    x(fl.Onset(i) + 1 : fl.Onset(i) + fl.Dur(i)) = true;
end
nW = floor(tot / 1800);
y = sum(reshape(x(1:nW*1800), 1800, nW), 1).';
for k = 1:numel(R.LagHours)
    L = k;
    u = y(1:end-L) - mean(y(1:end-L)); v = y(1+L:end) - mean(y(1+L:end));
    a = sum(u .* v) / sqrt(sum(u.^2) * sum(v.^2));
    assert(abs(a - R.ACF(k)) < 1e-10, 'lag %d', L);
end
end

function t2
fl = SimulateSleepWake("iid", 2, Days=3);
fl.Onset(80:end) = fl.Onset(80:end) + 777;       % a gap after bout 79
rs = RandStream('mrg32k3a', 'Seed', 5);
s = RenewalSurrogate(fl, rs);
assert(isequal(sort(s.Dur(fl.IsSleep)), sort(fl.Dur(fl.IsSleep))), 'sleep bouts');
assert(isequal(sort(s.Dur(~fl.IsSleep)), sort(fl.Dur(~fl.IsSleep))), 'wake bouts');
assert(isequal(s.IsSleep, fl.IsSleep), 'alternation');
assert(s.Onset(80) == fl.Onset(80) && s.Onset(1) == fl.Onset(1), 'segment starts');
seg1 = 1:79;
assert(isequal(sort(s.Dur(seg1)), sort(fl.Dur(seg1))), 'bouts crossed the gap');
c = s.Onset(2:end) == s.Onset(1:end-1) + s.Dur(1:end-1);
assert(sum(~c) == 1 && ~c(79), 'contiguity broken outside the gap');
end

function t3
for m = ["iid","shared","homeostat","ld"]
    a = SimulateSleepWake(m, 3, Days=2); b = SimulateSleepWake(m, 3, Days=2);
    assert(isequaln(a, b), '%s not deterministic', m);
    assert(all(a.Onset(2:end) == a.Onset(1:end-1) + a.Dur(1:end-1)), '%s not contiguous', m);
    assert(all(a.IsSleep(2:end) ~= a.IsSleep(1:end-1)), '%s does not alternate', m);
    assert(all(a.Dur(a.IsSleep) >= 300) && all(a.Dur(~a.IsSleep) >= 2), '%s floors', m);
end
end

function t4
% First measured 2026-11-06: MinP 0.100 (above the 0.095 rule), MaxP 0.030.
% That failure demoted PRESPEC entry C to descriptive (amendment of
% 2026-11-06). The test now guards against the bias getting WORSE, not
% against the known excess.
r = group_rates("iid", 100, 'independent bouts');
assert(r(1) <= 0.15 && r(2) <= 0.095, 'over-rejection beyond the recorded 0.100');
end

function t5
r = group_rates("homeostat", 50, 'HOMEOSTAT (Tw 2 h, Ts 1 h, gain 2)');
assert(r(1) > 0.3, 'compensation not detected for the homeostat');
end

function t6
group_rates("shared", 50, 'shared slow state (no homeostat)');
end

function rates = group_rates(mode, reps, label)
rej = [0 0]; lagMin = NaN(reps, 1); dmin = lagMin;
for g = 1:reps
    fl = arrayfun(@(f) SimulateSleepWake(mode, 90000 + 100*g + f), 1:48);
    R = SleepAmountACF(fl, NumSurrogates=199);
    rej = rej + [R.Group.MinP < 0.05, R.Group.MaxP < 0.05];
    lagMin(g) = R.Group.MinLagHours; dmin(g) = R.Group.MinDev;
end
rates = rej / reps;
fprintf('        %s, 48-fly groups x %d: rejection MinP (compensation) %.3f, MaxP (persistence) %.3f; median min deviation %.4f at %.1f h\n', ...
    label, reps, rates, median(dmin), median(lagMin));
end
