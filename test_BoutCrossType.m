function test_BoutCrossType(varargin)
%TEST_BOUTCROSSTYPE  Suite and calibration for BoutCrossType (PRESPEC entry A).
%
%   test_BoutCrossType           all tests
%   test_BoutCrossType('fast')   skip the calibrations (5-9)
%
% Ground truth: triples built by hand with corr(...,'Spearman') from the
% Statistics Toolbox (skipped if absent), the exact swap of r_pre and r_post
% under time reversal, and SimulateSleepWake flies whose structure is known.

fast = nargin > 0 && strcmpi(varargin{1}, 'fast');
assert(isequal(fileparts(which('BoutCrossType')), pwd), 'resolution check');
t0 = tic;
res = {};
res{end+1} = run_test('1: statistics equal hand-built Spearman correlations', @t1);
res{end+1} = run_test('2: non-adjacent bouts never form a triple', @t2);
res{end+1} = run_test('3: time reversal swaps r_pre and r_post, negates Delta', @t3);
res{end+1} = run_test('4: deterministic; global RNG untouched; strata need labels', @t4);
if ~fast
    res{end+1} = run_test('5: per-fly calibration, independent bouts (shift and strata)', @t5);
    res{end+1} = run_test('6: group calibration, 48-fly groups, independent bouts (shift)', @t6);
    res{end+1} = run_test('7: shared slow state, no homeostat', @t7);
    res{end+1} = run_test('8: LD-like profile, no coupling: strata null and count', @t8);
    res{end+1} = run_test('9: homeostat: power', @t9);
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

% ---------------------------------------------------------------------
function t1
if isempty(which('corr')), fprintf('        (corr absent: skipped)\n'); return, end
fl = SimulateSleepWake("homeostat", 1, Days=4);
R = BoutCrossType(fl, NumSurrogates=19);
d = fl.Dur; s = fl.IsSleep;
i = find(s); i = i(i > 1 & i < numel(d));
rp = corr(d(i), d(i - 1), 'Type', 'Spearman');
ro = corr(d(i), d(i + 1), 'Type', 'Spearman');
j = find(~s); j = j(j > 1 & j < numel(d));
qp = corr(d(j), d(j - 1), 'Type', 'Spearman');
i2 = i(i > 3);
r2 = corr(d(i2), d(i2 - 3), 'Type', 'Spearman');
T = R.PerFly;
assert(abs(T.r_pre - rp) < 1e-12 && abs(T.r_post - ro) < 1e-12, 'r_pre/r_post');
assert(abs(T.Delta - (rp - ro)) < 1e-12 && abs(T.q_pre - qp) < 1e-12, 'Delta/q_pre');
assert(abs(T.r_pre2 - r2) < 1e-12, 'r_pre2');
assert(T.NNonAdjacent == 0 && T.NTriples == numel(i), 'counts');
end

function t2
fl = SimulateSleepWake("iid", 2, Days=4);
k = 101;                                   % open a gap after bout k
fl.Onset(k+1:end) = fl.Onset(k+1:end) + 600;
R = BoutCrossType(fl, NumSurrogates=19);
full = BoutCrossType(SimulateSleepWake("iid", 2, Days=4), NumSurrogates=19);
assert(R.PerFly.NNonAdjacent == 1, 'gap not counted');
% the triples centred on k and k+1 need the k -> k+1 link: both are lost
assert(full.PerFly.NTriples + full.PerFly.NWakeTriples - ...
    (R.PerFly.NTriples + R.PerFly.NWakeTriples) == 2, 'a triple spanned the gap');
end

function t3
fl = SimulateSleepWake("homeostat", 3, Days=4);
rv = fl;
rv.Dur = flipud(fl.Dur); rv.IsSleep = flipud(fl.IsSleep);
rv.Onset = [0; cumsum(rv.Dur(1:end-1))];
A = BoutCrossType(fl, NumSurrogates=19); B = BoutCrossType(rv, NumSurrogates=19);
assert(abs(A.PerFly.r_pre - B.PerFly.r_post) < 1e-12, 'r_pre <-> r_post');
assert(abs(A.PerFly.Delta + B.PerFly.Delta) < 1e-12, 'Delta negates');
assert(abs(A.PerFly.DeltaW + B.PerFly.DeltaW) < 1e-12, 'DeltaW negates');
end

function t4
fl = SimulateSleepWake("iid", 4, Days=4);
rng(9); before = rand(1, 2); rng(9);
A = BoutCrossType(fl, NumSurrogates=49);
assert(isequal(before, rand(1, 2)), 'global stream used');
B = BoutCrossType(fl, NumSurrogates=49);
assert(isequaln(A.PerFly, B.PerFly) && isequaln(A.Group, B.Group), 'not deterministic');
fl2 = rmfield(fl, 'Strata');
try
    BoutCrossType(fl2, Null="strata", NumSurrogates=9);
    error('no error without Strata');
catch err
    assert(~contains(err.message, 'no error'), err.message);
end
end

% ---------------------------------------------------------------------
function t5
r1 = perfly(@(i) SimulateSleepWake("iid", 1000 + i), "shift", 200, 'independent, shift');
r2 = perfly(@(i) SimulateSleepWake("iid", 1000 + i), "strata", 200, 'independent, strata');
assert(all([r1 r2] <= 0.095), 'over-rejects on independent bouts');
end

function t6
rej = zeros(1, 7); reps = 100;
for g = 1:reps
    fl = arrayfun(@(f) SimulateSleepWake("iid", 20000 + 100*g + f), 1:48);
    R = BoutCrossType(fl, NumSurrogates=199);
    rej = rej + (R.Group.MeanP.' < 0.05);
end
show('group (48 flies, 100 groups), independent, shift: group-mean rejection', rej / reps);
assert(all(rej / reps <= 0.11), 'group mean over-rejects');
end

function t7
r = perfly(@(i) SimulateSleepWake("shared", 3000 + i), "shift", 200, 'shared slow state, shift');
assert(r(3) <= 0.095, 'Delta over-rejects under a shared slow state');
end

function t8
r = perfly(@(i) SimulateSleepWake("ld", 4000 + i), "strata", 200, 'LD-like profile, strata');
rIa = perfly(@(i) SimulateSleepWake("ld", 4000 + i), "shift", 200, 'LD-like profile, SHIFT (no profile in null)');
rej = zeros(1, 7); reps = 100;
for g = 1:reps
    fl = arrayfun(@(f) SimulateSleepWake("ld", 50000 + 100*g + f), 1:40);
    R = BoutCrossType(fl, Null="strata", NumSurrogates=199);
    rej = rej + (R.Group.CountP.' < 0.05);
end
show('count test (40 LD-like flies, 100 groups), strata', rej / reps);
fprintf('        (shift null on LD-like flies, r_pre per fly: %.3f -- the profile alone)\n', rIa(1));
assert(r(1) <= 0.095 && r(3) <= 0.095, 'strata null over-rejects r_pre or Delta per fly');
end

function t9
r = perfly(@(i) SimulateSleepWake("homeostat", 6000 + i), "shift", 200, 'HOMEOSTAT (Tw 2 h, Ts 1 h, gain 2), shift');
rej = zeros(1, 7); reps = 50;
for g = 1:reps
    fl = arrayfun(@(f) SimulateSleepWake("homeostat", 70000 + 100*g + f), 1:48);
    R = BoutCrossType(fl, NumSurrogates=199);
    rej = rej + (R.Group.MeanP.' < 0.05);
end
show('group (48 homeostat flies, 50 groups), shift: power', rej / reps);
assert(r(1) > 0.3, 'r_pre has little power against the homeostat');
end

% ---------------------------------------------------------------------
function rates = perfly(gen, nullName, reps, label)
rej = zeros(1, 7); vals = NaN(reps, 7);
for i = 1:reps
    R = BoutCrossType(gen(i), Null=nullName, NumSurrogates=199);
    T = R.PerFly;
    p = [T.r_pre_P, T.r_post_P, T.Delta_P, T.r_pre2_P, T.q_pre_P, T.q_post_P, T.DeltaW_P];
    rej = rej + (p < 0.05);
    vals(i, :) = [T.r_pre, T.r_post, T.Delta, T.r_pre2, T.q_pre, T.q_post, T.DeltaW];
end
rates = rej / reps;
show(sprintf('%s, %d flies: per-fly rejection', label, reps), rates);
fprintf('          medians: r_pre %.3f r_post %.3f Delta %.3f q_pre %.3f q_post %.3f DeltaW %.3f\n', ...
    median(vals(:, [1 2 3 5 6 7]), 'omitnan'));
end

function show(label, r)
fprintf('        %s\n          r_pre %.3f  r_post %.3f  Delta %.3f  r_pre2 %.3f  q_pre %.3f  q_post %.3f  DeltaW %.3f\n', ...
    label, r);
end
