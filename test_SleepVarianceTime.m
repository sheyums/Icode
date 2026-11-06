function test_SleepVarianceTime(varargin)
%TEST_SLEEPVARIANCETIME  Suite and calibration for SleepVarianceTime (PRESPEC entry E).
%
%   test_SleepVarianceTime           all tests
%   test_SleepVarianceTime('fast')   skip the calibrations (4-6)

fast = nargin > 0 && strcmpi(varargin{1}, 'fast');
assert(isequal(fileparts(which('SleepVarianceTime')), pwd), 'resolution check');
t0 = tic;
res = {};
res{end+1} = run_test('1: variance of the surplus equals a per-second brute force', @t1);
res{end+1} = run_test('2: runs counted by hand; gaps refused; deterministic', @t2);
res{end+1} = run_test('3: a fly identical to its surrogates gives logR = 0', @t3);
if ~fast
    res{end+1} = run_test('4: calibration, independent bouts (48-fly groups)', @t4);
    res{end+1} = run_test('5: homeostat: does the curve bend down?', @t5);
    res{end+1} = run_test('6: shared slow state: does it bend up (not down)?', @t6);
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
tot = fl.Onset(end) + fl.Dur(end);
s = false(tot, 1);
for i = find(fl.IsSleep).', s(fl.Onset(i) + 1 : fl.Onset(i) + fl.Dur(i)) = true; end
cs = [0; cumsum(s)]; pbar = cs(end) / tot;
src = fileread(which('SleepVarianceTime'));
d = fullfile(tempdir, 'svt_vt'); if ~isfolder(d), mkdir(d); end
k = strfind(src, 'function V = vt('); e = strfind(src, 'function z = runs(');
fid = fopen(fullfile(d, 'vt.m'), 'w'); fwrite(fid, src(k:e-1)); fclose(fid);
addpath(d); c = onCleanup(@() rmpath(d));
L = [3600 4*3600]; V = vt(fl, L, 900);
for j = 1:2
    ts = (0:900:tot - L(j)).';
    dd = cs(ts + L(j) + 1) - cs(ts + 1) - pbar * L(j);
    assert(abs(V(j) - var(dd)) < 1e-6 * var(dd), 'scale %d', j);
end
end

function t2
d = [300 10 400 20 500 30 350 40 600 50 320 60 700 70 310 80 800 90 330 95 900 99 305 15].';
n = numel(d); isS = mod((1:n).', 2) == 1; on = [0; cumsum(d(1:end-1))];
fl = struct('Dur', d, 'IsSleep', isS, 'Onset', on);
src = fileread(which('SleepVarianceTime'));
dd = fullfile(tempdir, 'svt_runs'); if ~isfolder(dd), mkdir(dd); end
k = strfind(src, 'function z = runs(');
fid = fopen(fullfile(dd, 'runs.m'), 'w'); fwrite(fid, src(k:end)); fclose(fid);
addpath(dd); c = onCleanup(@() rmpath(dd));
z = runs(fl);
v = d(isS); m = median(v); a = v(v ~= m) > m;
assert(z(1) == 1 + sum(a(2:end) ~= a(1:end-1)), 'sleep runs');
fl2 = fl; fl2.Onset(5:end) = fl2.Onset(5:end) + 10;
try
    SleepVarianceTime(fl2, NumSurrogates=9); error('no error');
catch err
    assert(strcmp(err.identifier, 'SleepVarianceTime:Gap'), err.message);
end
g = SimulateSleepWake("iid", 2, Days=3);
assert(isequaln(SleepVarianceTime(g, NumSurrogates=19).PerFly, SleepVarianceTime(g, NumSurrogates=19).PerFly));
end

function t3
% Constant bouts: every surrogate equals the data.
n = 400; d = repmat([600; 300], n/2, 1); isS = repmat([true; false], n/2, 1);
fl = struct('Dur', d, 'IsSleep', isS, 'Onset', [0; cumsum(d(1:end-1))]);
R = SleepVarianceTime(fl, NumSurrogates=9, ScalesHours=[1 2 4]);
assert(all(abs(R.LogR(isfinite(R.LogR))) < 1e-9), 'logR not 0');
end

function t4
r = rates("iid", 100, 'independent bouts');
assert(all(r <= 0.095), 'over-rejects on independent bouts');
end

function t5
rates("homeostat", 50, 'HOMEOSTAT (Tw 2 h, Ts 1 h, gain 2)');
end

function t6
rates("shared", 50, 'shared slow state (no homeostat)');
end

function r = rates(mode, reps, label)
rej = zeros(1, 4); ms = NaN(reps, 1); lr = NaN(reps, 6);
for g = 1:reps
    fl = arrayfun(@(f) SimulateSleepWake(mode, 150000 + 100*g + f), 1:48);
    R = SleepVarianceTime(fl, NumSurrogates=199);
    G = R.Group;
    rej = rej + [G.P_Reversion < 0.05, G.P_Persistence < 0.05, G.RunsP_Sleep < 0.05, G.RunsP_Wake < 0.05];
    ms(g) = G.MeanSlope; lr(g, :) = G.MeanLogR;
end
r = rej / reps;
fprintf(['        %s, 48-fly groups x %d: rejection reversion %.3f, persistence %.3f, runs sleep %.3f, wake %.3f\n' ...
    '          median group slope %.4f; median group logR at 1 2 4 8 16 24 h: %s\n'], ...
    label, reps, r, median(ms), mat2str(round(median(lr, 1), 3)));
end
