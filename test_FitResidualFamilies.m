function test_FitResidualFamilies(varargin)
%TEST_FITRESIDUALFAMILIES  Suite for FitResidualFamilies.
%
%   test_FitResidualFamilies           all tests
%   test_FitResidualFamilies('fast')   skip the null calibration (test 4)
%
% Ground truth: closed-form log-likelihoods summed directly from the
% densities; GH parameter counts from GeneralizedHyperbolic_MLE's own model
% names; and a calibration on residuals whose true family is known.

fast = nargin > 0 && strcmpi(varargin{1}, 'fast');
assert(isequal(fileparts(which('FitResidualFamilies')), pwd), 'resolution check');
t0 = tic;
res = {};
res{end+1} = run_test('1: closed-form LogLikes equal direct density sums', @t1);
res{end+1} = run_test('2: GH counts 5 parameters only for the full GH', @t2);
res{end+1} = run_test('3: each family recovers its own data at n = 400', @t3);
if ~fast
    res{end+1} = run_test('4: null calibration, ShiftedExp-true residuals (n = 180)', @t4);
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

function t1
rs = RandStream('mt19937ar', 'Seed', 1);
e = -500*log(rs.rand(200, 1)) - 300;
T = FitResidualFamilies(e);
loc = min(e); s = mean(e) - loc;
ll = sum(-log(s) - (e - loc)/s);
assert(abs(T.LogLike(1) - ll) < 1e-8, 'ShiftedExp');
mu = mean(e); v = mean((e - mu).^2);
ll = sum(-0.5*log(2*pi*v) - (e - mu).^2/(2*v));
assert(abs(T.LogLike(2) - ll) < 1e-8, 'Normal');
m = median(e); b = mean(abs(e - m));
ll = sum(-log(2*b) - abs(e - m)/b);
assert(abs(T.LogLike(3) - ll) < 1e-8, 'Laplace');
assert(all(abs(T.AIC - (2*T.K - 2*T.LogLike)) < 1e-9 | isnan(T.LogLike)), 'AIC');
assert(all(abs(T.BIC - (T.K*log(200) - 2*T.LogLike)) < 1e-9 | isnan(T.LogLike)), 'BIC');
end

function t2
rs = RandStream('mt19937ar', 'Seed', 2);
for trial = 1:3
    e = -500*log(rs.rand(150, 1)); e = e - mean(e);
    T = FitResidualFamilies(e);
    g = T(T.Family == "GH", :);
    assert(g.K == 4 + startsWith(g.Note, "Generalized"), 'K %d for "%s"', g.K, g.Note);
end
end

function t3
rs = RandStream('mt19937ar', 'Seed', 3);
n = 400;
cases = {'ShiftedExp', -500*log(rs.rand(n, 1)) - 200; ...
         'Normal', 100*rs.randn(n, 1); ...
         'Laplace', 100*sign(rs.rand(n, 1) - 0.5).*log(rs.rand(n, 1))};
for i = 1:size(cases, 1)
    T = FitResidualFamilies(cases{i, 2});
    rk = T.DeltaBIC(T.Family == cases{i, 1});
    fprintf('        %s data: own family DeltaBIC %.2f; best %s\n', cases{i, 1}, rk, ...
        T.Family(T.DeltaBIC == 0));
    assert(rk < 6, '%s data: own family %.2f BIC units from the best', cases{i, 1}, rk);
end
end

function t4
% Residuals of bouts with no dynamics: exponential minus its mean, integer
% seconds as in the data. How often does each family win by BIC, and how
% far is ShiftedExp from the best when it is the truth?
rs = RandStream('mt19937ar', 'Seed', 4);
reps = 100;
wins = strings(reps, 1);
dSE = NaN(reps, 1);
for i = 1:reps
    e = round(-600*log(rs.rand(180, 1)));
    e = e - mean(e);
    T = FitResidualFamilies(e);
    wins(i) = T.Family(find(T.DeltaBIC == 0, 1));
    dSE(i) = T.DeltaBIC(T.Family == "ShiftedExp");
end
[u, ~, j] = unique(wins);
c = accumarray(j, 1);
fprintf('        ShiftedExp-true, n = 180, %d sets; BIC winners: %s\n', reps, ...
    strjoin(u + " " + string(c), ', '));
fprintf('        ShiftedExp DeltaBIC: median %.2f, 90th pct %.2f, 95th %.2f, max %.2f\n', ...
    median(dSE), prctile(dSE, 90), prctile(dSE, 95), max(dSE));
assert(sum(wins == "ShiftedExp") / reps >= 0.5, 'ShiftedExp should usually win when true');
end
