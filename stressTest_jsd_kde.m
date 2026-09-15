function stressTest_jsd_kde
%STRESSTEST_JSD_KDE  Stress tests for jsd_kde.
%
%   Emphasis on the AUTO-LOG DECISION: whether the log-transform engages when
%   the data are skewed enough to warrant it, stays out of the way when they
%   are not, and always announces itself so the choice is never silent.
%
%   The decision rule is
%       skew_thresh = 0.8 + 10/sqrt(nP + nQ)
%       apply log  <=>  skewness(P) > skew_thresh  OR  skewness(Q) > skew_thresh
%
%   so the threshold TIGHTENS as sample size grows: a moderately skewed sample
%   left alone at n = 200 can take the log path at n = 12000, because more data
%   makes a given skewness more credible.  Group 2 pins that behaviour.
%
%   The announcement is delivered through warning('jsd_kde:logTransform', ...),
%   so these tests inspect the warning identifier AND its text rather than
%   assuming the message says what it should.
%
%   Groups:
%       1. Auto-log engagement and announcement
%       2. Sample-size dependence of the threshold
%       3. Explicit override in both directions
%       4. The transform does not corrupt the estimate
%       5. Core estimator sanity
%       6. Reproducibility
%       7. Validation and degenerate input
%
%   See also JSD_KDE

clc;
if isempty(which('jsd_kde'))
    error('stressTest_jsd_kde:FunctionNotFound', 'jsd_kde.m must be on the path.');
end

passed = 0; failed = 0; total = 0; failures = strings(0,1);

fprintf('====================================================\n');
fprintf('jsd_kde - stress test\n');
fprintf('====================================================\n\n');

skewed = @(n, s) exp(s * randn(n,1));    % lognormal: right-skewed
symm   = @(n)    randn(n,1) + 10;        % normal: no skew

%% Group 1: auto-log engagement and announcement
fprintf('-- Group 1: Auto-log engagement and announcement ----\n');
rng(11);

A = skewed(400, 1.0);  B = skewed(400, 1.0);
[wid, msg] = captureWarn(@() jsd_kde(A, B));
sA = skewness(A); sB = skewness(B); thr = 0.8 + 10/sqrt(800);
fprintf('       strongly skewed: skewness %.2f / %.2f, threshold %.2f\n', sA, sB, thr);
check('1A: fires on strongly right-skewed data', ...
    strcmp(wid, 'jsd_kde:logTransform'), sprintf('id "%s"', wid));
check('1B: announcement says it was automatic', ...
    contains(lower(msg), 'automatic'), msg);
check('1C: announcement reports both skewness values', ...
    contains(msg, sprintf('%.2f', sA)) && contains(msg, sprintf('%.2f', sB)), msg);
check('1D: announcement reports the threshold', ...
    contains(msg, sprintf('%.2f', thr)), msg);

C = symm(400); D = symm(400);
[wid2, ~] = captureWarn(@() jsd_kde(C, D));
fprintf('       symmetric: skewness %.2f / %.2f, threshold %.2f\n', ...
    skewness(C), skewness(D), thr);
check('1E: silent on symmetric data', ...
    ~strcmp(wid2, 'jsd_kde:logTransform'), wid2);

% The rule is OR, so one skewed sample is enough on its own.
[wid3, ~] = captureWarn(@() jsd_kde(skewed(400,1.0), symm(400)));
check('1F: fires when only P is skewed', ...
    strcmp(wid3, 'jsd_kde:logTransform'), wid3);
[wid4, ~] = captureWarn(@() jsd_kde(symm(400), skewed(400,1.0)));
check('1G: fires when only Q is skewed', ...
    strcmp(wid4, 'jsd_kde:logTransform'), wid4);

% A log transform cannot help left skew, so it must stay out.
L1 = -exp(randn(400,1)); L2 = -exp(randn(400,1));
[wid5, ~] = captureWarn(@() jsd_kde(L1, L2));
fprintf('       left-skewed: skewness %.2f / %.2f\n', skewness(L1), skewness(L2));
check('1H: does NOT fire on left-skewed data', ...
    ~strcmp(wid5, 'jsd_kde:logTransform'), wid5);
fprintf('\n');

%% Group 2: threshold tightens with sample size
fprintf('-- Group 2: Threshold depends on sample size --------\n');
rng(5);
nSmall = 200; nBig = 12000;
s1 = skewed(nSmall, 0.30); s2 = skewed(nSmall, 0.30);
b1 = skewed(nBig,   0.30); b2 = skewed(nBig,   0.30);
thrS = 0.8 + 10/sqrt(2*nSmall);
thrB = 0.8 + 10/sqrt(2*nBig);
fprintf('       n=%5d: skew %.2f / %.2f  threshold %.2f\n', ...
    nSmall, skewness(s1), skewness(s2), thrS);
fprintf('       n=%5d: skew %.2f / %.2f  threshold %.2f\n', ...
    nBig, skewness(b1), skewness(b2), thrB);
check('2A: threshold is looser for small n', thrS > thrB, ...
    sprintf('%.3f vs %.3f', thrS, thrB));

[wS, ~] = captureWarn(@() jsd_kde(s1, s2));
[wB, ~] = captureWarn(@() jsd_kde(b1, b2));
firedS = strcmp(wS, 'jsd_kde:logTransform');
firedB = strcmp(wB, 'jsd_kde:logTransform');
fprintf('       fired at n=%d: %d   |   at n=%d: %d\n', nSmall, firedS, nBig, firedB);
check('2B: decision matches the stated rule at small n', ...
    firedS == (max(skewness(s1), skewness(s2)) > thrS), ...
    sprintf('fired=%d', firedS));
check('2C: decision matches the stated rule at large n', ...
    firedB == (max(skewness(b1), skewness(b2)) > thrB), ...
    sprintf('fired=%d', firedB));
fprintf('\n');

%% Group 3: explicit override
fprintf('-- Group 3: Explicit override -----------------------\n');
rng(7);
S1 = skewed(400, 1.0); S2 = skewed(400, 1.0);
[wOff, ~] = captureWarn(@() jsd_kde(S1, S2, 512, NaN, false));
check('3A: false suppresses it on skewed data', ...
    ~strcmp(wOff, 'jsd_kde:logTransform'), wOff);

N1 = symm(400); N2 = symm(400);
[wOn, mOn] = captureWarn(@() jsd_kde(N1, N2, 512, NaN, true));
check('3B: true forces it on symmetric data', ...
    strcmp(wOn, 'jsd_kde:logTransform'), wOn);
check('3C: forced announcement says user-specified', ...
    contains(lower(mOn), 'user-specified'), mOn);
fprintf('\n');

%% Group 4: the transform does not corrupt the estimate
fprintf('-- Group 4: Transform preserves the estimand --------\n');
rng(3);
P = exp(randn(800,1)); Q = exp(randn(800,1) + 0.8);
dOn  = runQuiet(@() jsd_kde(P, Q, 4096, 0, true));
dOff = runQuiet(@() jsd_kde(P, Q, 4096, 0, false));
relDiff = abs(dOn - dOff) / max(dOff, eps);
fprintf('       log %.5f vs no-log %.5f  (relative %.2f%%)\n', dOn, dOff, 100*relDiff);
check('4A: log and no-log agree within 2 percent', relDiff < 0.02, ...
    sprintf('%.5f vs %.5f', dOn, dOff));
check('4B: identical samples give 0 under the log path', ...
    runQuiet(@() jsd_kde(P, P, 512, 0, true)) < 1e-6, '');
check('4C: far-separated samples give ~1 under the log path', ...
    runQuiet(@() jsd_kde(P, P*1e6, 512, 0, true)) > 0.95, '');
fprintf('\n');

%% Group 5: core estimator sanity
fprintf('-- Group 5: Core estimator --------------------------\n');
rng(2);
X = randn(500,1); Y = randn(500,1); Z = randn(500,1) + 8;
dXY = runQuiet(@() jsd_kde(X, Y));
dXZ = runQuiet(@() jsd_kde(X, Z));
dZX = runQuiet(@() jsd_kde(Z, X));
fprintf('       d(X,Y)=%.4f  d(X,Z)=%.4f\n', dXY, dXZ);
check('5A: d lies in [0,1]', all([dXY dXZ] >= 0 & [dXY dXZ] <= 1), ...
    sprintf('%.4f %.4f', dXY, dXZ));
check('5B: symmetric in its arguments', dXZ == dZX, ...
    sprintf('%.10f vs %.10f', dXZ, dZX));
check('5C: self-distance is zero', runQuiet(@() jsd_kde(X, X)) < 1e-9, '');
check('5D: same distribution scores below different', dXY < dXZ, '');
check('5E: row and column inputs agree', ...
    abs(runQuiet(@() jsd_kde(X.', Z.')) - dXZ) < 1e-12, '');
fprintf('\n');

%% Group 6: reproducibility
fprintf('-- Group 6: Reproducibility -------------------------\n');
rng(4);
U = randn(400,1); V = randn(400,1) + 1;
[d1, ci1, se1] = runQuiet3(@() jsd_kde(U, V, 512, NaN, [], 200, 0.05, 42));
[d2, ci2, se2] = runQuiet3(@() jsd_kde(U, V, 512, NaN, [], 200, 0.05, 42));
check('6A: same seed reproduces d, CI and SE', ...
    d1 == d2 && isequal(ci1, ci2) && se1 == se2, '');
[~, ci3, ~] = runQuiet3(@() jsd_kde(U, V, 512, NaN, [], 200, 0.05, 43));
check('6B: a different seed changes the CI', ~isequal(ci1, ci3), '');
check('6C: d does not depend on the seed', ...
    d1 == runQuiet(@() jsd_kde(U, V, 512, NaN, [], 200, 0.05, 999)), '');

rng(123, 'twister');
before = rng;
runQuiet(@() jsd_kde(U, V, 512, NaN, [], 50, 0.05, 7));
afterState = rng;
check('6D: caller RNG state is restored', ...
    isequal(before.State, afterState.State), '');
fprintf('\n');

%% Group 7: validation and degenerate input
fprintf('-- Group 7: Validation and degenerate input ---------\n');
expectErr('7A: empty P',               @() jsd_kde([], randn(10,1)));
expectErr('7B: NaN in data',           @() jsd_kde([1;NaN;3], randn(10,1)));
expectErr('7C: Inf in data',           @() jsd_kde([1;Inf;3], randn(10,1)));
expectErr('7D: value below cutoff',    @() jsd_kde([1;2;3], [4;5;6], 512, 2));
expectErr('7E: negative rng_seed',     @() jsd_kde(randn(50,1), randn(50,1), 512, NaN, [], 10, 0.05, -3));
expectErr('7F: ci_alpha out of range', @() jsd_kde(randn(50,1), randn(50,1), 512, NaN, [], 10, 1.5));

[wc, ~] = captureWarn(@() jsd_kde(ones(50,1), ones(50,1)));
check('7G: constant vs same constant warns degenerate', ...
    strcmp(wc, 'jsd_kde:degenerateInput'), wc);
check('7H: constant vs same constant gives d = 0', ...
    runQuiet(@() jsd_kde(ones(50,1), ones(50,1))) == 0, '');
check('7I: distinct constants give d = 1', ...
    runQuiet(@() jsd_kde(ones(50,1), 2*ones(50,1))) == 1, '');
check('7J: [] and NaN min_cutoff agree', ...
    runQuiet(@() jsd_kde(X, Z, 512, [])) == runQuiet(@() jsd_kde(X, Z, 512, NaN)), '');
fprintf('\n');

%% Summary
fprintf('====================================================\n');
fprintf('SUMMARY: %d passed | %d failed | %d total\n', passed, failed, total);
fprintf('====================================================\n');
if failed > 0
    fprintf('\nFailed tests:\n');
    for k = 1:numel(failures)
        fprintf('  %s\n', failures(k));
    end
    error('stressTest_jsd_kde:TestsFailed', '%d stress test(s) failed.', failed);
else
    fprintf('\nAll stress tests passed.\n');
end

%% ---------------- nested helpers (share the counters) ----------------
    function check(name, cond, note)
        total = total + 1;
        if isscalar(cond) && islogical(cond) && cond
            passed = passed + 1;
            fprintf('  [PASS] %s\n', name);
        else
            failed = failed + 1;
            failures(end+1,1) = name + " -- " + string(note);
            fprintf('  [FAIL] %s -- %s\n', name, note);
        end
    end

    function expectErr(name, fh)
        total = total + 1;
        try
            fh();
            failed = failed + 1;
            failures(end+1,1) = name + " -- did not error";
            fprintf('  [FAIL] %s -- did not error\n', name);
        catch
            passed = passed + 1;
            fprintf('  [PASS] %s\n', name);
        end
    end
end

%% ---------------- local helpers ----------------
function [wid, msg] = captureWarn(fh)
%CAPTUREWARN Run fh silently and return the identifier and text of the last
%   warning it raised.  jsd_kde announces the log-transform decision through
%   warning(), so this is how the announcement is inspected rather than assumed.
st = warning('off', 'all');
restoreW = onCleanup(@() warning(st));
lastwarn('', '');
try
    fh();
catch
    % An error still leaves whatever warning preceded it in lastwarn.
end
[msg, wid] = lastwarn;
end

function varargout = runQuiet(fh)
st = warning('off', 'all');
restoreW = onCleanup(@() warning(st));
[varargout{1:max(nargout,1)}] = fh();
end

function [a, b, c] = runQuiet3(fh)
st = warning('off', 'all');
restoreW = onCleanup(@() warning(st));
[a, b, c] = fh();
end
