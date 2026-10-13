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
%   Groups 8-10 cover the later additions and are less about the auto-log
%   decision: the label-permutation null, the exact scale invariance of JS
%   divergence, and the per-replicate bandwidth estimation that makes the
%   bootstrap noise floor honest at small n.
%
%   Two of those groups guard against specific regressions rather than generic
%   sanity, and are worth stating plainly so they are not "simplified" away:
%
%     Group 9 pins the fact that rescaling BOTH samples by a common factor
%     leaves d exactly unchanged.  This matters because it is tempting to read
%     a larger d as evidence that two groups are "more separated" when one
%     condition simply has narrower distributions -- that reading is wrong, and
%     9A/9B make the invariance a test rather than a claim.  9C/9D pin what
%     does move d: the ratio of between-group separation to within-group
%     spread in the transformed coordinate.
%
%     Group 10 pins the small-n noise floor.  Bandwidth used to be cached from
%     the observed samples and reused by every bootstrap replicate, so no
%     replicate paid the cost of choosing its own bandwidth and the noise floor
%     came out far too low.  Two samples drawn from the SAME distribution at
%     nP/nQ = 400/10 then returned d_corr = 0.198 against a true value of 0.
%     10C holds the new behaviour below that number explicitly.
%
%   Groups:
%       1. Auto-log engagement and announcement
%       2. Sample-size dependence of the threshold
%       3. Explicit override in both directions
%       4. The transform does not corrupt the estimate
%       5. Core estimator sanity
%       6. Reproducibility
%       7. Validation and degenerate input
%       8. Permutation null
%       9. Scale invariance
%      10. Per-replicate bandwidth and the noise floor
%
%   Runtime note: groups 8 and 10 run bootstrap and permutation loops, which
%   re-estimate bandwidth per replicate, so the suite takes tens of seconds
%   rather than a couple.  Replicate counts here are deliberately small -- they
%   are sized to catch regressions, not to produce publication-grade intervals.
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

%% Group 8: Permutation null
fprintf('-- Group 8: Permutation null ------------------------\n');
rng(21);
Pa = skewed(400, 0.6) * 400;        % right-skewed, scaled off 1 so cutoff 0 bites
Qa = skewed(120, 0.6) * 550;

[~,~,~,~,~,~,~,~,~, pOff, dOff] = runQuiet(@() jsd_kde(Pa, Qa, 512, 0, true, 0));
check('8A: n_perm omitted leaves p_perm NaN', isnan(pOff), sprintf('%g', pOff));
check('8B: n_perm omitted leaves d_perm empty', isempty(dOff), '');

[~,~,~,~,~,~,~,~,~, pOn, dOn] = runQuiet(@() jsd_kde(Pa, Qa, 512, 0, true, 0, 0.05, 5, 300));
check('8C: d_perm has one entry per replicate', numel(dOn) == 300, sprintf('%d', numel(dOn)));
check('8D: p_perm respects the add-one floor', ...
    pOn >= 1/301 - eps && pOn <= 1, sprintf('%g', pOn));
check('8E: d_perm lies in [0,1]', all(dOn >= 0 & dOn <= 1), '');

% The nine original outputs must not shift when a permutation test is added:
% the permutation loop draws only AFTER the bootstrap loop has finished.
c1 = cell(1,9); c2 = cell(1,11);
[c1{:}] = runQuiet(@() jsd_kde(Pa, Qa, 512, 0, true, 40, 0.05, 42));
[c2{:}] = runQuiet(@() jsd_kde(Pa, Qa, 512, 0, true, 40, 0.05, 42, 200));
check('8F: n_perm does not perturb the bootstrap outputs', ...
    isequal(c1, c2(1:9)), 'bootstrap outputs moved');

% Calibration: under H0 the p-value should be roughly uniform, so a 5% nominal
% level must not fire far more often than 5% of the time.
rng(23); nrep = 12; pv = zeros(nrep,1);
for k = 1:nrep
    A8 = skewed(250, 0.6); B8 = skewed(250, 0.6);
    [~,~,~,~,~,~,~,~,~, pv(k)] = runQuiet(@() jsd_kde(A8, B8, 512, 0, true, 0, 0.05, [], 200));
end
fprintf('       H0 calibration: %d runs, median p = %.2f, frac p<=0.05 = %.2f\n', ...
    nrep, median(pv), mean(pv <= 0.05));
check('8G: H0 p-values are not systematically small', ...
    median(pv) > 0.2 && mean(pv <= 0.05) <= 0.25, sprintf('median %.3f', median(pv)));

rng(27);
[dSep,~,~,~,~,~,~,~,~, pSep] = runQuiet(@() ...
    jsd_kde(skewed(300,0.5)*400, skewed(300,0.5)*1100, 512, 0, true, 0, 0.05, [], 300));
fprintf('       power: well-separated samples d = %.3f, p = %.4f\n', dSep, pSep);
check('8H: well-separated samples reach the p-value floor', pSep <= 2/301, sprintf('%g', pSep));

[wP, mP] = captureWarn(@() jsd_kde(ones(50,1), skewed(50,0.6), 512, 0, true, 0, 0.05, [], 50));
check('8I: degenerate input announces the skipped test', ...
    strcmp(wP, 'jsd_kde:permSkippedDegenerate'), wP);
check('8J: skipped-test message says p_perm is NaN', contains(mP, 'NaN'), mP);
[~,~,~,~,~,~,~,~,~, pDeg] = runQuiet(@() jsd_kde(ones(50,1), skewed(50,0.6), 512, 0, true, 0, 0.05, [], 50));
check('8K: degenerate input returns p_perm = NaN', isnan(pDeg), sprintf('%g', pDeg));
fprintf('\n');

%% Group 9: Scale invariance
% JS divergence is exactly invariant under a common monotone reparametrisation
% (see SCALE INVARIANCE in jsd_kde). Rescaling BOTH samples by one factor must
% therefore leave d untouched -- so "the distributions are narrower" does not by
% itself move d. What moves it is the ratio of between-group separation to
% within-group spread in the transformed coordinate, which 9D pins.
fprintf('-- Group 9: Scale invariance ------------------------\n');
rng(31);
A9 = skewed(300, 0.6) * 400;  B9 = skewed(300, 0.6) * 660;
d9ref  = runQuiet(@() jsd_kde(A9,       B9,       512, 0, true, 0));
d9x10  = runQuiet(@() jsd_kde(10*A9,    10*B9,    512, 0, true, 0));
d9x100 = runQuiet(@() jsd_kde(100*A9,   100*B9,   512, 0, true, 0));
fprintf('       d = %.10f / %.10f / %.10f at x1 / x10 / x100\n', d9ref, d9x10, d9x100);
check('9A: a common 10x rescaling leaves d unchanged', ...
    abs(d9x10 - d9ref) < 1e-9, sprintf('delta %.2e', abs(d9x10 - d9ref)));
check('9B: a common 100x rescaling leaves d unchanged', ...
    abs(d9x100 - d9ref) < 1e-9, sprintf('delta %.2e', abs(d9x100 - d9ref)));

rng(33);
dSigma = zeros(1,4); sigmas = [0.8 0.5 0.3 0.15];
for k = 1:4
    % Separation fixed at 0.5 in log units; only the within-group spread moves.
    dSigma(k) = runQuiet(@() jsd_kde(lognrnd(6, sigmas(k), 400, 1), ...
                                     lognrnd(6.5, sigmas(k), 400, 1), 512, 0, true, 0));
end
fprintf('       sigma %.2f -> d %.3f;  sigma %.2f -> d %.3f\n', ...
    sigmas(1), dSigma(1), sigmas(end), dSigma(end));
check('9C: shrinking within-group spread at fixed separation raises d', ...
    all(diff(dSigma) > 0), sprintf('%s', mat2str(round(dSigma,3))));
check('9D: that rise is substantial, not marginal', ...
    dSigma(end) - dSigma(1) > 0.3, sprintf('rise %.3f', dSigma(end)-dSigma(1)));
fprintf('\n');

%% Group 10: Per-replicate bandwidth and the noise floor
% Bandwidth is re-estimated inside every bootstrap and permutation replicate.
% The point of that is an honest noise floor at small n: a 10-point KDE is far
% noisier than a 400-point one, and d_corr must know it. With bandwidth cached
% from the observed sample this failed badly -- two samples from the SAME
% distribution at 400/10 returned d_corr = 0.198 against a truth of 0.
fprintf('-- Group 10: Per-replicate bandwidth ----------------\n');
rng(37);
P10 = skewed(400, 0.6) * 400;  Q10 = skewed(10, 0.6) * 400;   % SAME distribution
[d10,~,~,~, nf_P, nf_Q, dc10] = runQuiet(@() jsd_kde(P10, Q10, 512, 0, true, 80, 0.05, 3));
fprintf('       400 vs 10 from one distribution: d = %.3f, d_corr = %.3f\n', d10, dc10);
fprintf('       noise floor: P (n=400) %.3f, Q (n=10) %.3f\n', nf_P, nf_Q);
check('10A: the smaller sample carries the larger noise floor', ...
    nf_Q > nf_P, sprintf('%.3f vs %.3f', nf_Q, nf_P));
check('10B: d_corr removes most of the small-n inflation', ...
    dc10 < 0.5 * d10, sprintf('d_corr %.3f vs d %.3f', dc10, d10));
check('10C: d_corr stays below the old cached-bandwidth value of 0.198', ...
    dc10 < 0.198, sprintf('%.3f', dc10));

% Verbose gates the log-transform announcement and nothing else.
[wV, ~] = captureWarn(@() jsd_kde(Pa, Qa, 512, 0, true, 0, 0.05, [], 0, false));
check('10D: Verbose false suppresses the log-transform notice', isempty(wV), wV);
[wT, ~] = captureWarn(@() jsd_kde(Pa, Qa, 512, 0, true, 0, 0.05, [], 0, true));
check('10E: Verbose true still issues it', strcmp(wT, 'jsd_kde:logTransform'), wT);
[wD, ~] = captureWarn(@() jsd_kde(ones(40,1), ones(40,1), 512, NaN, [], 0, 0.05, [], 0, false));
check('10F: Verbose false does NOT mask degenerate input', ...
    strcmp(wD, 'jsd_kde:degenerateInput'), wD);
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
