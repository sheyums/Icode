function test_HyperexponentialLRT()
%TEST_HYPEREXPONENTIALLRT Smoke tests for the mixture-order bootstrap LRT.
%
%   Run from the directory containing HyperexponentialLRT.m:
%
%       test_HyperexponentialLRT
%
% The test that matters here is TWO-DIRECTIONAL (tests 3 and 4). A
% bootstrap LRT that fails in one direction is worse than no test at all:
% one that never rejects looks reassuringly conservative while being blind
% to a real component, and one that always rejects confirms the larger
% model whatever the data. Both failure modes emit perfectly plausible
% p-values, so neither is visible without simulating data of known order
% and checking the verdict against the truth.
%
% B and n are kept small deliberately. Each replicate refits BOTH orders,
% so the cost is B x (two mixture fits) and a realistic B would make the
% suite unrunnable in Octave, whose gammainc is interpreted. The
% consequence is that p-values here are coarse -- with B=49 the finest
% resolvable value is 0.02 -- so the assertions are on the side of the
% decision boundary, never on a precise value.
%
% Tests covered
% -------------
%  1.  SamplingInterval is required   - and non-positive is rejected
%  2.  K0 < K1 is enforced            - nesting is what makes LR meaningful
%  3.  NULL TRUE: does not reject     - 2-component data, K=2 vs K=3
%  4.  NULL FALSE: does reject        - 3-component data, K=2 vs K=3
%  5.  Atom at zero is present        - the feature chi-square cannot have
%  6.  LR matches the two fits        - LR = 2*(logL1 - logL0) exactly
%  7.  p-value is a valid bootstrap p - in [1/(B+1), 1], never 0
%  8.  Null LRs are usable            - few negatives, finite, right count
%  9.  Continuous mode runs           - and its simulator is exact too
% 10.  Simulator reproduces the fit   - draws match the fitted survival
% 11.  RandomSeed reproducibility     - same seed, same p
% 12.  Degenerate replicates are KEPT - not silently dropped
% 13.  All UseParallel modes equal    - the only check on the seeding

if exist('OCTAVE_VERSION', 'builtin')
    lrt = @HyperexponentialLRT_oct;
    hyp = @FitHyperexponentialMLE_oct;
else
    lrt = @HyperexponentialLRT;
    hyp = @FitHyperexponentialMLE;
end

DT = 10; XMIN = 100;
FAST = {'SamplingInterval', DT, 'B', 49, 'Verbose', false, ...
        'nStartsBase', 3, 'nStartsPerComponent', 3, 'maxStarts', 12};

nPassed = 0; nFailed = 0;
fprintf('\nRunning tests for HyperexponentialLRT...\n\n');

%% Test 1: SamplingInterval is required
try
    rng(1);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 300);
    e1 = ''; e2 = '';
    try
        lrt(d, XMIN, 1, 2, 'B', 19, 'Verbose', false);
    catch err
        e1 = err.identifier;
    end
    try
        lrt(d, XMIN, 1, 2, 'SamplingInterval', -1, 'B', 19, 'Verbose', false);
    catch err
        e2 = err.identifier;
    end
    ok = ~isempty(strfind(e1, 'SamplingIntervalRequired')) && ~isempty(e2);
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 1, ...
        'SamplingInterval is required and must be positive', ...
        sprintf('got "%s" and "%s"', e1, e2));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 1, '', err.message);
end

%% Test 2: K0 must be strictly less than K1
%  The LR is only interpretable for NESTED models. K0 >= K1 would compute a
%  difference of log-likelihoods that is not a likelihood ratio at all.
try
    rng(2);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 300);
    bad = {};
    for pair = {[3 2], [2 2]}
        eid = '';
        try
            lrt(d, XMIN, pair{1}(1), pair{1}(2), FAST{:});
        catch err
            eid = err.identifier;
        end
        if isempty(strfind(eid, 'NotNested'))
            bad{end+1} = sprintf('K0=%d,K1=%d gave "%s"', ...
                pair{1}(1), pair{1}(2), eid); %#ok<AGROW>
        end
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 2, ...
        'K0 >= K1 is rejected as not nested', strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 2, '', err.message);
end

%% Test 3: NULL TRUE -- two-component data must NOT yield a third
try
    rng(3);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 600);
    L = lrt(d, XMIN, 2, 3, FAST{:}, 'RandomSeed', 31);
    ok = L.pValue > 0.05;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 3, ...
        sprintf('2-component truth not rejected: LR=%.2f, p=%.3g', L.LR, L.pValue), ...
        sprintf('REJECTED a true null: LR=%.3f p=%.4g (test is anti-conservative)', ...
            L.LR, L.pValue));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 3, '', err.message);
end

%% Test 4: NULL FALSE -- three-component data must yield a third
%  Well-separated time constants, so the third component is unambiguous. A
%  test that passes test 3 but fails this one is simply blind.
try
    rng(4);
    d = simMix([60 500 4000], [0.35 0.45 0.20], XMIN, DT, 900);
    L = lrt(d, XMIN, 2, 3, FAST{:}, 'RandomSeed', 41);
    ok = L.pValue <= 0.05;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 4, ...
        sprintf('3-component truth rejected K=2: LR=%.2f, p=%.3g', L.LR, L.pValue), ...
        sprintf('FAILED to reject a false null: LR=%.3f p=%.4g (test is blind)', ...
            L.LR, L.pValue));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 4, '', err.message);
end

%% Test 5: the atom at zero exists
%  Under the null, K1 SHOULD frequently collapse onto K0, putting mass
%  exactly at LR=0. That atom is the whole reason chi-square is wrong here,
%  and its absence would mean the degenerate replicates are being dropped.
try
    rng(5);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 600);
    L = lrt(d, XMIN, 2, 3, FAST{:}, 'RandomSeed', 51);
    ok = L.AtomAtZero > 0.02 && L.AtomAtZero < 0.95;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 5, ...
        sprintf('%.0f%% of null LRs sit at zero (K=3 collapsed onto K=2)', ...
            100*L.AtomAtZero), ...
        sprintf('AtomAtZero = %.3f, outside (0.02, 0.95)', L.AtomAtZero));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 5, '', err.message);
end

%% Test 6: LR is exactly 2*(logL1 - logL0) from the reported fits
try
    rng(6);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 500);
    L = lrt(d, XMIN, 1, 2, FAST{:}, 'RandomSeed', 61);
    lrCheck = 2*(L.LogLik1 - L.LogLik0);
    % and the fits it reports are the ones the fitter gives on this data
    H = hyp(d, XMIN, 'SamplingInterval', DT, 'MaxComponents', 2, ...
        'ErrorOnNoValidFit', false, 'Verbose', false);
    ok = abs(L.LR - lrCheck) < 1e-12 && ...
         abs(L.LogLik0 - H.AllFits(1).LogLik) < 1e-6 && ...
         abs(L.LogLik1 - H.AllFits(2).LogLik) < 1e-6;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 6, ...
        sprintf('LR = 2*(logL1-logL0) exactly, and both match the fitter (LR=%.4f)', L.LR), ...
        sprintf('LR=%.6f vs %.6f; logL0 %.6f vs %.6f', L.LR, lrCheck, ...
            L.LogLik0, H.AllFits(1).LogLik));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 6, '', err.message);
end

%% Test 7: the p-value is a valid bootstrap p-value
%  p = (1 + #{LR_b >= LR})/(BValid + 1) can never be 0 and never exceed 1.
try
    rng(7);
    d = simMix([60 500 4000], [0.35 0.45 0.20], XMIN, DT, 700);
    L = lrt(d, XMIN, 2, 3, FAST{:}, 'RandomSeed', 71);
    floorP = 1/(L.BValid + 1);
    ok = L.pValue >= floorP - 1e-12 && L.pValue <= 1 && L.pValue > 0;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 7, ...
        sprintf('p=%.4g lies in [%.4g, 1] and is never 0', L.pValue, floorP), ...
        sprintf('p=%g outside [%g, 1]', L.pValue, floorP));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 7, '', err.message);
end

%% Test 8: the null sample is usable
%  Negative LRs are impossible in exact arithmetic (K0 nests in K1), so a
%  large fraction of them means the K1 replicate fits are under-optimized
%  and the p-value is biased DOWN. This is the suite's honesty check.
try
    rng(8);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 600);
    L = lrt(d, XMIN, 2, 3, FAST{:}, 'RandomSeed', 81);
    ok = numel(L.LRNull) == L.BValid && all(isfinite(L.LRNull)) && ...
         all(L.LRNull >= 0) && L.NegativeLRFraction < 0.35 && ...
         L.BValid >= 20;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 8, ...
        sprintf('%d null LRs, all finite and clamped >= 0, %.0f%% were negative', ...
            L.BValid, 100*L.NegativeLRFraction), ...
        sprintf('BValid=%d numel=%d negFrac=%.3f', L.BValid, numel(L.LRNull), ...
            L.NegativeLRFraction));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 8, '', err.message);
end

%% Test 9: continuous mode runs and behaves
try
    rng(9);
    d = simMixCont([150 1600], [0.6 0.4], XMIN, 600);
    L = lrt(d, XMIN, 2, 3, 'SamplingInterval', DT, ...
        'DistributionType', 'continuous', 'B', 49, 'Verbose', false, ...
        'nStartsBase', 3, 'nStartsPerComponent', 3, 'maxStarts', 12, ...
        'RandomSeed', 91);
    ok = isfinite(L.LR) && L.pValue > 0 && L.pValue <= 1 && ...
         strcmp(char(L.DistributionType), 'continuous');
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 9, ...
        sprintf('continuous mode runs: LR=%.3f, p=%.3g', L.LR, L.pValue), ...
        sprintf('LR=%g p=%g mode=%s', L.LR, L.pValue, char(L.DistributionType)));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 9, '', err.message);
end

%% Test 10: the replicate simulator reproduces the fitted model
%  The null is only the null if replicates really come from the fitted K0.
%  Draw a large sample the same way the bootstrap does and compare its
%  empirical survival against the fit's own SurvivalHandle.
try
    rng(10);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 800);
    H = hyp(d, XMIN, 'SamplingInterval', DT, 'MaxComponents', 2, ...
        'ErrorOnNoValidFit', false, 'Verbose', false);
    f = H.AllFits(2);
    q = f.WeightsObserved(:).'; r = f.Rates(:).';
    N = 20000;
    nmin = max(1, round(XMIN/DT));
    cw = cumsum(q/sum(q)); u = rand(N,1); c = ones(N,1);
    for j = 1:numel(cw)-1, c = c + (u > cw(j)); end
    lam = r(c); lam = lam(:);
    v = rand(N,1); v(v<=0) = eps;
    g = floor(log(v) ./ (-lam*DT)); g(~isfinite(g)) = 0;
    sim = (nmin + g) * DT;
    probe = XMIN * [1 2 4 8 16];
    worst = 0;
    for t = probe
        worst = max(worst, abs(mean(sim > t) - f.SurvivalHandle(t)));
    end
    % 2.5 binomial SE at the worst case is about 0.009 at N=20000
    ok = worst < 0.015;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 10, ...
        sprintf('simulated draws match the fitted survival (max dev %.4f on %d draws)', ...
            worst, N), ...
        sprintf('max deviation %.4f exceeds Monte Carlo error', worst));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 10, '', err.message);
end

%% Test 11: RandomSeed makes the test reproducible
try
    rng(11);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 500);
    A = lrt(d, XMIN, 2, 3, FAST{:}, 'RandomSeed', 111);
    B = lrt(d, XMIN, 2, 3, FAST{:}, 'RandomSeed', 111);
    ok = A.pValue == B.pValue && abs(A.LR - B.LR) < 1e-12 && ...
         isequal(size(A.LRNull), size(B.LRNull)) && ...
         max(abs(A.LRNull - B.LRNull)) < 1e-9;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 11, ...
        'the same RandomSeed reproduces LR, the null sample and p', ...
        sprintf('p %.6g vs %.6g', A.pValue, B.pValue));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 11, '', err.message);
end

%% Test 12: replicates where K1 collapses are KEPT, not discarded
%  Under the null, K1 collapsing is expected, and those near-zero LRs ARE
%  the atom. If they were dropped as "degenerate" the null would be biased
%  upward and every p-value too small. Check the bookkeeping adds up:
%  BValid + FailedReplicates == B, with failures rare.
try
    rng(12);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 600);
    L = lrt(d, XMIN, 2, 3, FAST{:}, 'RandomSeed', 121);
    ok = (L.BValid + L.FailedReplicates) == L.B && ...
         L.FailedReplicates < 0.5*L.B && ...
         nnz(L.LRNull < 1e-6) > 0;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 12, ...
        sprintf('%d valid + %d failed = %d requested; %d collapsed replicates retained', ...
            L.BValid, L.FailedReplicates, L.B, nnz(L.LRNull < 1e-6)), ...
        sprintf('BValid=%d failed=%d B=%d zeros=%d', L.BValid, ...
            L.FailedReplicates, L.B, nnz(L.LRNull < 1e-6)));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 12, '', err.message);
end

%% Test 13: every UseParallel setting gives IDENTICAL results
%  This is the only test that catches a seeding mistake, and a seeding
%  mistake is the whole risk of parallelising a bootstrap. A bare parfor
%  leaves every worker on its own stream and lets iterations finish in any
%  order, so the null sample silently becomes irreproducible -- and it
%  would still look perfectly plausible, because a p-value from the wrong
%  null is still a number in [0,1]. Nothing else here would notice.
%
%  So the assertion is not "all three run", nor "the p-values are close",
%  but EQUAL TO THE BIT across "auto", true and false, on the LR, the whole
%  null sample, the p-value and the replicate bookkeeping. That holds only
%  if each replicate carries its own seed and is therefore independent of
%  the order it happens to be executed in.
%
%  LIMIT, worth knowing when this passes. Where no pool is open and no
%  Parallel Computing Toolbox exists, all three settings run serially, and
%  the test then proves only that the seeding path reproduces the ordinary
%  path -- real order-independence is demonstrated only on a machine with
%  workers. L.RanInParallel says which case you are in, so a passing run
%  tells you how much it proved.
try
    rng(13);
    d = simMix([150 1600], [0.6 0.4], XMIN, DT, 500);
    modes = {false, true, 'auto'};
    names = {'false', 'true', 'auto'};
    out = cell(1, 3);
    for m = 1:3
        out{m} = lrt(d, XMIN, 2, 3, FAST{:}, 'RandomSeed', 131, ...
            'UseParallel', modes{m});
    end
    A = out{1};
    bad = {};
    for m = 2:3
        Bm = out{m};
        if A.pValue ~= Bm.pValue
            bad{end+1} = sprintf('%s: p %.12g vs %.12g', names{m}, ...
                Bm.pValue, A.pValue); %#ok<AGROW>
        end
        if abs(A.LR - Bm.LR) > 0
            bad{end+1} = sprintf('%s: LR differs by %.3g', names{m}, ...
                abs(A.LR - Bm.LR)); %#ok<AGROW>
        end
        if ~isequal(size(A.LRNull), size(Bm.LRNull))
            bad{end+1} = sprintf('%s: null sizes %d vs %d', names{m}, ...
                numel(Bm.LRNull), numel(A.LRNull)); %#ok<AGROW>
        elseif max(abs(A.LRNull - Bm.LRNull)) > 0
            bad{end+1} = sprintf('%s: null differs by up to %.3g', names{m}, ...
                max(abs(A.LRNull - Bm.LRNull))); %#ok<AGROW>
        end
        if A.BValid ~= Bm.BValid || A.FailedReplicates ~= Bm.FailedReplicates
            bad{end+1} = sprintf('%s: bookkeeping differs', names{m}); %#ok<AGROW>
        end
    end
    % false must never run in parallel, whatever the environment offers
    if out{1}.RanInParallel
        bad{end+1} = 'UseParallel=false ran in parallel';
    end
    % and a bad value must be rejected rather than silently treated as one
    % of the three
    eid = '';
    try
        lrt(d, XMIN, 2, 3, FAST{:}, 'UseParallel', 'yes please');
    catch err
        eid = err.identifier;
    end
    if isempty(strfind(eid, 'BadUseParallel'))
        bad{end+1} = sprintf('unknown UseParallel gave "%s"', eid);
    end
    % no ternary in MATLAB, and merge() is Octave-only
    ranWhere = 'serially';
    if out{3}.RanInParallel, ranWhere = 'in parallel'; end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 13, ...
        sprintf(['false, true and "auto" agree exactly (p=%.4g, %d null ' ...
            'LRs); auto ran %s'], A.pValue, A.BValid, ranWhere), ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 13, '', err.message);
end

fprintf('\nSummary: %d passed, %d failed, %d total\n\n', ...
    nPassed, nFailed, nPassed + nFailed);
end

% ======================================================================
function [nP, nF] = rep(ok, nP, nF, idx, passMsg, failMsg)
if ok
    fprintf('[PASS] Test %d: %s\n', idx, passMsg); nP = nP + 1;
else
    fprintf('[FAIL] Test %d: %s\n', idx, failMsg); nF = nF + 1;
end
end

function d = simMix(tau, w, xmin, dt, n)
% Mixture of exponentials on the dt grid, conditioned on clearing xmin.
% ceil, not round: the fitters map with round, so this deliberately does
% not assume the two agree.
nmin = max(1, round(xmin/dt));
w = w(:).' / sum(w); cw = cumsum(w);
out = zeros(n,1); filled = 0;
while filled < n
    m = 6*(n - filled) + 500;
    u = rand(m,1); comp = ones(m,1);
    for j = 1:numel(cw)-1, comp = comp + (u > cw(j)); end
    t = -tau(comp)' .* log(rand(m,1));
    nv = ceil(t(:)/dt); nv = nv(nv >= nmin);
    take = min(numel(nv), n - filled);
    out(filled+1:filled+take) = nv(1:take)*dt; filled = filled + take;
end
d = out;
end

function d = simMixCont(tau, w, xmin, n)
w = w(:).' / sum(w); cw = cumsum(w);
out = zeros(n,1); filled = 0;
while filled < n
    m = 6*(n - filled) + 500;
    u = rand(m,1); comp = ones(m,1);
    for j = 1:numel(cw)-1, comp = comp + (u > cw(j)); end
    t = -tau(comp)' .* log(rand(m,1));
    t = t(t >= xmin);
    take = min(numel(t), n - filled);
    out(filled+1:filled+take) = t(1:take); filled = filled + take;
end
d = out;
end
