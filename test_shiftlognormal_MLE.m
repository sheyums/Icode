%% test_shiftlognormal_MLE.m
%
% Smoke tests for shiftlognormal_MLE.
%
% Run this script from the directory containing shiftlognormal_MLE.m.
% Each test prints [PASS] or [FAIL] with a short description.
% A summary line is printed at the end.
%
% Tests covered
% -------------
%  1.  Parameter recovery          - estimates are close to known truth
%  2.  Required fields present     - all documented fields exist in output
%  3.  FitX / Fit consistency      - sorted, same length as N, non-negative
%  4.  Shift < min(x)              - hard constraint always holds
%  5.  AIC and BIC formulas        - exact arithmetic check
%  6.  NaN / Inf removal           - observations stripped, N reduced
%  7.  Minimum sample size (n=5)   - succeeds at the floor
%  8.  Insufficient data (n=4)     - correct error identifier
%  9.  Constant data               - correct error identifier
% 10.  Invalid mode argument       - correct error identifier
% 11.  Negative observations       - fits without error
% 12.  Near-lower-bound clipping   - NearLowerBound diagnostic is true
% 13.  VERBOSE mode                - does not error

nPassed = 0;
nFailed = 0;

% Suppress warnings expected by several tests so they do not clutter output.
% Individual tests restore the state after they run.
defaultWarnState = warning('query', 'all');

fprintf('\nRunning tests for shiftlognormal_MLE...\n\n');

%% Test 1: Parameter recovery

try
    rng(42);
    trueShift = 5;
    trueMu    = 1.5;
    trueSigma = 0.4;
    n         = 2000;

    x = trueShift + exp(trueMu + trueSigma * randn(n, 1));

    r = shiftlognormal_MLE(x);

    assert(r.Converged, ...
        'Optimizer did not converge.');
    assert(abs(r.Mu    - trueMu)    < 0.15, ...
        sprintf('mu error = %.4f', abs(r.Mu - trueMu)));
    assert(abs(r.Sigma - trueSigma) < 0.10, ...
        sprintf('sigma error = %.4f', abs(r.Sigma - trueSigma)));
    assert(abs(r.Shift - trueShift) < 0.50, ...
        sprintf('shift error = %.4f', abs(r.Shift - trueShift)));

    fprintf('[PASS] 1. Parameter recovery\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 1. Parameter recovery: %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Test 2: Required output fields present

requiredFields = { ...
    'Model', 'ParameterNames', 'Params', ...
    'Mu', 'Sigma', 'Shift', ...
    'LogLike', 'AIC', 'BIC', ...
    'FitX', 'Fit', ...
    'N', 'NParameters', ...
    'Converged', 'ExitFlag', 'OptimizerOutput', ...
    'Diagnostics'};

try
    rng(1);
    x = 3 + exp(1 + 0.5 * randn(200, 1));
    r = shiftlognormal_MLE(x);

    for k = 1:numel(requiredFields)
        assert(isfield(r, requiredFields{k}), ...
            sprintf('Missing field: %s', requiredFields{k}));
    end

    fprintf('[PASS] 2. Required output fields present\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 2. Required output fields present: %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Test 3: FitX / Fit consistency

try
    rng(2);
    x = 2 + exp(1 + 0.4 * randn(300, 1));
    r = shiftlognormal_MLE(x);

    assert(numel(r.FitX) == r.N, ...
        'FitX length does not match N.');
    assert(numel(r.Fit) == r.N, ...
        'Fit length does not match N.');
    assert(all(diff(r.FitX) >= 0), ...
        'FitX is not sorted.');
    assert(all(r.Fit >= 0), ...
        'Fit contains negative values.');
    assert(all(isfinite(r.Fit)), ...
        'Fit contains non-finite values.');

    fprintf('[PASS] 3. FitX / Fit consistency\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 3. FitX / Fit consistency: %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Test 4: Shift strictly less than min(x)

try
    rng(3);
    x = 7 + exp(2 + 0.6 * randn(500, 1));
    r = shiftlognormal_MLE(x);

    assert(r.Shift < min(x), ...
        sprintf('Shift (%.6g) is not less than min(x) (%.6g).', ...
        r.Shift, min(x)));

    fprintf('[PASS] 4. Shift < min(x)\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 4. Shift < min(x): %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Test 5: AIC and BIC formulas

try
    rng(4);
    x = 1 + exp(0.5 + 0.3 * randn(400, 1));
    r = shiftlognormal_MLE(x);

    k          = r.NParameters;   % 3
    expectedAIC = 2*k - 2*r.LogLike;
    expectedBIC = k*log(r.N) - 2*r.LogLike;

    assert(abs(r.AIC - expectedAIC) < 1e-8, ...
        sprintf('AIC mismatch: got %.10g, expected %.10g', r.AIC, expectedAIC));
    assert(abs(r.BIC - expectedBIC) < 1e-8, ...
        sprintf('BIC mismatch: got %.10g, expected %.10g', r.BIC, expectedBIC));

    fprintf('[PASS] 5. AIC and BIC formulas\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 5. AIC and BIC formulas: %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Test 6: NaN / Inf observations are removed

try
    rng(5);
    xClean  = 4 + exp(1.2 + 0.4 * randn(50, 1));
    xDirty  = [xClean; NaN; Inf; -Inf];

    ws = warning('off', 'all');
    r  = shiftlognormal_MLE(xDirty);
    warning(ws);

    assert(r.N == 50, ...
        sprintf('Expected N=50 after removal, got %d.', r.N));
    assert(r.Diagnostics.NonfiniteObservationsRemoved == 3, ...
        sprintf('Expected 3 removed, got %d.', ...
        r.Diagnostics.NonfiniteObservationsRemoved));
    assert(isstruct(r) && isfield(r, 'Shift'), ...
        'Fit did not complete after removal.');

    fprintf('[PASS] 6. NaN / Inf removal\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 6. NaN / Inf removal: %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Test 7: Minimum sample size n=5 succeeds

try
    rng(6);
    x = 2 + exp(1 + 0.5 * randn(5, 1));

    ws = warning('off', 'all');
    r  = shiftlognormal_MLE(x);
    warning(ws);

    assert(r.N == 5, 'Expected N=5.');
    assert(isstruct(r), 'Result is not a struct.');

    fprintf('[PASS] 7. Minimum sample size (n=5) succeeds\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 7. Minimum sample size (n=5) succeeds: %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Test 8: n=4 throws InsufficientData

try
    shiftlognormal_MLE([1; 2; 3; 4]);
    fprintf('[FAIL] 8. Insufficient data (n=4): no error was thrown\n');
    nFailed = nFailed + 1;
catch e
    if strcmp(e.identifier, 'shiftlognormal_MLE:InsufficientData')
        fprintf('[PASS] 8. Insufficient data (n=4) throws correct error\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] 8. Insufficient data (n=4): wrong error id: %s\n', ...
            e.identifier);
        nFailed = nFailed + 1;
    end
end

%% Test 9: Constant data throws ConstantData

try
    shiftlognormal_MLE([3; 3; 3; 3; 3]);
    fprintf('[FAIL] 9. Constant data: no error was thrown\n');
    nFailed = nFailed + 1;
catch e
    if strcmp(e.identifier, 'shiftlognormal_MLE:ConstantData')
        fprintf('[PASS] 9. Constant data throws correct error\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] 9. Constant data: wrong error id: %s\n', ...
            e.identifier);
        nFailed = nFailed + 1;
    end
end

%% Test 10: Invalid mode throws InvalidMode

try
    shiftlognormal_MLE([1; 2; 3; 4; 5], 99);
    fprintf('[FAIL] 10. Invalid mode: no error was thrown\n');
    nFailed = nFailed + 1;
catch e
    if strcmp(e.identifier, 'shiftlognormal_MLE:InvalidMode')
        fprintf('[PASS] 10. Invalid mode throws correct error\n');
        nPassed = nPassed + 1;
    else
        fprintf('[FAIL] 10. Invalid mode: wrong error id: %s\n', ...
            e.identifier);
        nFailed = nFailed + 1;
    end
end

%% Test 11: Negative observations

try
    rng(7);
    trueShift = -100;
    x         = trueShift + exp(1.5 + 0.5 * randn(300, 1));

    r = shiftlognormal_MLE(x);

    assert(isstruct(r),       'Result is not a struct.');
    assert(r.Shift < min(x),  'Shift is not less than min(x).');
    assert(isfinite(r.LogLike), 'Log-likelihood is not finite.');

    fprintf('[PASS] 11. Negative observations\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 11. Negative observations: %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Test 12: Near-lower-bound clipping
%
% A very tight distribution (small sigma) causes min(x) - range(x) to
% exceed the true shift, so the optimizer clips to the lower bound.
% Condition: exp(6*sigma) < 2, i.e., sigma < log(2)/6 ≈ 0.116.

try
    % sigma=0.02 guarantees exp(6*sigma)=1.13 << 2, so range << gap to
    % trueShift=0 regardless of seed.
    rng(8);
    trueShift = 0;
    x         = trueShift + exp(3 + 0.02 * randn(500, 1));

    ws = warning('off', 'all');
    r  = shiftlognormal_MLE(x);
    warning(ws);

    assert(r.Diagnostics.NearLowerBound, ...
        'Expected NearLowerBound to be true for a tight distribution.');
    assert(r.Shift < min(x), ...
        'Shift must still be less than min(x) even when clipped.');

    fprintf('[PASS] 12. Near-lower-bound clipping\n');
    nPassed = nPassed + 1;
catch e
    fprintf('[FAIL] 12. Near-lower-bound clipping: %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Test 13: VERBOSE mode does not error

try
    rng(9);
    x = 1 + exp(1 + 0.4 * randn(100, 1));

    % Redirect output to suppress the verbose lines during testing.
    diary_file = tempname;
    diary(diary_file);
    r = shiftlognormal_MLE(x, 'VERBOSE');
    diary off;
    if exist(diary_file, 'file'), delete(diary_file); end

    assert(isstruct(r), 'VERBOSE mode did not return a struct.');

    fprintf('[PASS] 13. VERBOSE mode does not error\n');
    nPassed = nPassed + 1;
catch e
    diary off;
    fprintf('[FAIL] 13. VERBOSE mode: %s\n', e.message);
    nFailed = nFailed + 1;
end

%% Summary

fprintf('\n----------------------------------------\n');
fprintf('Results: %d passed, %d failed (of %d)\n', ...
    nPassed, nFailed, nPassed + nFailed);

if nFailed == 0
    fprintf('All tests passed.\n');
else
    fprintf('%d test(s) FAILED.\n', nFailed);
end
fprintf('----------------------------------------\n\n');
