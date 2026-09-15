%% stress_test_GH_MLE.m
% Stress-test suite for GeneralizedHyperbolic_MLE.
%
% Run this script from the folder containing GeneralizedHyperbolic_MLE.m.
% Each test prints PASS or FAIL with a short explanation.
% A summary at the end counts failures.
%
% Random seed is fixed so results are reproducible.
rng(42);

nFail = 0;
nPass = 0;

% =========================================================================
%  HELPER: run a test case and catch errors / check results
% =========================================================================
function [passed, msg] = run_test(label, testFn)
    try
        [passed, msg] = testFn();
    catch ME
        passed = false;
        msg = ['Unexpected exception: ' ME.message];
    end
    if passed
        fprintf('  PASS  %s\n', label);
    else
        fprintf('  FAIL  %s  ->  %s\n', label, msg);
    end
end

% =========================================================================
%  SECTION 1: INPUT VALIDATION  (should error cleanly, not crash)
% =========================================================================
fprintf('\n=== SECTION 1: Input Validation ===\n');

tests1 = {
    'empty data',          @() expect_error(@() GeneralizedHyperbolic_MLE([]))
    'scalar data',         @() expect_error(@() GeneralizedHyperbolic_MLE(3.14))
    'single observation',  @() expect_error(@() GeneralizedHyperbolic_MLE([1]))
    'constant data',       @() expect_error(@() GeneralizedHyperbolic_MLE([5 5 5 5 5]))
    'NaN in data',         @() expect_error(@() GeneralizedHyperbolic_MLE([1 2 NaN 4]))
    'Inf in data',         @() expect_error(@() GeneralizedHyperbolic_MLE([1 2 Inf 4]))
    'bad modelFlag',       @() expect_error(@() GeneralizedHyperbolic_MLE(randn(100,1),'XYZ'))
    'numeric modelFlag',   @() expect_error(@() GeneralizedHyperbolic_MLE(randn(100,1), 99))
    'matrix data',         @() expect_error(@() GeneralizedHyperbolic_MLE(randn(10,10)))
    '2 observations ok',   @() expect_no_error(@() GeneralizedHyperbolic_MLE([0;1],'NIG'))
};

for i = 1:size(tests1,1)
    [p,m] = run_test(tests1{i,1}, tests1{i,2});
    nFail = nFail + ~p;  nPass = nPass + p;
end

% =========================================================================
%  SECTION 2: SYNTHETIC DATA — KNOWN NIG PARAMETERS
%  Generate NIG samples using the normal-variance-mean mixture.
%  X = mu + beta*V + sqrt(V)*Z,  V ~ InvGaussian(delta/gamma, delta^2)
% =========================================================================
fprintf('\n=== SECTION 2: NIG Recovery (known parameters) ===\n');

% True NIG params: mu=0, lambda=-0.5, alpha=3, beta=1, delta=1
% => gamma = sqrt(9-1)=sqrt(8)
trueParams_NIG = struct('mu',0,'lambda',-0.5,'alpha',3,'beta',1,'delta',1);
trueParams_NIG.gamma = sqrt(trueParams_NIG.alpha^2 - trueParams_NIG.beta^2);

for n = [50, 200, 1000, 5000, 50000]
    label = sprintf('NIG recovery n=%d', n);
    data  = sample_NIG(trueParams_NIG, n);
    testFn = @() check_fit(data, 'NIG', trueParams_NIG, n);
    [p,m] = run_test(label, testFn);
    nFail = nFail + ~p;  nPass = nPass + p;
end

% =========================================================================
%  SECTION 3: SYNTHETIC DATA — GH OPTIMIZER ROBUSTNESS
%  Note: exact GH sampling (lambda != -0.5) requires a GIG sampler not
%  included here.  sample_GH_via_NIG generates NIG-like data with matched
%  (alpha, beta, delta, mu) but lambda is not exactly reproduced.
%  These tests therefore check that the GH optimizer runs without error
%  and returns a finite, improved likelihood — not exact parameter recovery.
% =========================================================================
fprintf('\n=== SECTION 3: GH Optimizer Robustness (approximate GH data) ===\n');

% NIG data, then fit as GH: GH LogLike should be >= NIG LogLike
for n = [200, 1000, 5000]
    label = sprintf('GH fit of NIG data n=%d (LogLike >= NIG)', n);
    data_tmp = sample_NIG(trueParams_NIG, n);
    testFn = @() check_gh_improves_nig(data_tmp);
    [p,m] = run_test(label, testFn);
    nFail = nFail + ~p;  nPass = nPass + p;
end

% =========================================================================
%  SECTION 4: NEAR-BOUNDARY PARAMETER CONFIGURATIONS
% =========================================================================
fprintf('\n=== SECTION 4: Near-Boundary Configurations ===\n');

n = 500;

% Near-symmetric: beta very close to zero
p_sym = struct('mu',0,'lambda',-0.5,'alpha',2,'beta',1e-4,'delta',1);
p_sym.gamma = sqrt(p_sym.alpha^2 - p_sym.beta^2);
[pf,m] = run_test('NIG near-symmetric (beta~0)', ...
    @() check_finite_and_ordered(sample_NIG(p_sym,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Near-boundary skewness: |beta|/alpha = 0.95
p_skew = struct('mu',0,'lambda',-0.5,'alpha',4,'beta',3.8,'delta',1);
p_skew.gamma = sqrt(p_skew.alpha^2 - p_skew.beta^2);
[pf,m] = run_test('NIG high skew (|beta|/alpha=0.95)', ...
    @() check_finite_and_ordered(sample_NIG(p_skew,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Negative skew
p_negskew = struct('mu',0,'lambda',-0.5,'alpha',4,'beta',-3.8,'delta',1);
p_negskew.gamma = sqrt(p_negskew.alpha^2 - p_negskew.beta^2);
[pf,m] = run_test('NIG negative high skew', ...
    @() check_finite_and_ordered(sample_NIG(p_negskew,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Very small delta
p_sdelta = struct('mu',0,'lambda',-0.5,'alpha',5,'beta',1,'delta',0.01);
p_sdelta.gamma = sqrt(p_sdelta.alpha^2 - p_sdelta.beta^2);
[pf,m] = run_test('NIG tiny delta=0.01', ...
    @() check_finite_and_ordered(sample_NIG(p_sdelta,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Very large delta
p_ldelta = struct('mu',0,'lambda',-0.5,'alpha',2,'beta',0.5,'delta',50);
p_ldelta.gamma = sqrt(p_ldelta.alpha^2 - p_ldelta.beta^2);
[pf,m] = run_test('NIG large delta=50', ...
    @() check_finite_and_ordered(sample_NIG(p_ldelta,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Large location offset (mu=1e6)
p_largemu = struct('mu',1e6,'lambda',-0.5,'alpha',3,'beta',0.5,'delta',1);
p_largemu.gamma = sqrt(p_largemu.alpha^2 - p_largemu.beta^2);
[pf,m] = run_test('NIG large mu=1e6', ...
    @() check_finite_and_ordered(sample_NIG(p_largemu,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% =========================================================================
%  SECTION 5: GH LAMBDA EXTREMES
% =========================================================================
fprintf('\n=== SECTION 5: Extreme Lambda Values ===\n');

% Section 5 uses exact NIG data (lambda=-0.5) fitted as GH.
% We test that the optimizer runs and LogLike is finite.
% We do NOT require all CI bounds to be finite: when the GH lambda profile
% is flat (e.g., the data don't distinguish lambda values), an unbounded
% CI in one direction is the statistically correct result, not a bug.
% The check here is: no crash, finite LogLike, |beta_hat| < alpha_hat.
n = 500;
for lam = [-3, -1, -0.5, 0, 0.5, 1, 3, 5]
    label  = sprintf('GH fit, NIG data, sweep lambda grid lam=%.1f (no crash)', lam);
    data   = sample_NIG(trueParams_NIG, n);   % exact NIG data each time
    [pf,m] = run_test(label, @() check_runs_ok(data,'GH'));
    nFail  = nFail+~pf; nPass = nPass+pf;
end

% =========================================================================
%  SECTION 6: DATA PATHOLOGIES
% =========================================================================
fprintf('\n=== SECTION 6: Data Pathologies ===\n');

base_data = sample_NIG(trueParams_NIG, 500);

% A few extreme outliers (3 observations at 20 sigma)
data_out = [base_data; mean(base_data) + 20*std(base_data)*[1;-1;1]];
[pf,m] = run_test('NIG with 3 extreme outliers', ...
    @() check_finite_and_ordered(data_out,'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Near-normal data (very large alpha → GH approaches Normal)
data_norm = randn(500,1);  % pure Normal
[pf,m] = run_test('Near-normal data (NIG fit)', ...
    @() check_finite_and_ordered(data_norm,'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% GH fit of Normal data: lambda profile may be flat/unbounded — that is
% correct behavior, not a bug.  Just check it runs and LogLike is finite.
[pf,m] = run_test('Near-normal data (GH fit, no crash)', ...
    @() check_runs_ok(data_norm,'GH'));
nFail = nFail+~pf; nPass = nPass+pf;

% Heavily skewed data (log-normal, which GH cannot perfectly fit)
data_lnorm = exp(0.5*randn(500,1));
[pf,m] = run_test('Log-normal data (NIG fit, misspecified)', ...
    @() check_finite_and_ordered(data_lnorm,'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% =========================================================================
%  SECTION 7: AUTO MODEL SELECTION
% =========================================================================
fprintf('\n=== SECTION 7: AUTO Model Selection ===\n');

% True NIG data — AUTO should select NIG (or at least not crash)
data_NIG = sample_NIG(trueParams_NIG, 1000);
[pf,m] = run_test('AUTO on NIG data — expect NIG or GH, no crash', ...
    @() check_auto(data_NIG));
nFail = nFail+~pf; nPass = nPass+pf;

% True GH data with lambda far from -0.5
% NOTE: sample_GH_via_NIG ignores lambda (see its help), so this exercises the
% optimizer on GH-shaped input rather than testing lambda recovery.
trueParams_GH = struct('mu',0,'lambda',1.5,'alpha',3,'beta',1,'delta',1);
trueParams_GH.gamma = sqrt(trueParams_GH.alpha^2 - trueParams_GH.beta^2);
data_GH = sample_GH_via_NIG(trueParams_GH, 2000);
[pf,m] = run_test('AUTO on GH data (lambda=1.5, n=2000)', ...
    @() check_auto(data_GH));
nFail = nFail+~pf; nPass = nPass+pf;

% =========================================================================
%  SECTION 8: PROFILE CI SANITY CHECKS
% =========================================================================
fprintf('\n=== SECTION 8: Profile CI Sanity ===\n');

data_ci = sample_NIG(trueParams_NIG, 500);
res = GeneralizedHyperbolic_MLE(data_ci, 'NIG');
ci  = res.profile_ci;
params = res.Params;  % [mu; lambda; alpha; beta; delta]

[pf,m] = run_test('CI: mu bounds finite and ordered', ...
    @() check_ci_bounds(ci.mu, params(1)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('CI: alpha bounds finite, ordered, and positive', ...
    @() check_ci_bounds_positive(ci.alpha, params(3)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('CI: beta bounds finite and ordered', ...
    @() check_ci_bounds(ci.beta, params(4)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('CI: delta bounds finite, ordered, and positive', ...
    @() check_ci_bounds_positive(ci.delta, params(5)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('CI: alpha lower > 0 strictly', ...
    @() check_condition(ci.alpha(1) > 0, 'alpha lower CI not positive'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('CI: alpha CI contains true alpha=3', ...
    @() check_condition(ci.alpha(1)<=3 && ci.alpha(2)>=3, ...
    sprintf('true alpha=3 not in [%.3f, %.3f]',ci.alpha(1),ci.alpha(2))));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('CI: |beta| < alpha at MLE (constraint)', ...
    @() check_condition(abs(params(4)) < params(3), ...
    sprintf('|beta|=%.4f >= alpha=%.4f',abs(params(4)),params(3))));
nFail = nFail+~pf; nPass = nPass+pf;

% =========================================================================
%  SECTION 9: OUTPUT STRUCTURE COMPLETENESS
% =========================================================================
fprintf('\n=== SECTION 9: Output Structure Fields ===\n');

res_gh  = GeneralizedHyperbolic_MLE(sample_NIG(trueParams_NIG,200),'GH');
res_nig = GeneralizedHyperbolic_MLE(sample_NIG(trueParams_NIG,200),'NIG');
res_auto= GeneralizedHyperbolic_MLE(sample_NIG(trueParams_NIG,200),'AUTO');

required = {'Model','Params','LogLike','Fit','FitX','FitPDF','AIC','BIC', ...
    'theta_hat','theta_names','nll','hessian','hessian_coordinates', ...
    'profile_ci','exitflag','output','converged'};

for i = 1:numel(required)
    f = required{i};
    [pf,m] = run_test(['GH has field: ' f], ...
        @() check_condition(isfield(res_gh,f), ['missing field: ' f]));
    nFail = nFail+~pf; nPass = nPass+pf;
end

[pf,m] = run_test('AUTO has ModelSelection field', ...
    @() check_condition(isfield(res_auto,'ModelSelection'),'missing ModelSelection'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('AUTO has LRT field', ...
    @() check_condition(isfield(res_auto,'LRT'),'missing LRT'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('Params has 5 elements', ...
    @() check_condition(numel(res_gh.Params)==5,'Params not length 5'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('Fit length matches data length', ...
    @() check_condition(numel(res_gh.Fit)==200,'Fit length mismatch'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('FitX is sorted ascending', ...
    @() check_condition(all(diff(res_gh.FitX)>=0),'FitX not sorted'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('FitPDF all non-negative', ...
    @() check_condition(all(res_gh.FitPDF(isfinite(res_gh.FitPDF))>=0), ...
    'negative density values'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('LogLike is finite scalar', ...
    @() check_condition(isscalar(res_gh.LogLike) && isfinite(res_gh.LogLike), ...
    sprintf('LogLike = %g', res_gh.LogLike)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('NIG lambda fixed at -0.5', ...
    @() check_condition(res_nig.Params(2)==-0.5, ...
    sprintf('NIG lambda=%.4f',res_nig.Params(2))));
nFail = nFail+~pf; nPass = nPass+pf;

% =========================================================================
%  SECTION 10: DEGENERACY GUARDS
% =========================================================================
% This section exists because a real dataset produced LogLike > 0 with a large
% negative AIC and an obviously wrong fit, while every test above reported
% PASS.  The cause was the optimizer escaping to parameters where the
% log-density is dominated by floating-point cancellation, so the likelihood
% it was maximizing was noise rather than a likelihood.
%
% Each case below is a direction in which the GH likelihood is known to be
% badly behaved.  check_model_validity enforces that whatever comes back is
% still a probability model: density evaluable, likelihood self-consistent,
% and the fitted density integrating to one.
fprintf('\n=== SECTION 10: Degeneracy Guards ===\n');

degenCases = {
    'large n (cancellation regime)',   sample_NIG(trueParams_NIG, 5000)
    'very large n',                    sample_NIG(trueParams_NIG, 20000)
    'near-normal (alpha -> infinity)', randn(3000,1)
    'heavily tied integers',           round(sample_NIG(trueParams_NIG, 2000))
    'coarse discretization',           round(sample_NIG(trueParams_NIG, 2000)*2)/2
    'huge scale (x 1e6)',              sample_NIG(trueParams_NIG, 1000)*1e6
    'tiny scale (x 1e-6)',             sample_NIG(trueParams_NIG, 1000)*1e-6
    'shifted far from zero',           sample_NIG(trueParams_NIG, 1000)+1e5
};
for i = 1:size(degenCases,1)
    for flag = {'NIG','GH'}
        [pf,m] = run_test(sprintf('%s [%s]', degenCases{i,1}, flag{1}), ...
            @() check_no_degeneracy(degenCases{i,2}, flag{1}));
        nFail = nFail+~pf; nPass = nPass+pf;
    end
end

% =========================================================================
%  SECTION 11: REGRESSION GUARDS FROM THE SEP 2026 INDEPENDENT AUDIT
% =========================================================================
% Every test above passed while the estimator had the defects below. They
% were found by checking it against ground truth that shares no code with it
% (a Bessel-free mixture density, an independent Nelder-Mead profile, exact
% location-scale equivariance), and each test asserts a RULE rather than one
% realization, so a correct estimator cannot fail it by chance.
%
%   11A-B  An OPEN bound was reported as a tiny finite number. A non-finite
%          profile value (besselk overflow) was read as the profile crossing
%          the threshold, and bisection homed in on the overflow edge:
%          delta intervals like [2.4e-305, 0.91] on variance-gamma data.
%   11C    delta_at_bound stayed false for fits sitting on the delta floor:
%          interior-point stops strictly inside a bound, 1e-3 to 2e-3 above
%          it in log units, and the flag demanded 1e-8.
%   11D    Results depended on the units of the data: the mu interval moved
%          17% when one sample was multiplied by 1e4. Absolute tolerances,
%          and fzero's TolX, which it multiplies by |x| internally.
%   11E    single input ran the likelihood in single precision.
fprintf('\n=== SECTION 11: Audit Regression Guards ===\n');

% Variance-gamma-limit data: X = mu + beta*V + sqrt(V)*Z, V ~ Exp(1), i.e.
% GH with lambda = 1 and delta -> 0, where delta's lower bound is genuinely open.
rs11 = RandStream('twister', 'Seed', 1);
V11  = -log(rand(rs11, 600, 1));
xVG  = 2 + 0.5*V11 + sqrt(V11).*randn(rs11, 600, 1);
resVG = GeneralizedHyperbolic_MLE(xVG, 'GH');

[pf,m] = run_test('11A: open delta bound on VG-limit data is exactly 0', ...
    @() check_condition(resVG.profile_ci.delta(1) == 0, ...
    sprintf('delta lower bound = %.4g (an open bound must be 0, not a tiny number)', ...
    resVG.profile_ci.delta(1))));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,m] = run_test('11B: no CI bound is a spurious tiny or huge finite number', ...
    @() check_no_spurious_bounds(resVG));
nFail = nFail+~pf; nPass = nPass+pf;

% 11C: VG data with delta exactly 0 (gamma mixing), large n: the fit goes to
% the floor. The assertion is conditional, so it states the rule itself.
V11c = -log(rand(rs11, 1000, 1));
xVGc = 0.5*V11c + sqrt(V11c).*randn(rs11, 1000, 1);
resVGc = GeneralizedHyperbolic_MLE(xVGc, 'GH');
floor11 = 1e-3*std(xVGc);
[pf,m] = run_test('11C: delta within 1% of its floor => delta_at_bound', ...
    @() check_condition(resVGc.Params(5) > 1.01*floor11 || resVGc.delta_at_bound, ...
    sprintf('delta = %.6g, floor = %.6g, delta_at_bound = %d', ...
    resVGc.Params(5), floor11, resVGc.delta_at_bound)));
nFail = nFail+~pf; nPass = nPass+pf;

% 11D: location-scale equivariance, y = 1e4*x + 1e4.
x11  = sample_NIG(trueParams_NIG, 800);
r11a = GeneralizedHyperbolic_MLE(x11, 'NIG');
r11b = GeneralizedHyperbolic_MLE(1e4*x11 + 1e4, 'NIG');
Pm   = [1e4*r11a.Params(1) + 1e4; r11a.Params(2); r11a.Params(3)/1e4; ...
        r11a.Params(4)/1e4; r11a.Params(5)*1e4];
[pf,m] = run_test('11D: parameters are equivariant under 1e4*x + 1e4', ...
    @() check_condition(max(abs(r11b.Params - Pm) ./ abs(Pm)) < 1e-5, ...
    sprintf('max relative error %.3g', max(abs(r11b.Params - Pm) ./ abs(Pm)))));
nFail = nFail+~pf; nPass = nPass+pf;
dLL = r11b.LogLike - (r11a.LogLike - 800*log(1e4));
[pf,m] = run_test('11D: log-likelihood shifts by exactly -n*log(scale)', ...
    @() check_condition(abs(dLL) < 1e-6, sprintf('discrepancy %.3g', dLL)));
nFail = nFail+~pf; nPass = nPass+pf;
ciA = [1e4*r11a.profile_ci.mu + 1e4, r11a.profile_ci.alpha/1e4, ...
       r11a.profile_ci.beta/1e4, r11a.profile_ci.delta*1e4];
ciB = [r11b.profile_ci.mu, r11b.profile_ci.alpha, r11b.profile_ci.beta, r11b.profile_ci.delta];
ciErr = max(abs(ciB - ciA) ./ abs(ciA));
[pf,m] = run_test('11D: CI bounds are equivariant (to root-finder tolerance)', ...
    @() check_condition(ciErr < 5e-3, sprintf('max relative error %.3g (was 0.17 before the fix)', ciErr)));
nFail = nFail+~pf; nPass = nPass+pf;

% 11E: single precision input
r11s = GeneralizedHyperbolic_MLE(single(x11), 'NIG');
[pf,m] = run_test('11E: single input is fitted in double precision', ...
    @() check_condition(isa(r11s.LogLike, 'double') && ...
    abs(r11s.LogLike - GeneralizedHyperbolic_MLE(double(single(x11)), 'NIG').LogLike) < 1e-9, ...
    sprintf('LogLike class %s', class(r11s.LogLike))));
nFail = nFail+~pf; nPass = nPass+pf;

% 11F: ten points near the Normal limit. A sound fit must not fail its own
% validity check because of the log-density's floating-point noise: the fixed
% 1e-6 tolerance reported valid = false here (Bessel argument ~6e9).
x11f = [-0.92891230072865449 0.29770359610432284 0.19808604340722599 0.80962243337342477 0.37782371522682678 -0.44532913476610414 -0.19631927565757104 -0.057295620412519732 -0.92631446155746755 -0.29055067417114161]';
r11f = GeneralizedHyperbolic_MLE(x11f, 'AUTO');
[pf,m] = run_test('11F: tiny sample near the Normal limit stays valid', ...
    @() check_condition(r11f.valid, sprintf('valid = %d, LogLike = %.8g', r11f.valid, r11f.LogLike)));
nFail = nFail+~pf; nPass = nPass+pf;

% =========================================================================
%  SUMMARY
% =========================================================================
fprintf('\n=== SUMMARY ===\n');
fprintf('  Passed: %d\n', nPass);
fprintf('  Failed: %d\n', nFail);
fprintf('  Total : %d\n', nPass+nFail);
if nFail == 0
    fprintf('  All tests passed.\n');
else
    fprintf('  *** %d test(s) failed — review output above. ***\n', nFail);
end


% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function [passed, msg] = expect_error(fn)
% Pass if fn throws any error.
    try
        fn();
        passed = false;
        msg = 'Expected an error but none was thrown.';
    catch
        passed = true;
        msg = '';
    end
end

function [passed, msg] = expect_no_error(fn)
% Pass if fn runs without error.
    try
        fn();
        passed = true;
        msg = '';
    catch ME
        passed = false;
        msg = ME.message;
    end
end

function [passed, msg] = check_condition(cond, failMsg)
    if cond
        passed = true;  msg = '';
    else
        passed = false;  msg = failMsg;
    end
end

function [passed, msg] = check_finite_and_ordered(data, flag)
% Fit model, check that it ran without error and CIs are ordered.
%
% CI bounds are checked by MEANING, not merely by finiteness.
%
% Section 5 already states the principle: when a profile is flat, an unbounded
% CI is the statistically correct answer rather than a bug.  This helper used
% to demand finite bounds unconditionally, which contradicted that and produced
% a spurious failure on 'NIG large delta=50'.  On that data the beta profile is
% flat to 1e-4 across the whole plausible range: beta genuinely is not
% identified, however well alpha and delta are pinned down.
%
% The estimator now distinguishes the two cases that used to look identical:
%   +/-Inf  the profile was evaluated throughout and never crossed - OPEN,
%           a legitimate result.
%   NaN     the profile could not be evaluated - UNKNOWN, a real failure.
% So Inf is accepted and NaN is not.  Note this is a per-parameter question:
% whole-model identifiability (beating a Normal on AIC) does not imply every
% individual parameter is identified, which is why that is not the criterion.
    res = GeneralizedHyperbolic_MLE(data, flag);
    [passed, msg] = check_model_validity(res, data);
    if ~passed, return; end

    ci  = res.profile_ci;
    fns = fieldnames(ci);
    for i = 1:numel(fns)
        b = ci.(fns{i});
        if any(isnan(b))
            passed = false;
            msg = sprintf(['CI for %s is NaN: the profile could not be ' ...
                'evaluated, so the bound is unknown (an open interval would ' ...
                'have been reported as +/-Inf).'], fns{i});
            return
        end
        if b(1) > b(2)
            passed = false;
            msg = sprintf('CI for %s is reversed: [%.4g, %.4g].', fns{i}, b(1), b(2));
            return
        end
    end
    if ~isfinite(res.LogLike)
        passed = false;  msg = 'LogLike is non-finite.';  return
    end
    passed = true;  msg = '';
end

function [passed, msg] = check_fit(data, flag, trueP, n)
% Fit model, check CI contains true params (loose check — may fail for small n).
% For large n the estimates should be close to the truth.
    res = GeneralizedHyperbolic_MLE(data, flag);
    [passed, msg] = check_model_validity(res, data);
    if ~passed, return; end
    if ~isfinite(res.LogLike)
        passed = false;  msg = 'LogLike not finite.';  return
    end
    % For n >= 1000 check that MLE is within 20% of true values.
    if n >= 1000
        est = res.Params;  % [mu; lambda; alpha; beta; delta]
        fields = {'mu','lambda','alpha','beta','delta'};
        trueVals = [trueP.mu, trueP.lambda, trueP.alpha, trueP.beta, trueP.delta];
        for i = 1:5
            tol = max(0.2*abs(trueVals(i)), 0.5);  % 20% or 0.5 absolute
            if abs(est(i) - trueVals(i)) > tol
                passed = false;
                msg = sprintf('%s estimate %.4f far from true %.4f (n=%d).', ...
                    fields{i}, est(i), trueVals(i), n);
                return
            end
        end
    end
    passed = true;  msg = '';
end

function [passed, msg] = check_runs_ok(data, flag)
% Pass if the fit completes without error, LogLike is finite, and
% the basic physical constraint |beta| < alpha holds at the MLE.
    try
        res = GeneralizedHyperbolic_MLE(data, flag);
    catch ME
        passed = false;  msg = ['Crashed: ' ME.message];  return
    end
    if ~isfinite(res.LogLike)
        passed = false;  msg = sprintf('LogLike = %g', res.LogLike);  return
    end
    [passed, msg] = check_model_validity(res, data);
    if ~passed, return; end
    mu_=res.Params(1); lam_=res.Params(2); alp_=res.Params(3);
    bet_=res.Params(4); del_=res.Params(5);
    if abs(bet_) >= alp_
        passed = false;
        msg = sprintf('Constraint violated: |beta|=%.4g >= alpha=%.4g', abs(bet_), alp_);
        return
    end
    if del_ <= 0
        passed = false;
        msg = sprintf('delta=%.4g not positive', del_);  return
    end
    passed = true;  msg = '';
end

function [passed, msg] = check_gh_improves_nig(data)
% GH (5 params) must achieve LogLike >= NIG (4 params) on the same data,
% since NIG is a restricted submodel of GH.
    try
        resGH  = GeneralizedHyperbolic_MLE(data, 'GH');
        resNIG = GeneralizedHyperbolic_MLE(data, 'NIG');
    catch ME
        passed = false;  msg = ['Crashed: ' ME.message];  return
    end
    if ~isfinite(resGH.LogLike) || ~isfinite(resNIG.LogLike)
        passed = false;
        msg = sprintf('Non-finite LogLike: GH=%.4g, NIG=%.4g', resGH.LogLike, resNIG.LogLike);
        return
    end
    [passed, msg] = check_model_validity(resGH, data);
    if ~passed, msg = ['GH: ' msg];  return; end
    [passed, msg] = check_model_validity(resNIG, data);
    if ~passed, msg = ['NIG: ' msg];  return; end
    % Allow a tiny numerical tolerance.
    if resGH.LogLike < resNIG.LogLike - 1e-4
        passed = false;
        msg = sprintf('GH LogLike (%.4f) < NIG LogLike (%.4f)', resGH.LogLike, resNIG.LogLike);
        return
    end
    passed = true;  msg = '';
end

function [passed, msg] = check_auto(data)
    res = GeneralizedHyperbolic_MLE(data, 'AUTO');
    if ~isfield(res,'ModelSelection')
        passed = false;  msg = 'Missing ModelSelection field.';  return
    end
    if ~isfinite(res.LogLike)
        passed = false;  msg = 'LogLike not finite.';  return
    end
    passed = true;  msg = '';
end

function [passed, msg] = check_ci_bounds(ci, est)
    if ~all(isfinite(ci))
        passed = false;  msg = sprintf('CI not finite: [%g, %g]', ci(1), ci(2));
    elseif ci(1) > ci(2)
        passed = false;  msg = sprintf('CI reversed: [%g, %g]', ci(1), ci(2));
    elseif est < ci(1) || est > ci(2)
        passed = false;  msg = sprintf('MLE %.4g not inside CI [%.4g, %.4g]', est, ci(1), ci(2));
    else
        passed = true;  msg = '';
    end
end

function [passed, msg] = check_ci_bounds_positive(ci, est)
    [passed, msg] = check_ci_bounds(ci, est);
    if passed && ci(1) <= 0
        passed = false;  msg = sprintf('Lower bound %.4g not positive', ci(1));
    end
end

% =========================================================================
%  RANDOM VARIATE GENERATORS
% =========================================================================

function x = sample_NIG(p, n)
% Generate n samples from NIG(mu, alpha, beta, delta) using the
% normal-variance-mean mixture representation:
%   X = mu + beta*V + sqrt(V)*Z
% where V ~ InverseGaussian(delta/gamma, delta^2) and Z ~ N(0,1).
    gam = p.gamma;
    mu_ig  = p.delta / gam;          % mean of the InvGaussian
    lam_ig = p.delta^2;              % shape of the InvGaussian
    V = sample_inverse_gaussian(mu_ig, lam_ig, n);
    Z = randn(n, 1);
    x = p.mu + p.beta*V + sqrt(V).*Z;
end

function v = sample_inverse_gaussian(mu, lam, n)
% Generate n samples from InvGaussian(mu, lambda) using the
% Michael-Schucany-Haas Wacker algorithm.
    y = randn(n,1).^2;
    x = mu + (mu^2 * y)/(2*lam) - (mu/(2*lam)) * sqrt(4*mu*lam*y + mu^2*y.^2);
    u = rand(n,1);
    v = x;
    idx = u > mu ./ (mu + x);
    v(idx) = mu^2 ./ x(idx);
end

function x = sample_GH_via_NIG(p, n)
% Generate approximate GH samples by sampling from NIG with the same
% (mu, alpha, beta, delta), ignoring the lambda parameter.
% This is sufficient for checking optimizer robustness but does NOT
% produce data with the target lambda — do not use for parameter recovery
% tests.  For exact GH sampling, a GIG variate generator is required.
    p_approx        = struct('mu', p.mu, 'lambda', -0.5, ...
                             'alpha', p.alpha, 'beta', p.beta, 'delta', p.delta);
    p_approx.gamma  = sqrt(p_approx.alpha^2 - p_approx.beta^2);
    x = sample_NIG(p_approx, n);
end

% =========================================================================
%  UNIVERSAL MODEL-VALIDITY INVARIANTS
% =========================================================================
function [passed, msg] = check_model_validity(res, data)
% Assert that a returned fit is a valid probability model.
%
% WHY THIS EXISTS.  Every check in this suite used to test either parameter
% closeness (recovery cases only) or "did not crash".  Neither catches a fit
% that converged cleanly onto numerical garbage.  A real failure looked like
% this: LogLike = +4.5e19, so AIC = 2k - 2*LogLike came out large and
% NEGATIVE, results.Fit was entirely NaN, and results.converged was true.
% isfinite(LogLike) is satisfied by 4.5e19, so the suite reported PASS.
%
% The invariants below are what "valid" actually means. They are cheap, they
% apply to every fit regardless of the test, and any one of them would have
% caught that failure.
    passed = false;
    P = res.Params;                 % [mu; lambda; alpha; beta; delta]
    mu = P(1); lam = P(2); al = P(3); be = P(4); de = P(5);

    if ~all(isfinite(P))
        msg = 'Params contain non-finite values.';  return
    end
    if ~(al > abs(be))
        msg = sprintf('Constraint violated: alpha (%.4g) must exceed |beta| (%.4g).', ...
            al, abs(be));  return
    end
    if ~(de > 0) || ~(al > 0)
        msg = sprintf('alpha (%.4g) and delta (%.4g) must be positive.', al, de);  return
    end
    if isfield(res, 'valid') && ~res.valid
        msg = 'results.valid is false (estimator flagged the fit as unusable).';  return
    end

    % The density must be evaluable at every observation.
    F = res.Fit(:);
    nBad = nnz(~isfinite(F) | F <= 0);
    if nBad > 0
        msg = sprintf('%d of %d fitted density values are NaN/Inf/non-positive.', ...
            nBad, numel(F));  return
    end

    % The reported likelihood must be the one the density implies.
    llFromFit = sum(log(F));
    if abs(llFromFit - res.LogLike) > max(1e-6, 1e-8*abs(res.LogLike))
        msg = sprintf('LogLike %.6g disagrees with sum(log(Fit)) = %.6g.', ...
            res.LogLike, llFromFit);  return
    end

    % AIC/BIC must follow from LogLike by their definitions.  A negative AIC
    % is not itself an error, but an AIC inconsistent with LogLike is.
    k = numel(res.theta_hat);  n = numel(data);

    % The Hessian must be k-by-k.  A k-by-1 here means an optimizer output was
    % read off the wrong position: fminunc returns 6 outputs ending
    % [grad, hessian], while fmincon returns 7 ending [lambda, grad, hessian].
    % Keeping fminunc's 6-output pattern after switching to fmincon silently
    % stores the GRADIENT as results.hessian.  Nothing else notices, because a
    % gradient at a converged solution is a small plausible-looking vector.
    if isfield(res, 'hessian') && ~isempty(res.hessian)
        sz = size(res.hessian);
        if ~isequal(sz, [k k])
            msg = sprintf(['hessian is %s but theta has %d elements; ' ...
                'expected %dx%d. A kx1 usually means the gradient was ' ...
                'captured instead.'], mat2str(sz), k, k, k);
            return
        end
    end
    if abs(res.AIC - (2*k - 2*res.LogLike)) > 1e-6
        msg = sprintf('AIC %.6g inconsistent with 2k-2LL = %.6g.', ...
            res.AIC, 2*k - 2*res.LogLike);  return
    end
    if abs(res.BIC - (k*log(n) - 2*res.LogLike)) > 1e-6
        msg = sprintf('BIC %.6g inconsistent with k*log(n)-2LL = %.6g.', ...
            res.BIC, k*log(n) - 2*res.LogLike);  return
    end

    % Cross-check against an INDEPENDENT density implementation.
    Fref = gh_pdf_ref(data(:), mu, lam, al, be, de);
    if ~all(isfinite(Fref))
        msg = 'Independent density evaluation produced non-finite values.';  return
    end
    relErr = max(abs(Fref - F) ./ max(F, realmin));
    if relErr > 1e-6
        msg = sprintf('Fit disagrees with independent density (max rel err %.3g).', ...
            relErr);  return
    end

    % The decisive test: a probability density integrates to one.  Degenerate
    % parameter vectors fail this no matter how good their reported LogLike.
    %
    % Integrate in STANDARDIZED coordinates u = (x - mu)/sc.  A naive
    % integral(f, -Inf, Inf) silently returns 0 whenever the density is narrow
    % relative to the axis it lives on: adaptive quadrature samples an infinite
    % domain coarsely, never lands on the peak, and concludes there is nothing
    % there.  That is a property of the quadrature, not of the fit, and it made
    % this check fail on perfectly good fits to data at mu = 1e6 or with a
    % scale of 1e-6.  After the substitution the integrand always has width of
    % order one, so the same routine resolves it.
    sc = max([de, 1/max(al, realmin), std(data), realmin]);
    try
        Z = integral(@(u) sc * gh_pdf_ref(mu + sc*u, mu, lam, al, be, de), ...
            -Inf, Inf, 'AbsTol', 1e-12, 'RelTol', 1e-10);
    catch ME
        msg = ['Density integration failed: ' ME.message];  return
    end
    if ~isfinite(Z) || abs(Z - 1) > 1e-3
        msg = sprintf(['Fitted density integrates to %.6g, not 1 ' ...
            '(mu=%.4g, alpha=%.4g, delta=%.4g, scale used %.4g).'], ...
            Z, mu, al, de, sc);  return
    end

    passed = true;  msg = '';
end

function f = gh_pdf_ref(x, mu, lambda, alpha, beta, delta)
% Independent GH density, written directly from Barndorff-Nielsen (1977):
%
%   f(x) = a * (delta^2+(x-mu)^2)^((lambda-1/2)/2)
%            * K_{lambda-1/2}(alpha*sqrt(delta^2+(x-mu)^2)) * exp(beta*(x-mu))
%   a    = (gamma/delta)^lambda / (sqrt(2*pi) * alpha^(lambda-1/2) * K_lambda(delta*gamma))
%
% Written separately from the estimator's internal routine on purpose: a fit
% cross-checked only against the code that produced it cannot reveal an error
% in that code.
    gam = sqrt(max(alpha^2 - beta^2, realmin));
    xm  = x - mu;
    z   = sqrt(delta^2 + xm.^2);

    % Scaled Bessel evaluations keep this usable for large arguments:
    % K(nu,y) = besselk(nu,y,1)*exp(-y).
    logA = lambda*log(gam/delta) - 0.5*log(2*pi) - (lambda - 0.5)*log(alpha) ...
         - (log(besselk(lambda, delta*gam, 1)) - delta*gam);
    logf = logA + (lambda - 0.5)*log(z) ...
         + (log(besselk(lambda - 0.5, alpha*z, 1)) - alpha*z) + beta*xm;
    f = exp(logf);
end

function [passed, msg] = check_no_degeneracy(data, flag)
% A fit must be a valid probability model AND its likelihood must be within
% the range a real density can produce.
    data = data(:);
    try
        res = GeneralizedHyperbolic_MLE(data, flag);
    catch ME
        passed = false;  msg = ['Crashed: ' ME.message];  return
    end

    [passed, msg] = check_model_validity(res, data);
    if ~passed, return; end

    % Guard against the specific symptom that motivated this section: a huge
    % positive LogLike, which drags AIC and BIC far negative.  A positive
    % LogLike is legitimate on its own (a sharply peaked density can exceed 1),
    % so the test is not "LogLike <= 0" but "LogLike is attainable": it can
    % never exceed n*log(max density observed).
    n  = numel(data);
    ub = n * log(max(res.Fit));
    if res.LogLike > ub + 1e-6
        passed = false;
        msg = sprintf('LogLike %.6g exceeds its upper bound n*log(max pdf) = %.6g.', ...
            res.LogLike, ub);
        return
    end

    % Parameters that have run away are unusable even if the arithmetic
    % happens to stay finite.  Scale-relative, so this holds for any data.
    s = std(data);
    if res.Params(3) > 1e8/max(s,realmin) || res.Params(5) > 1e8*s
        passed = false;
        msg = sprintf('Parameters ran away: alpha=%.4g, delta=%.4g for std(data)=%.4g.', ...
            res.Params(3), res.Params(5), s);
        return
    end

    passed = true;  msg = '';
end

function [passed, msg] = check_no_spurious_bounds(res)
% Open bounds are 0 or +/-Inf by contract. A finite bound of magnitude below
% 1e-100 or above 1e100 is an evaluation artifact presented as a result.
    passed = true;  msg = '';
    fns = fieldnames(res.profile_ci);
    for i = 1:numel(fns)
        b = res.profile_ci.(fns{i});
        bad = isfinite(b) & b ~= 0 & (abs(b) < 1e-100 | abs(b) > 1e100);
        if any(bad)
            passed = false;
            msg = sprintf('CI for %s is %s', fns{i}, mat2str(b, 4));
            return
        end
    end
end
