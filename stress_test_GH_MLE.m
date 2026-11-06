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
    'single observation',  @() expect_error(@() GeneralizedHyperbolic_MLE(1))
    'constant data',       @() expect_error_id(@() GeneralizedHyperbolic_MLE(5*ones(20,1)), 'GeneralizedHyperbolic_MLE:ConstantData')
    'NaN in data',         @() expect_error(@() GeneralizedHyperbolic_MLE([1 2 NaN 4]))
    'Inf in data',         @() expect_error(@() GeneralizedHyperbolic_MLE([1 2 Inf 4]))
    'bad modelFlag',       @() expect_error(@() GeneralizedHyperbolic_MLE(randn(100,1),'XYZ'))
    'numeric modelFlag',   @() expect_error(@() GeneralizedHyperbolic_MLE(randn(100,1), 99))
    'matrix data',         @() expect_error(@() GeneralizedHyperbolic_MLE(randn(10,10)))
    % A 5-parameter model needs at least 10 observations (Nov 2026; 2 used
    % to be accepted). 9 must error, 10 must run.
    '2 observations rejected', @() expect_error_id(@() GeneralizedHyperbolic_MLE([0;1],'NIG'), 'GeneralizedHyperbolic_MLE:TooFewObservations')
    '9 observations rejected', @() expect_error_id(@() GeneralizedHyperbolic_MLE((1:9)','NIG'), 'GeneralizedHyperbolic_MLE:TooFewObservations')
    '10 observations allowed', @() expect_no_error(@() GeneralizedHyperbolic_MLE([-0.93 0.30 0.20 0.81 0.38 -0.45 -0.20 -0.06 -0.93 -0.29]','NIG'))
};

for i = 1:size(tests1,1)
    [p,~] = run_test(tests1{i,1}, tests1{i,2});
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
    [p,~] = run_test(label, testFn);
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
    [p,~] = run_test(label, testFn);
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
[pf,~] = run_test('NIG near-symmetric (beta~0)', ...
    @() check_finite_and_ordered(sample_NIG(p_sym,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Near-boundary skewness: |beta|/alpha = 0.95
p_skew = struct('mu',0,'lambda',-0.5,'alpha',4,'beta',3.8,'delta',1);
p_skew.gamma = sqrt(p_skew.alpha^2 - p_skew.beta^2);
[pf,~] = run_test('NIG high skew (|beta|/alpha=0.95)', ...
    @() check_finite_and_ordered(sample_NIG(p_skew,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Negative skew
p_negskew = struct('mu',0,'lambda',-0.5,'alpha',4,'beta',-3.8,'delta',1);
p_negskew.gamma = sqrt(p_negskew.alpha^2 - p_negskew.beta^2);
[pf,~] = run_test('NIG negative high skew', ...
    @() check_finite_and_ordered(sample_NIG(p_negskew,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Very small delta
p_sdelta = struct('mu',0,'lambda',-0.5,'alpha',5,'beta',1,'delta',0.01);
p_sdelta.gamma = sqrt(p_sdelta.alpha^2 - p_sdelta.beta^2);
[pf,~] = run_test('NIG tiny delta=0.01', ...
    @() check_finite_and_ordered(sample_NIG(p_sdelta,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Very large delta
p_ldelta = struct('mu',0,'lambda',-0.5,'alpha',2,'beta',0.5,'delta',50);
p_ldelta.gamma = sqrt(p_ldelta.alpha^2 - p_ldelta.beta^2);
[pf,~] = run_test('NIG large delta=50', ...
    @() check_finite_and_ordered(sample_NIG(p_ldelta,n),'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Large location offset (mu=1e6)
p_largemu = struct('mu',1e6,'lambda',-0.5,'alpha',3,'beta',0.5,'delta',1);
p_largemu.gamma = sqrt(p_largemu.alpha^2 - p_largemu.beta^2);
[pf,~] = run_test('NIG large mu=1e6', ...
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
    [pf,~] = run_test(label, @() check_runs_ok(data,'GH'));
    nFail  = nFail+~pf; nPass = nPass+pf;
end

% =========================================================================
%  SECTION 6: DATA PATHOLOGIES
% =========================================================================
fprintf('\n=== SECTION 6: Data Pathologies ===\n');

base_data = sample_NIG(trueParams_NIG, 500);

% A few extreme outliers (3 observations at 20 sigma)
data_out = [base_data; mean(base_data) + 20*std(base_data)*[1;-1;1]];
[pf,~] = run_test('NIG with 3 extreme outliers', ...
    @() check_finite_and_ordered(data_out,'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% Near-normal data (very large alpha → GH approaches Normal)
data_norm = randn(500,1);  % pure Normal
[pf,~] = run_test('Near-normal data (NIG fit)', ...
    @() check_finite_and_ordered(data_norm,'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% GH fit of Normal data: lambda profile may be flat/unbounded — that is
% correct behavior, not a bug.  Just check it runs and LogLike is finite.
[pf,~] = run_test('Near-normal data (GH fit, no crash)', ...
    @() check_runs_ok(data_norm,'GH'));
nFail = nFail+~pf; nPass = nPass+pf;

% Heavily skewed data (log-normal, which GH cannot perfectly fit)
data_lnorm = exp(0.5*randn(500,1));
[pf,~] = run_test('Log-normal data (NIG fit, misspecified)', ...
    @() check_finite_and_ordered(data_lnorm,'NIG'));
nFail = nFail+~pf; nPass = nPass+pf;

% =========================================================================
%  SECTION 7: AUTO MODEL SELECTION
% =========================================================================
fprintf('\n=== SECTION 7: AUTO Model Selection ===\n');

% True NIG data — AUTO should select NIG (or at least not crash)
data_NIG = sample_NIG(trueParams_NIG, 1000);
[pf,~] = run_test('AUTO on NIG data — expect NIG or GH, no crash', ...
    @() check_auto(data_NIG));
nFail = nFail+~pf; nPass = nPass+pf;

% True GH data with lambda far from -0.5
% NOTE: sample_GH_via_NIG ignores lambda (see its help), so this exercises the
% optimizer on GH-shaped input rather than testing lambda recovery.
trueParams_GH = struct('mu',0,'lambda',1.5,'alpha',3,'beta',1,'delta',1);
trueParams_GH.gamma = sqrt(trueParams_GH.alpha^2 - trueParams_GH.beta^2);
data_GH = sample_GH_via_NIG(trueParams_GH, 2000);
[pf,~] = run_test('AUTO on GH data (lambda=1.5, n=2000)', ...
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

[pf,~] = run_test('CI: mu bounds finite and ordered', ...
    @() check_ci_bounds(ci.mu, params(1)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('CI: alpha bounds finite, ordered, and positive', ...
    @() check_ci_bounds_positive(ci.alpha, params(3)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('CI: beta bounds finite and ordered', ...
    @() check_ci_bounds(ci.beta, params(4)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('CI: delta bounds finite, ordered, and positive', ...
    @() check_ci_bounds_positive(ci.delta, params(5)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('CI: alpha lower > 0 strictly', ...
    @() check_condition(ci.alpha(1) > 0, 'alpha lower CI not positive'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('CI: alpha CI contains true alpha=3', ...
    @() check_condition(ci.alpha(1)<=3 && ci.alpha(2)>=3, ...
    sprintf('true alpha=3 not in [%.3f, %.3f]',ci.alpha(1),ci.alpha(2))));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('CI: |beta| < alpha at MLE (constraint)', ...
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
    [pf,~] = run_test(['GH has field: ' f], ...
        @() check_condition(isfield(res_gh,f), ['missing field: ' f]));
    nFail = nFail+~pf; nPass = nPass+pf;
end

[pf,~] = run_test('AUTO has ModelSelection field', ...
    @() check_condition(isfield(res_auto,'ModelSelection'),'missing ModelSelection'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('AUTO has LRT field', ...
    @() check_condition(isfield(res_auto,'LRT'),'missing LRT'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('Params has 5 elements', ...
    @() check_condition(numel(res_gh.Params)==5,'Params not length 5'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('Fit length matches data length', ...
    @() check_condition(numel(res_gh.Fit)==200,'Fit length mismatch'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('FitX is sorted ascending', ...
    @() check_condition(all(diff(res_gh.FitX)>=0),'FitX not sorted'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('FitPDF all non-negative', ...
    @() check_condition(all(res_gh.FitPDF(isfinite(res_gh.FitPDF))>=0), ...
    'negative density values'));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('LogLike is finite scalar', ...
    @() check_condition(isscalar(res_gh.LogLike) && isfinite(res_gh.LogLike), ...
    sprintf('LogLike = %g', res_gh.LogLike)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('NIG lambda fixed at -0.5', ...
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
        [pf,~] = run_test(sprintf('%s [%s]', degenCases{i,1}, flag{1}), ...
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

[pf,~] = run_test('11A: open delta bound on VG-limit data is exactly 0', ...
    @() check_condition(resVG.profile_ci.delta(1) == 0, ...
    sprintf('delta lower bound = %.4g (an open bound must be 0, not a tiny number)', ...
    resVG.profile_ci.delta(1))));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('11B: no CI bound is a spurious tiny or huge finite number', ...
    @() check_no_spurious_bounds(resVG));
nFail = nFail+~pf; nPass = nPass+pf;

% 11C: VG data with delta exactly 0 (gamma mixing), large n: the fit goes to
% the floor. The assertion is conditional, so it states the rule itself.
V11c = -log(rand(rs11, 1000, 1));
xVGc = 0.5*V11c + sqrt(V11c).*randn(rs11, 1000, 1);
resVGc = GeneralizedHyperbolic_MLE(xVGc, 'GH');
floor11 = 1e-3*std(xVGc);
[pf,~] = run_test('11C: delta within 1% of its floor => delta_at_bound', ...
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
[pf,~] = run_test('11D: parameters are equivariant under 1e4*x + 1e4', ...
    @() check_condition(max(abs(r11b.Params - Pm) ./ abs(Pm)) < 1e-5, ...
    sprintf('max relative error %.3g', max(abs(r11b.Params - Pm) ./ abs(Pm)))));
nFail = nFail+~pf; nPass = nPass+pf;
dLL = r11b.LogLike - (r11a.LogLike - 800*log(1e4));
[pf,~] = run_test('11D: log-likelihood shifts by exactly -n*log(scale)', ...
    @() check_condition(abs(dLL) < 1e-6, sprintf('discrepancy %.3g', dLL)));
nFail = nFail+~pf; nPass = nPass+pf;
ciA = [1e4*r11a.profile_ci.mu + 1e4, r11a.profile_ci.alpha/1e4, ...
       r11a.profile_ci.beta/1e4, r11a.profile_ci.delta*1e4];
ciB = [r11b.profile_ci.mu, r11b.profile_ci.alpha, r11b.profile_ci.beta, r11b.profile_ci.delta];
ciErr = max(abs(ciB - ciA) ./ abs(ciA));
[pf,~] = run_test('11D: CI bounds are equivariant (to root-finder tolerance)', ...
    @() check_condition(ciErr < 5e-3, sprintf('max relative error %.3g (was 0.17 before the fix)', ciErr)));
nFail = nFail+~pf; nPass = nPass+pf;

% 11E: single precision input
r11s = GeneralizedHyperbolic_MLE(single(x11), 'NIG');
[pf,~] = run_test('11E: single input is fitted in double precision', ...
    @() check_condition(isa(r11s.LogLike, 'double') && ...
    abs(r11s.LogLike - GeneralizedHyperbolic_MLE(double(single(x11)), 'NIG').LogLike) < 1e-9, ...
    sprintf('LogLike class %s', class(r11s.LogLike))));
nFail = nFail+~pf; nPass = nPass+pf;

% 11F: ten points near the Normal limit. A sound fit must not fail its own
% validity check because of the log-density's floating-point noise: the fixed
% 1e-6 tolerance reported valid = false here (Bessel argument ~6e9).
x11f = [-0.92891230072865449 0.29770359610432284 0.19808604340722599 0.80962243337342477 0.37782371522682678 -0.44532913476610414 -0.19631927565757104 -0.057295620412519732 -0.92631446155746755 -0.29055067417114161]';
r11f = GeneralizedHyperbolic_MLE(x11f, 'AUTO');
[pf,~] = run_test('11F: tiny sample near the Normal limit stays valid', ...
    @() check_condition(r11f.valid, sprintf('valid = %d, LogLike = %.8g', r11f.valid, r11f.LogLike)));
nFail = nFail+~pf; nPass = nPass+pf;

% =========================================================================
%  SECTION 12: REGRESSION GUARDS FROM THE NOV 2026 AUDIT (8 findings)
% =========================================================================
% One guard per finding, each asserting a RULE that a correct estimator
% satisfies on any data, not where one optimizer run happens to land.
%
%   12A  AUTO preferred an invalid fit if its LogLike/BIC won.
%   12B  The profile CIs never checked that the point estimate was optimal
%        (new field profile_refits).
%   12C  AUTO computed full CIs for both models and discarded one, and the GH
%        fit re-ran the NIG multistart.
%   12D  The skewed Student-t limit (lambda < 0, gamma -> 0) went unreported.
%   12E  alpha_core's gamma lost all accuracy from xi = 19 (local function):
%        covered ONLY BY INSPECTION -- no public output isolates it.
%   12F  MAX_BESSEL_ARG ignored n (local function): tested through behaviour,
%        the noise bound n*argMax*eps of a near-Normal large-n fit.
%   12G  n = 2 was accepted (Section 1 holds the n = 9 / n = 10 boundary).
%   12H  Validity was judged on exp(logpdf), which underflows (LogFit).
fprintf('\n=== SECTION 12: Nov 2026 Audit Regression Guards ===\n');
rs12 = RandStream('twister', 'Seed', 12);

% ---- 12C first: it produces the fits several later guards reuse ----------
% AUTO must return CIs (computed for the selected model only), must reach the
% same GH and NIG point estimates as the stand-alone fits, and must cost less
% than running both stand-alone fits.
x12 = sample_NIG_rs(rs12, trueParams_NIG, 600);
t0 = tic; r12gh  = GeneralizedHyperbolic_MLE(x12, 'GH');   tGH   = toc(t0);
t0 = tic; r12nig = GeneralizedHyperbolic_MLE(x12, 'NIG');  tNIG  = toc(t0);
t0 = tic; r12au  = GeneralizedHyperbolic_MLE(x12, 'AUTO'); tAUTO = toc(t0);
fprintf('  [timing] n=600: GH %.1f s, NIG %.1f s, GH+NIG %.1f s, AUTO %.1f s\n', ...
    tGH, tNIG, tGH + tNIG, tAUTO);

[pf,~] = run_test('12C: AUTO result carries computed profile CIs', ...
    @() check_condition(isfield(r12au, 'profile_ci') && ...
    ~all(isnan(cell2mat(struct2cell(r12au.profile_ci))), 'all'), 'profile_ci missing or all NaN'));
nFail = nFail+~pf; nPass = nPass+pf;

% Point estimates are deterministic, so AUTO's LRT log-likelihoods equal the
% stand-alone fits' unless a profile refit moved one of them.
okGH  = abs(r12au.LRT.LogLike_GH  - r12gh.LogLike)  < 1e-9 || r12gh.profile_refits  > 0 || r12au.profile_refits > 0;
okNIG = abs(r12au.LRT.LogLike_NIG - r12nig.LogLike) < 1e-9 || r12nig.profile_refits > 0 || r12au.profile_refits > 0;
[pf,~] = run_test('12C: AUTO reaches the stand-alone GH and NIG optima', ...
    @() check_condition(okGH && okNIG, sprintf('dGH = %.3g, dNIG = %.3g', ...
    r12au.LRT.LogLike_GH - r12gh.LogLike, r12au.LRT.LogLike_NIG - r12nig.LogLike)));
nFail = nFail+~pf; nPass = nPass+pf;

[pf,~] = run_test('12C: AUTO costs less than GH + NIG run separately', ...
    @() check_condition(tAUTO < tGH + tNIG, ...
    sprintf('AUTO %.1f s versus GH+NIG %.1f s', tAUTO, tGH + tNIG)));
nFail = nFail+~pf; nPass = nPass+pf;

% ---- 12A: AUTO's selection rule -------------------------------------------
% The returned model may be invalid only when both were. The branch where
% exactly one fit is invalid is covered by inspection: no input is known
% that makes one model's fit invalid and not the other's.
autoRes = {r12au, r11f};
okA = true;  msgA = '';
for i = 1:numel(autoRes)
    r = autoRes{i};
    if ~(r.valid || contains(r.ModelSelection, 'no fit is valid'))
        okA = false;  msgA = r.ModelSelection;
    end
    % The reported LRT must agree with the returned model's LogLike.
    llSel = r.LRT.(['LogLike_' strtok(r.ModelSelection)]);
    if abs(llSel - r.LogLike) > 1e-9
        okA = false;  msgA = sprintf('LRT LogLike %.10g versus result %.10g', llSel, r.LogLike);
    end
end
[pf,~] = run_test('12A: AUTO returns an invalid fit only if all are invalid', ...
    @() check_condition(okA, msgA));
nFail = nFail+~pf; nPass = nPass+pf;

% ---- 12B: profile refits ----------------------------------------------------
fitsB = {r12gh, r12nig, r12au};
okB = all(cellfun(@(r) isfield(r, 'profile_refits') && ismember(r.profile_refits, 0:2), fitsB));
[pf,~] = run_test('12B: profile_refits reported, an integer in 0..2', ...
    @() check_condition(okB, 'profile_refits missing or out of range'));
nFail = nFail+~pf; nPass = nPass+pf;

% ---- 12D: Student-t limit -------------------------------------------------
% The flag must agree with its definition, checked with the suite's OWN
% density: with gamma pushed down by exp(-20) at fixed (lambda, beta, delta,
% mu), the log-likelihood drops by less than 1e-3 nats, and lambda < 0.
%
% NOT asserted: that Student-t data RAISE the flag. Measured before writing
% this (Nov 2026), it does not hold: t(3) samples of n = 1000 fitted as GH
% (4 seeds), t(1.5) as GH, and exact Cauchy samples fitted as NIG (whose
% t-limit IS the Cauchy, 3 seeds) all landed at an interior gamma of order
% 1/max|x| (0.11-0.19 for t(3), 7e-4 to 3e-3 for Cauchy) and none was
% flagged. A finite sample has no data beyond max|x|, so a little exponential
% damping there is free in-sample and buys normalizing mass: the optimum
% usually sits at gamma ~ 1/max|x|. The likelihood IS flat at tiny gamma; the
% MLE is usually just not there. In the first full run of this section ONE of
% the three fits below (a t(3) GH fit, lambda = -1.55) was flagged, at a cost
% of 0.00085 nats -- so whether it fires turns on a fraction of a millinat,
% and the positive branch IS exercised here on this RandStream. The count is
% printed, not asserted.
tFits = cell(1, 3);
for s = 1:2
    rsT = RandStream('twister', 'Seed', 100 + s);
    xt  = trnd_rs(rsT, 3, 1000);
    tFits{s} = struct('x', xt, 'r', GeneralizedHyperbolic_MLE(xt, 'GH'));
end
rsT = RandStream('twister', 'Seed', 103);
xt  = tan(pi*(rand(rsT, 1000, 1) - 0.5));     % standard Cauchy
tFits{3} = struct('x', xt, 'r', GeneralizedHyperbolic_MLE(xt, 'NIG'));
nFlag = 0;  okD = true;  msgD = '';
for s = 1:3
    r = tFits{s}.r;  xt = tFits{s}.x;
    [ok, why] = student_t_flag_consistent(r, xt);
    if ~ok, okD = false; msgD = why; end
    nFlag = nFlag + r.student_t_limit;
end
for c = {r12gh, r12nig}
    [ok, why] = student_t_flag_consistent(c{1}, x12);
    if ~ok, okD = false; msgD = why; end
end
fprintf('  [info] student_t_limit raised on %d of 3 heavy-tailed fits (t(3) GH x2, Cauchy NIG)\n', nFlag);
[pf,~] = run_test('12D: student_t_limit agrees with an independent re-evaluation', ...
    @() check_condition(okD, msgD));
nFail = nFail+~pf; nPass = nPass+pf;

% ---- 12F: the Bessel-argument cap scales with n ---------------------------
% Near-Normal data push the Bessel arguments toward the cap. Whatever the fit
% does, its NLL noise n*argMax*eps must stay at or under 1e-2 nats (with a
% fixed cap of 1e10 it could reach 0.044 here). argMax from the parameters.
xN = randn(rs12, 20000, 1);
rN = GeneralizedHyperbolic_MLE(xN, 'NIG');
P  = rN.Params;  gN = sqrt(max(P(3)^2 - P(4)^2, 0));
argMax = max([P(5)*gN; P(3)*sqrt(P(5)^2 + (xN - P(1)).^2)]);
noiseN = numel(xN) * argMax * eps;
[pf,~] = run_test('12F: near-Normal n=20000 fit keeps n*argMax*eps <= 1e-2', ...
    @() check_condition(noiseN <= 1e-2*(1 + 1e-6), ...
    sprintf('n*argMax*eps = %.3g (argMax %.3g)', noiseN, argMax)));
nFail = nFail+~pf; nPass = nPass+pf;
[pf,~] = run_test('12F: that fit is a valid probability model', ...
    @() check_model_validity(rN, xN));
nFail = nFail+~pf; nPass = nPass+pf;

% ---- 12G: small samples warn ----------------------------------------------
% The warning is promoted to an error so it is caught by identifier, before
% any fitting, whatever other warnings the fit would raise afterwards.
wstate = warning('error', 'GeneralizedHyperbolic_MLE:SmallSample');
xG = randn(rs12, 30, 1);
[pf,~] = run_test('12G: n=30 raises SmallSample', ...
    @() expect_error_id(@() GeneralizedHyperbolic_MLE(xG, 'NIG'), ...
    'GeneralizedHyperbolic_MLE:SmallSample'));
warning(wstate);
nFail = nFail+~pf; nPass = nPass+pf;

% ---- 12H: validity is judged on the log-density ----------------------------
% Rule on every fit above: LogFit exists, LogLike = sum(LogFit), Fit =
% exp(LogFit), and valid is never false merely because Fit underflowed.
allH = [{r12gh, r12nig, r12au, rN}, cellfun(@(c) c.r, tFits, 'UniformOutput', false)];
okH = true;  msgH = '';
for i = 1:numel(allH)
    r = allH{i};
    if ~isfield(r, 'LogFit')
        okH = false;  msgH = 'no LogFit field';  break
    end
    if abs(sum(r.LogFit) - r.LogLike) > max(1e-6, 1e-8*abs(r.LogLike)) || ...
            any(abs(exp(r.LogFit) - r.Fit) > 0)
        okH = false;  msgH = sprintf('fit %d: LogFit inconsistent', i);  break
    end
end
[pf,~] = run_test('12H: LogFit reported and consistent with LogLike and Fit', ...
    @() check_condition(okH, msgH));
nFail = nFail+~pf; nPass = nPass+pf;

% ---- 12S: stress pass -------------------------------------------------------
% Many seeds and sizes, heavy tails, near-Normal data, ties and outliers.
% Asserts invariants only: a valid probability model, ordered CIs with no
% UNKNOWN (NaN) bound, the flags consistent with their definitions, and, on
% two cases, location-scale equivariance.
fprintf('  -- stress pass --\n');
stress = {};
for n = [60 250 800]
    for s = 1:2
        rsS = RandStream('twister', 'Seed', 1000*n + s);
        stress(end+1, :) = {sprintf('NIG data n=%d seed %d', n, s), sample_NIG_rs(rsS, trueParams_NIG, n), 'NIG'}; %#ok<SAGROW>
    end
end
for nu = [1.5 4]
    rsS = RandStream('twister', 'Seed', round(10*nu));
    stress(end+1, :) = {sprintf('t(%.1f) data n=500', nu), trnd_rs(rsS, nu, 500), 'NIG'}; %#ok<SAGROW>
end
rsS = RandStream('twister', 'Seed', 77);
stress(end+1, :) = {'near-Normal n=400 [GH]', randn(rsS, 400, 1), 'GH'};
stress(end+1, :) = {'heavily tied n=600', round(2*sample_NIG_rs(rsS, trueParams_NIG, 600)), 'NIG'};
stress(end+1, :) = {'outliers at 50 sd n=400', [sample_NIG_rs(rsS, trueParams_NIG, 397); 50*[1; -1; 1]], 'NIG'};
stress(end+1, :) = {'t(2.5) data n=600 [GH]', trnd_rs(rsS, 2.5, 600), 'GH'};
for i = 1:size(stress, 1)
    [pf,~] = run_test(['12S: ' stress{i,1}], @() check_stress_case(stress{i,2}, stress{i,3}));
    nFail = nFail+~pf; nPass = nPass+pf;
end
for i = [3 7]   % NIG n=250 seed 1, t(1.5)
    xs = stress{i,2};
    [pf,~] = run_test(['12S: equivariance under 1e3*x - 7, ' stress{i,1}], ...
        @() check_equivariance(xs, 1e3, -7, stress{i,3}));
    nFail = nFail+~pf; nPass = nPass+pf;
end

% =========================================================================
%  SECTION 13: HYPERBOLIC MODEL, STANDARD ERRORS, WALD CIs (Nov 2026)
% =========================================================================
% Ground truth that shares no code with the estimator: the hyperbolic
% (lambda = 1) density in closed form,
%   f(x) = gamma/(2*alpha*delta*K1(delta*gamma))
%          * exp(-alpha*sqrt(delta^2 + (x-mu)^2) + beta*(x-mu)),
% which needs only besselk(1, .), an independent Nelder-Mead multistart on
% it, and an independent finite-difference Hessian in PHYSICAL coordinates.
%   13A  'HYP' and 'HYPERBOLIC' fit lambda = 1 and agree.
%   13B  HYP LogLike equals the closed-form density sum.
%   13C  HYP reaches the global optimum of the closed form.
%   13D  GH nests both: LogLike_GH >= max(LogLike_NIG, LogLike_HYP).
%   13E  AUTO's restricted model is the better valid one; LRT p-values agree.
%   13F  SEs equal an independent physical-coordinate Hessian's.
%   13G  Wald CIs: formula, positivity, contain the estimate.
%   13H  'NoProfile' returns the same fit and SEs, no profiles, faster.
%   13I  SEs equivariant under c*x + b.
%   13J  At the delta floor, no Wald SE is reported (all NaN, SEValid false).
%   13K  A crossing made by the Bessel-argument guard is reported open.
fprintf('\n=== SECTION 13: Hyperbolic model, SEs, Wald CIs ===\n');
rs13 = RandStream('twister', 'Seed', 13);
pHYP = struct('mu', 0.3, 'lambda', 1, 'alpha', 2, 'beta', 0.6, 'delta', 0.8);
x13  = sample_HYP_rs(rs13, pHYP, 800);

t0 = tic; r13  = GeneralizedHyperbolic_MLE(x13, 'HYP');               tHYPfull = toc(t0);
r13b = GeneralizedHyperbolic_MLE(x13, 'hyperbolic');
t0 = tic; r13w = GeneralizedHyperbolic_MLE(x13, 'HYP', 'NoProfile');  tHYPwald = toc(t0);

[pf,~] = run_test('13A: HYP fits lambda = 1; HYPERBOLIC alias identical', ...
    @() check_condition(r13.Params(2) == 1 && isequal(r13.profile_ci.lambda, [1 1]) && ...
    startsWith(r13.Model, 'Hyperbolic') && isequal(r13.Params, r13b.Params) && ...
    r13.LogLike == r13b.LogLike, 'lambda, Model or alias mismatch'));
nFail = nFail+~pf; nPass = nPass+pf;

llRef = hyp_ll_ref(x13, r13.Params);
[pf,~] = run_test('13B: HYP LogLike equals the closed-form density sum', ...
    @() check_condition(abs(llRef - r13.LogLike) < 1e-8*abs(llRef), ...
    sprintf('closed form %.12g versus %.12g', llRef, r13.LogLike)));
nFail = nFail+~pf; nPass = nPass+pf;

llInd = hyp_independent_max(x13, rs13);
[pf,~] = run_test('13C: HYP reaches the independent closed-form optimum', ...
    @() check_condition(llInd <= r13.LogLike + 1e-4, ...
    sprintf('independent search higher by %.3g', llInd - r13.LogLike)));
nFail = nFail+~pf; nPass = nPass+pf;

xN = sample_NIG_rs(rs13, trueParams_NIG, 600);
okD = true;  msgD = '';  autoR = {};
for d = {x13, xN}
    ra = GeneralizedHyperbolic_MLE(d{1}, 'AUTO', 'NoProfile');
    autoR{end+1} = ra; %#ok<SAGROW>
    gap = max(ra.LRT.LogLike_NIG, ra.LRT.LogLike_HYP) - ra.LRT.LogLike_GH;
    if gap > 1e-6
        okD = false;  msgD = sprintf('a nested model beats GH by %.3g', gap);
    end
end
[pf,~] = run_test('13D: GH LogLike >= NIG and hyperbolic LogLike', ...
    @() check_condition(okD, msgD));
nFail = nFail+~pf; nPass = nPass+pf;

okE = true;  msgE = '';
for i = 1:numel(autoR)
    L = autoR{i}.LRT;
    want = 'NIG';  if L.LogLike_HYP > L.LogLike_NIG, want = 'HYP'; end
    pN = chi2cdf(max(2*(L.LogLike_GH - L.LogLike_NIG), 0), 1, 'upper');
    pH = chi2cdf(max(2*(L.LogLike_GH - L.LogLike_HYP), 0), 1, 'upper');
    if ~strcmp(L.Restricted, want) || abs(L.pValue_NIG - pN) > 1e-12 || ...
            abs(L.pValue_HYP - pH) > 1e-12 || ...
            abs(L.pValue - L.(['pValue_' want])) > 1e-12
        okE = false;  msgE = sprintf('Restricted %s, expected %s', L.Restricted, want);
    end
end
[pf,~] = run_test('13E: AUTO tests GH against the better restricted model', ...
    @() check_condition(okE, msgE));
nFail = nFail+~pf; nPass = nPass+pf;

seRef = hyp_se_ref(x13, r13.Params);
relF  = max(abs(r13.SE([1 3 4 5]) - seRef) ./ seRef);
[pf,~] = run_test('13F: SEs match an independent physical-coordinate Hessian', ...
    @() check_condition(r13.SEValid && relF < 1e-3 && r13.SE(2) == 0, ...
    sprintf('max relative SE difference %.3g', relF)));
nFail = nFail+~pf; nPass = nPass+pf;

zc = 1.959963984540054;  W = r13.CIWald;  P = r13.Params;  S = r13.SE;
okG = max(abs(W.mu - (P(1) + [-1 1]*zc*S(1)))) < 1e-12*max(1, abs(P(1))) && ...
    max(abs(W.beta - (P(4) + [-1 1]*zc*S(4)))) < 1e-12 && ...
    max(abs(log(W.alpha) - (log(P(3)) + [-1 1]*zc*S(3)/P(3)))) < 1e-12 && ...
    max(abs(log(W.delta) - (log(P(5)) + [-1 1]*zc*S(5)/P(5)))) < 1e-12 && ...
    W.alpha(1) > 0 && W.delta(1) > 0 && isequal(W.lambda, [1 1]);
names13 = {'mu','alpha','beta','delta'};  iP13 = [1 3 4 5];
for j = 1:4
    okG = okG && W.(names13{j})(1) <= P(iP13(j)) && P(iP13(j)) <= W.(names13{j})(2);
end
[pf,~] = run_test('13G: Wald CIs follow their formulas and contain the estimate', ...
    @() check_condition(okG, 'Wald interval formula or ordering violated'));
nFail = nFail+~pf; nPass = nPass+pf;

okH = isequal(r13w.Params, r13.Params) || r13.profile_refits > 0;
okH = okH && max(abs(r13w.SE - r13.SE)) <= 1e-12 + 1e-9*max(abs(r13.SE)) || r13.profile_refits > 0;
okH = okH && strcmp(r13w.CIMethod, 'wald') && strcmp(r13.CIMethod, 'profile') && ...
    all(isnan([r13w.profile_ci.mu r13w.profile_ci.alpha r13w.profile_ci.beta r13w.profile_ci.delta])) && ...
    tHYPwald < tHYPfull;
[pf,~] = run_test('13H: NoProfile gives the same fit and SEs, no profiles, faster', ...
    @() check_condition(okH, sprintf('NoProfile %.1f s versus %.1f s', tHYPwald, tHYPfull)));
nFail = nFail+~pf; nPass = nPass+pf;

okI = true;  msgI = '';
for c = [1e-3 -250]
    rc = GeneralizedHyperbolic_MLE(c*x13 + 40, 'HYP', 'NoProfile');
    want = r13w.SE .* [abs(c); 1; 1/abs(c); 1/abs(c); abs(c)];
    rel = max(abs(rc.SE - want) ./ max(want, eps));
    if ~(rel < 1e-4)
        okI = false;  msgI = sprintf('c = %g: relative SE difference %.3g', c, rel);
    end
end
[pf,~] = run_test('13I: SEs equivariant under c*x + 40 (c = 1e-3, -250)', ...
    @() check_condition(okI, msgI));
nFail = nFail+~pf; nPass = nPass+pf;

% Laplace data: the hyperbolic delta -> 0 limit IS the (asymmetric) Laplace,
% so the fit sits on the delta floor.
xL = -log(rand(rs13, 3000, 1)) .* sign(rand(rs13, 3000, 1) - 0.5);
prevW = warning('off', 'GeneralizedHyperbolic_MLE:DeltaAtBound');
rL = GeneralizedHyperbolic_MLE(xL, 'HYP', 'NoProfile');
warning(prevW);
[pf,~] = run_test('13J: at the delta floor, no Wald SEs (NaN, SEValid false)', ...
    @() check_condition(rL.delta_at_bound && ~rL.SEValid && all(isnan(rL.SE([1 3 4 5]))) && ...
    rL.SE(2) == 0 && all(isnan([rL.CIWald.alpha rL.CIWald.delta rL.CIWald.mu])), ...
    sprintf('delta_at_bound %d, SE = %s', rL.delta_at_bound, mat2str(rL.SE', 3))));
nFail = nFail+~pf; nPass = nPass+pf;

% 13K: a crossing made by the Bessel-argument guard is an OPEN bound. On this
% hyperbolic sample (coverage replicate 2038) alpha and beta upper bounds
% used to read ~8.5e8 -- where the guard starts, not where the likelihood
% rises.
rsK = RandStream('mt19937ar', 'Seed', 2038);
xK  = sample_HYP_rs(rsK, struct('mu', 0, 'lambda', 1, 'alpha', 2, 'beta', 0.5, 'delta', 1), 500);
rK  = GeneralizedHyperbolic_MLE(xK, 'HYP');
hiA = rK.profile_ci.alpha(2);  hiB = rK.profile_ci.beta(2);
okK = ~(isfinite(hiA) && hiA > 1e6*rK.Params(3)) && ~(isfinite(hiB) && hiB > 1e6*abs(rK.Params(4))) && ...
    isfinite(rK.profile_ci.alpha(1)) && rK.profile_ci.alpha(1) < rK.Params(3);
[pf,~] = run_test('13K: guard-made crossings reported open, not as ~1e9 bounds', ...
    @() check_condition(okK, sprintf('alpha upper %.3g, beta upper %.3g', hiA, hiB)));
nFail = nFail+~pf; nPass = nPass+pf;

% =========================================================================
%  SECTION 14: SKEWNESS OUTPUTS (Nov 2026)
% =========================================================================
fprintf('\n--- Section 14: skewness outputs ---\n');
prevW14 = warning('off', 'GeneralizedHyperbolic_MLE:NoGainOverNormal');
rs14 = RandStream('mt19937ar', 'Seed', 1401);
P14  = struct('mu', 0, 'alpha', 2, 'beta', 0.8, 'delta', 1);  P14.gamma = sqrt(P14.alpha^2 - P14.beta^2);
x14  = sample_NIG_rs(rs14, P14, 800);
rN14 = GeneralizedHyperbolic_MLE(x14, 'NIG', 'NoProfile');
rH14 = GeneralizedHyperbolic_MLE(x14, 'HYP', 'NoProfile');

% 14A: NIG moment skewness has the closed form 3*beta/(alpha*sqrt(delta*gamma)).
q = rN14.Params;  gA = sqrt(q(3)^2 - q(4)^2);
gRef = 3*q(4) / (q(3)*sqrt(q(5)*gA));
[pf,~] = run_test('14A: NIG Skewness equals 3*beta/(alpha*sqrt(delta*gamma))', ...
    @() check_condition(abs(rN14.Skewness - gRef) < 1e-9*max(1, abs(gRef)), ...
    sprintf('%.12g versus %.12g', rN14.Skewness, gRef)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14B: hyperbolic moment skewness against GIG moments integrated numerically
% (no Bessel functions, so an independent route).
gRef = skew_by_integration(rH14.Params);
[pf,~] = run_test('14B: HYP Skewness matches GIG moments by numerical integration', ...
    @() check_condition(abs(rH14.Skewness - gRef) < 1e-7*max(1, abs(gRef)), ...
    sprintf('%.12g versus %.12g', rH14.Skewness, gRef)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14C: model Bowley against quartiles of a brute-force CDF on a fine grid.
bRef = bowley_by_grid(rH14.Params);
[pf,~] = run_test('14C: SkewnessBowleyModel matches a brute-force CDF', ...
    @() check_condition(abs(rH14.SkewnessBowleyModel - bRef) < 1e-6, ...
    sprintf('%.10g versus %.10g', rH14.SkewnessBowleyModel, bRef)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14D: data Bowley equals quantile()'s on tie-free data, and the hand-worked
% mid-quantile value on tied data (Fmid = .1 .4 .7 .9 at 1 2 3 10, quartiles
% 1.5, 7/3, 4.75).
qd  = quantile(x14, [0.25 0.5 0.75]);  bQ = (qd(3) + qd(1) - 2*qd(2)) / (qd(3) - qd(1));
prevT = warning('off', 'all');
rT  = GeneralizedHyperbolic_MLE(repelem([1 2 2 3 10], 2), 'NIG', 'NoProfile');
warning(prevT);
bHand = (1.5 + 4.75 - 2*7/3) / 3.25;
[pf,~] = run_test('14D: SkewnessBowleyData = quantile() without ties, mid-quantiles with ties', ...
    @() check_condition(abs(rN14.SkewnessBowleyData - bQ) < 1e-12 && ...
    abs(rT.SkewnessBowleyData - bHand) < 1e-12, ...
    sprintf('%.12g vs %.12g; tied %.12g vs %.12g', rN14.SkewnessBowleyData, bQ, rT.SkewnessBowleyData, bHand)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14E: SkewnessSE equals the delta method through Cov (an independent path:
% Cov is in the physical parameters, the code works in theta).
gfun = @(v) skew_closed(v);
gr = zeros(5, 1);  v0 = rH14.Params;
for i = [1 3 4 5]
    e = zeros(5, 1);  e(i) = 1e-5 * max(1, abs(v0(i)));
    gr(i) = (gfun(v0 + e) - gfun(v0 - e)) / (2*e(i));
end
seRef = sqrt(gr.' * rH14.Cov * gr);
[pf,~] = run_test('14E: SkewnessSE equals the delta method through Cov', ...
    @() check_condition(abs(rH14.SkewnessSE - seRef) < 1e-4*seRef, ...
    sprintf('%.8g versus %.8g', rH14.SkewnessSE, seRef)));
nFail = nFail+~pf; nPass = nPass+pf;

% 14F: equivariance. A shift and a rescaling by 1e-6 leave every skewness
% output unchanged; negating the data flips their signs. Also the regression
% for the Wald fix: at scale 1e-6 the SEs used to come back NaN (the
% data-unit information failed the rcond test); they must scale by 1e-6.
rS = GeneralizedHyperbolic_MLE(1e-6*x14 - 7, 'NIG', 'NoProfile');
rM = GeneralizedHyperbolic_MLE(-x14, 'NIG', 'NoProfile');
f1 = [rN14.Skewness rN14.SkewnessSE rN14.SkewnessBowleyModel rN14.SkewnessBowleyData];
fS = [rS.Skewness   rS.SkewnessSE   rS.SkewnessBowleyModel   rS.SkewnessBowleyData];
fM = [rM.Skewness   rM.SkewnessSE   rM.SkewnessBowleyModel   rM.SkewnessBowleyData];
okF = all(abs(fS - f1) < 1e-5*max(1, abs(f1))) && ...
    all(abs(fM - f1.*[-1 1 -1 -1]) < 1e-5*max(1, abs(f1))) && rS.SEValid && ...
    abs(rS.SE(1)/rN14.SE(1) - 1e-6) < 1e-9 && abs(rS.SE(5)/rN14.SE(5) - 1e-6) < 1e-9;
[pf,~] = run_test('14F: skewness outputs equivariant; Wald SEs valid at scale 1e-6', ...
    @() check_condition(okF, sprintf('orig %s scaled %s negated %s', ...
    mat2str(f1, 6), mat2str(fS, 6), mat2str(fM, 6))));
nFail = nFail+~pf; nPass = nPass+pf;
warning(prevW14);

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

function [passed, msg] = expect_error_id(fn, id)
% Pass if fn throws an error with identifier ID.
    try
        fn();
        passed = false;
        msg = sprintf('Expected error %s but none was thrown.', id);
    catch ME
        passed = strcmp(ME.identifier, id);
        msg = '';
        if ~passed
            msg = sprintf('Expected error %s, got %s: %s', id, ME.identifier, ME.message);
        end
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
    alp_=res.Params(3);
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

    % The density must be evaluable at every observation. Judged on the
    % LOG-density when the estimator reports it: a legitimate log-density
    % below about -745 underflows to Fit = 0 without anything being wrong.
    F = res.Fit(:);
    if isfield(res, 'LogFit')
        LF = res.LogFit(:);
    else
        LF = log(F);
    end
    nBad = nnz(~isfinite(LF));
    if nBad > 0
        msg = sprintf('%d of %d fitted log-density values are NaN/Inf.', ...
            nBad, numel(F));  return
    end

    % The reported likelihood must be the one the density implies.
    llFromFit = sum(LF);
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
    %
    % Centre on the data median, not on mu: a strongly skewed fit puts its
    % mass far from mu (a near-Normal NIG fit at n = 20000 had mu = -212 with
    % its mass near 0, about 50 scale units away, and the mu-centred integral
    % returned 0 for a density that integrates to 1). Any fit that explains
    % the data has mass at the data, so the median is always under the peak.
    sc = max([de, 1/max(al, realmin), std(data), realmin]);
    c  = median(data);
    try
        Z = integral(@(u) sc * gh_pdf_ref(c + sc*u, mu, lam, al, be, de), ...
            -Inf, Inf, 'AbsTol', 1e-12, 'RelTol', 1e-10);
    catch ME
        msg = ['Density integration failed: ' ME.message];  return
    end
    if ~isfinite(Z) || abs(Z - 1) > 1e-3
        msg = sprintf(['Fitted density integrates to %.6g, not 1 ' ...
            '(mu=%.4g, alpha=%.4g, delta=%.4g, centre %.4g, scale %.4g).'], ...
            Z, mu, al, de, c, sc);  return
    end

    passed = true;  msg = '';
end

function [f, logf] = gh_pdf_ref(x, mu, lambda, alpha, beta, delta, gamOverride)
% Independent GH density, written directly from Barndorff-Nielsen (1977):
%
%   f(x) = a * (delta^2+(x-mu)^2)^((lambda-1/2)/2)
%            * K_{lambda-1/2}(alpha*sqrt(delta^2+(x-mu)^2)) * exp(beta*(x-mu))
%   a    = (gamma/delta)^lambda / (sqrt(2*pi) * alpha^(lambda-1/2) * K_lambda(delta*gamma))
%
% Written separately from the estimator's internal routine on purpose: a fit
% cross-checked only against the code that produced it cannot reveal an error
% in that code.
    % Optional gamOverride: give gamma directly (alpha is then hypot(beta,
    % gamma)), because alpha^2 - beta^2 cancels to nothing once gamma is
    % ~1e-8 of |beta|, which is exactly where the Student-t check probes.
    if nargin >= 7
        gam = gamOverride;  alpha = hypot(beta, gam);
    else
        gam = sqrt(max(alpha^2 - beta^2, realmin));
    end
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

% =========================================================================
%  SECTION 12 HELPERS
% =========================================================================
function x = sample_NIG_rs(rs, p, n)
% sample_NIG drawing from RandStream rs, so Section 12 does not depend on
% how many global draws the sections before it consumed.
    gam = p.gamma;
    m   = p.delta / gam;   lam = p.delta^2;
    y   = randn(rs, n, 1).^2;
    v   = m + (m^2*y)/(2*lam) - (m/(2*lam))*sqrt(4*m*lam*y + m^2*y.^2);
    u   = rand(rs, n, 1);
    flip = u > m ./ (m + v);
    v(flip) = m^2 ./ v(flip);
    x = p.mu + p.beta*v + sqrt(v).*randn(rs, n, 1);
end

function t = trnd_rs(rs, nu, n)
% Student-t variates by Bailey's (1994) polar method, exact for any nu > 0,
% with no toolbox call: (U, V) uniform on the unit disc, W = U^2 + V^2,
% T = U*sqrt(nu*(W^(-2/nu) - 1)/W).
    t = zeros(n, 1);  k = 0;
    while k < n
        U = 2*rand(rs, 2*n, 1) - 1;  V = 2*rand(rs, 2*n, 1) - 1;
        W = U.^2 + V.^2;
        ok = W > 0 & W <= 1;
        T = U(ok) .* sqrt(nu*(W(ok).^(-2/nu) - 1) ./ W(ok));
        take = min(numel(T), n - k);
        t(k+1:k+take) = T(1:take);  k = k + take;
    end
end

function [ok, why] = student_t_flag_consistent(r, x)
% student_t_limit must equal its definition, re-evaluated with this suite's
% own density: lambda < 0 and the log-likelihood drops by less than 1e-3
% when gamma is multiplied by exp(-20) at fixed (lambda, beta, delta, mu).
% Within 1e-6 of the 1e-3 threshold either answer is accepted.
    ok = true;  why = '';
    P = r.Params;  lam = P(2);
    if ~isfield(r, 'student_t_limit')
        ok = false;  why = 'no student_t_limit field';  return
    end
    if lam >= 0
        ok = ~r.student_t_limit;
        if ~ok, why = sprintf('flagged with lambda = %.3g >= 0', lam); end
        return
    end
    iPsi = numel(r.theta_hat) - 2;       % psi is 3rd from the end in both models
    gam  = exp(r.theta_hat(iPsi));
    [~, l0] = gh_pdf_ref(x(:), P(1), lam, [], P(4), P(5), gam);
    [~, l1] = gh_pdf_ref(x(:), P(1), lam, [], P(4), P(5), gam*exp(-20));
    cost = sum(l0) - sum(l1);
    if ~isfinite(cost)
        expected = false;
    else
        expected = cost < 1e-3;
        if abs(cost - 1e-3) < 1e-6, return; end
    end
    ok = (expected == r.student_t_limit);
    if ~ok
        why = sprintf('student_t_limit = %d but the independent cost is %.4g nats (lambda %.3g)', ...
            r.student_t_limit, cost, lam);
    end
end

function [passed, msg] = check_stress_case(x, flag)
% Invariants every fit must satisfy, whatever the data.
    t0 = tic;
    res = GeneralizedHyperbolic_MLE(x, flag);
    fprintf('        [%s n=%d] lambda %.3g, refits %d, t-limit %d, delta_at_bound %d, %.0f s\n', ...
        flag, numel(x), res.Params(2), res.profile_refits, res.student_t_limit, ...
        res.delta_at_bound, toc(t0));
    [passed, msg] = check_model_validity(res, x);
    if ~passed, return; end
    passed = false;
    fns = fieldnames(res.profile_ci);
    for i = 1:numel(fns)
        b = res.profile_ci.(fns{i});
        if any(isnan(b))
            msg = sprintf('CI for %s has an UNKNOWN (NaN) bound: %s', fns{i}, mat2str(b, 4));  return
        end
        if b(1) > b(2)
            msg = sprintf('CI for %s is reversed: %s', fns{i}, mat2str(b, 4));  return
        end
    end
    if ~ismember(res.profile_refits, 0:2)
        msg = sprintf('profile_refits = %g', res.profile_refits);  return
    end
    [ok, why] = student_t_flag_consistent(res, x);
    if ~ok, msg = why; return; end
    if ~(res.Params(5) > 1.01*1e-3*std(x) || res.delta_at_bound)
        msg = 'delta within 1% of its floor but delta_at_bound is false';  return
    end
    passed = true;  msg = '';
end

function [passed, msg] = check_equivariance(x, c, sh, flag)
% Location-scale equivariance under y = c*x + sh. On the Student-t limit,
% alpha (hence gamma) is not identified and is excluded from the comparison.
    ra = GeneralizedHyperbolic_MLE(x, flag);
    rb = GeneralizedHyperbolic_MLE(c*x + sh, flag);
    passed = false;
    dLL = rb.LogLike - (ra.LogLike - numel(x)*log(c));
    if abs(dLL) > 1e-6
        msg = sprintf('LogLike discrepancy %.3g', dLL);  return
    end
    Pm = [c*ra.Params(1) + sh; ra.Params(2); ra.Params(3)/c; ra.Params(4)/c; ra.Params(5)*c];
    use = true(5, 1);
    if ra.student_t_limit || rb.student_t_limit, use(3) = false; end
    rel = abs(rb.Params - Pm) ./ max(abs(Pm), realmin);
    if max(rel(use)) > 1e-5
        msg = sprintf('parameter relative error %.3g', max(rel(use)));  return
    end
    ciA = [c*ra.profile_ci.mu + sh, ra.profile_ci.beta/c, ra.profile_ci.delta*c];
    ciB = [rb.profile_ci.mu, rb.profile_ci.beta, rb.profile_ci.delta];
    if ~(ra.student_t_limit || rb.student_t_limit)
        ciA = [ciA, ra.profile_ci.alpha/c];  ciB = [ciB, rb.profile_ci.alpha];
    end
    fin = isfinite(ciA) & ciA ~= 0;
    if ~isequal(isfinite(ciA), isfinite(ciB))
        msg = 'an open bound in one fit is finite in the other';  return
    end
    err = max(abs(ciB(fin) - ciA(fin)) ./ abs(ciA(fin)));
    if ~isempty(err) && err > 5e-3
        msg = sprintf('CI relative error %.3g', err);  return
    end
    passed = true;  msg = '';
end

% =========================================================================
%  SECTION 13 HELPERS: closed-form hyperbolic ground truth
% =========================================================================
function x = sample_HYP_rs(rs, p, n)
% Hyperbolic variates, X = mu + beta*W + sqrt(W)*Z with W ~ GIG(1, delta^2,
% gamma^2), W drawn by inverting its CDF on a fine log-grid (error ~1e-7).
    gam = sqrt(p.alpha^2 - p.beta^2);
    lf  = @(u) u - 0.5*(p.delta^2*exp(-u) + gam^2*exp(u));   % log-density of log W
    u   = linspace(-60, 60, 200001);  L = lf(u);  m = max(L);
    k   = find(L > m - 60);
    u   = linspace(u(max(k(1)-1, 1)), u(min(k(end)+1, numel(u))), 400001);
    F   = cumtrapz(u, exp(lf(u) - m));  F = F / F(end);
    [Fu, iu] = unique(F);
    W   = exp(interp1(Fu, u(iu), rand(rs, n, 1), 'pchip'));
    x   = p.mu + p.beta*W + sqrt(W).*randn(rs, n, 1);
end

function ll = hyp_ll_ref(x, P)
% Closed-form hyperbolic log-likelihood; P = [mu; lambda; alpha; beta; delta].
    mu = P(1);  a = P(3);  b = P(4);  d = P(5);  g = sqrt(a^2 - b^2);
    ll = numel(x)*(log(g) - log(2*a*d) - log(besselk(1, d*g, 1)) + d*g) ...
        + sum(-a*sqrt(d^2 + (x - mu).^2) + b*(x - mu));
end

function nll = hyp_nll_phys(x, q)
% NLL in physical coordinates q = [mu; alpha; beta; delta]; Inf if infeasible.
    if q(2) <= abs(q(3)) || q(4) <= 0, nll = Inf; return; end
    nll = -hyp_ll_ref(x, [q(1); 1; q(2); q(3); q(4)]);
end

function llMax = hyp_independent_max(x, rs)
% Independent 30-start Nelder-Mead on the closed form, in its own
% unconstrained coordinates [mu, log(alpha-|beta|)... ] -- here
% [mu, log(gamma), beta, log(delta)].
    s = std(x);  m = median(x);  best = Inf;
    f = @(v) hyp_nll_phys(x, [v(1); hypot(v(3), exp(v(2))); v(3); exp(v(4))]);
    o = optimset('Display', 'off', 'MaxFunEvals', 8000, 'MaxIter', 8000, ...
        'TolX', 1e-10, 'TolFun', 1e-10);
    for k = 1:30
        v0 = [m + 0.5*s*randn(rs); log((0.3 + 2*rand(rs))/s); (rand(rs) - 0.5)/s; ...
            log((0.1 + 2*rand(rs))*s)];
        v = fminsearch(f, v0, o);  v = fminsearch(f, v, o);
        best = min(best, f(v));
    end
    llMax = -best;
end

function se = hyp_se_ref(x, P)
% SEs of [mu; alpha; beta; delta] from a central-difference Hessian of the
% closed-form NLL in PHYSICAL coordinates (relative steps 1e-4), inverted.
    q = [P(1); P(3); P(4); P(5)];
    s = std(x);
    h = 1e-4 * max(abs(q), 0.1*[s; 1/s; 1/s; s]);
    f = @(v) hyp_nll_phys(x, v);
    H = zeros(4);
    for i = 1:4
        for j = i:4
            ei = zeros(4,1); ei(i) = h(i);  ej = zeros(4,1); ej(j) = h(j);
            H(i,j) = (f(q+ei+ej) - f(q+ei-ej) - f(q-ei+ej) + f(q-ei-ej)) / (4*h(i)*h(j));
            H(j,i) = H(i,j);
        end
    end
    se = sqrt(diag(inv(H)));
end

function g = skew_closed(v)
% Moment skewness from Params v = [mu; lambda; alpha; beta; delta] by the
% Bessel closed form (used only as the function differentiated in 14E).
    gam = sqrt(v(3)^2 - v(4)^2);  w = v(5)*gam;  s = v(5)/gam;  L = v(2);
    K = @(nu) besselk(nu, w);
    m1 = s*K(L+1)/K(L);  m2 = s^2*K(L+2)/K(L);  m3 = s^3*K(L+3)/K(L);
    k2 = m2 - m1^2;  k3 = m3 - 3*m1*m2 + 2*m1^3;
    g  = (3*v(4)*k2 + v(4)^3*k3) / (m1 + v(4)^2*k2)^1.5;
end

function g = skew_by_integration(v)
% Moment skewness of GH with Params v, from central moments of the GIG mixing
% variable integrated numerically in u = log V (no Bessel functions).
    gam = sqrt(v(3)^2 - v(4)^2);  L = v(2);  d = v(5);  b = v(4);
    lf  = @(u) L*u - 0.5*(d^2*exp(-u) + gam^2*exp(u));
    u0  = fminbnd(@(u) -lf(u), -50, 50);  c = lf(u0);
    % Finite limits where the log-density has fallen 700 below its peak; an
    % infinite range evaluates (exp(u) - m1)^3 * 0 = Inf*0 far in the tail.
    a   = fzero(@(u) lf(u) - c + 700, [u0 - 50, u0]);
    z   = fzero(@(u) lf(u) - c + 700, [u0, u0 + 50]);
    f   = @(u, r) exp(r*u + lf(u) - c);
    o   = {'AbsTol', 0, 'RelTol', 1e-13};
    Z   = integral(@(u) f(u, 0), a, z, o{:});
    m1  = integral(@(u) f(u, 1), a, z, o{:}) / Z;
    k2  = integral(@(u) (exp(u) - m1).^2 .* f(u, 0), a, z, o{:}) / Z;
    k3  = integral(@(u) (exp(u) - m1).^3 .* f(u, 0), a, z, o{:}) / Z;
    g   = (3*b*k2 + b^3*k3) / (m1 + b^2*k2)^1.5;
end

function b = bowley_by_grid(v)
% Bowley skewness of GH with Params v from a trapezoid CDF on a fine grid.
    mu = v(1);  L = v(2);  al = v(3);  be = v(4);  de = v(5);
    gam = sqrt(al^2 - be^2);  w = 60/(al - abs(be)) + 60*de;
    xs = linspace(mu - w, mu + w, 4e6);  q = sqrt(de^2 + (xs - mu).^2);
    l  = L*log(gam/de) - (L - 0.5)*log(al) - log(besselk(L, de*gam)) + ...
         (L - 0.5)*log(q) + log(besselk(L - 0.5, al*q, 1)) - al*q + be*(xs - mu);
    F  = cumtrapz(xs, exp(l));  F = F / F(end);
    [Fu, iu] = unique(F);
    qq = interp1(Fu, xs(iu), [0.25 0.5 0.75]);
    b  = (qq(3) + qq(1) - 2*qq(2)) / (qq(3) - qq(1));
end
