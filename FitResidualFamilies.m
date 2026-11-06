function T = FitResidualFamilies(e)
%FITRESIDUALFAMILIES  Fit eight candidate noise families to one residual vector.
%
%   T = FitResidualFamilies(eps)
%
% For the question "what is the distribution of eps?" (e.g. the residuals of
% BoutAR1Noise), fits every family below to the SAME data by maximum
% likelihood and returns comparable log-likelihoods. All are densities on
% the real line with respect to Lebesgue measure, so their LogLike, AIC and
% BIC compare directly.
%
%   Family         K   Fitter
%   ShiftedExp     2   closed form: loc = min(eps), scale = mean - min. The
%                      shape eps has when the bouts are exponential and
%                      there are no dynamics (or EAR(1) noise).
%   Normal         2   closed form.
%   Laplace        2   closed form: loc = median, scale = mean |eps - median|.
%   AsymLaplace    3   AsymLaplace_MLE(eps, 'NOPROFILE')
%   ShiftLogNormal 3   shiftlognormal_MLE(eps)
%   ShiftWeibull   3   shiftweibull_MLE(eps, 'NOPROFILE')
%   GenNorm        3   GenNorm_MLE(eps, 1)  (one component)
%   GH             4-5 GeneralizedHyperbolic_MLE(eps, 'AUTO', 'NoProfile'):
%                      NIG or hyperbolic (4) or full GH (5), whichever AUTO
%                      picks; the pick is in Note.
%
% OUTPUT  table, one row per family: Family, K, LogLike, AIC (true AIC,
% 2K - 2 LogLike), BIC (K log n - 2 LogLike), DeltaBIC (from the best row),
% Note (the error message when a fit failed; then LogLike etc. are NaN).
%
% CAVEATS
%   - ShiftedExp's location is the sample minimum, a boundary MLE: its
%     likelihood is not regular, so information criteria involving it are
%     conventional, not exact. Calibrate any decision on simulation (see
%     test_FitResidualFamilies) rather than on a |dBIC| rule of thumb.
%   - MEASURED (test_FitResidualFamilies, 2026-11-05): on 100 simulated
%     residual sets that ARE shifted exponential (integer seconds, n = 180),
%     ShiftedExp won by BIC in 100 of 100 (every other family scored worse),
%     so a flexible family beating it on real data is not chance
%     flexibility. The boundary MLE favours ShiftedExp somewhat: this
%     shows the comparison is not biased AGAINST it, not that it is fair.
%   - eps from a regression are not independent and not identically
%     distributed when the model is wrong; the likelihoods treat them as iid.
%
% Needs the toolboxes the six fitters need (Optimization, Statistics).
%
% CHANGE LOG
%   2026-11-05  First version.

e = double(e(:));
e = e(isfinite(e));
n = numel(e);
fam = ["ShiftedExp";"Normal";"Laplace";"AsymLaplace";"ShiftLogNormal";"ShiftWeibull";"GenNorm";"GH"];
K = [2;2;2;3;3;3;3;NaN];
LL = NaN(8, 1);
note = strings(8, 1);

ws = warning('off', 'all');
c = onCleanup(@() warning(ws));

% closed forms
s = mean(e) - min(e);
if s > 0, LL(1) = -n*log(s) - n; end
v = mean((e - mean(e)).^2);
if v > 0, LL(2) = -n/2 * (log(2*pi*v) + 1); end
b = mean(abs(e - median(e)));
if b > 0, LL(3) = -n*log(2*b) - n; end

calls = {@() AsymLaplace_MLE(e, 'NOPROFILE'), @() shiftlognormal_MLE(e), ...
         @() shiftweibull_MLE(e, 'NOPROFILE'), @() GenNorm_MLE(e, 1)};
for i = 1:4
    try
        r = calls{i}();
        LL(3 + i) = r(1).LogLike;
        if isfield(r, 'NParameters') && ~isempty(r(1).NParameters)
            K(3 + i) = r(1).NParameters;
        end
    catch err
        note(3 + i) = string(err.message);
    end
end
try
    r = GeneralizedHyperbolic_MLE(e, 'AUTO', 'NoProfile');
    LL(8) = r.LogLike;
    m = string(r.Model);
    % Model strings: "Normal-Inverse Gaussian (...)", "Hyperbolic (mu,
    % lambda = 1, ...)", "Generalized Hyperbolic (...)". Only the last
    % estimates lambda. (Matching "hyperbolic (" would also catch it.)
    if startsWith(m, "Generalized")
        K(8) = 5;
    elseif startsWith(m, "Normal-Inverse") || startsWith(m, "Hyperbolic")
        K(8) = 4;
    else
        error('FitResidualFamilies:GHModel', 'Unrecognised GH model "%s".', m);
    end
    note(8) = m;
catch err
    note(8) = string(err.message);
end

AIC = 2*K - 2*LL;
BIC = K*log(n) - 2*LL;
T = table(fam, K, LL, AIC, BIC, BIC - min(BIC), note, ...
    'VariableNames', {'Family','K','LogLike','AIC','BIC','DeltaBIC','Note'});
end
