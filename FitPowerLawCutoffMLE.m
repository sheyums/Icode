function H = FitPowerLawCutoffMLE(eventseries, xmin, varargin)
%FITPOWERLAWCUTOFFMLE Left-truncated power law with an exponential cutoff.
%
%   H = FITPOWERLAWCUTOFFMLE(eventseries, xmin, SamplingInterval=dt, ...)
%
%   Fits f(t) ~ t^(-alpha) * exp(-t/tau) and reports [alpha, tau].
%
%   THIS IS THE GAMMA, RE-PARAMETRIZED, NOT A SEPARATE MODEL. Gamma's
%   density t^(k-1)*exp(-t/theta) is identical to the above with
%   alpha = 1-k and tau = theta, and under left truncation the two are the
%   same distribution: truncation deletes (0,xmin), the only region their
%   normalisations differ over. Verified on real data to 3e-10 in
%   log-likelihood. So this calls FITGAMMAMLE and relabels; the fit, the
%   log-likelihood, the parameter count and every information criterion are
%   identical.
%
%   It exists because (alpha, tau) is the interpretable parametrization for
%   bout data -- "scale-free with a cutoff at tau" says more than "gamma
%   with shape 0.134" -- and because a model sweep usually wants that row
%   named. DO NOT enter both this and the gamma in one comparison table and
%   treat them as independent evidence: they are one model with two names,
%   and H.EquivalentTo records that.
%
%   Note alpha < 1 here means the untruncated density is not integrable at
%   0; that is untroubling under left truncation, which never evaluates it
%   there, but it does mean alpha should not be read as a power-law
%   exponent valid down to zero duration.
%
%   See also FITGAMMAMLE, FITPOWERLAWMLE, FITTRUNCATEDDISCRETEMLE.

H = FitTruncatedDiscreteMLE(eventseries, xmin, "gamma", varargin{:});
H.Model = 'powerlaw_cutoff';
H.ParamNames = {'alpha','cutoff'};
% alpha = 1 - shape, cutoff = scale. d(alpha)/d(shape) = -1, so the
% standard error carries over unchanged.
H.Params = [1 - H.Params(1), H.Params(2)];
H.ParamSE = [H.ParamSE(1), H.ParamSE(2)];
H.EquivalentTo = 'gamma (alpha = 1 - shape, cutoff = scale)';
end
