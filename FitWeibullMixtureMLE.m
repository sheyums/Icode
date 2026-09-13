function H = FitWeibullMixtureMLE(eventseries, xmin, varargin)
%FITWEIBULLMIXTUREMLE Left-truncated two-component Weibull mixture.
%
%   H = FITWEIBULLMIXTUREMLE(eventseries, xmin, SamplingInterval=dt, ...)
%
%   Fits w*Weibull(scale1,shape1) + (1-w)*Weibull(scale2,shape2), reporting
%   [w1, scale1, shape1, scale2, shape2]. Components are sorted by scale.
%
%   WHY THIS FAMILY EXISTS. Every other model in this package has a
%   MONOTONE hazard, and a hyperexponential has a strictly decreasing one
%   at any order K -- a sum of decreasing exponential terms cannot be
%   anything else. So data whose hazard falls, rises, then falls again lie
%   outside all of them by construction, and no number of exponential
%   components will help: extra components refine a monotone shape, they do
%   not create a hump. On per0 DD WAKE bouts, K=3, K=4 and K=5 together
%   bought 0.33 nats over K=2 while the goodness-of-fit test rejected every
%   candidate in the library, which is what that looks like from inside a
%   model comparison.
%
%   This reaches the three regimes with two components: shape1 < 1 gives a
%   steeply falling hazard at short durations, shape2 > 1 contributes a
%   rising one in the middle, and whichever component is heavier-tailed
%   dominates again at the top. An exponential first component will not
%   substitute -- its hazard is flat, so the mixture starts flat instead of
%   falling.
%
%   WHAT IT NESTS, and why that is useful:
%     shape1 = shape2 = 1   the 2-component hyperexponential
%     shape1 = 1            exponential + Weibull
%     identical components  a single Weibull
%   shape = 1 is an INTERIOR point of shape > 0 and the model remains
%   identified there, so "is the hazard monotone?" is a REGULAR hypothesis:
%   test shape2 = 1 by an ordinary likelihood ratio against chi2(1). That
%   is unlike the mixture-ORDER question, where the null sits on a boundary
%   and a ridge at once and needs HYPEREXPONENTIALLRT.
%
%   FIVE PARAMETERS ON A TRUNCATED SAMPLE IS A LOT. Check the guard before
%   reading the estimates: H.Diagnostics.GuardOK is false when the weight
%   has reached a boundary (one component holding fewer than
%   MinMixtureCount observations, default 5) or when the two components
%   have converged on the same law. Either way the extra parameters are
%   unidentified, the standard errors are meaningless, and the 1-component
%   weibull is the honest fit.
%
%   See also FITTRUNCATEDDISCRETEMLE, FITWEIBULLMLE, FITHYPEREXPONENTIALMLE.

H = FitTruncatedDiscreteMLE(eventseries, xmin, "weibull_mix", varargin{:});
end
