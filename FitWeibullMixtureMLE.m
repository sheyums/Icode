function H = FitWeibullMixtureMLE(eventseries, xmin, varargin)
%FITWEIBULLMIXTUREMLE Left-truncated two-component Weibull mixture.
%
%   H = FITWEIBULLMIXTUREMLE(eventseries, xmin, SamplingInterval=dt, ...)
%
%   THE MODEL. Two Weibull components in the SCALE-SHAPE convention, with
%   a_i = scale_i > 0 and b_i = shape_i > 0:
%
%       S_i(t) = exp( -(t/a_i)^b_i )
%       f_i(t) = (b_i/a_i) (t/a_i)^(b_i - 1) S_i(t)
%       h_i(t) = (b_i/a_i) (t/a_i)^(b_i - 1)
%
%   so b_i = 1 is the exponential with constant hazard 1/a_i, b_i < 1 gives
%   a DECREASING component hazard and b_i > 1 an INCREASING one. That pair
%   is the whole point of the family. Mixed over weights w, w_2 = 1 - w_1:
%
%       S(t) = w_1 S_1(t) + w_2 S_2(t)
%       f(t) = w_1 f_1(t) + w_2 f_2(t)
%       F(t) = 1 - S(t)
%
%   The mixture hazard f/S is NOT the average of h_1 and h_2: it is their
%   average weighted by which component is still alive, and that is what
%   lets it fall, rise and fall again while each component's own hazard is
%   monotone.
%
%   WHAT IS MAXIMIZED is neither f nor F but the LEFT-TRUNCATED likelihood
%   of DISCRETE durations. With bin width dt and n_min = max(1,
%   round(xmin/dt)), the contribution of a bout of n bins is
%
%       P(N = n | N >= n_min) = [F(n*dt) - F((n-1)*dt)] / S((n_min-1)*dt)
%
%   and in DistributionType="continuous" mode, f(t)/S(xmin) instead. The
%   denominator is the only place xmin enters; see FITTRUNCATEDDISCRETEMLE.
%
%   REPORTED PARAMETERS: [q1, scale1, shape1, q2, scale2, shape2] -- SIX
%   numbers while k = 5, because sum(q) = 1 leaves one weight determined
%   but both worth printing. Components are sorted by scale.
%
%   The q are the OBSERVED weights, q_i proportional to w_i * S_i(xmin) --
%   component i's share of the bouts that SURVIVED truncation -- not the
%   mixing weights w above, which are what the likelihood is written in.
%   The two differ whenever a component's mass sits near or below xmin.
%   H.Theta holds w; H.Params holds q. Never read one as the other.
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
