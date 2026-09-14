function H = FitPowerLawMLE(eventseries, xmin, varargin)
%FITPOWERLAWMLE Left-truncated power law (Pareto) fit, discrete by default.
%
%   H = FITPOWERLAWMLE(eventseries, xmin, SamplingInterval=dt, ...)
%
%   Parameter [alpha], the SURVIVAL exponent: S(t) ~ t^(-alpha). The
%   density exponent is alpha+1, so Clauset/Shalizi/Newman's alpha, which
%   is the density exponent, corresponds to alpha+1 here. Compare like with
%   like before quoting a number against the literature.
%
%   The Pareto scale is PINNED, not fitted. Under left truncation
%   S(t)/S(xmin) = (xmin/t)^alpha, in which a free scale cancels exactly,
%   so it is unidentified -- the same cancellation that un-identifies the
%   exponentiated Weibull's alpha when xmin sits deep in the tail.
%
%   THIS IS THE BINNED CONTINUOUS POWER LAW, NOT THE ZETA/ZIPF
%   DISTRIBUTION. Clauset et al. recommend zeta, P(n) ~ n^(-a)/zeta(a,nmin),
%   for discrete data; that is a different model with a different
%   log-likelihood. Binned Pareto is used here because every other family
%   in the registry is discretized the same way, and an AIC table may only
%   mix fits made under one convention. If you have a separate zeta-based
%   power-law fitter, its log-likelihood does NOT belong in the same table
%   as these.
%
%   See also FITTRUNCATEDDISCRETEMLE, FITPOWERLAWCUTOFFMLE.
H = FitTruncatedDiscreteMLE(eventseries, xmin, "powerlaw", varargin{:});
end
