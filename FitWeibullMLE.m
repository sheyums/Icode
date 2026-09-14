function H = FitWeibullMLE(eventseries, xmin, varargin)
%FITWEIBULLMLE Left-truncated Weibull fit (discrete by default).
%
%   H = FITWEIBULLMLE(eventseries, xmin, SamplingInterval=dt, ...)
%
%   Parameters [scale lambda, shape k], f(t) ~ t^(k-1)*exp(-(t/lambda)^k).
%   k < 1 gives a decreasing hazard and a stretched-exponential tail --
%   heavier than exponential, lighter than a power law.
%
%   This is scipy's weibull_min. It is also FITEXPONENTIATEDWEIBULLMLE with
%   FixAlpha=true; that route additionally gives the likelihood-ratio test
%   of alpha=1 against the exponentiated Weibull, which this one cannot.
%   The two agree on the fit itself.
%
%   See also FITTRUNCATEDDISCRETEMLE, FITEXPONENTIATEDWEIBULLMLE.
H = FitTruncatedDiscreteMLE(eventseries, xmin, "weibull", varargin{:});
end
