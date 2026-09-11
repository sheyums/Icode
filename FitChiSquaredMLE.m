function H = FitChiSquaredMLE(eventseries, xmin, varargin)
%FITCHISQUAREDMLE Left-truncated chi-squared fit (discrete by default).
%
%   H = FITCHISQUAREDMLE(eventseries, xmin, SamplingInterval=dt, ...)
%
%   Parameter [nu]. chi2(nu) IS Gamma(nu/2, scale 2): the scale is FIXED at
%   2, which for durations measured in seconds is an arbitrary constant
%   with no physical meaning. This is therefore FITGAMMAMLE with one
%   parameter deleted for no reason, and it is provided only so a model
%   sweep can be reported as complete.
%
%   Expect it to lose badly rather than merely lose. On a real 3415-bout
%   sleep dataset it fell about 88000 nats below the free gamma, because no
%   value of nu can compensate for a scale pinned at 2 s. If you want a
%   scaled chi-squared, that is just the gamma -- use FITGAMMAMLE.
%
%   See also FITTRUNCATEDDISCRETEMLE, FITGAMMAMLE.
H = FitTruncatedDiscreteMLE(eventseries, xmin, "chisquared", varargin{:});
end
