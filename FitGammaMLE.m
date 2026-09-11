function H = FitGammaMLE(eventseries, xmin, varargin)
%FITGAMMAMLE Left-truncated gamma fit (discrete by default).
%
%   H = FITGAMMAMLE(eventseries, xmin, SamplingInterval=dt, ...)
%
%   Parameters [shape k, scale theta], density f(t) ~ t^(k-1)*exp(-t/theta).
%   k < 1 gives a decreasing hazard, which bout durations generally want;
%   k = 1 is the exponential; k > 1 an increasing hazard.
%
%   THIS IS ALSO THE POWER LAW WITH AN EXPONENTIAL CUTOFF. Writing
%   alpha = 1-k and tau = theta, f(t) ~ t^(-alpha)*exp(-t/tau), and under
%   left truncation the two are the SAME distribution -- truncation deletes
%   (0,xmin), the only region their normalisations differ over. Verified
%   numerically to 3e-10. FITPOWERLAWCUTOFFMLE fits exactly this model and
%   reports (alpha, tau) instead; do not enter both in a model comparison
%   and treat them as independent evidence.
%
%   Not to be confused with the Weibull, where the exponential acts on a
%   POWER of t, f(t) ~ t^(k-1)*exp(-(t/lambda)^k), giving a
%   stretched-exponential rather than exponential tail.
%
%   All options are passed to FITTRUNCATEDDISCRETEMLE; see it for
%   SamplingInterval (required), DistributionType, the guards and the
%   output layout.
%
%   See also FITTRUNCATEDDISCRETEMLE, FITPOWERLAWCUTOFFMLE, FITERLANGMLE,
%   FITCHISQUAREDMLE, FITPEARSON3MLE.
H = FitTruncatedDiscreteMLE(eventseries, xmin, "gamma", varargin{:});
end
