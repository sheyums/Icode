function H = FitPearson3MLE(eventseries, xmin, varargin)
%FITPEARSON3MLE Left-truncated Pearson type III fit (discrete by default).
%
%   H = FITPEARSON3MLE(eventseries, xmin, SamplingInterval=dt, ...)
%
%   Parameters [shape k, scale theta, location]: a gamma shifted along the
%   time axis. Three parameters, and on left-truncated bout data it is
%   usually the strongest non-mixture model -- but the location is exactly
%   where it can go wrong.
%
%   THE LOCATION IS GUARDED. The likelihood diverges as the location
%   approaches the smallest retained duration, the same pathology
%   SHIFTLOGNORMAL_MLE documents for its shift. The engine parametrizes
%   location = xmin - exp(z) so it can never cross the threshold, reports
%   the gap as H.Diagnostics.LocationGap, and warns
%   (FitTruncatedDiscreteMLE:LocationAtBoundary) when the gap falls below
%   MinLocationGap. A warned fit is a stopping place on a divergence, not a
%   maximum: prefer the two-parameter gamma in that case.
%
%   See also FITTRUNCATEDDISCRETEMLE, FITGAMMAMLE, SHIFTLOGNORMAL_MLE.
H = FitTruncatedDiscreteMLE(eventseries, xmin, "pearson3", varargin{:});
end
