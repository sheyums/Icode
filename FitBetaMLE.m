function H = FitBetaMLE(eventseries, xmin, varargin)
%FITBETAMLE Left-truncated beta fit on [0, UpperBound], discrete by default.
%
%   H = FITBETAMLE(eventseries, xmin, SamplingInterval=dt, UpperBound=B, ...)
%
%   Parameters [p, q]. Beta has BOUNDED support, so UpperBound is required.
%
%   UpperBound is a DURATION IN THE DATA'S UNITS, exactly like xmin -- not
%   a number of samples. For a 24 h recording scored in seconds it is
%   86400; for the same data rescaled to minutes it is 1440. Dividing a
%   duration by SamplingInterval converts it into a BIN COUNT, which is a
%   different quantity and coincides only when SamplingInterval happens to
%   be 1.
%
%   It must be FIXED, not fitted. Estimating the upper bound makes the
%   likelihood unbounded as it approaches max(eventseries), the same
%   pathology SHIFTLOGNORMAL_MLE documents for its shift. Use a bound the
%   process physically cannot exceed, such as the recording length.
%
%   Be aware of how much of the support the data actually occupy. q governs
%   the density near UpperBound, so if the longest bout reaches only a
%   fraction of the way there, q is extrapolation into a region holding no
%   observations. H.Diagnostics.SupportCoverage reports max(data)/UpperBound.
%
%   See also FITTRUNCATEDDISCRETEMLE.
H = FitTruncatedDiscreteMLE(eventseries, xmin, "beta", varargin{:});
end
