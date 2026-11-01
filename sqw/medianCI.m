function [m, lo, hi] = medianCI(v, options)
%MEDIANCI  Median over flies with a fly-bootstrap percentile interval (sqw17_report.ci).
%
%   [m, lo, hi] = medianCI(v)
%   [m, lo, hi] = medianCI(v, NBoot=2000, Level=0.95, Stream=rs, MinN=3)
%
%   v      one value per fly; NaN and +-Inf are DROPPED first (as in the
%          Python report). Use MEDIANCI on values with Inf only if you mean it:
%          for EPR_chain, dropping Inf removes the most dissipative flies.
%   m      median of the kept values; lo, hi type-7 percentiles of NBoot
%          bootstrap medians (flies drawn with replacement). All NaN if
%          fewer than MinN values.
%   Stream a RandStream (default: a fresh twister seeded 20). Pass one stream
%          through a whole report to get one reproducible sequence.

arguments
    v (:,1) double
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.Level (1,1) double {mustBeInRange(options.Level, 0, 1, "exclusive")} = 0.95
    options.Stream = []
    options.MinN (1,1) double {mustBePositive, mustBeInteger} = 3
end
v = v(isfinite(v));
if numel(v) < options.MinN
    m = NaN; lo = NaN; hi = NaN; return
end
rs = options.Stream;
if isempty(rs), rs = RandStream('twister', 'Seed', 20); end
n = numel(v);
b = median(v(randi(rs, n, n, options.NBoot)), 1)';
a = (1 - options.Level) / 2;
q = quantile7(b, [a, 1 - a]);
m = median(v); lo = q(1); hi = q(2);
end
