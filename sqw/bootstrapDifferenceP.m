function [ci, p] = bootstrapDifferenceP(bootA, bootB, options)
%BOOTSTRAPDIFFERENCEP  Interval and two-sided p for the difference of two bootstrap distributions (sqw_ref.bootstrap_difference_p).
%
%   [ci, p] = bootstrapDifferenceP(bootA, bootB)
%   [ci, p] = bootstrapDifferenceP(bootA, bootB, Level=0.95)
%
%   d = bootA - bootB, paired by replicate (two INDEPENDENT bootstraps of
%   equal length). ci = type-7 percentiles of d (as numpy.percentile);
%   p = min(1, 2 min(mean(d <= 0), mean(d >= 0))).

arguments
    bootA (:,1) double
    bootB (:,1) double
    options.Level (1,1) double {mustBeInRange(options.Level, 0, 1, "exclusive")} = 0.95
end
if numel(bootA) ~= numel(bootB)
    error('bootstrapDifferenceP:length', 'The two bootstraps must have the same length.');
end
d = bootA - bootB;
a = (1 - options.Level) / 2;
ci = quantile7(d, [a, 1 - a]);
p = min(1, 2 * min(mean(d <= 0), mean(d >= 0)));
end
