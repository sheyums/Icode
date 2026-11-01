function q = quantile7(x, p)
%QUANTILE7  Sample quantiles by linear interpolation (Hyndman-Fan type 7).
%
%   q = quantile7(x, p)      p in [0, 1] (NOT percent)
%
%   The definition of numpy.percentile's default, used for every interval
%   in sqw_ref. MATLAB's prctile / quantile use a different rule (type 5)
%   and give slightly different values on the same sample. NaNs in x are an
%   error: remove them first, deliberately.

arguments
    x (:,1) double
    p double {mustBeInRange(p, 0, 1)}
end
if any(isnan(x)), error('quantile7:nan', 'x holds NaN.'); end
if isempty(x), q = NaN(size(p)); return, end
x = sort(x); n = numel(x);
h = (n - 1) * p(:) + 1;                % columns throughout: x(lo) takes x's shape
lo = floor(h); hi = min(lo + 1, n);
q = x(lo) + (h - lo) .* (x(hi) - x(lo));
q = reshape(q, size(p));
end
