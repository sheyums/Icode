function p = wilcoxonP(x)
%WILCOXONP  Two-sided Wilcoxon signed-rank p of x against 0, as scipy.stats.wilcoxon's defaults.
%
%   p = wilcoxonP(x)
%
%   x   vector; NaN and zeros are dropped (scipy's zero_method 'wilcox')
%   p   two-sided p. As scipy's method 'auto': the EXACT null distribution when n <= 50 and the
%       |x| have no ties (SIGNRANK 'exact'); otherwise the normal approximation of scipy, with the
%       tie correction sum(t^3 - t) / 48 in the variance and NO continuity correction:
%           z = (R+ - n(n+1)/4) / sqrt(n(n+1)(2n+1)/24 - sum(t^3 - t)/48),  p = 2 P(Z > |z|)
%       (R+ = sum of the average ranks of |x| over x > 0). NaN with fewer than 2 values.
%   Values that are equal only up to round-off are NOT ties here: round them first (e.g.
%   round(J, 10) for the quantised cycle flux), as the 2026-11-04 convention does.
%   (Until 2026-11-04 this used MATLAB's exact method with mid-ranks when ties were present, which
%   differs from scipy in the second digit: sss dark 0.675 against scipy's 0.656.)
%
%   Example:
%       p = wilcoxonP(round(J(cycle), 10));
%
%   See also BINOMGREATER, HOLMADJUST, SIGNRANK, SPEARMANPERMTEST.

x = x(:); x = x(~isnan(x) & x ~= 0); n = numel(x);
if n < 2, p = NaN; return, end
a = abs(x); [s, i] = sort(a); r = zeros(n, 1); t = []; k = 1;
while k <= n
    j = k; while j < n && s(j + 1) == s(k), j = j + 1; end
    r(i(k:j)) = (k + j) / 2; t(end + 1) = j - k + 1; k = j + 1; %#ok<AGROW>
end
if n <= 50 && all(t == 1)
    p = signrank(x, 0, 'method', 'exact');
else
    se = sqrt(n * (n + 1) * (2 * n + 1) / 24 - sum(t .^ 3 - t) / 48);
    z = (sum(r(x > 0)) - n * (n + 1) / 4) / se; p = min(1, erfc(abs(z) / sqrt(2)));
end
end
