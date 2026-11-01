function p = wilcoxonP(x)
%WILCOXONP  Two-sided Wilcoxon signed-rank p of x against 0 (the role of scipy.stats.wilcoxon in the 20p/20q reports).
%
%   p = wilcoxonP(x)
%
%   NaN are dropped. Zeros are dropped (scipy's default zero_method 'wilcox').
%   Exact distribution for n <= 50, normal approximation above, as scipy's
%   'auto' (version >= 1.9); with ties MATLAB's exact method uses mid-ranks
%   where scipy switches to the normal approximation, so p can differ in the
%   second digit there. NaN with fewer than 2 values.
%
%   See also BINOMGREATER, HOLMADJUST, SIGNRANK.

x = x(:); x = x(~isnan(x) & x ~= 0);
if numel(x) < 2, p = NaN; return, end
if numel(x) <= 50, p = signrank(x, 0, 'method', 'exact');
else, p = signrank(x, 0, 'method', 'approximate'); end
end
