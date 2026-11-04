function [rho, p] = spearmanPermTest(x, y, nPerm, rs, options)
%SPEARMANPERMTEST  Spearman rank correlation with a two-sided permutation p (sqw22.py perm()).
%
%   [rho, p] = spearmanPermTest(x, y)
%   [rho, p] = spearmanPermTest(x, y, nPerm, rs)
%   [rho, p] = spearmanPermTest(x, y, nPerm, rs, RoundDigits=Inf)   % no rounding
%
%   x, y    vectors of the same length (rows with a non-finite value must be removed by the caller)
%   nPerm   number of permutations (default 10000)
%   rs      RandStream for the permutations (default twister, seed 22)
%   RoundDigits  x and y are rounded to this many decimals before ranking (default 10), so values
%           equal up to round-off (quantised J and c_obs: 1/24, 1/12 ... cycles/h) rank as ties, as
%           sqw22.py since 2026-11-04; Inf = no rounding
%   rho     Spearman rho: Pearson correlation of average ranks (ties get their mean rank, as
%           scipy rankdata / MATLAB tiedrank)
%   p       (1 + #{|rho_null| >= |rho| - 1e-12}) / (nPerm + 1), rho_null from permuting y's ranks
%           against x's. Monte Carlo: equal to the Python's only within permutation error.
%
%   Example:
%       [rho, p] = spearmanPermTest(sleep, J, 10000, RandStream('twister', 'Seed', 22));
%
%   See also SPEARMANRHO, STRATIFIEDSPEARMANTEST, SLEEPVSCYCLETESTS.

arguments
    x (:,1) double
    y (:,1) double
    nPerm (1,1) double {mustBeNonnegative, mustBeInteger} = 10000
    rs = RandStream('twister', 'Seed', 22)
    options.RoundDigits (1,1) double = 10
end
if numel(x) ~= numel(y), error('sqw:size', 'x and y must have the same length.'); end
if isfinite(options.RoundDigits), x = round(x, options.RoundDigits); y = round(y, options.RoundDigits); end
rx = avgRank(x); ry = avgRank(y); rho = pc(rx, ry);
if nPerm == 0 || ~isfinite(rho), p = NaN; return, end
n = numel(ry); [~, idx] = sort(rand(rs, n, nPerm), 1);
R = ry(idx); R = R - mean(R, 1); a = rx - mean(rx);
null = (a' * R) ./ sqrt((a' * a) * sum(R .^ 2, 1));
p = (1 + sum(abs(null) >= abs(rho) - 1e-12)) / (nPerm + 1);
end

function r = pc(a, b)
a = a - mean(a); b = b - mean(b); r = (a' * b) / sqrt((a' * a) * (b' * b));
end

function r = avgRank(v)
[s, i] = sort(v); n = numel(v); r = zeros(n, 1); k = 1;
while k <= n
    j = k; while j < n && s(j + 1) == s(k), j = j + 1; end
    r(i(k:j)) = (k + j) / 2; k = j + 1;
end
end
