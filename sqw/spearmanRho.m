function rho = spearmanRho(a, b)
%SPEARMANRHO  Spearman rank correlation with average ranks for ties (no toolbox).
%
%   rho = spearmanRho(a, b)     pairs with a non-finite value are dropped;
%                               NaN if fewer than 3 pairs remain.

arguments
    a (:,1) double
    b (:,1) double
end
k = isfinite(a) & isfinite(b);
if nnz(k) < 3, rho = NaN; return, end
ra = avgRank(a(k)); rb = avgRank(b(k));
ra = ra - mean(ra); rb = rb - mean(rb);
rho = (ra' * rb) / sqrt((ra' * ra) * (rb' * rb));
end

function r = avgRank(v)
[s, i] = sort(v);
n = numel(v); r = zeros(n, 1); k = 1;
while k <= n
    j = k; while j < n && s(j + 1) == s(k), j = j + 1; end
    r(i(k:j)) = (k + j) / 2;
    k = j + 1;
end
end
