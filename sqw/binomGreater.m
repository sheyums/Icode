function p = binomGreater(k, n, p0)
%BINOMGREATER  One-sided binomial test P(X >= k), X ~ Binomial(n, p0) (scipy binomtest, alternative 'greater').
%
%   p = binomGreater(k, n, p0)
%
%   Used for "more flies with p < 0.05 than chance": k of n flies, p0 = 0.05.

if n == 0, p = NaN; return, end
p = 1 - binocdf(k - 1, n, p0);
end
