function p = binom_upper(k, n, q)
% P(X >= k) for X ~ Binomial(n, q), toolbox-free.
if n == 0, p = NaN; return, end
j = k:n;
lp = gammaln(n + 1) - gammaln(j + 1) - gammaln(n - j + 1) + j*log(q) + (n - j)*log1p(-q);
p = min(1, sum(exp(lp)));
end
