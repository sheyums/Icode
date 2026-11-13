function y = qsim(x, sc, rs)
y = x; s = x.s; d = x.d; n = numel(d);
switch sc
    case 1, g = randn(rs, n, 1);
    case 2, z = ar1(0.8, n, rs); g = sqrt(0.075) * z + sqrt(0.925) * randn(rs, n, 1);
    case 3, g = ar1(0.06, n, rs);
end
for t = [true false]
    ix = find(s == t); v = d(ix); v = sort(v(randi(rs, numel(v), numel(v), 1)));
    [~, ord] = sort(g(ix)); d(ix(ord)) = v;
end
y.d = d;
end
function z = ar1(rho, n, rs)
e = randn(rs, n, 1); z = zeros(n, 1); z(1) = e(1);
for i = 2:n, z(i) = rho * z(i - 1) + sqrt(1 - rho^2) * e(i); end
end
