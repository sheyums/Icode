function q = tQuantile(p, nu)
%TQUANTILE  Quantile of Student's t distribution without the Statistics Toolbox.
%
%   q = tQuantile(p, nu)      same as tinv(p, nu)
%
%   Uses the regularized incomplete beta: for p > 1/2,
%   x = betaincinv(2(1 - p), nu/2, 1/2) and q = sqrt(nu (1/x - 1)); the lower
%   tail by symmetry.

arguments
    p double {mustBeInRange(p, 0, 1)}
    nu double {mustBePositive}
end
q = zeros(size(p + nu));
p = p + zeros(size(q)); nu = nu + zeros(size(q));
up = p > 0.5; lo = p < 0.5;
x = betaincinv(2 * (1 - p(up)), nu(up) / 2, 0.5);
q(up) = sqrt(nu(up) .* (1 ./ x - 1));
x = betaincinv(2 * p(lo), nu(lo) / 2, 0.5);
q(lo) = -sqrt(nu(lo) .* (1 ./ x - 1));
q(p == 0) = -Inf; q(p == 1) = Inf;
end
