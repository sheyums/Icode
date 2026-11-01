function [J, D, B] = turExactCT(T, d, options)
%TUREXACTCT  TUREXACT for the continuous-time generator T - I (sqw_ref.tur_exact_ct; PRESPEC 20q amendment C).
%
%   [J, D, B] = turExactCT(T)
%   [J, D, B] = turExactCT(T, d, H=1e-3)
%
%   Same inputs and outputs as TUREXACT, for the jump process with generator
%   T - I (rates per minute = the off-diagonal of T). The tilt acts on the
%   off-diagonal only: lambda(s) = largest real eigenvalue of
%   (T - I) with off-diagonal entries x exp(s d), J = lambda'(0),
%   D = lambda''(0) / 2 by central differences, x 60 per hour.
%   Here the steady-state TUR B <= sigma holds exactly; kappa = B_discrete / B_ct
%   is the 1-minute chain's inflation (<= 1.01 on the fitted fly chains).
%
%   sqw_ref writes D as (lambda(h) + lambda(-h)) / h^2 / 2, using lambda(0) = 0
%   for a generator. This port subtracts 2 lambda(0) explicitly; for a
%   stochastic T (checked) the two agree to round-off.
%
%   See also TUREXACT.

arguments
    T (:,:) double {mustBeNonnegative}
    d (:,:) double = cycleIncrements()
    options.H (1,1) double {mustBePositive} = 1e-3
end
if size(T, 1) ~= size(T, 2) || ~isequal(size(d), size(T))
    error('turExact:size', 'T must be square and d the same size.');
end
if any(abs(sum(T, 2) - 1) > 1e-9)
    error('turExact:rows', 'Every row of T must sum to 1 (FILLEMPTYROWS for empty rows).');
end
h = options.H;
lp = lam(T, d, h); lm = lam(T, d, -h); l0 = lam(T, d, 0);
J = (lp - lm) / (2 * h) * 60;
D = (lp - 2 * l0 + lm) / h^2 / 2 * 60;
if D > 0, B = J^2 / D; else, B = NaN; end
end

function l = lam(T, d, s)
n = size(T, 1);
M = T .* exp(s * d);
M(1:n + 1:end) = diag(T);                   % the diagonal is not tilted
l = max(real(eig(M))) - 1;
end
