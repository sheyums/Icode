function [J, D, B] = turExact(T, d, options)
%TUREXACT  Exact current, diffusion and TUR bound of a per-minute Markov chain (sqw_ref.tur_exact; PRESPEC 20q).
%
%   [J, D, B] = turExact(T)
%   [J, D, B] = turExact(T, d, H=1e-3)
%
%   T  (n x n) per-minute transition matrix (rows sum to 1), at its stationary
%      state; e.g. a phase of CHAINFITPHASED's P after FILLEMPTYROWS, or the
%      9-state chain of HMMOBSERVEDCHAIN
%   d  (n x n) increments of the current on each transition (default
%      CYCLEINCREMENTS, the forward S/Q/W cycle)
%
%   J  mean current per hour, D diffusion per hour, B = J^2 / D (k_B per hour;
%      NaN unless D > 0), from the tilted matrix T .* exp(s d):
%      lambda(s) = log of its largest real eigenvalue (Perron root),
%      J = lambda'(0), D = lambda''(0) / 2, by central differences with step H,
%      x 60 per hour.
%
%   The steady-state TUR (B <= sigma) is proven in continuous time; a 1-minute
%   discrete chain can exceed it slightly, so 20q deflates by
%   kappa = B / B_ct (TUREXACTCT).
%
%   See also TUREXACTCT, TURESTIMATE, MODELENERGETICS.

arguments
    T (:,:) double {mustBeNonnegative}
    d (:,:) double = cycleIncrements()
    options.H (1,1) double {mustBePositive} = 1e-3
end
checkChain(T, d);
h = options.H;
lp = logPerron(T .* exp(h * d)); lm = logPerron(T .* exp(-h * d)); l0 = logPerron(T);
J = (lp - lm) / (2 * h) * 60;
D = (lp - 2 * l0 + lm) / h^2 / 2 * 60;
if D > 0, B = J^2 / D; else, B = NaN; end
end

function l = logPerron(M)
l = log(max(real(eig(M))));
end

function checkChain(T, d)
if size(T, 1) ~= size(T, 2) || ~isequal(size(d), size(T))
    error('turExact:size', 'T must be square and d the same size.');
end
if any(abs(sum(T, 2) - 1) > 1e-9)
    error('turExact:rows', 'Every row of T must sum to 1 (FILLEMPTYROWS for empty rows).');
end
end
