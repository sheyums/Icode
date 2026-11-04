function [z, p, rhoEquiv] = stratifiedSpearmanTest(X, Y, nPerm, rs, options)
%STRATIFIEDSPEARMANTEST  Fisher-z combination of within-block Spearman rhos, with a within-block permutation p (sqw22.py strat()).
%
%   [z, p, rhoEquiv] = stratifiedSpearmanTest(X, Y)
%   [z, p, rhoEquiv] = stratifiedSpearmanTest(X, Y, nPerm, rs)
%
%   X, Y      cell arrays, one cell per block: the x and y values of that block (finite only)
%   nPerm     number of permutations (default 2000)
%   rs        RandStream (default twister, seed 22)
%   RoundDigits  values rounded to this many decimals before ranking (default 10; Inf = none), as
%             SPEARMANPERMTEST
%   z         weighted mean of atanh(rho_b) over blocks, weights n_b - 3 (rho clipped to +-0.999999)
%   p         (1 + #{|z_null| >= |z|}) / (nPerm + 1), z_null with x permuted WITHIN each block
%   rhoEquiv  tanh(z)
%
%   Example:
%       [z, p] = stratifiedSpearmanTest({sleepCS, sleepFumin}, {cCS, cFumin}, 2000, RandStream('twister', 'Seed', 22));
%
%   See also SPEARMANPERMTEST, SLEEPVSCYCLETESTS.

arguments
    X cell
    Y cell
    nPerm (1,1) double {mustBeNonnegative, mustBeInteger} = 2000
    rs = RandStream('twister', 'Seed', 22)
    options.RoundDigits (1,1) double = 10
end
if isfinite(options.RoundDigits)
    X = cellfun(@(v) round(v, options.RoundDigits), X, 'UniformOutput', false); Y = cellfun(@(v) round(v, options.RoundDigits), Y, 'UniformOutput', false);
end
w = cellfun(@numel, X) - 3;
zf = @(XX) sum(w .* cellfun(@(x, y) atanh(min(max(spearmanRho(x, y), -0.999999), 0.999999)), XX, Y)) / sum(w);
z = zf(X); rhoEquiv = tanh(z); zn = zeros(nPerm, 1);
for k = 1:nPerm
    zn(k) = zf(cellfun(@(x) x(randperm(rs, numel(x))), X, 'UniformOutput', false));
end
p = (1 + sum(abs(zn) >= abs(z))) / (nPerm + 1);
end
