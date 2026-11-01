function Lb = onsagerBootstrap(fly, J, EPR, options)
%ONSAGERBOOTSTRAP  Fly bootstrap of the Onsager conductance L (sqw_ref.onsager_bootstrap).
%
%   Lb = onsagerBootstrap(fly, J, EPR)
%   Lb = onsagerBootstrap(fly, J, EPR, NBoot=2000, RandomSeed=20)
%
%   fly      (n x 1) fly identifier of each point (numeric or string); a fly
%            usually has two points (light, dark)
%   J, EPR   the points, as ONSAGERFIT
%   Lb       (NBoot x 1) bootstrap L values: flies drawn with replacement,
%            each with ALL its points
%
%   The draws use a private RandStream, so the global generator is left
%   alone. Seeds do not carry across MATLAB and numpy: compare with Python
%   statistically (median within Monte Carlo error), never draw for draw.
%
%   See also ONSAGERFIT, BOOTSTRAPDIFFERENCEP.

arguments
    fly (:,1)
    J (:,1) double
    EPR (:,1) double
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
end
if numel(fly) ~= numel(J) || numel(J) ~= numel(EPR)
    error('onsagerBootstrap:length', 'fly, J and EPR must have the same length.');
end
[~, ~, g] = unique(fly);
nF = max(g);
rows = accumarray(g, (1:numel(g))', [nF 1], @(v) {v});
rs = RandStream('twister', 'Seed', options.RandomSeed);
Lb = zeros(options.NBoot, 1);
for b = 1:options.NBoot
    pick = randi(rs, nF, nF, 1);
    r = vertcat(rows{pick});
    x = J(r) .^ 2; e = EPR(r);
    Lb(b) = (x' * x) / (x' * e);          % L = 1/s, s = x'e / x'x
end
end
