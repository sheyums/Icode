function M = eprJModels(J, EPR, options)
%EPRJMODELS  EPR against J across flies: constant affinity vs linear response (sqw_ref.epr_j_models; PRESPEC 20j).
%
%   M = eprJModels(J, EPR)
%   M = eprJModels(J, EPR, Offset=0.02, JMin=0.05, EPRMin=0.05, MinExponentPoints=5)
%
%   M1 (constant affinity)   EPR = b + a |J|
%   M2 (linear response)     EPR = b + c J^2       (c = 1/L)
%   Both by least squares with a free intercept b (it absorbs the estimator
%   offset; it is not physics), both 2 parameters, compared by residual sum
%   of squares. Exponent: slope of ln(EPR - Offset) on ln|J| over points with
%   |J| >= JMin and EPR >= EPRMin; NaN if fewer than MinExponentPoints.
%
%   M  struct: a, b_linear, rss_linear, c, b_quadratic, rss_quadratic,
%      exponent, n_points, n_exponent
%
%   See also ONSAGERFIT.

arguments
    J (:,1) double
    EPR (:,1) double
    options.Offset (1,1) double = 0.02
    options.JMin (1,1) double {mustBeNonnegative} = 0.05
    options.EPRMin (1,1) double = 0.05
    options.MinExponentPoints (1,1) double {mustBePositive, mustBeInteger} = 5
end
if numel(J) ~= numel(EPR), error('eprJModels:length', 'J and EPR must have the same length.'); end
[a, b1, r1] = ls2(abs(J), EPR);
[c, b2, r2] = ls2(J .^ 2, EPR);
m = abs(J) >= options.JMin & EPR >= options.EPRMin;
if nnz(m) >= options.MinExponentPoints
    k = polyfit(log(abs(J(m))), log(EPR(m) - options.Offset), 1);
    ex = k(1);
else
    ex = NaN;
end
M = struct('a', a, 'b_linear', b1, 'rss_linear', r1, 'c', c, 'b_quadratic', b2, ...
    'rss_quadratic', r2, 'exponent', ex, 'n_points', numel(J), 'n_exponent', nnz(m));
end

function [slope, icpt, rss] = ls2(x, e)
A = [ones(size(x)), x];
coef = A \ e;
icpt = coef(1); slope = coef(2);
rss = sum((A * coef - e) .^ 2);
end
