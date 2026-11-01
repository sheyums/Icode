function F = onsagerFit(J, EPR, options)
%ONSAGERFIT  Onsager conductance L from EPR = J^2/L through the origin (sqw_ref.onsager_fit; Methods Eq. 25).
%
%   F = onsagerFit(J, EPR)
%   F = onsagerFit(J, EPR, Level=0.95)
%
%   J, EPR  per fly x phase cycle current (cycles/h) and EPR (nats/h)
%
%   F  struct: s (slope = 1/L), se_s (OLS through the origin, df = n - 1),
%      L = 1/s, se_L = se_s/s^2, ci_L (t interval of s inverted; upper limit
%      Inf if the lower slope limit <= 0), n, r2 = 1 - RSS/sum(EPR.^2)
%      (R^2 through the origin).
%
%   Needs no toolbox: the t quantile comes from BETAINCINV (TQUANTILE).
%
%   See also ONSAGERBOOTSTRAP, EPRJMODELS.

arguments
    J (:,1) double
    EPR (:,1) double
    options.Level (1,1) double {mustBeInRange(options.Level, 0, 1, "exclusive")} = 0.95
end
if numel(J) ~= numel(EPR), error('onsagerFit:length', 'J and EPR must have the same length.'); end
if any(~isfinite(J)) || any(~isfinite(EPR))
    error('onsagerFit:finite', 'J and EPR must be finite: exclude rows before fitting.');
end
x = J .^ 2; n = numel(x);
if n < 2, error('onsagerFit:n', 'Need at least 2 points.'); end
s = (x' * EPR) / (x' * x);
res = EPR - s * x;
se = sqrt((res' * res) / (n - 1) / (x' * x));
q = tQuantile(1 - (1 - options.Level) / 2, n - 1);
lo = s - q * se; hi = s + q * se;
if lo > 0, up = 1 / lo; else, up = Inf; end
F = struct('s', s, 'se_s', se, 'L', 1 / s, 'se_L', se / s ^ 2, 'ci_L', [1 / hi, up], ...
    'n', n, 'r2', 1 - (res' * res) / (EPR' * EPR));
end
