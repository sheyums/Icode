function y = labCodeToSqw(s)
%LABCODETOSQW  Lab-coded symbols (S = 0, Q = 1, W = 2, missing = NaN) to the port's internal 1..3.
%
%   y = labCodeToSqw(s)
%
%   s   0 = S, 1 = Q, 2 = W, NaN = missing, -1 = wake of unknown type
%   y   same size: 1 = Q, 2 = W, 3 = S, 0 = missing, 4 = wake of unknown type
%
%   Anything else (e.g. 3, 0.5, Inf) is an error rather than a guess.
%
%   See also SQWTOLABCODE.

s = double(s); y = zeros(size(s));
ok = isnan(s) | ismember(s, [-1 0 1 2]);
if ~all(ok(:)), error('sqw:labCode', 'Lab symbols must be 0 (S), 1 (Q), 2 (W), NaN (missing) or -1 (unknown wake).'); end
y(s == 1) = 1; y(s == 2) = 2; y(s == 0) = 3; y(s == -1) = 4;
end
