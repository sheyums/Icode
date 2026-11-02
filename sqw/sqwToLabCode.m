function s = sqwToLabCode(y)
%SQWTOLABCODE  Port symbols to the lab's coding: S = 0, Q = 1, W = 2, missing = NaN.
%
%   s = sqwToLabCode(y)
%
%   y   symbols as the port uses them internally (SQWSYMBOLS): 1 = Q, 2 = W, 3 = S,
%       0 = missing, 4 = wake of unknown type
%   s   same size: 0 = S, 1 = Q, 2 = W, NaN = missing, -1 = wake of unknown type
%
%   Sleep stays 0 and wake stays > 0, as in SLEEPTRACE (0 = sleep, 1 = wake). Every
%   file this folder writes uses this coding; the internal 1..3 is what MATLAB needs to
%   index matrices by symbol. LABCODETOSQW is the inverse.
%
%   See also LABCODETOSQW, SQWSYMBOLS, SQWBATCHSYMBOLS.

if ~all(ismember(y(:), 0:4)), error('sqw:symbols', 'Port symbols must be 0..4.'); end
map = [NaN 1 2 0 -1];                       % y = 0 1 2 3 4
s = reshape(map(y + 1), size(y));
end
