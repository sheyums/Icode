function S = blockKLStats(y, ph, nP)
%BLOCKKLSTATS  Block-irreversibility statistics per phase (sqw_ref.block_kl_stats; PRESPEC 20p).
%
%   S = blockKLStats(y, ph, nP)
%
%   S  (nP x 1) struct array, nats/h:
%        d2, d3, d4   d_n = (D_n - D_{n-1}) x 60 (BLOCKKL, all symbols)
%        D3 = d3 - d2 (PRIMARY), D4 = d4 - d3
%        D3QW         D3 on wake-only blocks (alphabet {Q, W})
%
%   No statistic here is interpretable on its own. The excess
%   E = D3 - median(D3 of records simulated from the fly's own first-order
%   chain, same mask, same estimator) is what PRESPEC 20p reports; raw d3/d2
%   ratios are estimator bias.
%
%   See also BLOCKKL, CHAINFITPHASED, FILLEMPTYROWS, HMMSIMULATEPHASED.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
end
D = blockKL(y, ph, nP); Dq = blockKL(y, ph, nP, QW=true);
d = diff(D, 1, 2) * 60; dq = diff(Dq, 1, 2) * 60;    % columns d2, d3, d4
S = struct('d2', num2cell(d(:, 1)), 'd3', num2cell(d(:, 2)), 'd4', num2cell(d(:, 3)), ...
    'D3', num2cell(d(:, 2) - d(:, 1)), 'D4', num2cell(d(:, 3) - d(:, 2)), 'D3QW', num2cell(dq(:, 2) - dq(:, 1)));
end
