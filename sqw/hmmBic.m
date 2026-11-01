function [b, k, nEff] = hmmBic(ll, K, nP, n, y)
%HMMBIC  BIC and parameter count of the phase-switched HMM (sqw_ref.bic, n_params, n_eff).
%
%   [b, k] = hmmBic(ll, K, nP, n)
%   [b, k, nEff] = hmmBic(ll, K, nP, [], y)      % n = number of symbol changes
%
%   k = (K - 1) + nP K (K - 1) + 2K   free parameters (3 symbols per state)
%   b = -2 ll + k ln n
%   Primary n = number of minutes; sensitivity n_eff = max(2, number of
%   symbol changes, missing minutes included, as sqw_ref.n_eff). Flag a fly
%   if K = 2 or K = 4 beats K = 3 by more than 10; routing needs K = 3.

arguments
    ll double
    K (1,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    n double = []
    y (:,1) double = []
end
k = (K - 1) + nP * K * (K - 1) + 2 * K;
nEff = NaN;
if ~isempty(y)
    nEff = max(2, nnz(diff(y) ~= 0));
    if isempty(n), n = nEff; end
end
b = -2 * ll + k * log(n);
end
