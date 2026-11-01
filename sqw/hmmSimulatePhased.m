function [z, y] = hmmSimulatePhased(pi0, P, B, ph, rs)
%HMMSIMULATEPHASED  Hidden path and symbols from the phase-switched HMM (sqw_ref.simulate).
%
%   [z, y] = hmmSimulatePhased(pi0, P, B, ph)
%   [z, y] = hmmSimulatePhased(pi0, P, B, ph, rs)     % rs: a RandStream
%
%   z (T x 1) states 1..K, y (T x 1) symbols 1..3, T = numel(ph); the step
%   INTO minute t uses P(:,:,ph(t)). A row of P that sums to zero cannot be
%   simulated from: fill it first (FILLEMPTYROWS).

arguments
    pi0 (:,1) double
    P (:,:,:) double
    B (:,3) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    rs = RandStream('twister', 'Seed', 0)
end
T = numel(ph); K = numel(pi0);
if any(abs(sum(P, 2) - 1) > 1e-8, 'all')
    error('hmmSimulatePhased:rows', 'Every row of P must sum to 1 (FILLEMPTYROWS for empty rows).');
end
cP = cumsum(P, 2); cB = cumsum(B, 2); cpi = cumsum(pi0);
u = rand(rs, T, 1); v = rand(rs, T, 1);
z = zeros(T, 1);
z(1) = min(find(rand(rs) <= cpi, 1), K);
for t = 2:T
    k = find(u(t) <= cP(z(t - 1), :, ph(t)), 1);
    if isempty(k), k = K; end
    z(t) = k;
end
y = sum(v > cB(z, :), 2) + 1;
y = min(y, 3);
end
