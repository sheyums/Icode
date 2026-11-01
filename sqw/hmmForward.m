function [A, c, ll] = hmmForward(pi0, P, B, y, ph)
%HMMFORWARD  Scaled forward pass of the phase-switched HMM (sqw_ref.forward, loglik).
%
%   [A, c, ll] = hmmForward(pi0, P, B, y, ph)
%
%   pi0  (K x 1) initial distribution
%   P    (K x K x nP) per-phase transition matrices; the step INTO minute t
%        uses P(:,:,ph(t))
%   B    (K x 3) emissions shared over phases
%   y    (T x 1) symbols (HMMEMISSION), ph (T x 1) phases 1..nP
%
%   A    (T x K) filtered probabilities alpha_hat(t,:) = P(z_t | y_1..t)
%   c    (T x 1) scale factors, ll = sum(log(c)) = log-likelihood

arguments
    pi0 (:,1) double
    P (:,:,:) double
    B (:,3) double
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
end
E = hmmEmission(B, y);
[T, K] = size(E);
if numel(ph) ~= T, error('sqw:length', 'y and ph must have the same length.'); end
A = zeros(T, K); c = zeros(T, 1);
a = pi0' .* E(1, :);
c(1) = sum(a); A(1, :) = a / c(1);
for t = 2:T
    a = (A(t - 1, :) * P(:, :, ph(t))) .* E(t, :);
    c(t) = sum(a); A(t, :) = a / c(t);
end
ll = sum(log(c));
end
