function [pi1, P1, B1, ll] = hmmEmStep(pi0, P, B, y, ph)
%HMMEMSTEP  One Baum-Welch update of the phase-switched HMM (sqw_ref.em_step).
%
%   [pi1, P1, B1, ll] = hmmEmStep(pi0, P, B, y, ph)
%
%   Inputs as HMMFORWARD. Returns the updated parameters and the
%   log-likelihood of the INPUT parameters (deterministic).
%
%   Transition expected counts are split by the phase of the ARRIVAL minute;
%   emissions are pooled over phases. An unknown-wake minute (4) has its
%   posterior weight gamma(k) split between Q and W in proportion to B(k,Q)
%   and B(k,W), so it counts toward wake but toward neither type. A row with
%   no expected transitions (a phase never entered, or a state never left
%   within it) comes back as zeros, exactly as in sqw_ref; callers that
%   simulate must fill it (FILLEMPTYROWS).

arguments
    pi0 (:,1) double
    P (:,:,:) double
    B (:,3) double
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
end
nP = size(P, 3);
E = hmmEmission(B, y);
[A, c, ll] = hmmForward(pi0, P, B, y, ph);
[T, K] = size(E);
Bk = ones(T, K);
for t = T - 1:-1:1
    Bk(t, :) = ((P(:, :, ph(t + 1)) * (E(t + 1, :) .* Bk(t + 1, :))') / c(t + 1))';
end
G = A .* Bk; G = G ./ sum(G, 2);
R = E .* Bk ./ c;
P1 = zeros(size(P));
for p = 1:nP
    idx = find(ph(2:end) == p) + 1;
    Xi = (A(idx - 1, :)' * R(idx, :)) .* P(:, :, p);
    P1(:, :, p) = Xi ./ max(sum(Xi, 2), 1e-300);
end
onehot = double([y == 1, y == 2, y == 3]);
Bn = G' * onehot;
u = y == 4;
if any(u)
    g = sum(G(u, :), 1)';
    qw = B(:, 1) + B(:, 2); qw(qw <= 0) = 1;
    Bn(:, 1) = Bn(:, 1) + g .* B(:, 1) ./ qw;
    Bn(:, 2) = Bn(:, 2) + g .* B(:, 2) ./ qw;
end
B1 = Bn ./ max(sum(Bn, 2), 1e-300);
pi1 = G(1, :)';
end
