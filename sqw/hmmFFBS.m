function Z = hmmFFBS(pi0, P, B, y, ph, options)
%HMMFFBS  Posterior hidden paths by forward filtering, backward sampling (sqw_ref.ffbs).
%
%   Z = hmmFFBS(pi0, P, B, y, ph)
%   Z = hmmFFBS(pi0, P, B, y, ph, NPaths=200, RandomSeed=0)
%
%   Inputs as HMMFORWARD. Z (NPaths x T), states 1..K, each row one path
%   drawn from P(z | y). Forward pass first; the last state is drawn from the
%   final filtered distribution, then backward
%   p(z_t = i | z_{t+1} = j, y_1..t) ~ alpha_t(i) P_{ph(t+1)}(i, j).
%   Draws come from a private RandStream; numpy's cannot be reproduced, so
%   compare with Python statistically (G6), never path for path.
%
%   See also POSTERIORSUMMARIES, HMMVITERBIPHASED.

arguments
    pi0 (:,1) double
    P (:,:,:) double
    B (:,3) double
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    options.NPaths (1,1) double {mustBePositive, mustBeInteger} = 200
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 0
end
A = hmmForward(pi0, P, B, y, ph);
[T, K] = size(A);
n = options.NPaths;
rs = RandStream('twister', 'Seed', options.RandomSeed);
Z = zeros(n, T);
cdf = cumsum(A(T, :)); cdf = cdf / cdf(end);
Z(:, T) = sum(rand(rs, n, 1) > cdf, 2) + 1;
for t = T - 1:-1:1
    w = A(t, :)' .* P(:, Z(:, t + 1), ph(t + 1));      % K x n
    cw = cumsum(w, 1); cw = cw ./ cw(end, :);
    Z(:, t) = sum(rand(rs, 1, n) > cw, 1)' + 1;
end
Z = min(Z, K);                                          % guard against cdf(end) < u by rounding
end
