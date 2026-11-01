function z = hmmViterbiPhased(pi0, P, B, y, ph)
%HMMVITERBIPHASED  Most probable hidden path of the phase-switched HMM (sqw_ref.viterbi).
%
%   z = hmmViterbiPhased(pi0, P, B, y, ph)
%
%   Inputs as HMMFORWARD. z (T x 1) in 1..K. Log domain; ties go to the
%   lowest state index (numpy.argmax). Kept for comparison only: the
%   analysis decodes with HMMFFBS, because Viterbi is biased when emissions
%   are degenerate.

arguments
    pi0 (:,1) double
    P (:,:,:) double
    B (:,3) double
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
end
lE = log(hmmEmission(B, y));
[T, K] = size(lE);
lP = log(P);
d = log(pi0)' + lE(1, :);
bp = zeros(T, K);
for t = 2:T
    m = d' + lP(:, :, ph(t));            % m(i,j): from i to j
    [mx, bp(t, :)] = max(m, [], 1);      % first maximum = lowest index
    d = mx + lE(t, :);
end
z = zeros(T, 1);
[~, z(T)] = max(d);
for t = T:-1:2
    z(t - 1) = bp(t, z(t));
end
end
