function d = cycleIncrements()
%CYCLEINCREMENTS  Increments of the counted cycle current on symbol transitions (sqw_ref D_CYCLE; PRESPEC 20q).
%
%   d = cycleIncrements()
%
%   d  (3 x 3), symbols 1 = Q, 2 = W, 3 = S: d(from, to) = +1 for the forward
%      cycle W -> Q, Q -> S, S -> W, -1 for the reverses, 0 on the diagonal.
%
%   See also TUREXACT, TURESTIMATE, HMMOBSERVEDCHAIN.

d = zeros(3);
fw = [2 1; 1 3; 3 2];                       % W->Q, Q->S, S->W
for r = 1:3
    d(fw(r, 1), fw(r, 2)) = 1; d(fw(r, 2), fw(r, 1)) = -1;
end
end
