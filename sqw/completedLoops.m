function L = completedLoops(y, ph)
%COMPLETEDLOOPS  Completed Q/W/S cycles on one fly's minute symbols, found by loop erasure (cycle_times.py loops()).
%
%   L = completedLoops(y, ph)
%
%   y    (T x 1) minute symbols in the port's coding (SQWSYMBOLS / SQWFLYINPUT): 1 = Q, 2 = W, 3 = S,
%        0 = missing, 4 = wake of unknown type (treated as missing, as every negative code in Python)
%   ph   (T x 1) phase of each minute (MINUTEPHASES)
%
%   The run-compressed path is walked minute by minute. A missing minute resets the path. When a
%   state already on the path is entered again, the loop since its last entry is erased; if that
%   loop holds all three states it is a COMPLETED CYCLE: forward if it is a rotation of
%   W -> Q -> S -> W, backward otherwise. After the erasure the returning state's entry minute is
%   reset to the current minute.
%
%   L  table, one row per completed cycle: t (the completion minute, 1-based), phase (ph(t), the
%      phase it is assigned to), dir (+1 forward, -1 backward), dur (minutes from entering the
%      loop's first state to returning to it).
%
%   Example (from the Icode\sqw folder):
%       [y, ph] = sqwFlyInput("D:\flies\CS\fly1.txt", true);
%       L = completedLoops(y, ph);
%       netPerHour = sum(L.dir(L.phase == 2)) / (nnz(ph == 2 & y >= 1 & y <= 3) / 60)   % = chain J (dark)
%
%   See also CYCLEWALK, COMPLETEDLOOPRATES, JCYCLECOUNT.

arguments
    y (:,1) double
    ph (:,1) double
end
st = zeros(3, 1); en = zeros(3, 1); nS = 0; prev = -1;      % path as a stack of (state, entry minute)
t0 = zeros(0, 1); p0 = t0; d0 = t0; u0 = t0;
for t = 1:numel(y)
    x = y(t);
    if x < 1 || x > 3, nS = 0; prev = -1; continue, end
    if x == prev, continue, end
    prev = x;
    k = find(st(1:nS) == x, 1);
    if isempty(k), nS = nS + 1; st(nS) = x; en(nS) = t; continue, end
    if nS - k + 1 == 3
        lp = st(k:nS)';
        t0(end + 1, 1) = t; p0(end + 1, 1) = ph(t); u0(end + 1, 1) = t - en(k); %#ok<AGROW>
        d0(end + 1, 1) = 2 * isForward(lp) - 1; %#ok<AGROW>
    end
    nS = k; en(k) = t;
end
L = table(t0, p0, d0, u0, 'VariableNames', {'t', 'phase', 'dir', 'dur'});
end

function tf = isForward(lp)
% rotations of (W, Q, S) = (2, 1, 3) in the port's coding
tf = isequal(lp, [2 1 3]) || isequal(lp, [1 3 2]) || isequal(lp, [3 2 1]);
end
