function w = cycleWalk(y)
%CYCLEWALK  Cumulative net-cycle walk of one fly: +1 at each S -> W, -1 at each W -> S (fig_act1_walk.py).
%
%   w = cycleWalk(y)
%
%   y   (T x 1) minute symbols in the port's coding: 1 = Q, 2 = W, 3 = S, 0 = missing (4 = unknown wake
%       counts as missing)
%   w   (T x 1) w(1) = 0; w(t) = number of sleep -> active-wake steps (S -> W) minus active-wake ->
%       sleep steps (W -> S) up to minute t, counting only steps between two consecutive SCORED minutes.
%       Its slope is the cycle flux; it equals the net number of completed cycles (COMPLETEDLOOPS)
%       to within one, the loop in progress.
%
%   Example (from the Icode\sqw folder):
%       [y, ~] = sqwFlyInput("D:\flies\CS\fly1.txt", true);
%       plot((0:numel(y) - 1) / 60, cycleWalk(y)); xlabel('time (h)')
%
%   See also COMPLETEDLOOPS, EXPORTCUMULATIVEWALKFIGURE.

arguments
    y (:,1) double
end
a = y(1:end - 1); b = y(2:end);
st = double(a == 3 & b == 2) - double(a == 2 & b == 3);
w = [0; cumsum(st)];
end
