function [tf, tv] = isDegenerateH1(B, options)
%ISDEGENERATEH1  Is H1 unidentified? (sqw_ref.is_degenerate)
%
%   [tf, tv] = isDegenerateH1(B)
%   [tf, tv] = isDegenerateH1(B, Threshold=0.05)
%
%   B   (3 x 3) LABELLED emissions (rows H1, H2, H3)
%   tf  true if H1's emission row lies within total variation Threshold of
%       H2's or H3's; such flies are excluded from routing (report how many)
%   tv  [TV(H1,H2), TV(H1,H3)], TV = half the L1 distance

arguments
    B (3,3) double
    options.Threshold (1,1) double {mustBeNonnegative} = 0.05
end
tv = [0.5 * sum(abs(B(1, :) - B(2, :))), 0.5 * sum(abs(B(1, :) - B(3, :)))];
tf = min(tv) < options.Threshold;
end
