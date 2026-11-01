function [P, B, pi0, order] = labelStates(P, B, pi0)
%LABELSTATES  Order a K = 3 fit as (H1, H2, H3) (sqw_ref.label_states).
%
%   [P, B, pi0, order] = labelStates(P, B, pi0)
%
%   H3 = the state with the largest P(S | state); of the other two, H2 = the
%   larger P(W | state); H1 = the remaining one. The same permutation is
%   applied to every phase's matrix (P is 3 x 3 x nP). Ties go to the lower
%   index, as numpy.argmax. order(k) = original index of the state now in
%   position k.
%
%   See also ISDEGENERATEH1, HMMFITPHASED.

arguments
    P (3,3,:) double
    B (3,3) double
    pi0 (3,1) double
end
[~, h3] = max(B(:, 3));
rest = setdiff(1:3, h3);
[~, k] = max(B(rest, 2)); h2 = rest(k);
h1 = setdiff(rest, h2);
order = [h1 h2 h3];
P = P(order, order, :);
B = B(order, :);
pi0 = pi0(order);
end
