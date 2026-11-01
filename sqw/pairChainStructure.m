function [B2, Allow] = pairChainStructure()
%PAIRCHAINSTRUCTURE  Fixed emissions and structural zeros of the second-order chain as 9 pair states (sqw_ref B_PAIR, ALLOW_PAIR).
%
%   [B2, Allow] = pairChainStructure()
%
%   Pair state k = 3 (y_{t-1} - 1) + y_t, k = 1..9, symbols 1 = Q, 2 = W, 3 = S
%   (sqw_ref's 0-based 3 y_{t-1} + y_t, plus one).
%   B2     (9 x 3) state k emits its second symbol y_t with probability 1
%   Allow  (9 x 9) logical: (a,b) -> (b,c) only, i.e. the second symbol of k
%          equals the first symbol of l
%
%   See also CHAIN2FITPHASED, FILLEMPTYROWS.

k = (1:9)';
second = mod(k - 1, 3) + 1;                 % y_t of pair k
first = floor((k - 1) / 3) + 1;             % y_{t-1} of pair k
B2 = zeros(9, 3);
B2(sub2ind([9 3], k, second)) = 1;
Allow = second == first';                   % Allow(k,l): second(k) == first(l)
end
