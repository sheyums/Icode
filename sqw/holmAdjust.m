function adj = holmAdjust(p)
%HOLMADJUST  Holm step-down adjusted p-values (the 20p/20q reports' holm()).
%
%   adj = holmAdjust(p)
%
%   Sorted ascending, adj_(k) = max over j <= k of (m - j + 1) p_(j), capped at 1,
%   returned in the input order. The family is whatever p holds: the reports
%   take all pre-specified blocks (20p: 10 blocks; 20q: the 2 primary DD blocks).

p = p(:); m = numel(p); [ps, o] = sort(p);
adj = zeros(m, 1); run = 0;
for k = 1:m
    run = max(run, (m - k + 1) * ps(k)); adj(o(k)) = min(1, run);
end
end
