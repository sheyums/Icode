function [kap, ci, kb] = qkappa(fast, ds, nB, seed)
kap = mean(fast(:, 2), 'omitnan') / mean(fast(:, 1), 'omitnan');
rs = RandStream('mrg32k3a', 'Seed', seed); kb = NaN(nB, 1); u = unique(ds);
for b = 1:nB
    ix = [];
    for k = u(:).', m = find(ds(:) == k); ix = [ix; m(randi(rs, numel(m), numel(m), 1))]; end %#ok<AGROW>
    kb(b) = mean(fast(ix, 2), 'omitnan') / mean(fast(ix, 1), 'omitnan');
end
ci = quantile(kb, [.025 .975]);
end
