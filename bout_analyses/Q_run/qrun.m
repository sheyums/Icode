function R = qrun(X, nR, seed)
nF = numel(X); obs = NaN(nF, 5); nm = NaN(nF, 5); g2 = zeros(nR, 1); cnt = zeros(nR, 1);
for f = 1:nF
    x = X{f}; s = x.s; rk0 = NaN(size(x.d)); rk0(s) = tiedrank(x.d(s)); rk0(~s) = tiedrank(x.d(~s));
    obs(f, :) = qstats(x.d, x, rk0);
    rs = RandStream('mrg32k3a', 'Seed', seed); rs.Substream = f; acc = NaN(nR, 5);
    for r = 1:nR
        idx = (1:numel(x.d)).';
        for q = 1:numel(x.W2), m = x.W2{q}; idx(m) = m(randperm(rs, numel(m))); end
        acc(r, :) = qstats([], x, rk0(idx));
    end
    nm(f, :) = mean(acc, 1, 'omitnan');
    if isfinite(obs(f, 2)), ok = isfinite(acc(:, 2)); g2(ok) = g2(ok) + acc(ok, 2); cnt(ok) = cnt(ok) + 1; end
end
R.obs = obs; R.nullmean = nm; R.fast = obs - nm;
R.F = mean(R.fast, 1, 'omitnan');
R.g2null = g2 ./ cnt; R.pT1 = tp(mean(obs(isfinite(obs(:, 2)), 2)), R.g2null);
end
