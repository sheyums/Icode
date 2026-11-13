[X, ds] = qload(); x = X{1}; rs = RandStream('mrg32k3a','Seed',3);
idx = (1:numel(x.d)).'; for q = 1:numel(x.W2), m = x.W2{q}; idx(m) = m(randperm(rs, numel(m))); end
s = x.s; rk0 = NaN(size(x.d)); rk0(s) = tiedrank(x.d(s)); rk0(~s) = tiedrank(x.d(~s));
a = qstats(x.d(idx), x); b = qstats([], x, rk0(idx)); fprintf('max abs diff %.2e\n', max(abs(a - b)));
tic; R = qrun(X, 20, 1); fprintf('113 flies x 20 perms %.1f s\n', toc);
