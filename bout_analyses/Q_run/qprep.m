function x = qprep(fl, cutH)
o = fl.Onset(:); d = fl.Dur(:); s = logical(fl.IsSleep(:));
k = o >= cutH * 3600; o = o(k); d = d(k); s = s(k);
seg = floor(o / 43200); n = numel(o);
a = [o(2:end) == o(1:end-1) + d(1:end-1) & seg(2:end) == seg(1:end-1); false];
x.o = o; x.d = d; x.s = s;
for lag = 1:3
    ok = true(n, 1); for j = 0:lag-1, ok = ok & [a(1+j:end); false(j, 1)]; end
    i = find(ok); i = i(i + lag <= n);
    x.P{lag} = {i(s(i)), i(~s(i))};          % kinds: starting with sleep / with wake
end
w = floor(o / 7200); sets = {};
for t = [true false]
    ix = find(s == t); [~, ~, gi] = unique(w(ix));
    for u = 1:max(gi), m = ix(gi == u); if numel(m) > 1, sets{end+1} = m; end, end %#ok<AGROW>
end
x.W2 = sets;
end
