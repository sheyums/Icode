function v = qstats(d, x, rk)
% v = [r1 r2 r3 rSS rWW]; rk optional precomputed within-type ranks of d
if nargin < 3
    s = x.s; rk = NaN(size(d)); rk(s) = tiedrank(d(s)); rk(~s) = tiedrank(d(~s));
end
v = NaN(1, 5); rr = NaN(3, 2);
for lag = 1:3
    for kd = 1:2
        i = x.P{lag}{kd};
        if numel(i) >= 5
            a = rk(i) - mean(rk(i)); b = rk(i + lag) - mean(rk(i + lag));
            rr(lag, kd) = (a.' * b) / sqrt((a.' * a) * (b.' * b));
        end
    end
    v(lag) = mean(rr(lag, :), 'omitnan');
end
v(4) = rr(2, 1); v(5) = rr(2, 2);
end
