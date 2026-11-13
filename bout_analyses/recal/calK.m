L = load('C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\homeo_run\homeo_results.mat', 'flies');
P = L.flies.per0DD; nF = numel(P);
on = cell(1, nF); d = on; T = on; S12 = on;
for f = 1:nF
    on{f} = P(f).Onset(P(f).IsSleep); d{f} = P(f).Dur(P(f).IsSleep);
    T{f} = on{f} / 86400; S12{f} = floor(on{f} / 43200);
end
% ---------------- calibration ----------------
sc = ["a stationary","b common trend +10%/d","c fly trends N(0,.1)","d first-day x1.5","e 12-h state SD .3"];
rej = zeros(5, 2);
for k = 1:5
    for g = 1:100
        rs = RandStream('mrg32k3a', 'Seed', 1000 * k + g);
        sl = 0.1 * randn(rs, nF, 1); Y = cell(1, nF);
        for f = 1:nF
            mu = mean(d{f}) - 300; t = T{f}; tm = mean(t);
            switch k
                case 1, m = mu * ones(size(t));
                case 2, m = mu * exp(log(1.1) * (t - tm));
                case 3, m = mu * exp(sl(f) * (t - tm));
                case 4, m = mu * (1 + 0.5 * (t < 1));
                case 5, z = randn(rs, max(S12{f}) + 1, 1); m = mu * exp(0.3 * z(S12{f} + 1) - 0.045);
            end
            Y{f} = log(300 + m .* -log(rand(rs, size(t))));
        end
        o = fluct(Y, T, S12, 199, 40000 + 100 * k + g, true); tr = trend(Y, T, 2000, 40000 + 100 * k + g);
        rej(k, :) = rej(k, :) + [o.p < 0.05, tr.p < 0.05];
    end
    fprintf('CAL %-24s fluctuation rejection %.2f | trend rejection %.2f\n', sc(k), rej(k, 1) / 100, rej(k, 2) / 100);
end
okCal = all(rej(1:4, 1) / 100 <= 0.095);
fprintf('CALIBRATION fluctuation test size ok in a-d: %d\n', okCal);
save('calK.mat','rej','okCal');
% ================= functions =================
function o = trend(Y, T, nR, seed)
nF = numel(Y); b = NaN(nF, 1);
for f = 1:nF, X = [ones(numel(T{f}), 1) T{f}(:)]; c = X \ Y{f}(:); b(f) = c(2); end
rs = RandStream('mrg32k3a', 'Seed', seed); s = 2 * (rand(rs, nF, nR) > 0.5) - 1;
n = mean(b .* s, 1); o.b = b; o.mean = mean(b); o.p = (1 + sum(abs(n) >= abs(o.mean))) / (nR + 1);
end

function o = fluct(Y, T, S, nR, seed, detr)
nF = numel(Y); nS = max(cellfun(@max, S)) + 1;
Q = cell(1, nF); si = cell(1, nF); cnt = cell(1, nF);
for f = 1:nF
    [Q{f}, ~] = qr([ones(numel(T{f}), 1) T{f}(:)], 0); si{f} = S{f}(:) + 1;
    cnt{f} = accumarray(si{f}, 1, [nS 1]);
end
[Vo, Mo, res] = vset(Y, Q, si, cnt, nS, detr);
Vn = NaN(nF, nR); l1n = NaN(1, nR);
rs = cell(1, nF);
for f = 1:nF, rs{f} = RandStream('mrg32k3a', 'Seed', seed); rs{f}.Substream = f; end
for r = 1:nR
    Yp = cell(1, nF);
    for f = 1:nF, y = Y{f}(:); Yp{f} = y(randperm(rs{f}, numel(y))); end
    [Vn(:, r), Mn] = vset(Yp, Q, si, cnt, nS, detr); l1n(r) = lag1(Mn);
end
mV = mean(Vn, 2, 'omitnan');
o.lr = log(Vo ./ mV); o.lrn = log(Vn ./ mV); o.ratio = Vo ./ mV;
o.G = mean(o.lr, 'omitnan'); o.Gn = mean(o.lrn, 1, 'omitnan');
o.p = (1 + sum(o.Gn >= o.G)) / (nR + 1); o.lag1 = lag1(Mo); o.lag1n = l1n; o.res = res; o.M = Mo;
end

function [V, M, res] = vset(Y, Q, si, cnt, nS, detr)
nF = numel(Y); M = NaN(nF, nS); R = cell(1, nF);
for f = 1:nF
    y = Y{f}(:);
    if detr, r = y - Q{f} * (Q{f}' * y); else, r = y - mean(y); end
    R{f} = r; s = accumarray(si{f}, r, [nS 1]) ./ cnt{f}; s(cnt{f} < 3) = NaN; M(f, :) = s.';
end
if detr, c = mean(M, 1, 'omitnan'); c(~isfinite(c)) = 0; M = M - c; else, c = zeros(1, nS); end
V = NaN(nF, 1); res = cell(1, nF);
for f = 1:nF
    m = M(f, isfinite(M(f, :))); if numel(m) >= 3, V(f) = var(m); end
    res{f} = R{f} - c(si{f}).';
end
end

function r = lag1(M)
a = M(:, 1:end-1); b = M(:, 2:end); k = isfinite(a) & isfinite(b);
r = NaN(size(M, 1), 1);
for f = 1:size(M, 1)
    if sum(k(f, :)) >= 4, cc = corrcoef(a(f, k(f, :)), b(f, k(f, :))); r(f) = cc(2); end
end
r = mean(r, 'omitnan');
end

function [c0, cn] = clus(R, S, nR, seed)
nF = numel(R); co = NaN(nF, 1); cnn = NaN(nF, nR);
for f = 1:nF
    rs = RandStream('mrg32k3a', 'Seed', seed); rs.Substream = f; x = R{f}(:); s = S{f}(:);
    co(f) = cs(x, s);
    for r = 1:nR, cnn(f, r) = cs(x(randperm(rs, numel(x))), s); end
end
c0 = mean(co, 'omitnan'); cn = mean(cnn, 1, 'omitnan');
end
function c = cs(d, seg)
c = NaN; n = numel(d); if n < 10, return, end
sd = sort(d); pr = find(seg(2:end) == seg(1:end-1)); e = d >= sd(ceil(0.9 * n));
s = pr(e(pr)); if isempty(s), return, end
c = mean(e(s + 1)) / mean(e);
end
function h = holm(p)
[ps, ix] = sort(p); m = numel(p); h = zeros(1, m); r = 0;
for k = 1:m, r = max(r, min(1, (m - k + 1) * ps(k))); h(ix(k)) = r; end
end
function p = tp(o, n), n = n(isfinite(n)); p = min(1, 2 * min(1 + sum(n >= o), 1 + sum(n <= o)) / (numel(n) + 1)); end
function v = prc(x, q), x = sort(x(isfinite(x))); v = x(max(1, round(q / 100 * numel(x)))); end
