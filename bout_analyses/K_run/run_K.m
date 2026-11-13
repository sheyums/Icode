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
        o = fluct(Y, T, S12, 199, 7, true); tr = trend(Y, T, 2000, 7);
        rej(k, :) = rej(k, :) + [o.p < 0.05, tr.p < 0.05];
    end
    fprintf('CAL %-24s fluctuation rejection %.2f | trend rejection %.2f\n', sc(k), rej(k, 1) / 100, rej(k, 2) / 100);
end
okCal = all(rej(1:4, 1) / 100 <= 0.095);
fprintf('CALIBRATION fluctuation test size ok in a-d: %d\n', okCal);
% ---------------- data: per0 DD ----------------
Y = cellfun(@log, d, 'UniformOutput', false);
tr = trend(Y, T, 10000, 3);
fl = fluct(Y, T, S12, 1000, 3, true);
raw = fluct(Y, T, S12, 1000, 3, false);
fprintf('\nper0 DD TREND: mean slope %.4f /day (%+.1f%%/day); flies with positive slope %d/%d; p %.4f\n', ...
    tr.mean, 100 * (exp(tr.mean) - 1), sum(tr.b > 0), nF, tr.p);
fprintf('per0 DD FLUCTUATION: group mean logRatio %.3f (null 95%% upper %.3f) p %.4f; flies logRatio>0 %d/%d\n', ...
    fl.G, prc(fl.Gn, 97.5), fl.p, sum(fl.lr > 0), nF);
fprintf('  no detrending: logRatio %.3f p %.4f; excess variance (V/nullV - 1): raw %.3f, detrended %.3f -> share removed by trends %.2f\n', ...
    raw.G, raw.p, mean(raw.ratio, 'omitnan') - 1, mean(fl.ratio, 'omitnan') - 1, 1 - (mean(fl.ratio, 'omitnan') - 1) / (mean(raw.ratio, 'omitnan') - 1));
fprintf('  lag-1 correlation of successive 12-h levels (detrended): obs %.3f, null mean %.3f, p(greater) %.4f\n', ...
    fl.lag1, mean(fl.lag1n), (1 + sum(fl.lag1n >= fl.lag1)) / (numel(fl.lag1n) + 1));
h = holm([tr.p, fl.p]); fprintf('PRIMARY Holm: trend %.4f, fluctuation %.4f (calibration ok %d)\n', h, okCal);
recD = cellfun(@max, T);
lbl = ["< 4 days", ">= 5 days"]; grps = {recD < 4, recD >= 5};
for q = 1:2
    gi = grps{q}; Gs = mean(fl.lr(gi), 'omitnan'); Gn = mean(fl.lrn(gi, :), 1, 'omitnan');
    fprintf('  record %s: %d flies, logRatio %.3f p %.4f\n', lbl(q), sum(gi), Gs, (1 + sum(Gn >= Gs)) / (numel(Gn) + 1));
end
k1 = cellfun(@(t) t >= 1, T, 'UniformOutput', false);
sub = @(A) cellfun(@(y, k) y(k), A, k1, 'UniformOutput', false);
tk = trend(sub(Y), sub(T), 10000, 3); fk = fluct(sub(Y), sub(T), sub(S12), 1000, 3, true);
fprintf('  first 24 h dropped: trend %+.1f%%/day p %.4f; fluctuation logRatio %.3f p %.4f\n', 100 * (exp(tk.mean) - 1), tk.p, fk.G, fk.p);
S4 = cellfun(@(o) floor(o / 14400), on, 'UniformOutput', false);
f4 = fluct(Y, T, S4, 1000, 3, true); fprintf('  4-h segments: fluctuation logRatio %.3f p %.4f\n', f4.G, f4.p);
[c0, cn] = clus(fl.res, S12, 1000, 3);
fprintf('  long-bout clustering C on detrended durations: %.3f vs shuffle %.3f p %.4f (raw C was 1.340)\n', c0, mean(cn), tp(c0, cn));
% ---------------- CS LD descriptive ----------------
C = L.flies.CSLD; nm = ["day","night"];
for p = 1:2
    Y2 = {}; T2 = {}; S2 = {};
    for f = 1:numel(C)
        o = C(f).Onset(C(f).IsSleep); dd = C(f).Dur(C(f).IsSleep); zt = mod(o, 86400) / 3600;
        k = (zt >= 12) + 1 == p; Y2{end+1} = log(dd(k)); T2{end+1} = o(k) / 86400; S2{end+1} = floor(o(k) / 86400); %#ok<SAGROW>
    end
    t2 = trend(Y2, T2, 10000, 3); f2 = fluct(Y2, T2, S2, 1000, 3, true); r2 = fluct(Y2, T2, S2, 1000, 3, false);
    [c2, c2n] = clus(f2.res, S2, 1000, 3);
    fprintf('CS LD %s: trend %+.1f%%/day p %.4f; fluctuation (day-to-day) logRatio %.3f p %.4f; raw %.3f p %.4f; C detrended %.3f vs %.3f p %.4f\n', ...
        nm(p), 100 * (exp(t2.mean) - 1), t2.p, f2.G, f2.p, r2.G, r2.p, c2, mean(c2n), tp(c2, c2n));
end
save('K_results.mat', 'tr', 'fl', 'raw', 'rej', 'tk', 'fk', 'f4', 'c0', 'cn', 'T', 'S12', '-v7.3');

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
