L = load('C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\homeo_run\homeo_results.mat', 'flies');
P = L.flies.per0DD; nF = numel(P); K = 4;
nul = ["P","S12","W2"];
% ---------------- calibration ----------------
[on, d] = getb(P, 24);
sc = ["a stationary","b trend -2%/d","c acclimation tau 18h","d fly trends","e 3-h state (power)"];
rej = zeros(5, 3);
for k = 1:5
    for g = 1:100
        rs = RandStream('mrg32k3a', 'Seed', 5000 + 100 * k + g); sl = 0.1 * randn(rs, nF, 1); D = cell(1, nF);
        for f = 1:nF
            t = on{f} / 86400; tm = mean(t); mu = mean(d{f}) - 300;
            switch k
                case 1, m = mu * ones(size(t));
                case 2, m = mu * exp(log(0.98) * (t - tm));
                case 3, m = mu * (1 + 0.5 * exp(-t / 0.75));
                case 4, m = mu * exp(sl(f) * (t - tm));
                case 5, b = floor(on{f} / 10800); z = randn(rs, max(b) + 1, 1); m = mu * exp(0.3 * z(b + 1) - 0.045);
            end
            D{f} = 300 + m .* -log(rand(rs, size(t)));
        end
        R = runset(on, D, 199, 30000 + 100 * k + g, 1);
        rej(k, :) = rej(k, :) + (R.p < 0.05);
    end
    fprintf('CAL %-22s rejection P %.2f | S12 %.2f | W2 %.2f\n', sc(k), rej(k, :) / 100);
end
okP = rej(1, 1) / 100 <= 0.095; okS = all(rej(1:4, 2) / 100 <= 0.095); okW = all(rej(1:4, 3) / 100 <= 0.095);
fprintf('CALIBRATION size ok: P %d, S12 %d, W2 %d\n', okP, okS, okW);
save('calL.mat','rej','okP','okS','okW');

% ================= functions =================
function [on, d] = getb(P, cut)
on = cell(1, numel(P)); d = on;
for f = 1:numel(P)
    o = P(f).Onset(P(f).IsSleep); x = P(f).Dur(P(f).IsSleep); k = o >= cut * 3600;
    on{f} = o(k); d{f} = x(k);
end
end

function R = runset(on, D, nR, seed, K)
nF = numel(on); obs = NaN(nF, K, 2); n3 = NaN(nF, nR, 2, 3); nk = NaN(nF, nR, K); pf = NaN(nF, 3);
for f = 1:nF
    o = on{f}(:); d = D{f}(:); seg = floor(o / 43200);
    for lag = 1:K, obs(f, lag, :) = permute(cs(d, seg, lag), [1 3 2]); end
    G = {ones(size(d)), seg, floor(o / 7200)};
    rs = RandStream('mrg32k3a', 'Seed', seed); rs.Substream = f;
    for j = 1:3
        for r = 1:nR
            dp = pw(d, G{j}, rs);
            n3(f, r, :, j) = permute(cs(dp, seg, 1), [1 3 2]);
            if j == 2 && K > 1
                for lag = 1:K, x = cs(dp, seg, lag); nk(f, r, lag) = x(1); end
            end
        end
        if isfinite(obs(f, 1, 1)), pf(f, j) = tp(obs(f, 1, 1), n3(f, :, 1, j)); end
    end
end
R.obs = obs; R.G = squeeze(mean(obs(:, 1, :), 1, 'omitnan')).';
R.Gn = squeeze(mean(n3, 1, 'omitnan'));                % nR x 2 x 3
R.pt = NaN(2, 3);
for t = 1:2, for j = 1:3, R.pt(t, j) = tp(R.G(t), R.Gn(:, t, j)); end, end
R.p = R.pt(1, :);
R.nsig = sum(pf < 0.05, 1); R.nf = sum(isfinite(pf(:, 1)));
if K > 1
    R.gk = squeeze(mean(obs(:, :, 1), 1, 'omitnan')); gn = squeeze(mean(nk, 1, 'omitnan'));
    R.mk = mean(gn, 1); ex = R.gk - R.mk; R.D = ex(1) - mean(ex(2:K));
    en = gn - R.mk; R.Dn = en(:, 1) - mean(en(:, 2:K), 2); R.pD = tp(R.D, R.Dn);
end
end

function x = pw(x, g, rs)
[~, ~, gi] = unique(g);
for u = 1:max(gi)
    i = find(gi == u); if numel(i) > 1, x(i) = x(i(randperm(rs, numel(i)))); end
end
end
function c = cs(d, seg, lag)
c = NaN(1, 2); n = numel(d); if n < 10, return, end
sd = sort(d); pr = find(seg(1+lag:end) == seg(1:end-lag));
E = {d >= sd(ceil(0.9 * n)), d <= sd(ceil(0.1 * n))};
for k = 1:2
    e = E{k}; s = pr(e(pr)); if isempty(s), continue, end
    c(k) = mean(e(s + lag)) / mean(e);
end
end
function h = holm(p)
[ps, ix] = sort(p); m = numel(p); h = zeros(1, m); r = 0;
for k = 1:m, r = max(r, min(1, (m - k + 1) * ps(k))); h(ix(k)) = r; end
end
function p = tp(o, n), n = n(isfinite(n)); p = min(1, 2 * min(1 + sum(n >= o), 1 + sum(n <= o)) / (numel(n) + 1)); end
function v = prc(x, q), x = sort(x(isfinite(x))); v = x(max(1, round(q / 100 * numel(x)))); end
