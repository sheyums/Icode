T = load(fullfile('..', 'N_run', 'tim0_flies.mat'), 'tim0');
P = T.tim0; nF = numel(P); K = 4;
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
% ---------------- data ----------------
for cut = 24
    [on, d] = getb(P, cut);
    R = runset(on, d, 1000, 4, K);
    fprintf('\n#### tim0 DD from %d h: %d flies, median bouts %g\n', cut, sum(isfinite(R.obs(:, 1, 1))), median(cellfun(@numel, d)));
    tn = ["upper","lower"];
    for t = 1:2
        fprintf('  %s 10%%: observed %.3f |', tn(t), R.G(t));
        for j = 1:3, fprintf(' %s %.3f p %.4f |', nul(j), mean(R.Gn(:, t, j)), R.pt(t, j)); end
        fprintf('\n');
    end
    fprintf('  upper lag profile C_k %s; S12 mean %s; excess %s; Dlag %.3f (S12 95%% %.3f..%.3f) p %.4f\n', ...
        mat2str(round(R.gk, 3)), mat2str(round(R.mk, 3)), mat2str(round(R.gk - R.mk, 3)), R.D, prc(R.Dn, 2.5), prc(R.Dn, 97.5), R.pD);
    fprintf('  per fly upper: outside S12 range %d/%d, outside W2 range %d/%d\n', R.nsig(2), R.nf, R.nsig(3), R.nf);
    if cut == 24
        h = holm(R.p);
        fprintf('PRIMARY (from 24 h, upper): P p %.4f Holm %.4f | S12 p %.4f Holm %.4f | W2 p %.4f Holm %.4f (calibration ok P %d S12 %d W2 %d)\n', ...
            R.p(1), h(1), R.p(2), h(2), R.p(3), h(3), okP, okS, okW);
        R24 = R;
    end
end
save('OL_results.mat', 'R24', 'rej', '-v7.3');

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
