T = load(fullfile('..', 'N_run', 'tim0_flies.mat'), 'tim0'); L.flies.tim0DD = T.tim0;
sets = {"tim0DD", 24};
G = struct();
for k = 1:size(sets, 1)
    F = L.flies.(sets{k, 1}); X = cell(1, numel(F));
    for f = 1:numel(F), X{f} = prep(F(f), sets{k, 2}); end
    G.(sets{k, 1}) = X;
end
% ---------------- calibration ----------------
scn = ["a independent","b shared 2-h state","c sleep decline 2%/d","d homeostat (power)"];
for k = 1:size(sets, 1)
    X = G.(sets{k, 1}); nF = numel(X); rej = zeros(4, 3);
    for sc = 1:4
        for g = 1:100
            rs = RandStream('mrg32k3a', 'Seed', 9000 + 100 * sc + g);
            Y = cellfun(@(x) simfly(x, sc, rs), X, 'UniformOutput', false);
            R = runset(Y, 199, 20000 + 100 * sc + g, ["P","W2"]);
            pD = flipp(R.obs(:, 4), 2000, 20000 + 100 * sc + g);
            rej(sc, :) = rej(sc, :) + ([R.p(1), R.p(2), pD] < 0.05);
        end
        fprintf('CAL %s %-22s rejection R vs P %.2f | R vs W2 %.2f | Dtail %.2f\n', sets{k, 1}, scn(sc), rej(sc, :) / 100);
    end
    ok.(sets{k, 1}) = [rej(1, 1), max(rej(1:3, 2)), max(rej(1:3, 3))] / 100 <= 0.095;
    fprintf('CALIBRATION %s size ok: R vs P %d, R vs W2 %d, Dtail %d\n', sets{k, 1}, ok.(sets{k, 1}));
end
% ---------------- data ----------------
prim = []; lab = strings(0); out = struct();
for k = 1:size(sets, 1)
    nm = sets{k, 1}; X = G.(nm);
    nulls = ["P","S12","W2"]; if nm == "CSLD", nulls = [nulls "Z"]; end
    R = runset(X, 1000, 5, nulls);
    npre = cellfun(@(x) numel(x.pre), X); npost = cellfun(@(x) numel(x.post), X);
    fprintf('\n#### %s: %d flies; pairs per fly median pre %g, post %g\n', nm, numel(X), median(npre), median(npost));
    fprintf('  coupling R observed %.4f |', R.G(1));
    for j = 1:numel(nulls), fprintf(' %s %.4f p %.4f |', nulls(j), mean(R.Gn(:, 1, j)), R.p(j)); end
    fprintf('\n  r_pre %.4f, r_post %.4f; Delta (pre - post) %.4f sign-flip p %.4f\n', R.G(2), R.G(3), R.G(2) - R.G(3), flipp(R.obs(:, 2) - R.obs(:, 3), 10000, 5));
    pD = flipp(R.obs(:, 4), 10000, 5); pD2 = flipp(R.obs(:, 5), 10000, 5);
    fprintf('  Dtail (long wake: sleep after vs before) %.4f [Cpre %.3f, Cpost %.3f] sign-flip p %.4f; flies > 0: %d/%d\n', ...
        R.G(4), mean(R.obs(:, 6), 'omitnan'), mean(R.obs(:, 7), 'omitnan'), pD, sum(R.obs(:, 4) > 0), sum(isfinite(R.obs(:, 4))));
    fprintf('  reverse tail (long sleep: wake after vs before) %.4f sign-flip p %.4f\n', R.G(5), pD2);
    fprintf('  r_pre per null:'); for j = 1:numel(nulls), fprintf(' %s %.4f p %.4f;', nulls(j), mean(R.Gn(:, 2, j)), tp(R.G(2), R.Gn(:, 2, j))); end; fprintf('\n');
    jW = find(nulls == "W2");
    if nm == "tim0DD"
        prim = [prim, R.p(1), R.p(jW), pD]; lab = [lab, "tim0 R vs P", "tim0 R vs W2", "tim0 Dtail"]; %#ok<AGROW>
    else
        prim = [prim, R.p(1), R.p(jW)]; lab = [lab, "CSLD R vs P", "CSLD R vs W2"]; %#ok<AGROW>
    end
    out.(nm) = R;
end
h = holm(prim);
fprintf('\nPRIMARY Holm: '); for i = 1:numel(prim), fprintf('%s p %.4f Holm %.4f; ', lab(i), prim(i), h(i)); end; fprintf('\n');
fprintf('Calibration flags (R vs P, R vs W2, Dtail): tim0 %s\n', mat2str(ok.tim0DD));
save('OM_results.mat', 'out', 'prim', 'h', 'lab', 'ok', '-v7.3');

% ================= functions =================
function x = prep(fl, cutH)
o = fl.Onset(:); d = fl.Dur(:); s = logical(fl.IsSleep(:));
k = o >= cutH * 3600; o = o(k); d = d(k); s = s(k);
seg = floor(o / 43200);
adj = o(2:end) == o(1:end-1) + d(1:end-1) & seg(2:end) == seg(1:end-1);
x.o = o; x.d = d; x.s = s; x.seg = seg;
x.pre = find(adj & ~s(1:end-1) & s(2:end));
x.post = find(adj & s(1:end-1) & ~s(2:end));
grp.P = zeros(size(o)); grp.S12 = seg; grp.W2 = floor(o / 7200); grp.Z = floor(mod(o, 86400) / 3600);
for nm = string(fieldnames(grp)).'
    sets = {};
    for t = [true false]
        [~, ~, gi] = unique(grp.(nm)(s == t));
        ix = find(s == t);
        for u = 1:max(gi), m = ix(gi == u); if numel(m) > 1, sets{end+1} = m; end, end %#ok<AGROW>
    end
    x.g.(nm) = sets;
end
end

function v = stats(d, x)
s = x.s; v = NaN(1, 7);
rk = NaN(size(d)); rk(s) = rankavg(d(s)); rk(~s) = rankavg(d(~s));
a = x.pre; b = x.post;
if numel(a) >= 5, c = corrcoef(rk(a), rk(a + 1)); v(2) = c(2); end
if numel(b) >= 5, c = corrcoef(rk(b), rk(b + 1)); v(3) = c(2); end
v(1) = mean(v(2:3));
ds = sort(d(s)); dw = sort(d(~s));
if numel(ds) < 10 || numel(dw) < 10, return, end
lS = s & d >= ds(ceil(0.9 * numel(ds))); lW = ~s & d >= dw(ceil(0.9 * numel(dw)));
pS = mean(lS(s)); pW = mean(lW(~s));
i1 = a(lW(a)); i2 = b(lW(b + 1));
if ~isempty(i1) && ~isempty(i2)
    v(6) = mean(lS(i1 + 1)) / pS; v(7) = mean(lS(i2)) / pS; v(4) = v(6) - v(7);
end
j1 = b(lS(b)); j2 = a(lS(a + 1));                    % long sleep: wake after (post) vs before (pre)
if ~isempty(j1) && ~isempty(j2), v(5) = mean(lW(j1 + 1)) / pW - mean(lW(j2)) / pW; end
end

function R = runset(X, nR, seed, nulls)
nF = numel(X); obs = NaN(nF, 7); nn = NaN(nF, nR, 3, numel(nulls));
for f = 1:nF
    x = X{f}; obs(f, :) = stats(x.d, x);
    rs = RandStream('mrg32k3a', 'Seed', seed); rs.Substream = f;
    for j = 1:numel(nulls)
        sets = x.g.(nulls(j));
        for r = 1:nR
            idx = (1:numel(x.d)).';
            for q = 1:numel(sets), m = sets{q}; idx(m) = m(randperm(rs, numel(m))); end
            v = stats(x.d(idx), x); nn(f, r, :, j) = permute(v(1:3), [1 3 2]);
        end
    end
end
R.obs = obs; R.G = mean(obs, 1, 'omitnan');
R.Gn = squeeze(mean(nn, 1, 'omitnan'));
if numel(nulls) == 1, R.Gn = reshape(R.Gn, nR, 3, 1); end
R.p = arrayfun(@(j) tp(R.G(1), R.Gn(:, 1, j)), 1:numel(nulls));
R.nulls = nulls;
end

function y = simfly(x, sc, rs)
y = x; s = x.s; d = x.d;
ds = d(s); dw = d(~s);
d(s) = ds(randi(rs, numel(ds), sum(s), 1)); d(~s) = dw(randi(rs, numel(dw), sum(~s), 1));
switch sc
    case 2
        w = floor(x.o / 7200); [~, ~, wi] = unique(w); z = 0.3 * randn(rs, max(wi), 1); d = d .* exp(z(wi));
    case 3
        t = x.o / 86400; d(s) = d(s) .* 0.98 .^ (t(s) - mean(t));
    case 4
        mw = median(d(~s)); a = x.pre; d(a + 1) = d(a + 1) .* (d(a) / mw) .^ 0.15;
end
y.d = d;
end

function r = rankavg(v)
[u, ~, j] = unique(v); c = accumarray(j, 1); cs = cumsum(c);
avg = cs - (c - 1) / 2; r = avg(j); %#ok<NASGU>
r = reshape(avg(j), size(v));
end
function p = flipp(v, nR, seed)
v = v(isfinite(v)); rs = RandStream('mrg32k3a', 'Seed', seed);
s = 2 * (rand(rs, numel(v), nR) > 0.5) - 1; n = mean(v .* s, 1); o = mean(v);
p = (1 + sum(abs(n) >= abs(o))) / (nR + 1);
end
function h = holm(p)
[ps, ix] = sort(p); m = numel(p); h = zeros(1, m); r = 0;
for k = 1:m, r = max(r, min(1, (m - k + 1) * ps(k))); h(ix(k)) = r; end
end
function p = tp(o, n), n = n(isfinite(n)); p = min(1, 2 * min(1 + sum(n >= o), 1 + sum(n <= o)) / (numel(n) + 1)); end
