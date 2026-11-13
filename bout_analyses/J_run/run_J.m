L = load('C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\homeo_run\homeo_results.mat', 'flies');
nR = 1000; K = 4; bws = [2 1 3];
cases = {{"CSLD", 1, "CSLD_day"}, {"CSLD", 2, "CSLD_night"}, {"per0DD", 1, "per0DD"}};
res = struct(); pPrim = NaN(1, 2);
for c = 1:3
    nm = cases{c}{1}; p = cases{c}{2}; pn = cases{c}{3}; F = L.flies.(nm); nF = numel(F);
    obs = NaN(nF, K, 2); n4 = NaN(nF, nR, 2, 3); n1 = NaN(nF, nR, K); occ = [];
    for f = 1:nF
        fl = F(f); on = fl.Onset(fl.IsSleep); d = fl.Dur(fl.IsSleep);
        zt = mod(on, 86400) / 3600; dy = floor(on / 86400);
        if nm == "CSLD", ph = double(zt >= 12) + 1; seg = 2 * dy + ph; else, ph = ones(size(on)); seg = floor(on / 43200); end
        k = ph == p; d = d(k); seg = seg(k); zt = zt(k);
        for lag = 1:K, obs(f, lag, :) = permute(cstat(d, seg, lag), [1 3 2]); end
        G4 = arrayfun(@(b) {seg * 100 + floor(zt / b)}, bws);
        [~, ~, gi] = unique(G4{1}); cn = accumarray(gi, 1); occ = [occ; cn(gi)]; %#ok<AGROW>
        rs = RandStream('mrg32k3a', 'Seed', 2); rs.Substream = f;
        for j = 1:3
            for r = 1:nR, n4(f, r, :, j) = permute(cstat(pw(d, G4{j}, rs), seg, 1), [1 3 2]); end
        end
        for r = 1:nR
            dp = pw(d, floor(zt), rs);
            for lag = 1:K, x = cstat(dp, seg, lag); n1(f, r, lag) = x(1); end
        end
    end
    fprintf('\n#### %s: %d flies; day x 2-h bins: bouts alone %.2f, in bins >= 3 %.2f\n', pn, nF, mean(occ == 1), mean(occ >= 3));
    tn = ["upper","lower"];
    for t = 1:2
        go = mean(obs(:, 1, t), 'omitnan'); fprintf('  %s 10%%: observed %.3f |', tn(t), go);
        for j = 1:3
            g = squeeze(mean(n4(:, :, t, j), 1, 'omitnan')); pj = tp(go, g);
            fprintf(' N4 %dh %.3f p %.4f |', bws(j), mean(g), pj);
            if c == 1 && t == 1 && j == 1, pPrim(1) = pj; end
        end
        fprintf('\n');
    end
    gobs = squeeze(mean(obs(:, :, 1), 1, 'omitnan'));            % 1 x K
    gn1 = squeeze(mean(n1, 1, 'omitnan'));                         % nR x K
    m1 = mean(gn1, 1); ex = gobs - m1; D = ex(1) - mean(ex(2:K));
    exn = gn1 - m1; Dn = exn(:, 1) - mean(exn(:, 2:K), 2); pD = tp(D, Dn);
    if c == 1, pPrim(2) = pD; end
    fprintf('  upper lag profile: observed C_k %s; N1 mean %s; excess %s; Dlag %.3f (N1 95%% %.3f..%.3f) p %.4f\n', ...
        mat2str(round(gobs, 3)), mat2str(round(m1, 3)), mat2str(round(ex, 3)), D, prc(Dn, 2.5), prc(Dn, 97.5), pD);
    res.(pn) = struct('obs', obs, 'n4', n4, 'n1', n1, 'occ', occ, 'D', D, 'Dn', Dn);
end
h = holm(pPrim);
fprintf('\nPRIMARY CS LD day upper: vs N4 2h p %.4f Holm %.4f; Dlag vs N1 p %.4f Holm %.4f\n', pPrim(1), h(1), pPrim(2), h(2));
save(fullfile(fileparts(mfilename('fullpath')), 'J_results.mat'), 'res', 'pPrim', 'h', '-v7.3');

function h = holm(p)
[ps, ix] = sort(p); m = numel(p); h = zeros(1, m); r = 0;
for k = 1:m, r = max(r, min(1, (m - k + 1) * ps(k))); h(ix(k)) = r; end
end
function x = pw(x, g, rs)
[~, ~, gi] = unique(g);
for u = 1:max(gi)
    i = find(gi == u); if numel(i) > 1, x(i) = x(i(randperm(rs, numel(i)))); end
end
end
function c = cstat(d, seg, lag)
c = NaN(1, 2); n = numel(d); if n < 10, return, end
sd = sort(d); pr = find(seg(1+lag:end) == seg(1:end-lag));
E = {d >= sd(ceil(0.9 * n)), d <= sd(ceil(0.1 * n))};
for k = 1:2
    e = E{k}; s = pr(e(pr)); if isempty(s), continue, end
    c(k) = mean(e(s + lag)) / mean(e);
end
end
function p = tp(o, n)
n = n(isfinite(n)); p = min(1, 2 * min(1 + sum(n >= o), 1 + sum(n <= o)) / (numel(n) + 1));
end
function v = prc(x, q), x = sort(x(isfinite(x))); v = x(max(1, round(q / 100 * numel(x)))); end
